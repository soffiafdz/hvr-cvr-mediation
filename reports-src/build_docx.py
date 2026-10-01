#!/usr/bin/env python3
"""Build a single Word file of the manuscript for Editorial Manager.

Alzheimer's & Dementia requires one .docx at initial submission (it is
parsed to auto-populate title/authors/abstract). Quarto's docx target
cannot render the TikZ CONSORT diagram or the gt LaTeX tables, so this
converts the already-rendered LaTeX instead. Nothing here touches the
.qmd, so the PDF build cannot regress.

Prerequisite: render the PDF first, which produces the .tex and figures.

    quarto render manuscript_ad_submission.qmd --to pdf
    python3 build_docx.py

Output: ../outputs/reports/manuscript_ad_submission.docx

Four pandoc-LaTeX incompatibilities are worked around below; each is
load-bearing, so read the comments before editing.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

SRC = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(SRC)
TEX = os.path.join(SRC, "manuscript_ad_submission.tex")
FIGDIR = os.path.join(SRC, "manuscript_ad_submission_files", "figure-pdf")
OUTDIR = os.path.join(ROOT, "outputs", "reports")
OUT = os.path.join(OUTDIR, "manuscript_ad_submission.docx")


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


# ---------------------------------------------------------------
# 1. Figures: Word cannot embed PDF images, so rasterise to PNG.
# ---------------------------------------------------------------
def figures_to_png():
    n = 0
    for f in sorted(os.listdir(FIGDIR)):
        if not f.endswith(".pdf"):
            continue
        base = f[:-4]
        run(["pdftocairo", "-png", "-r", "300", "-singlefile",
             os.path.join(FIGDIR, f), os.path.join(FIGDIR, base)])
        n += 1
    return n


# ---------------------------------------------------------------
# 2. CONSORT diagram exists only as TikZ; compile it standalone.
# ---------------------------------------------------------------
def build_consort(tex):
    m = re.search(r"\\begin\{tikzpicture\}.*?\\end\{tikzpicture\}", tex, re.S)
    if not m:
        sys.exit("CONSORT tikzpicture not found in .tex")
    tmp = tempfile.mkdtemp()
    stand = os.path.join(tmp, "consort.tex")
    with open(stand, "w", encoding="utf-8") as fh:
        fh.write(
            "\\documentclass[border=6pt]{standalone}\n"
            "\\usepackage{tikz}\n"
            "\\usetikzlibrary{shapes.geometric, arrows.meta, positioning}\n"
            "\\usepackage{fontspec}\n\\setmainfont{TeX Gyre Termes}\n"
            "\\begin{document}\n" + m.group(0) + "\n\\end{document}\n")
    run(["xelatex", "-interaction=nonstopmode", "consort.tex"], cwd=tmp)
    pdf = os.path.join(tmp, "consort.pdf")
    if not os.path.exists(pdf):
        sys.exit("CONSORT standalone failed to compile")
    run(["pdftocairo", "-png", "-r", "300", "-singlefile", pdf,
         os.path.join(tmp, "consort")])
    shutil.copy(os.path.join(tmp, "consort.png"),
                os.path.join(FIGDIR, "fig-consort.png"))
    shutil.rmtree(tmp, ignore_errors=True)


# ---------------------------------------------------------------
# 3. Reference doc: journal wants 12-pt Times, double spaced.
# ---------------------------------------------------------------
def build_reference_doc(path):
    base = path + ".base"
    with open(base, "wb") as fh:
        subprocess.run(["pandoc", "--print-default-data-file",
                        "reference.docx"], stdout=fh, check=True)
    zin = zipfile.ZipFile(base)
    st = zin.read("word/styles.xml").decode("utf-8")
    st = re.sub(r"<w:rFonts[^/]*/>",
                '<w:rFonts w:ascii="Times New Roman" '
                'w:hAnsi="Times New Roman" w:eastAsia="Times New Roman" '
                'w:cs="Times New Roman"/>', st, count=1)
    st = re.sub(r'<w:sz w:val="\d+"/>', '<w:sz w:val="24"/>', st)
    st = re.sub(r'<w:szCs w:val="\d+"/>', '<w:szCs w:val="24"/>', st)
    st = re.sub(r'<w:spacing([^/]*?)w:line="\d+"([^/]*?)/>',
                r'<w:spacing\1w:line="480"\2/>', st)
    if 'w:line="480"' not in st:
        st = st.replace("<w:pPrDefault>",
                        '<w:pPrDefault><w:pPr><w:spacing w:line="480" '
                        'w:lineRule="auto"/></w:pPr>', 1)
    zout = zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED)
    for item in zin.infolist():
        data = zin.read(item.filename)
        if item.filename == "word/styles.xml":
            data = st.encode("utf-8")
        zout.writestr(item, data)
    zout.close()
    zin.close()
    os.remove(base)


# ---------------------------------------------------------------
# 4. LaTeX fixups. Pandoc's reader silently DROPS a braced group after
#    \centering and \noindent - which loses the image in every figure
#    and the entire tabular in every table. These are not cosmetic.
# ---------------------------------------------------------------
def flatten_group(s, macro):
    """Turn `\\macro{ ... }` into `\\macro ...`, brace-matched."""
    key = "\\" + macro + "{"
    out, i, n = [], 0, 0
    while True:
        j = s.find(key, i)
        if j < 0:
            out.append(s[i:])
            return "".join(out), n
        out.append(s[i:j])
        out.append("\\" + macro + "\n")
        k, depth = j + len(key), 1
        while k < len(s) and depth:
            if s[k] == "{":
                depth += 1
            elif s[k] == "}":
                depth -= 1
            k += 1
        out.append(s[j + len(key):k - 1])
        i, n = k, n + 1


def preprocess(tex):
    log = {}
    tex, log["tikz->png"] = re.subn(
        r"\\begin\{tikzpicture\}.*?\\end\{tikzpicture\}",
        r"\\includegraphics[width=\\linewidth]"
        r"{manuscript_ad_submission_files/figure-pdf/fig-consort.png}",
        tex, flags=re.S)
    tex, log["pdf->png"] = re.subn(
        r"(figure-pdf/[A-Za-z0-9_-]+)\.pdf", r"\1.png", tex)
    # \pandocbounded is a macro pandoc WRITES but cannot READ.
    tex, log["pandocbounded"] = re.subn(
        r"\\pandocbounded\{(\\includegraphics\[[^\]]*\]\{[^}]*\})\}",
        r"\1", tex)
    tex, log["centering"] = flatten_group(tex, "centering")
    # \cmidrule leaks its argument as literal text into header rows.
    tex, log["cmidrule"] = re.subn(
        r"\\cmidrule(\([a-z]+\))?\{[0-9]+-[0-9]+\}\s*", "", tex)
    # Quarto wraps display math in \phantomsection\label{}{...}, which
    # pandoc's math parser rejects.
    tex, log["eq-open"] = re.subn(
        r"\\begin\{equation\}\\phantomsection\\label\{[^}]*\}\{",
        r"\\begin{equation}", tex)
    tex, log["eq-close"] = re.subn(
        r"\}\\end\{equation\}", r"\\end{equation}", tex)
    for pat in [r"\\usepackage\[mathlines\]\{lineno\}\n?", r"\\linenumbers\n?",
                r"\\usepackage\{tikz\}\n?", r"\\usetikzlibrary\{[^}]*\}\n?",
                r"\\newunicodechar\{[^}]*\}\{[^}]*\}\n?",
                r"\\usepackage\{newunicodechar\}\n?", r"\\doublespacing\n?"]:
        tex = re.sub(pat, "", tex)
    # Unwrap the first-page abbreviation footnote into a paragraph.
    tex = re.sub(r"\\newcommand\\blfootnote\[1\]\{%.*?\n\}\n", "", tex,
                 flags=re.S)
    key = "\\blfootnote{"
    i = tex.find(key)
    if i >= 0:
        k, depth = i + len(key), 1
        while k < len(tex) and depth:
            if tex[k] == "{":
                depth += 1
            elif tex[k] == "}":
                depth -= 1
            k += 1
        inner = tex[i + len(key):k - 1].replace("\\footnotesize", "").strip()
        tex = tex[:i] + "\n\n" + inner + "\n\n" + tex[k:]
        log["abbrev-footnote"] = 1
    return tex, log


def main():
    if not os.path.exists(TEX):
        sys.exit("%s not found - render the PDF first" % TEX)
    os.makedirs(OUTDIR, exist_ok=True)
    tex = open(TEX, encoding="utf-8").read()

    print("figures rasterised :", figures_to_png())
    build_consort(tex)
    tex, log = preprocess(tex)
    for k, v in log.items():
        print("  %-16s %s" % (k, v))

    tmpdir = tempfile.mkdtemp()
    ref = os.path.join(tmpdir, "reference.docx")
    build_reference_doc(ref)
    pre = os.path.join(tmpdir, "for_docx.tex")
    open(pre, "w", encoding="utf-8").write(tex)

    r = run(["pandoc", pre, "-f", "latex", "-o", OUT,
             "--reference-doc=" + ref, "--resource-path=" + SRC])
    shutil.rmtree(tmpdir, ignore_errors=True)
    if r.returncode or not os.path.exists(OUT):
        sys.exit("pandoc failed:\n" + r.stderr[:2000])

    z = zipfile.ZipFile(OUT)
    doc = z.read("word/document.xml").decode("utf-8", "replace")
    runs = "".join(re.findall(r"<w:t[^>]*>(.*?)</w:t>", doc, re.S))
    leak = len(re.findall(r"\(lr\)|cmidrule|toprule|\\[a-zA-Z]{3,}", runs))
    print("\n%s (%.1f KB)" % (OUT, os.path.getsize(OUT) / 1024))
    print("  images %d | tables %d | equations %d | latex leakage %d"
          % (len([n for n in z.namelist() if n.startswith("word/media/")]),
             doc.count("<w:tbl>"), doc.count("<m:oMath"), leak))
    if leak:
        print("  WARNING: LaTeX leaked into the Word text - inspect before use")


if __name__ == "__main__":
    main()
