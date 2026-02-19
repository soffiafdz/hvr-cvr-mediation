#!/usr/bin/env Rscript

# =============================================================================
# 01_preprocess_volumes.R
# =============================================================================
# Preprocess ADNI hippocampal and ventricular volumes with head-size adjustment.
#
# This script loads raw segmentation volumes and applies three head-size
# adjustment methods (NON, PRP, RES) to create the input data for z-score
# calculation in script 02.
#
# Inputs (from data/raw/adni/):
#   - vols_hcvc_adni.csv: HC and VC volumes in stereotaxic (stx) space (mm³)
#   - adni_icc.csv: Intracranial cavity volume with session metadata
#   - adni_scalefactor.csv: SCALEFACTOR to convert stx → native space
#
# Outputs:
#   - data/derivatives/adni_hc-hvr_adj.rds: Head-size adjusted volumes (long format)
#
# Head-size adjustment methods:
#   - NON: No adjustment (raw volumes in cc)
#   - PRP: Proportional (volume / ICC)
#   - RES: Residual (volume - B*(ICC - mean_ICC_CN)), where B is the
#          regression coefficient from lm(volume ~ ICC) and mean_ICC_CN
#          is the mean ICC of cognitively normal subjects
#
# HVR (Hippocampus-to-Ventricle Ratio):
#   - Computed as HC / (HC + VC) for NON adjustment only
#   - PRP would be identical (head-size cancels in ratio), RES unused
#   - GAMLSS normative models use NON
# =============================================================================

# --- Setup ---
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(stringr)
})

# Source utilities
source(here("R/utils/logging.R"))
source(here("R/utils/config.R"))
source(here("R/utils/data_io.R"))
source(here("R/utils/validation.R"))

# Initialize
log_script_start("01_preprocess_volumes.R")
config <- load_config()

# --- Configuration ---
FORCE_REGENERATE <- get_script_setting(

"force_regenerate", "preprocess_volumes",
  default = FALSE
)

output.path <- get_data_path("derivatives", "adni_hchvr")

# Check if output exists and skip if not forcing
if (!FORCE_REGENERATE && file.exists(output.path)) {
  log_info("Output file exists and force_regenerate=FALSE: %s", output.path)
  log_info("Skipping volume preprocessing. Set force_regenerate=TRUE to rerun.")
  log_script_end("01_preprocess_volumes.R", success = TRUE)
  quit(status = 0)
}

# --- Load Input Data ---
log_section("Loading Input Data")

# Raw volumes
vols.path <- here("data/raw/adni/vols_hcvc_adni.csv")
check_files_exist(vols.path)
vols_raw.dt <- fread(vols.path)
log_info("Raw volumes loaded: %d observations", nrow(vols_raw.dt))

# ICC and session metadata
icc.path <- here("data/raw/adni/adni_icc.csv")
check_files_exist(icc.path)
icc.dt <- fread(icc.path)
icc.dt[, EXAMDATE := as.Date(EXAMDATE)]
setkey(icc.dt, PTID, EXAMDATE)
log_info("ICC data loaded: %d sessions", nrow(icc.dt))

# SCALEFACTOR for stx → native conversion
sf.path <- here("data/raw/adni/adni_scalefactor.csv")
check_files_exist(sf.path)
sf.dt <- fread(sf.path)
# Handle both VISIT (YYYYMMDD) and EXAMDATE (YYYY-MM-DD) formats
if ("VISIT" %in% names(sf.dt)) {
  sf.dt[, EXAMDATE := as.Date(as.character(VISIT), format = "%Y%m%d")]
  sf.dt[, VISIT := NULL]
} else {
 sf.dt[, EXAMDATE := as.Date(EXAMDATE)]
}
setkey(sf.dt, PTID, EXAMDATE)
log_info("SCALEFACTOR data loaded: %d sessions", nrow(sf.dt))

# Log derived SCALEFACTOR sessions (these were computed from backup, not original pipeline)
if ("DERIVED" %in% names(sf.dt)) {
  n_derived <- sum(sf.dt$DERIVED, na.rm = TRUE)
  if (n_derived > 0) {
    log_warn(
      "%d sessions have SCALEFACTOR derived from backup (not original pipeline)",
      n_derived
    )
    derived_ptids <- unique(sf.dt[DERIVED == TRUE, PTID])
    log_info("Derived SCALEFACTOR subjects: %s", paste(derived_ptids, collapse = ", "))
  }
  # Remove DERIVED column - not propagated downstream
  sf.dt[, DERIVED := NULL]
}

# --- Parse Volume IDs ---
log_section("Parsing Volume Identifiers")

# ID format: stx2_XXX_S_XXXX_YYYYMMDD_t1_hcvc or stx2_XXX_S_XXXX_YYYY-MM-DD_t1_hcvc
# Extract PTID (XXX_S_XXXX) and EXAMDATE (YYYYMMDD or YYYY-MM-DD)
vols_raw.dt[, PTID := str_extract(ID, "[0-9]{3}_S_[0-9]{4}")]

# Handle both date formats
vols_raw.dt[, date_str := str_extract(ID, "[0-9]{8}|[0-9]{4}-[0-9]{2}-[0-9]{2}")]
vols_raw.dt[grepl("-", date_str), EXAMDATE := as.Date(date_str, format = "%Y-%m-%d")]
vols_raw.dt[!grepl("-", date_str) & !is.na(date_str),
            EXAMDATE := as.Date(date_str, format = "%Y%m%d")]
vols_raw.dt[, date_str := NULL]

# Validate parsing
n_parse_fail <- vols_raw.dt[is.na(PTID) | is.na(EXAMDATE), .N]
if (n_parse_fail > 0) {
  log_warn("Failed to parse %d IDs", n_parse_fail)
  vols_raw.dt <- vols_raw.dt[!is.na(PTID) & !is.na(EXAMDATE)]
}

log_info("Parsed %d volume observations", nrow(vols_raw.dt))

# --- Filter Zero Values ---
log_section("Filtering Invalid Observations")

n_before <- nrow(vols_raw.dt)

# Remove observations with 0 in any HC or LV measurement
# LHC/RHC = Left/Right Hippocampus, LCSF/RCSF = Left/Right Ventricles
vols_raw.dt <- vols_raw.dt[
  LHC > 0 & RHC > 0 & LCSF > 0 & RCSF > 0
]

n_after <- nrow(vols_raw.dt)
n_removed <- n_before - n_after

log_info(
  "Removed %d observations with zero values (%.1f%%)",
  n_removed, 100 * n_removed / n_before
)
log_info("Remaining observations: %d", n_after)

# --- Merge with ICC and SCALEFACTOR Data ---
log_section("Merging with ICC and SCALEFACTOR Data")

setkey(vols_raw.dt, PTID, EXAMDATE)

# First merge with ICC (metadata)
vols.dt <- merge(
  vols_raw.dt[, .(PTID, EXAMDATE, LHC, RHC, LCSF, RCSF)],
  icc.dt,
  by = c("PTID", "EXAMDATE")
)

n_merged_icc <- nrow(vols.dt)
n_no_icc <- n_after - n_merged_icc

if (n_no_icc > 0) {
  log_warn("%d observations dropped (no matching ICC data)", n_no_icc)
}

# Then merge with SCALEFACTOR (for stx → native conversion)
vols.dt <- merge(
  vols.dt,
  sf.dt[, .(PTID, EXAMDATE, SCALEFACTOR)],
  by = c("PTID", "EXAMDATE")
)

n_merged <- nrow(vols.dt)
n_no_sf <- n_merged_icc - n_merged

if (n_no_sf > 0) {
  log_warn("%d observations dropped (no matching SCALEFACTOR)", n_no_sf)
}

# Filter out zero SCALEFACTOR (invalid)
n_zero_sf <- vols.dt[SCALEFACTOR == 0, .N]
if (n_zero_sf > 0) {
  log_warn("Removing %d observations with SCALEFACTOR=0", n_zero_sf)
  vols.dt <- vols.dt[SCALEFACTOR != 0]
}

log_info("Merged dataset: %d observations", nrow(vols.dt))
log_info("Unique subjects: %d", uniqueN(vols.dt$PTID))

# --- Convert to Long Format ---
log_section("Converting to Long Format")
# Volumes are in stx mm³, convert to native cc:
#   native_cc = stx_mm3 / (SCALEFACTOR * 1000)
# Create separate rows for Left (L) and Right (R) sides

vols_long.dt <- rbind(
  vols.dt[, .(
    PTID, EXAMDATE, VISCODE, DX, ICC,
    SIDE = "L",
    HC = LHC / (SCALEFACTOR * 1000),
    VC = LCSF / (SCALEFACTOR * 1000)
  )],
  vols.dt[, .(
    PTID, EXAMDATE, VISCODE, DX, ICC,
    SIDE = "R",
    HC = RHC / (SCALEFACTOR * 1000),
    VC = RCSF / (SCALEFACTOR * 1000)
  )]
)

setkey(vols_long.dt, PTID, EXAMDATE, SIDE)
log_info("Long format: %d rows (2 sides per session)", nrow(vols_long.dt))

# --- Head-Size Adjustment ---
log_section("Applying Head-Size Adjustments")

# Compute regression coefficients for residual adjustment
# Using cognitively normal (CN) subjects as reference
mean_icc_cn <- vols_long.dt[DX == "CN", mean(ICC)]
log_info("Mean ICC for CN subjects: %.2f cc", mean_icc_cn)

# Compute B coefficients (slope of volume ~ ICC) for each ROI and side
b_coefs.dt <- vols_long.dt[, .(
  B_HC = coef(lm(HC ~ ICC))[2],
  B_VC = coef(lm(VC ~ ICC))[2]
), by = SIDE]

log_info("Regression coefficients (B) for residual adjustment:")
log_info("  Left HC: %.6f, Left VC: %.6f",
         b_coefs.dt[SIDE == "L", B_HC], b_coefs.dt[SIDE == "L", B_VC])
log_info("  Right HC: %.6f, Right VC: %.6f",
         b_coefs.dt[SIDE == "R", B_HC], b_coefs.dt[SIDE == "R", B_VC])

# Merge B coefficients
vols_long.dt <- merge(vols_long.dt, b_coefs.dt, by = "SIDE")

# Apply adjustments
vols_long.dt[, `:=`(
  # NON: No adjustment (raw cc)
  HC_NON = HC,
  VC_NON = VC,
  # PRP: Proportional (divide by ICC)
  HC_PRP = HC / ICC,
  VC_PRP = VC / ICC,
  # RES: Residual (subtract ICC effect relative to CN mean)
  HC_RES = HC - B_HC * (ICC - mean_icc_cn),
  VC_RES = VC - B_VC * (ICC - mean_icc_cn)
)]

# Remove intermediate columns
vols_long.dt[, c("HC", "VC", "B_HC", "B_VC") := NULL]

# --- Reshape to Final Format ---
log_section("Reshaping to Final Format")

# Melt to long format with ADJ column
hc_hvr.dt <- melt(
  vols_long.dt,
  id.vars = c("PTID", "EXAMDATE", "VISCODE", "DX", "ICC", "SIDE"),
  measure.vars = patterns("^HC_", "^VC_"),
  variable.name = "ADJ",
  value.name = c("HC", "VC")
)

# Fix ADJ labels
hc_hvr.dt[, ADJ := c("NON", "PRP", "RES")[ADJ]]

# Compute HVR for NON only (PRP is mathematically identical, RES unused)
# GAMLSS models use NON since ratio inherently normalizes head-size
hc_hvr.dt[ADJ == "NON", HVR := HC / (HC + VC)]

# Set column order to match expected output
setcolorder(hc_hvr.dt, c("PTID", "EXAMDATE", "VISCODE", "DX", "ICC",
                          "SIDE", "ADJ", "HC", "VC", "HVR"))
setkey(hc_hvr.dt, PTID, EXAMDATE)

# --- Validation ---
log_section("Validating Output")

# Check dimensions
n_subjects <- uniqueN(hc_hvr.dt$PTID)
n_sessions <- uniqueN(hc_hvr.dt[, .(PTID, EXAMDATE)])
n_rows <- nrow(hc_hvr.dt)
expected_rows <- n_sessions * 2 * 3  # 2 sides × 3 adjustments

log_info("Output dimensions:")
log_info("  Subjects: %d", n_subjects)
log_info("  Sessions: %d", n_sessions)
log_info("  Rows: %d (expected: %d)", n_rows, expected_rows)

if (n_rows != expected_rows) {
  log_warn("Row count mismatch! Check for missing data.")
}

# Check value ranges
log_info("Value ranges:")
log_info("  ICC: [%.1f, %.1f] cc", min(hc_hvr.dt$ICC), max(hc_hvr.dt$ICC))
log_info("  HC (NON): [%.3f, %.3f] cc",
         hc_hvr.dt[ADJ == "NON", min(HC)], hc_hvr.dt[ADJ == "NON", max(HC)])
log_info("  HVR (NON): [%.3f, %.3f]",
         hc_hvr.dt[ADJ == "NON", min(HVR, na.rm = TRUE)],
         hc_hvr.dt[ADJ == "NON", max(HVR, na.rm = TRUE)])

# --- Save Output ---
log_section("Saving Output")

ensure_directory(dirname(output.path))
write_rds_safe(
  hc_hvr.dt, output.path,
  description = "ADNI head-size adjusted HC/VC/HVR volumes"
)

# --- Cleanup ---
log_script_end("01_preprocess_volumes.R", success = TRUE)
