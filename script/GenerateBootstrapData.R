
# GenerateBootstrapData.R
# R equivalent of GenerateBootstrapData.do
# Runs 4 bootstrap models (1001 iterations: run=0 is baseline, runs 1-1000 are resamples)
#
# Output files (in data/output/bootstrap/):
#   bootstrap_noLag.csv         - pooled, no lags
#   bootstrap_richpoor.csv      - rich/poor split, no lags
#   bootstrap_5lag.csv          - pooled, 5 lags
#   bootstrap_richpoor_5lag.csv - rich/poor split, 5 lags
#
# Note: country-specific linear and quadratic time trends are absorbed via fixest's
#   new_id[time, time2] notation, equivalent to Stata's _yi_* / _y2_* columns.

rm(list = ls())

library(fixest)
library(data.table)

setwd("/Users/yl3699/Downloads/BurkeHsiangMiguel2015_Replication")
dir.create("data/output/bootstrap", recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# Load and prepare data
# ----------------------------------------------------------------------------
dta <- as.data.table(read.csv("data/input/GrowthClimateDataset.csv"))
dta[, temp  := UDel_temp_popweight]
dta[, temp2 := temp^2]
dta[, prec  := UDel_precip_popweight]
dta[, prec2 := UDel_precip_popweight_2]
dta[, time  := year - 1985]
dta[, time2 := time^2]
dta[, poor  := fifelse(is.na(GDPpctile_WDIppp), NA_integer_,
                       as.integer(GDPpctile_WDIppp < 50))]
# Baseline uses original iso_id as new_id
dta[, new_id := iso_id]

set.seed(8675309)

# ----------------------------------------------------------------------------
# Helper: cluster-resample countries with replacement.
# Each copy of a country gets a new unique integer new_id so that FEs and
# country-specific trends are estimated separately for each copy.
# ----------------------------------------------------------------------------
cluster_resample <- function(dt) {
  countries <- unique(dt$iso_id)
  n         <- length(countries)
  sampled   <- sample(countries, n, replace = TRUE)
  parts     <- lapply(seq_along(sampled), function(i) {
    rows <- dt[iso_id == sampled[i]]
    rows[, new_id := i]
    rows
  })
  rbindlist(parts)
}

# ----------------------------------------------------------------------------
# Helper: add lag variables L1-L5 for temp, temp2, prec, prec2.
# Operates on a copy; sorted by new_id then year before lagging.
# ----------------------------------------------------------------------------
add_lags <- function(dt) {
  d <- copy(dt)
  setorderv(d, c("new_id", "year"))
  for (v in c("temp", "temp2", "prec", "prec2")) {
    for (l in 1:5) {
      d[, (paste0(v, "_L", l)) := data.table::shift(get(v), l), by = "new_id"]
    }
  }
  d
}

# RHS variable names for the 5-lag model (L0 through L5 for each of 4 vars)
lag_rhs <- c("temp", "temp2", "prec", "prec2",
             paste0("temp_L",  1:5),
             paste0("temp2_L", 1:5),
             paste0("prec_L",  1:5),
             paste0("prec2_L", 1:5))


# ============================================================================
# MODEL 1: Pooled, no lags  ->  bootstrap_noLag.csv
# Columns: run, temp, temp2, prec, prec2
# ============================================================================
cat("Bootstrap model 1: pooled no-lag (1001 iterations)...\n")

run_nolag <- function(d) {
  fit <- feols(growthWDI ~ temp + temp2 + prec + prec2 |
                 new_id[time, time2] + year,
               data = d, warn = FALSE, notes = FALSE)
  b <- coef(fit)
  data.frame(temp  = unname(b["temp"]),  temp2 = unname(b["temp2"]),
             prec  = unname(b["prec"]),  prec2 = unname(b["prec2"]))
}

na_nolag <- data.frame(temp = NA_real_, temp2 = NA_real_,
                       prec = NA_real_, prec2 = NA_real_)

results <- vector("list", 1001)
results[[1]] <- cbind(run = 0, run_nolag(dta))
for (nn in seq_len(1000)) {
  if (nn %% 100 == 0) cat("  iteration", nn, "\n")
  results[[nn + 1]] <- tryCatch(
    cbind(run = nn, run_nolag(cluster_resample(dta))),
    error = function(e) cbind(run = nn, na_nolag)
  )
}
write.csv(rbindlist(results), "data/output/bootstrap/bootstrap_noLag.csv",
          row.names = FALSE)
cat("Wrote bootstrap_noLag.csv\n\n")


# ============================================================================
# MODEL 2: Rich/poor, no lags  ->  bootstrap_richpoor.csv
# Columns: run, temp, temppoor, temp2, temp2poor, prec, precpoor, prec2, prec2poor
#
# Note: temppoor = TOTAL poor coefficient (rich coef + differential), matching
#   Stata's 0.poor#c.var / 1.poor#c.var parameterisation.
# ============================================================================
cat("Bootstrap model 2: rich/poor no-lag (1001 iterations)...\n")

run_richpoor <- function(d) {
  tmp <- d[!is.na(poor)]
  tmp[, temp_poor  := temp  * poor]
  tmp[, temp2_poor := temp2 * poor]
  tmp[, prec_poor  := prec  * poor]
  tmp[, prec2_poor := prec2 * poor]
  fit <- feols(growthWDI ~ temp + temp2 + temp_poor + temp2_poor +
                 prec + prec2 + prec_poor + prec2_poor |
                 new_id[time, time2] + year,
               data = tmp, warn = FALSE, notes = FALSE)
  b <- coef(fit)
  data.frame(
    temp      = unname(b["temp"]),
    temppoor  = unname(b["temp"]  + b["temp_poor"]),   # total poor coefficient
    temp2     = unname(b["temp2"]),
    temp2poor = unname(b["temp2"] + b["temp2_poor"]),
    prec      = unname(b["prec"]),
    precpoor  = unname(b["prec"]  + b["prec_poor"]),
    prec2     = unname(b["prec2"]),
    prec2poor = unname(b["prec2"] + b["prec2_poor"])
  )
}

na_richpoor <- data.frame(
  temp = NA_real_, temppoor = NA_real_, temp2 = NA_real_, temp2poor = NA_real_,
  prec = NA_real_, precpoor = NA_real_, prec2 = NA_real_, prec2poor = NA_real_)

results <- vector("list", 1001)
results[[1]] <- cbind(run = 0, run_richpoor(dta))
for (nn in seq_len(1000)) {
  if (nn %% 100 == 0) cat("  iteration", nn, "\n")
  results[[nn + 1]] <- tryCatch(
    cbind(run = nn, run_richpoor(cluster_resample(dta))),
    error = function(e) cbind(run = nn, na_richpoor)
  )
}
write.csv(rbindlist(results), "data/output/bootstrap/bootstrap_richpoor.csv",
          row.names = FALSE)
cat("Wrote bootstrap_richpoor.csv\n\n")


# ============================================================================
# MODEL 3: Pooled, 5 lags  ->  bootstrap_5lag.csv
# Columns: run, temp, L1temp..L5temp, temp2, L1temp2..L5temp2, tlin, tsq
#   tlin = sum of L0-L5 linear temp coefficients
#   tsq  = sum of L0-L5 quadratic temp coefficients
# ============================================================================
cat("Bootstrap model 3: pooled 5-lag (1001 iterations)...\n")

fml_5lag <- as.formula(
  paste("growthWDI ~", paste(lag_rhs, collapse = " + "),
        "| new_id[time, time2] + year"))

run_5lag <- function(d) {
  d   <- add_lags(d)
  fit <- feols(fml_5lag, data = d, warn = FALSE, notes = FALSE)
  b   <- coef(fit)
  tlin <- sum(b[c("temp",  "temp_L1",  "temp_L2",  "temp_L3",  "temp_L4",  "temp_L5")],  na.rm = TRUE)
  tsq  <- sum(b[c("temp2", "temp2_L1", "temp2_L2", "temp2_L3", "temp2_L4", "temp2_L5")], na.rm = TRUE)
  data.frame(
    temp    = unname(b["temp"]),    L1temp  = unname(b["temp_L1"]),
    L2temp  = unname(b["temp_L2"]), L3temp  = unname(b["temp_L3"]),
    L4temp  = unname(b["temp_L4"]), L5temp  = unname(b["temp_L5"]),
    temp2   = unname(b["temp2"]),   L1temp2 = unname(b["temp2_L1"]),
    L2temp2 = unname(b["temp2_L2"]),L3temp2 = unname(b["temp2_L3"]),
    L4temp2 = unname(b["temp2_L4"]),L5temp2 = unname(b["temp2_L5"]),
    tlin = tlin, tsq = tsq
  )
}

na_5lag <- data.frame(
  temp = NA_real_,  L1temp = NA_real_,  L2temp = NA_real_,
  L3temp = NA_real_, L4temp = NA_real_, L5temp = NA_real_,
  temp2 = NA_real_, L1temp2 = NA_real_, L2temp2 = NA_real_,
  L3temp2 = NA_real_, L4temp2 = NA_real_, L5temp2 = NA_real_,
  tlin = NA_real_, tsq = NA_real_)

results <- vector("list", 1001)
results[[1]] <- cbind(run = 0, run_5lag(dta))
for (nn in seq_len(1000)) {
  if (nn %% 100 == 0) cat("  iteration", nn, "\n")
  results[[nn + 1]] <- tryCatch(
    cbind(run = nn, run_5lag(cluster_resample(dta))),
    error = function(e) cbind(run = nn, na_5lag)
  )
}
write.csv(rbindlist(results), "data/output/bootstrap/bootstrap_5lag.csv",
          row.names = FALSE)
cat("Wrote bootstrap_5lag.csv\n\n")


# ============================================================================
# MODEL 4: Rich/poor, 5 lags  ->  bootstrap_richpoor_5lag.csv
# Key outputs used by ComputeMainProjections.R:
#   tlin     = sum of rich linear temp coefficients (L0-L5)
#   tlinpoor = sum of TOTAL poor linear temp coefficients (L0-L5)
#   tsq      = sum of rich quadratic temp coefficients (L0-L5)
#   tsqpoor  = sum of TOTAL poor quadratic temp coefficients (L0-L5)
# ============================================================================
cat("Bootstrap model 4: rich/poor 5-lag (1001 iterations)...\n")

run_richpoor_5lag <- function(d) {
  tmp <- d[!is.na(poor)]
  tmp <- add_lags(tmp)

  # Create poor-interaction columns for all lag variables
  for (v in c("temp", "temp2", "prec", "prec2")) {
    tmp[, paste0(v, "_poor") := get(v) * poor]
    for (l in 1:5) {
      base_col <- paste0(v, "_L", l)
      tmp[, paste0(base_col, "_poor") := get(base_col) * poor]
    }
  }

  poor_vars <- c(
    paste0(c("temp", "temp2", "prec", "prec2"), "_poor"),
    paste0(rep(c("temp", "temp2", "prec", "prec2"), each = 5),
           "_L", rep(1:5, 4), "_poor"))

  fml <- as.formula(
    paste("growthWDI ~", paste(c(lag_rhs, poor_vars), collapse = " + "),
          "| new_id[time, time2] + year"))

  fit <- feols(fml, data = tmp, warn = FALSE, notes = FALSE)
  b   <- coef(fit)

  # Sums for rich countries
  tlin <- sum(b[c("temp",  "temp_L1",  "temp_L2",  "temp_L3",  "temp_L4",  "temp_L5")],  na.rm = TRUE)
  tsq  <- sum(b[c("temp2", "temp2_L1", "temp2_L2", "temp2_L3", "temp2_L4", "temp2_L5")], na.rm = TRUE)
  # Sums of poor differentials (to add to rich sums for total poor effect)
  tlin_diff <- sum(b[c("temp_poor", "temp_L1_poor", "temp_L2_poor",
                        "temp_L3_poor", "temp_L4_poor", "temp_L5_poor")],  na.rm = TRUE)
  tsq_diff  <- sum(b[c("temp2_poor", "temp2_L1_poor", "temp2_L2_poor",
                        "temp2_L3_poor", "temp2_L4_poor", "temp2_L5_poor")], na.rm = TRUE)

  data.frame(
    temp        = unname(b["temp"]),
    temppoor    = unname(b["temp"]    + b["temp_poor"]),
    L1temp      = unname(b["temp_L1"]),
    L1temppoor  = unname(b["temp_L1"] + b["temp_L1_poor"]),
    L2temp      = unname(b["temp_L2"]),
    L2temppoor  = unname(b["temp_L2"] + b["temp_L2_poor"]),
    L3temp      = unname(b["temp_L3"]),
    L3temppoor  = unname(b["temp_L3"] + b["temp_L3_poor"]),
    L4temp      = unname(b["temp_L4"]),
    L4temppoor  = unname(b["temp_L4"] + b["temp_L4_poor"]),
    L5temp      = unname(b["temp_L5"]),
    L5temppoor  = unname(b["temp_L5"] + b["temp_L5_poor"]),
    temp2       = unname(b["temp2"]),
    temp2poor   = unname(b["temp2"]    + b["temp2_poor"]),
    L1temp2     = unname(b["temp2_L1"]),
    L1temp2poor = unname(b["temp2_L1"] + b["temp2_L1_poor"]),
    L2temp2     = unname(b["temp2_L2"]),
    L2temp2poor = unname(b["temp2_L2"] + b["temp2_L2_poor"]),
    L3temp2     = unname(b["temp2_L3"]),
    L3temp2poor = unname(b["temp2_L3"] + b["temp2_L3_poor"]),
    L4temp2     = unname(b["temp2_L4"]),
    L4temp2poor = unname(b["temp2_L4"] + b["temp2_L4_poor"]),
    L5temp2     = unname(b["temp2_L5"]),
    L5temp2poor = unname(b["temp2_L5"] + b["temp2_L5_poor"]),
    tlin     = tlin,
    tlinpoor = tlin + tlin_diff,   # total poor effect = rich + differential sums
    tsq      = tsq,
    tsqpoor  = tsq  + tsq_diff
  )
}

na_rp5lag <- data.frame(
  temp = NA_real_,     temppoor = NA_real_,
  L1temp = NA_real_,   L1temppoor = NA_real_,
  L2temp = NA_real_,   L2temppoor = NA_real_,
  L3temp = NA_real_,   L3temppoor = NA_real_,
  L4temp = NA_real_,   L4temppoor = NA_real_,
  L5temp = NA_real_,   L5temppoor = NA_real_,
  temp2 = NA_real_,    temp2poor = NA_real_,
  L1temp2 = NA_real_,  L1temp2poor = NA_real_,
  L2temp2 = NA_real_,  L2temp2poor = NA_real_,
  L3temp2 = NA_real_,  L3temp2poor = NA_real_,
  L4temp2 = NA_real_,  L4temp2poor = NA_real_,
  L5temp2 = NA_real_,  L5temp2poor = NA_real_,
  tlin = NA_real_,    tlinpoor = NA_real_,
  tsq  = NA_real_,    tsqpoor  = NA_real_)

results <- vector("list", 1001)
results[[1]] <- cbind(run = 0, run_richpoor_5lag(dta))
for (nn in seq_len(1000)) {
  if (nn %% 100 == 0) cat("  iteration", nn, "\n")
  results[[nn + 1]] <- tryCatch(
    cbind(run = nn, run_richpoor_5lag(cluster_resample(dta))),
    error = function(e) cbind(run = nn, na_rp5lag)
  )
}
write.csv(rbindlist(results), "data/output/bootstrap/bootstrap_richpoor_5lag.csv",
          row.names = FALSE)
cat("Wrote bootstrap_richpoor_5lag.csv\n")

cat("\nGenerateBootstrapData.R complete.\n")
cat("Run ComputeMainProjections.R next.\n")
