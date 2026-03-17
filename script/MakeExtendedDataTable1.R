
# MakeExtendedDataTable1.R
# R equivalent of MakeExtendedDataTable1.do
# Runs 11 robustness specifications for the baseline regression.
#
# Outputs:
#   data/output/Coefficients_robustness.csv  (b1, b2 for each spec)
#   figures/ExtendedDataFigs_Input/Table1.tex (LaTeX regression table)
#
# Requires: install.packages("modelsummary")

rm(list = ls())

library(fixest)
library(data.table)
library(modelsummary)

setwd("/Users/yl3699/Downloads/BurkeHsiangMiguel2015_Replication")
dir.create("data/output",                      recursive = TRUE, showWarnings = FALSE)
dir.create("figures/ExtendedDataFigs_Input",   recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# Load and prepare data
# ----------------------------------------------------------------------------
dta <- as.data.table(read.csv("data/input/GrowthClimateDataset.csv"))

dta[, temp    := UDel_temp_popweight]
dta[, temp2   := temp^2]
dta[, precip  := UDel_precip_popweight / 1000]  # rescale so coefficients are legible
dta[, precip2 := precip^2]
dta[, time    := year - 1985]
dta[, time2   := time^2]

# Lagged dependent variables for specs 9 and 10 (sorted by iso_id, year)
setorderv(dta, c("iso_id", "year"))
dta[, lag1_growth := shift(growthWDI, 1), by = iso_id]
dta[, lag2_growth := shift(growthWDI, 2), by = iso_id]
dta[, lag3_growth := shift(growthWDI, 3), by = iso_id]

# ----------------------------------------------------------------------------
# Run 11 robustness specifications
# Stata -> R translation:
#   reg ... i.year _yi_* _y2_* i.iso_id, cl(iso_id)
#     -> feols(... | iso_id[time, time2] + year, cluster = ~iso_id)
#   Continent-year FE (_cy_*):  year^continent
#   L.growthWDI:                lag1_growth (pre-computed above)
# ----------------------------------------------------------------------------

base_rhs <- "temp + temp2 + precip + precip2"

specs <- list(
  # 1. Full sample baseline
  feols(as.formula(paste("growthWDI ~", base_rhs,
                         "| iso_id[time, time2] + year")),
        cluster = ~iso_id, data = dta),

  # 2. Drop countries with fewer than 20 WDI observations
  feols(as.formula(paste("growthWDI ~", base_rhs,
                         "| iso_id[time, time2] + year")),
        cluster = ~iso_id, data = dta[wdinomiss >= 20]),

  # 3. Drop oil-exporting countries
  feols(as.formula(paste("growthWDI ~", base_rhs,
                         "| iso_id[time, time2] + year")),
        cluster = ~iso_id, data = dta[oil == 0]),

  # 4. Drop US and China
  feols(as.formula(paste("growthWDI ~", base_rhs,
                         "| iso_id[time, time2] + year")),
        cluster = ~iso_id, data = dta[iso != "USA" & iso != "CHN"]),

  # 5. Continent-year FE + country time trends (replaces i.year)
  feols(as.formula(paste("growthWDI ~", base_rhs,
                         "| iso_id[time, time2] + year^continent")),
        cluster = ~iso_id, data = dta),

  # 6. Continent-year FE, no country time trends
  feols(as.formula(paste("growthWDI ~", base_rhs,
                         "| iso_id + year^continent")),
        cluster = ~iso_id, data = dta),

  # 7. No year FE (country trends only)
  feols(as.formula(paste("growthWDI ~", base_rhs,
                         "| iso_id[time, time2]")),
        cluster = ~iso_id, data = dta),

  # 8. Linear country time trends only (no quadratic trend)
  feols(as.formula(paste("growthWDI ~", base_rhs,
                         "| iso_id[time] + year")),
        cluster = ~iso_id, data = dta),

  # 9. One lag of GDP growth as control (lagged dependent variable), no year FE
  feols(as.formula(paste("growthWDI ~ lag1_growth +", base_rhs,
                         "| iso_id[time, time2]")),
        cluster = ~iso_id, data = dta),

  # 10. Three lags of GDP growth, no year FE
  feols(as.formula(paste("growthWDI ~ lag1_growth + lag2_growth + lag3_growth +",
                         base_rhs, "| iso_id[time, time2]")),
        cluster = ~iso_id, data = dta),

  # 11. Penn World Tables data as outcome, no year FE
  feols(as.formula(paste("rgdpCAPgr ~", base_rhs,
                         "| iso_id[time, time2]")),
        cluster = ~iso_id, data = dta)
)

# ----------------------------------------------------------------------------
# Save coefficients CSV (mirrors Stata's postfile: mod, b1, b2)
# ----------------------------------------------------------------------------
coef_df <- data.frame(
  mod = seq_along(specs),
  b1  = sapply(specs, function(f) unname(coef(f)["temp"])),
  b2  = sapply(specs, function(f) unname(coef(f)["temp2"]))
)
write.csv(coef_df, "data/output/Coefficients_robustness.csv", row.names = FALSE)
cat("Wrote data/output/Coefficients_robustness.csv\n")

# ----------------------------------------------------------------------------
# LaTeX table via modelsummary
# ----------------------------------------------------------------------------
col_labels <- c("Base", ">20yrs", "No oil", "No US/CHN",
                "ContYr FE", "ContYr+noTrend", "No YrFE",
                "LinearTime", "LDV 1lag", "LDV 3lags", "PWT")
names(specs) <- col_labels

# Optimum temperature for each spec: -b1 / (2 * b2)
opts <- sapply(specs, function(f) {
  b <- coef(f)
  round(unname(b["temp"] / (-2 * b["temp2"])), 1)
})

# Build add_rows data.frame (extra GOF row: Optimum temperature)
add_rows_df <- data.frame(term = "Optimum (C)", stringsAsFactors = FALSE)
for (i in seq_along(col_labels)) add_rows_df[[col_labels[i]]] <- opts[i]

modelsummary(
  specs,
  output   = "figures/ExtendedDataFigs_Input/Table1.tex",
  coef_map = c(
    "temp"    = "Temp.",
    "temp2"   = "Temp. sq.",
    "precip"  = "Precip.",
    "precip2" = "Precip. sq."
  ),
  stars    = c("*" = 0.10, "**" = 0.05, "***" = 0.01),
  gof_map  = list(
    list(raw = "nobs",      clean = "Observations", fmt = 0),
    list(raw = "r.squared", clean = "R squared",    fmt = 3)
  ),
  add_rows = add_rows_df,
  title    = "Extended Data Table 1: Robustness checks",
  escape   = FALSE,
  notes    = paste("Clustered standard errors by country in parentheses.",
                   "Precip. rescaled by 1/1000 relative to main text.")
)
cat("Wrote figures/ExtendedDataFigs_Input/Table1.tex\n")
cat("MakeExtendedDataTable1.R complete.\n")
