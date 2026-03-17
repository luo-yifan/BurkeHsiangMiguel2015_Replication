
# MakeExtendedDataTable2.R
# R equivalent of MakeExtendedDataTable2.do
# Runs 6 robustness specifications for the rich/poor income-heterogeneity model.
#
# Output:
#   figures/ExtendedDataFigs_Input/Table2.tex (LaTeX regression table)
#
# Requires: install.packages("modelsummary")

rm(list = ls())

library(fixest)
library(data.table)
library(modelsummary)

setwd("/Users/yl3699/Downloads/BurkeHsiangMiguel2015_Replication")
dir.create("figures/ExtendedDataFigs_Input", recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# Load and prepare data
# Drop countries without PPP data (Stata: "drop if GDPpctile_WDIppp==.")
# ----------------------------------------------------------------------------
dta <- as.data.table(read.csv("data/input/GrowthClimateDataset.csv"))
dta <- dta[!is.na(GDPpctile_WDIppp)]

dta[, temp  := UDel_temp_popweight]
dta[, prec  := UDel_precip_popweight]
dta[, poor  := as.integer(GDPpctile_WDIppp < 50)]
dta[, temp2 := temp^2]
dta[, prec2 := prec^2]

# Interaction variables (poor differential coding: rich = base, poor = base + interaction)
dta[, temp_poor  := temp  * poor]
dta[, temp2_poor := temp2 * poor]
dta[, prec_poor  := prec  * poor]
dta[, prec2_poor := prec2 * poor]

dta[, time  := year - 1985]
dta[, time2 := time^2]

# ----------------------------------------------------------------------------
# 6 specifications
# Stata -> R:
#   areg ... a(iso_id) cl(iso_id)  ->  feols(... | iso_id, cluster = ~iso_id)
#   _yi_* _y2_* i.year             ->  iso_id[time, time2] + year
#   i.poorWDIppp#i.year            ->  year^poor  (year FE interacted with poor)
#   _cy_*  (continent-year FE)     ->  year^continent
# ----------------------------------------------------------------------------
rhs <- paste("temp + temp2 + temp_poor + temp2_poor",
             "+ prec + prec2 + prec_poor + prec2_poor")

specs <- list(
  # 1. Baseline: country trends + year FE
  feols(as.formula(paste("growthWDI ~", rhs,
                         "| iso_id[time, time2] + year")),
        cluster = ~iso_id, data = dta),

  # 2. Poor-specific year FE (i.poorWDIppp#i.year)
  feols(as.formula(paste("growthWDI ~", rhs,
                         "| iso_id[time, time2] + year^poor")),
        cluster = ~iso_id, data = dta),

  # 3. Restrict to countries with >= 20 WDI observations
  feols(as.formula(paste("growthWDI ~", rhs,
                         "| iso_id[time, time2] + year")),
        cluster = ~iso_id, data = dta[wdinomiss >= 20]),

  # 4. >= 20 obs + poor-specific year FE
  feols(as.formula(paste("growthWDI ~", rhs,
                         "| iso_id[time, time2] + year^poor")),
        cluster = ~iso_id, data = dta[wdinomiss >= 20]),

  # 5. Continent-year FE + country trends (no year FE)
  feols(as.formula(paste("growthWDI ~", rhs,
                         "| iso_id[time, time2] + year^continent")),
        cluster = ~iso_id, data = dta),

  # 6. Continent-year FE only (no country time trends)
  feols(as.formula(paste("growthWDI ~", rhs,
                         "| iso_id + year^continent")),
        cluster = ~iso_id, data = dta)
)

# ----------------------------------------------------------------------------
# Compute lincom-equivalent statistics (mirrors Stata's lincom command)
#   linb  = b[temp] + b[temp_poor]          (total poor linear coef)
#   linse = SE of that sum                  (delta method)
#   qb    = b[temp2] + b[temp2_poor]        (total poor quadratic coef)
#   qse   = SE of that sum
#   optr  = -b[temp]  / (2 * b[temp2])      (optimum for rich countries)
#   optp  = -(b[temp]+b[temp_poor]) / (2*(b[temp2]+b[temp2_poor]))  (optimum, poor)
# ----------------------------------------------------------------------------
compute_stats <- function(fit) {
  b <- coef(fit)
  V <- vcov(fit)

  se_sum <- function(n1, n2) sqrt(V[n1, n1] + V[n2, n2] + 2 * V[n1, n2])

  list(
    linb  = round(b["temp"]  + b["temp_poor"],  4),
    linse = round(se_sum("temp",  "temp_poor"),  4),
    qb    = round(b["temp2"] + b["temp2_poor"], 4),
    qse   = round(se_sum("temp2", "temp2_poor"), 4),
    optr  = round(b["temp"] / (-2 * b["temp2"]), 1),
    optp  = round((b["temp"] + b["temp_poor"]) /
                    (-2 * (b["temp2"] + b["temp2_poor"])), 1)
  )
}

stats_list <- lapply(specs, compute_stats)

# ----------------------------------------------------------------------------
# LaTeX table via modelsummary
# ----------------------------------------------------------------------------
col_labels <- c("Base", "poor-yr FE", ">20Yrs", ">20Yrs+poor-yr FE",
                "ContYr FE", "ContYr+noTrend")
names(specs) <- col_labels

# Build add_rows: one row per extra statistic
stat_names <- c(
  "beta1 + beta3",     "se(beta1 + beta3)",
  "beta2 + beta4",     "se(beta2 + beta4)",
  "Rich optimum (C)",  "Poor optimum (C)"
)
stat_keys <- c("linb", "linse", "qb", "qse", "optr", "optp")

add_rows_df <- data.frame(term = stat_names, stringsAsFactors = FALSE)
for (j in seq_along(col_labels)) {
  s <- stats_list[[j]]
  add_rows_df[[col_labels[j]]] <- sapply(stat_keys, function(k) s[[k]])
}

modelsummary(
  specs,
  output   = "figures/ExtendedDataFigs_Input/Table2.tex",
  coef_map = c(
    "temp"       = "Temperature (b1)",
    "temp2"      = "Temperature sq. (b2)",
    "temp_poor"  = "Temperature x poor (b3)",
    "temp2_poor" = "Temperature sq. x poor (b4)"
  ),
  stars    = c("*" = 0.10, "**" = 0.05, "***" = 0.01),
  gof_map  = list(
    list(raw = "nobs",      clean = "Observations", fmt = 0),
    list(raw = "r.squared", clean = "R squared",    fmt = 3)
  ),
  add_rows = add_rows_df,
  title    = "Extended Data Table 2: Income heterogeneity robustness",
  escape   = FALSE,
  notes    = "Clustered standard errors by country in parentheses. Countries without PPP data excluded."
)
cat("Wrote figures/ExtendedDataFigs_Input/Table2.tex\n")
cat("MakeExtendedDataTable2.R complete.\n")
