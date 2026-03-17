
# GenerateFigure2Data.R
# R equivalent of GenerateFigure2Data.do
# Generates all output files needed for Figure 2 and downstream R scripts

rm(list = ls())

library(fixest)

setwd("/Users/yl3699/Downloads/BurkeHsiangMiguel2015_Replication")
dir.create("data/output", recursive = TRUE, showWarnings = FALSE)

# Read data
dta <- read.csv("data/input/GrowthClimateDataset.csv")
dta$temp  <- dta$UDel_temp_popweight
dta$temp2 <- dta$temp^2   # explicit column avoids I() naming issues in coef()

# Get pre-computed country time trend variable names (created by Stata xi command).
# read.csv() prepends "X" to column names starting with "_", so match both forms.
yi_cols <- grep("^X?_yi_", names(dta), value = TRUE)
y2_cols <- grep("^X?_y2_", names(dta), value = TRUE)

# Build trend formula string (backtick-quote to handle special characters)
trend_vars <- c(yi_cols, y2_cols)
if (length(trend_vars) > 0) {
  trend_formula <- paste(paste0("`", trend_vars, "`"), collapse = " + ")
} else {
  trend_formula <- NULL
}

z90 <- qnorm(0.95)  # 90% two-sided CI: 1.645

# Helper: build formula string safely (handles NULL trend_formula)
make_fml <- function(lhs, rhs_base, fe = "iso_id + year") {
  rhs <- if (!is.null(trend_formula)) paste(rhs_base, "+", trend_formula) else rhs_base
  as.formula(paste(lhs, "~", rhs, "|", fe))
}

# -----------------------------------------------------------------------
# 1. SAVE MAIN DATASET
# -----------------------------------------------------------------------
mainDataset <- dta[, c("iso", "year", "UDel_temp_popweight", "Pop", "TotGDP",
                        "growthWDI", "GDPpctile_WDIppp", "continent", "countryname")]
write.csv(mainDataset, "data/output/mainDataset.csv", row.names = FALSE)
cat("Wrote data/output/mainDataset.csv\n")

# -----------------------------------------------------------------------
# 2. BASELINE REGRESSION (Panel A of Figure 2)
# -----------------------------------------------------------------------
cat("Running baseline regression...\n")

fml_base <- make_fml("growthWDI",
                     "temp + temp2 + UDel_precip_popweight + UDel_precip_popweight_2")
fit_base <- feols(fml_base, cluster = ~iso_id, data = dta)

b1 <- coef(fit_base)["temp"]
b2 <- coef(fit_base)["temp2"]
V  <- vcov(fit_base)[c("temp", "temp2"), c("temp", "temp2")]

# Save coefficients
coef_df <- data.frame(b1 = b1, b2 = b2)
write.csv(coef_df, "data/output/estimatedCoefficients.csv", row.names = FALSE)

# Predicted response + 90% CI at temp = -5 to 35
temps    <- (-5):35
estimate <- b1 * temps + b2 * temps^2
se_est   <- sqrt(pmax(0, temps^2 * V[1,1] + temps^4 * V[2,2] + 2 * temps^3 * V[1,2]))
resp_df  <- data.frame(x = temps, estimate = estimate,
                       min90 = estimate - z90 * se_est,
                       max90 = estimate + z90 * se_est)
write.csv(resp_df, "data/output/estimatedGlobalResponse.csv", row.names = FALSE)
cat("Wrote estimatedGlobalResponse.csv and estimatedCoefficients.csv\n")

# -----------------------------------------------------------------------
# 3. INCOME HETEROGENEITY (Panels B, D, E)
# -----------------------------------------------------------------------
cat("Running income heterogeneity regressions...\n")

outcome_vars <- c("growthWDI", "AgrGDPgrowthCap", "NonAgrGDPgrowthCap")
all_het <- list()

for (v in outcome_vars) {
  cat(" ", v, "\n")
  tmp <- dta
  tmp$poor <- ifelse(is.na(tmp$GDPpctile_WDIppp), NA,
                     as.integer(tmp$GDPpctile_WDIppp < 50))
  tmp <- tmp[!is.na(tmp$poor), ]

  # Separate regressors: base = rich (poor=0), differential = poor-specific
  tmp$temp_poor  <- tmp$temp  * tmp$poor
  tmp$temp2_poor <- tmp$temp2 * tmp$poor
  tmp$prec_poor  <- tmp$UDel_precip_popweight   * tmp$poor
  tmp$prec2_poor <- tmp$UDel_precip_popweight_2 * tmp$poor

  fml_het <- make_fml(v,
    "temp + temp2 + temp_poor + temp2_poor + UDel_precip_popweight + UDel_precip_popweight_2 + prec_poor + prec2_poor")
  fit_het <- feols(fml_het, cluster = ~iso_id, data = tmp)

  b    <- coef(fit_het)
  Vmat <- vcov(fit_het)

  # Rich (poor=0): coefficients on temp, temp2
  b1r <- b["temp"];   b2r <- b["temp2"]
  V_r <- Vmat[c("temp","temp2"), c("temp","temp2")]

  # Poor (poor=1): rich + differential
  b1p <- b["temp"] + b["temp_poor"]
  b2p <- b["temp2"] + b["temp2_poor"]
  idx4   <- c("temp", "temp2", "temp_poor", "temp2_poor")
  V_all4 <- Vmat[idx4, idx4]

  rows <- list()
  for (tt in 0:30) {
    cv_r  <- c(tt, tt^2)
    est_r <- b1r * tt + b2r * tt^2
    se_r  <- sqrt(pmax(0, as.numeric(t(cv_r) %*% V_r %*% cv_r)))

    cv_p  <- c(tt, tt^2, tt, tt^2)
    est_p <- b1p * tt + b2p * tt^2
    se_p  <- sqrt(pmax(0, as.numeric(t(cv_p) %*% V_all4 %*% cv_p)))

    rows <- c(rows, list(
      data.frame(x=tt, interact=0, estimate=est_r,
                 min90=est_r - z90*se_r, max90=est_r + z90*se_r,
                 model=v, stringsAsFactors=FALSE),
      data.frame(x=tt, interact=1, estimate=est_p,
                 min90=est_p - z90*se_p, max90=est_p + z90*se_p,
                 model=v, stringsAsFactors=FALSE)
    ))
  }
  all_het[[v]] <- do.call(rbind, rows)
}

het_df <- do.call(rbind, all_het)
write.csv(het_df, "data/output/EffectHeterogeneity.csv", row.names = FALSE)
cat("Wrote data/output/EffectHeterogeneity.csv\n")

# -----------------------------------------------------------------------
# 4. TIME HETEROGENEITY (Panel C)
# -----------------------------------------------------------------------
cat("Running time heterogeneity regression...\n")

tmp <- dta
tmp$early       <- as.integer(tmp$year < 1990)
tmp$temp_early  <- tmp$temp  * tmp$early
tmp$temp2_early <- tmp$temp2 * tmp$early
tmp$prec_early  <- tmp$UDel_precip_popweight   * tmp$early
tmp$prec2_early <- tmp$UDel_precip_popweight_2 * tmp$early

fml_time <- make_fml("growthWDI",
  "temp + temp2 + temp_early + temp2_early + UDel_precip_popweight + UDel_precip_popweight_2 + prec_early + prec2_early")
fit_time <- feols(fml_time, cluster = ~iso_id, data = tmp)

b    <- coef(fit_time)
Vmat <- vcov(fit_time)

# Late (base, early=0)
b1_late <- b["temp"];   b2_late <- b["temp2"]
V_late  <- Vmat[c("temp","temp2"), c("temp","temp2")]

# Early = late + differential
b1_early <- b["temp"] + b["temp_early"]
b2_early <- b["temp2"] + b["temp2_early"]
idx4   <- c("temp", "temp2", "temp_early", "temp2_early")
V_all4 <- Vmat[idx4, idx4]

rows <- list()
for (tt in 0:30) {
  cv_late <- c(tt, tt^2)
  est_l   <- b1_late * tt + b2_late * tt^2
  se_l    <- sqrt(pmax(0, as.numeric(t(cv_late) %*% V_late %*% cv_late)))

  cv_e  <- c(tt, tt^2, tt, tt^2)
  est_e <- b1_early * tt + b2_early * tt^2
  se_e  <- sqrt(pmax(0, as.numeric(t(cv_e) %*% V_all4 %*% cv_e)))

  rows <- c(rows, list(
    data.frame(x=tt, interact=0, estimate=est_l,
               min90=est_l - z90*se_l, max90=est_l + z90*se_l),
    data.frame(x=tt, interact=1, estimate=est_e,
               min90=est_e - z90*se_e, max90=est_e + z90*se_e)
  ))
}
time_df <- do.call(rbind, rows)
write.csv(time_df, "data/output/EffectHeterogeneityOverTime.csv", row.names = FALSE)
cat("Wrote data/output/EffectHeterogeneityOverTime.csv\n")

cat("\nGenerateFigure2Data.R complete.\n")
