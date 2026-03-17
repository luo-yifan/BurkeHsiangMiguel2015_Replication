
# NRK_MS_Analysis.R
# Alternative temperature specifications for BHM baseline:
#
#   (a) NRK shock estimator (simple and sophisticated):
#       Replace temp  → shock_it
#       Replace temp² → shock_it × T̄_i
#       where T̄_i = country mean temperature over first 30 years of sample
#
#       Simple shock:       shock_it = T_it - T̄_i^{30yr}
#       Sophisticated shock: innovation to temp from nonlinear AR model
#                            (Nath, Ramey, Klenow 2024, eq. 4):
#           T_it = Σ_j γ_j T_{i,t-j} + Σ_j θ_j T_{i,t-j}·T̄_i + μ_i + μ_t + τ_it
#           p = 3 lags, year FE, country FE
#
#   (b) McIntosh–Schlenker within estimator:
#       Replace temp² with (temp - T̄_i^{full})²
#
#   (c) McIntosh–Schlenker hybrid estimator:
#       Keep temp + temp², add (temp - T̄_i^{full})²
#
# FE structure kept identical to BHM: iso_id[time, time2] + year
# Clustering: by country (iso_id)

rm(list = ls())

library(fixest)
library(data.table)

setwd("/Users/yl3699/Downloads/BurkeHsiangMiguel2015_Replication")

# ============================================================
# 1. LOAD AND PREPARE DATA
# ============================================================
dta <- as.data.table(read.csv("data/input/GrowthClimateDataset.csv"))
setorder(dta, iso, year)

dta[, temp  := UDel_temp_popweight]
dta[, temp2 := temp^2]
dta[, prec  := UDel_precip_popweight]
dta[, prec2 := prec^2]
dta[, time  := year - 1985]
dta[, time2 := time^2]

# ============================================================
# 2. COUNTRY MEAN TEMPERATURES
# ============================================================
# First 30 years of sample: 1960-1989 (used by NRK approaches)
# Full-sample mean (used by McIntosh-Schlenker approaches)
dta[, T_bar_30 := mean(temp[year <= 1989], na.rm = TRUE), by = iso]
dta[, T_bar_full := mean(temp, na.rm = TRUE), by = iso]

cat("Countries with T_bar_30 computed:", dta[!is.na(T_bar_30), uniqueN(iso)], "\n")
cat("  (countries with no 1960-1989 obs use NA for T_bar_30)\n")

# For countries with no 1960-1989 observations, fall back to full-sample mean
dta[is.na(T_bar_30), T_bar_30 := T_bar_full]

# ============================================================
# 3. SIMPLE NRK SHOCK: deviation from 30-year mean
# ============================================================
dta[, shock_simple := temp - T_bar_30]

# ============================================================
# 4. SOPHISTICATED NRK SHOCK: innovation to AR model
#    T_it = Σ_j γ_j T_{i,t-j} + Σ_j θ_j T_{i,t-j}·T̄_i + μ_i + μ_t + τ_it
#    p = 3 lags (as in NRK Table 3 notes)
# ============================================================
# Create lag variables
for (j in 1:3) {
  dta[, paste0("temp_L", j) := shift(temp, j), by = iso]
}

# Fit AR model: temp on 3 lags + lag×T̄_i interactions, with country + year FE
ar_fit <- feols(
  temp ~ temp_L1 + temp_L2 + temp_L3 +
         I(temp_L1 * T_bar_30) + I(temp_L2 * T_bar_30) + I(temp_L3 * T_bar_30) |
         iso + year,
  data = dta
)
cat(sprintf("\nNRK AR model: N = %d, R² = %.4f\n",
            nobs(ar_fit), r2(ar_fit, "r2")))

# NRK shock = AR residual (only available where all 3 lags exist)
dta[, shock_nrk := NA_real_]
obs_ar <- obs(ar_fit)
dta[obs_ar, shock_nrk := resid(ar_fit)]
cat("NRK shock obs:", sum(!is.na(dta$shock_nrk)), "\n")

# Interaction terms for NRK regressions
dta[, shock_simple_x_Tbar := shock_simple * T_bar_30]
dta[, shock_nrk_x_Tbar    := shock_nrk   * T_bar_30]

# ============================================================
# 5. McINTOSH-SCHLENKER WITHIN TERM
#    within_sq = (temp - T̄_i^{full})²
# ============================================================
dta[, within_dev := temp - T_bar_full]
dta[, within_sq  := within_dev^2]

# ============================================================
# 6. RUN ALL FIVE MODELS
# ============================================================
fe_str <- "iso_id[time, time2] + year"

# (0) BHM baseline
fit0 <- feols(
  as.formula(paste("growthWDI ~ temp + temp2 + prec + prec2 |", fe_str)),
  cluster = ~iso_id, data = dta
)

# (a1) NRK simple shock
fit_a1 <- feols(
  as.formula(paste("growthWDI ~ shock_simple + shock_simple_x_Tbar + prec + prec2 |", fe_str)),
  cluster = ~iso_id, data = dta
)

# (a2) NRK sophisticated shock
fit_a2 <- feols(
  as.formula(paste("growthWDI ~ shock_nrk + shock_nrk_x_Tbar + prec + prec2 |", fe_str)),
  cluster = ~iso_id, data = dta
)

# (b) MS within estimator: replace temp² with (temp - T̄)²
fit_b <- feols(
  as.formula(paste("growthWDI ~ temp + within_sq + prec + prec2 |", fe_str)),
  cluster = ~iso_id, data = dta
)

# (c) MS hybrid estimator: keep temp + temp², add (temp - T̄)²
fit_c <- feols(
  as.formula(paste("growthWDI ~ temp + temp2 + within_sq + prec + prec2 |", fe_str)),
  cluster = ~iso_id, data = dta
)

# ============================================================
# 7. SUMMARIZE RESULTS
# ============================================================
extract_t <- function(fit, var) {
  b <- coef(fit)[var]
  se <- sqrt(vcov(fit)[var, var])
  c(est = b, se = se, t = b / se)
}

stars <- function(t) {
  a <- abs(t)
  if (a > 2.576) "***" else if (a > 1.960) "**" else if (a > 1.645) "*" else ""
}

cat("\n")
cat("==========================================================================\n")
cat("  MODEL COMPARISON: BHM baseline vs NRK and McIntosh-Schlenker variants\n")
cat("==========================================================================\n\n")

# --- (0) BHM baseline ---
b0 <- coef(fit0)
cat("(0) BHM BASELINE: growthWDI ~ temp + temp² + prec + prec² | FE\n")
cat(sprintf("    temp:   %.6f  (SE %.6f, t = %.3f%s)\n",
            b0["temp"],  sqrt(vcov(fit0)["temp","temp"]),
            b0["temp"] / sqrt(vcov(fit0)["temp","temp"]),
            stars(b0["temp"] / sqrt(vcov(fit0)["temp","temp"]))))
cat(sprintf("    temp²:  %.7f  (SE %.7f, t = %.3f%s)\n",
            b0["temp2"], sqrt(vcov(fit0)["temp2","temp2"]),
            b0["temp2"] / sqrt(vcov(fit0)["temp2","temp2"]),
            stars(b0["temp2"] / sqrt(vcov(fit0)["temp2","temp2"]))))
opt0 <- -b0["temp"] / (2 * b0["temp2"])
cat(sprintf("    Optimal T: %.1f°C  |  ME @ 20°C: %.4f pp\n", opt0,
            b0["temp"] + 2 * b0["temp2"] * 20))
cat(sprintf("    N = %d\n\n", nobs(fit0)))

# --- (a1) NRK simple shock ---
ba1 <- coef(fit_a1)
V_a1 <- vcov(fit_a1)
cat("(a1) NRK SIMPLE SHOCK: growthWDI ~ shock_simple + shock_simple×T̄ + prec | FE\n")
cat(sprintf("     δ₁ (shock):       %.6f  (SE %.6f, t = %.3f%s)\n",
            ba1["shock_simple"], sqrt(V_a1["shock_simple","shock_simple"]),
            ba1["shock_simple"] / sqrt(V_a1["shock_simple","shock_simple"]),
            stars(ba1["shock_simple"] / sqrt(V_a1["shock_simple","shock_simple"]))))
cat(sprintf("     δ₂ (shock×T̄):    %.6f  (SE %.6f, t = %.3f%s)\n",
            ba1["shock_simple_x_Tbar"], sqrt(V_a1["shock_simple_x_Tbar","shock_simple_x_Tbar"]),
            ba1["shock_simple_x_Tbar"] / sqrt(V_a1["shock_simple_x_Tbar","shock_simple_x_Tbar"]),
            stars(ba1["shock_simple_x_Tbar"] / sqrt(V_a1["shock_simple_x_Tbar","shock_simple_x_Tbar"]))))
bliss_a1 <- -ba1["shock_simple"] / ba1["shock_simple_x_Tbar"]
cat(sprintf("     Bliss T̄: %.1f°C  |  ME @ T̄=20°C: %.4f pp\n", bliss_a1,
            ba1["shock_simple"] + ba1["shock_simple_x_Tbar"] * 20))
cat(sprintf("     N = %d\n\n", nobs(fit_a1)))

# --- (a2) NRK sophisticated shock ---
ba2 <- coef(fit_a2)
V_a2 <- vcov(fit_a2)
cat("(a2) NRK SOPHISTICATED SHOCK (AR innovation): growthWDI ~ τ + τ×T̄ + prec | FE\n")
cat(sprintf("     δ₁ (shock):       %.6f  (SE %.6f, t = %.3f%s)\n",
            ba2["shock_nrk"], sqrt(V_a2["shock_nrk","shock_nrk"]),
            ba2["shock_nrk"] / sqrt(V_a2["shock_nrk","shock_nrk"]),
            stars(ba2["shock_nrk"] / sqrt(V_a2["shock_nrk","shock_nrk"]))))
cat(sprintf("     δ₂ (shock×T̄):    %.6f  (SE %.6f, t = %.3f%s)\n",
            ba2["shock_nrk_x_Tbar"], sqrt(V_a2["shock_nrk_x_Tbar","shock_nrk_x_Tbar"]),
            ba2["shock_nrk_x_Tbar"] / sqrt(V_a2["shock_nrk_x_Tbar","shock_nrk_x_Tbar"]),
            stars(ba2["shock_nrk_x_Tbar"] / sqrt(V_a2["shock_nrk_x_Tbar","shock_nrk_x_Tbar"]))))
bliss_a2 <- -ba2["shock_nrk"] / ba2["shock_nrk_x_Tbar"]
cat(sprintf("     Bliss T̄: %.1f°C  |  ME @ T̄=20°C: %.4f pp\n", bliss_a2,
            ba2["shock_nrk"] + ba2["shock_nrk_x_Tbar"] * 20))
cat(sprintf("     N = %d\n\n", nobs(fit_a2)))

# --- (b) MS within ---
bb <- coef(fit_b)
V_b <- vcov(fit_b)
cat("(b) MS WITHIN: growthWDI ~ temp + (temp-T̄_full)² + prec | FE\n")
cat(sprintf("    β₁ (temp):        %.6f  (SE %.6f, t = %.3f%s)\n",
            bb["temp"], sqrt(V_b["temp","temp"]),
            bb["temp"] / sqrt(V_b["temp","temp"]),
            stars(bb["temp"] / sqrt(V_b["temp","temp"]))))
cat(sprintf("    β₃ (within_sq):   %.7f  (SE %.7f, t = %.3f%s)\n",
            bb["within_sq"], sqrt(V_b["within_sq","within_sq"]),
            bb["within_sq"] / sqrt(V_b["within_sq","within_sq"]),
            stars(bb["within_sq"] / sqrt(V_b["within_sq","within_sq"]))))
# ME = β₁ + 2β₃·(T - T̄): at T = T̄, ME = β₁; bliss where T = T̄ - β₁/(2β₃)
# Note: for within estimator, ME varies by deviation from COUNTRY mean, not absolute T
cat(sprintf("    ME at dev=0:      %.4f pp  (ME = β₁ when T = T̄)\n", bb["temp"]))
cat(sprintf("    N = %d\n\n", nobs(fit_b)))

# --- (c) MS hybrid ---
bc <- coef(fit_c)
V_c <- vcov(fit_c)
cat("(c) MS HYBRID: growthWDI ~ temp + temp² + (temp-T̄_full)² + prec | FE\n")
cat(sprintf("    β₁ (temp):        %.6f  (SE %.6f, t = %.3f%s)\n",
            bc["temp"], sqrt(V_c["temp","temp"]),
            bc["temp"] / sqrt(V_c["temp","temp"]),
            stars(bc["temp"] / sqrt(V_c["temp","temp"]))))
cat(sprintf("    β₂ (temp²):       %.7f  (SE %.7f, t = %.3f%s)\n",
            bc["temp2"], sqrt(V_c["temp2","temp2"]),
            bc["temp2"] / sqrt(V_c["temp2","temp2"]),
            stars(bc["temp2"] / sqrt(V_c["temp2","temp2"]))))
cat(sprintf("    β₃ (within_sq):   %.7f  (SE %.7f, t = %.3f%s)\n",
            bc["within_sq"], sqrt(V_c["within_sq","within_sq"]),
            bc["within_sq"] / sqrt(V_c["within_sq","within_sq"]),
            stars(bc["within_sq"] / sqrt(V_c["within_sq","within_sq"]))))
cat(sprintf("    N = %d\n", nobs(fit_c)))

# ============================================================
# 8. RESPONSE-FUNCTION COMPARISON PLOT
#    Plot marginal effects by temperature (or T̄ for NRK)
# ============================================================
cat("\n=== Marginal effect at various temperatures ===\n")
cat(sprintf("%-8s  %-10s  %-10s  %-10s  %-10s  %-10s\n",
            "Temp", "BHM", "NRK-simp", "NRK-AR", "MS-within", "MS-hybrid"))
cat(strrep("-", 65), "\n")
for (T in c(-5, 0, 5, 10, 13, 15, 20, 25, 30)) {
  me0   <- b0["temp"]  + 2 * b0["temp2"]  * T
  # For NRK: ME depends on T̄_i not current T; show at T̄_i = T
  me_a1 <- ba1["shock_simple"] + ba1["shock_simple_x_Tbar"] * T
  me_a2 <- ba2["shock_nrk"]    + ba2["shock_nrk_x_Tbar"]    * T
  # For MS within: dev = T - T̄; at country mean T̄ = T, dev = 0, ME = β₁ + 2β₃·0 = β₁
  # But to show curvature: assume dev = 0, so ME_within = β₁ (same for all T, not very informative)
  # Better: show ME relative to own mean (dev = T - T̄_avg) where T̄_avg = world avg ~14°C
  T_world_avg <- 13.9
  me_b  <- bb["temp"] + 2 * bb["within_sq"] * (T - T_world_avg)
  me_c  <- bc["temp"] + 2 * bc["temp2"] * T + 2 * bc["within_sq"] * (T - T_world_avg)
  cat(sprintf("%-8.0f  %-10.4f  %-10.4f  %-10.4f  %-10.4f  %-10.4f\n",
              T, me0, me_a1, me_a2, me_b, me_c))
}
cat("\nNote: NRK ME at T means 'marginal effect for a country with T̄_i = T degrees'\n")
cat("Note: MS within ME uses world average T̄ ≈ 13.9°C for reference\n")

# ============================================================
# 9. OPTIONAL: SAVE RESULTS FOR PLOTTING
# ============================================================
T_seq <- seq(-5, 30, by = 0.5)
resp_df <- data.frame(
  temp        = T_seq,
  bhm         = b0["temp"]  + 2 * b0["temp2"]  * T_seq,
  nrk_simple  = ba1["shock_simple"] + ba1["shock_simple_x_Tbar"] * T_seq,
  nrk_ar      = ba2["shock_nrk"]    + ba2["shock_nrk_x_Tbar"]    * T_seq,
  ms_within   = bb["temp"]  + 2 * bb["within_sq"]  * (T_seq - mean(dta$T_bar_full, na.rm = TRUE)),
  ms_hybrid   = bc["temp"]  + 2 * bc["temp2"] * T_seq +
                  2 * bc["within_sq"] * (T_seq - mean(dta$T_bar_full, na.rm = TRUE))
)
dir.create("data/output", recursive = TRUE, showWarnings = FALSE)
write.csv(resp_df, "data/output/NRK_MS_MarginalEffects.csv", row.names = FALSE)
cat("\nSaved data/output/NRK_MS_MarginalEffects.csv\n")

cat("\nNRK_MS_Analysis.R complete.\n")
