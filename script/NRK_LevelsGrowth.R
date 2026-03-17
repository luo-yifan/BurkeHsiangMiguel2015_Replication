
# NRK_LevelsGrowth.R
# Q5: Estimate models incorporating NRK insights
#
# NRK's key critiques of BHM:
#   1. Temperature is autocorrelated → OLS on levels conflates persistent trends with shocks
#   2. BHM's quadratic in T_it combines cross-sectional (between) and within variation
#   3. BHM's growth-rate specification tests for growth effects, not level effects
#   4. Omitting lags of temperature (and GDP growth) may bias estimates
#
# This script addresses these critiques via:
#   (A) Distributed-lag model: lags of temp + lags of GDP growth (addresses #1, #4)
#   (B) Local projections (Jordà 2005): cumulative IRF at h = 0..10 horizons
#       → If cumulative IRF grows permanently ⟹ growth effect
#       → If cumulative IRF plateaus or reverses ⟹ level effect
#   (C) Levels regression: log(GDPpc) as outcome (directly tests level effect)
#   (D) State-dependent local projections: interact shock with T̄_i
#       (mimics NRK Table 3: Δy_it = δ₁τ_it + δ₂(τ_it × T̄_i) + controls)
#
# Output: figures/ExtendedDataFigs_Input/NRK_IRF.pdf
#         data/output/nrk_levels_growth_results.csv

rm(list = ls())

library(fixest)
library(data.table)

setwd("/Users/yl3699/Downloads/BurkeHsiangMiguel2015_Replication")
dir.create("figures/ExtendedDataFigs_Input", recursive = TRUE, showWarnings = FALSE)
dir.create("data/output",                   recursive = TRUE, showWarnings = FALSE)

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

# Log GDP per capita — need to recover level from growth rates
# growthWDI = (GDPpc_t - GDPpc_{t-1}) / GDPpc_{t-1}  [not in %]
# To get log-levels: cumulate growth rates from a base
# Strategy: set log(GDPpc) at the first observed year = 0 per country,
#           then cumsum the approximate log-differences (≈ growth rates for small growth)
dta[, lnGDP := {
  gr <- growthWDI    # approximately Δln(GDP)
  gr[is.na(gr)] <- 0
  cumsum(gr)
}, by = iso]

# Country mean temperature (30-yr and full sample)
dta[, T_bar_30   := mean(temp[year <= 1989], na.rm = TRUE), by = iso]
dta[, T_bar_full := mean(temp, na.rm = TRUE), by = iso]
dta[is.na(T_bar_30), T_bar_30 := T_bar_full]

# ============================================================
# 2. NRK AR SHOCK (same as NRK_MS_Analysis.R)
#    T_it = Σ γ_j T_{i,t-j} + Σ θ_j T_{i,t-j}·T̄_i + μ_i + μ_t + τ_it
# ============================================================
for (j in 1:3) {
  dta[, paste0("temp_L", j) := shift(temp, j), by = iso]
}
ar_fit <- feols(
  temp ~ temp_L1 + temp_L2 + temp_L3 +
         I(temp_L1 * T_bar_30) + I(temp_L2 * T_bar_30) + I(temp_L3 * T_bar_30) |
         iso + year,
  data = dta
)
dta[, shock := NA_real_]
dta[obs(ar_fit), shock := resid(ar_fit)]
dta[, shock_x_Tbar := shock * T_bar_30]
cat(sprintf("NRK AR shock: N=%d, AR R²=%.4f, shock obs=%d\n",
            nobs(ar_fit), r2(ar_fit,"r2"), sum(!is.na(dta$shock))))

# Lags of growthWDI (for AR-in-outcome controls)
for (j in 1:3) {
  dta[, paste0("g_L", j) := shift(growthWDI, j), by = iso]
}

# Lags of temp for distributed-lag model
for (j in 1:6) {
  dta[, paste0("temp_lag", j)  := shift(temp,  j), by = iso]
  dta[, paste0("temp2_lag", j) := shift(temp2, j), by = iso]
}

# ============================================================
# 3. MODEL A: BHM WITH OUTCOME AND TEMP LAGS
#    Addresses NRK concerns about serial correlation and omitted lags.
#    Control for 3 lags of GDP growth and 3 lags of temperature.
# ============================================================
cat("\n=== Model A: BHM + lags of growth and temperature ===\n")

fit_bhm <- feols(
  growthWDI ~ temp + temp2 + prec + prec2 | iso_id[time, time2] + year,
  cluster = ~iso_id, data = dta
)

fit_A <- feols(
  growthWDI ~ temp + temp2 + prec + prec2 +
              g_L1 + g_L2 + g_L3 +
              temp_lag1 + temp2_lag1 + temp_lag2 + temp2_lag2 + temp_lag3 + temp2_lag3 |
              iso_id[time, time2] + year,
  cluster = ~iso_id, data = dta
)

bhm_b  <- coef(fit_bhm)
A_b    <- coef(fit_A)
A_vcov <- vcov(fit_A)

stars <- function(t) {
  a <- abs(t)
  if (a > 2.576) "***" else if (a > 1.96) "**" else if (a > 1.645) "*" else ""
}

cat(sprintf("BHM baseline: temp=%.5f (t=%.2f%s), temp²=%.6f (t=%.2f%s), N=%d\n",
    bhm_b["temp"],  bhm_b["temp"] / sqrt(vcov(fit_bhm)["temp","temp"]),
    stars(bhm_b["temp"] / sqrt(vcov(fit_bhm)["temp","temp"])),
    bhm_b["temp2"], bhm_b["temp2"] / sqrt(vcov(fit_bhm)["temp2","temp2"]),
    stars(bhm_b["temp2"] / sqrt(vcov(fit_bhm)["temp2","temp2"])),
    nobs(fit_bhm)))

cat(sprintf("Model A (+lags): temp=%.5f (t=%.2f%s), temp²=%.6f (t=%.2f%s), N=%d\n",
    A_b["temp"],  A_b["temp"] / sqrt(A_vcov["temp","temp"]),
    stars(A_b["temp"] / sqrt(A_vcov["temp","temp"])),
    A_b["temp2"], A_b["temp2"] / sqrt(A_vcov["temp2","temp2"]),
    stars(A_b["temp2"] / sqrt(A_vcov["temp2","temp2"])),
    nobs(fit_A)))

cat(sprintf("  OptT BHM: %.1f°C  |  OptT A: %.1f°C\n",
    -bhm_b["temp"] / (2 * bhm_b["temp2"]),
    -A_b["temp"]   / (2 * A_b["temp2"])))

# Cumulative temp effect (contemporaneous + 3 lags) at a reference T = 20°C
cum_temp_A <- A_b["temp"] + A_b["temp_lag1"] + A_b["temp_lag2"] + A_b["temp_lag3"]
cum_temp2_A <- A_b["temp2"] + A_b["temp2_lag1"] + A_b["temp2_lag2"] + A_b["temp2_lag3"]
cat(sprintf("  Cumulative temp effect (lag 0-3): temp=%.5f, temp²=%.6f\n",
            cum_temp_A, cum_temp2_A))
cat(sprintf("  Cumulative ME @ 20°C: %.5f  (BHM contemporaneous: %.5f)\n",
            cum_temp_A + 2 * cum_temp2_A * 20,
            bhm_b["temp"] + 2 * bhm_b["temp2"] * 20))

# ============================================================
# 4. MODEL B: LOCAL PROJECTIONS (Jordà 2005)
#    Levels vs. growth test:
#      Regress cumulative log-growth [y_{t+h} - y_t] on shock τ_t,
#      controlling for FE and lagged outcomes/temperature.
#      h = 0, 1, ..., 10
#
#    If IRF (β_h) keeps growing with h → persistent growth-rate effect
#    If IRF plateaus or reverses → temporary (level) effect
#
#    We implement state-dependent: separate IRFs for hot (T̄ > 20°C) and cold (T̄ ≤ 20°C)
# ============================================================
cat("\n=== Model B: Local Projections (levels vs. growth) ===\n")

H <- 10  # max horizon
lp_results <- data.frame(
  h      = 0:H,
  beta   = NA_real_,   # unconditional shock response
  se     = NA_real_,
  beta_cold = NA_real_, # response for T̄ ≤ 20°C countries
  se_cold   = NA_real_,
  beta_hot  = NA_real_, # response for T̄ > 20°C countries
  se_hot    = NA_real_
)

# Median T̄ as hot/cold threshold
tbar_med <- median(dta$T_bar_30, na.rm = TRUE)
dta[, hot := as.integer(T_bar_30 > tbar_med)]
dta[, shock_hot  := shock * hot]
dta[, shock_cold := shock * (1 - hot)]
cat(sprintf("Hot/cold threshold: median T̄ = %.1f°C\n", tbar_med))

for (h in 0:H) {
  # Lead the outcome h periods: cumulative growth from t to t+h
  dta[, cum_growth := {
    fut <- shift(lnGDP, -h)  # lnGDP at t+h
    fut - lnGDP               # cumulative log-growth
  }, by = iso]

  # Regression: cum_growth ~ shock + shock_hot + controls | FE
  # Controls: 3 lags of outcome + 3 lags of temp shock
  # FE: iso_id[time, time2] + year (same as BHM)
  tryCatch({
    fit_lp <- feols(
      cum_growth ~ shock + shock_hot +
                   g_L1 + g_L2 + g_L3 |
                   iso_id[time, time2] + year,
      cluster = ~iso_id, data = dta
    )
    bb <- coef(fit_lp)
    vv <- vcov(fit_lp)

    # Cold response = coef on shock (shock_hot=0 for cold)
    # Hot  response = coef on shock + coef on shock_hot
    lp_results$beta[h+1]      <- bb["shock"]
    lp_results$se[h+1]        <- sqrt(vv["shock", "shock"])
    lp_results$beta_cold[h+1] <- bb["shock"]
    lp_results$se_cold[h+1]   <- sqrt(vv["shock", "shock"])
    lp_results$beta_hot[h+1]  <- bb["shock"] + bb["shock_hot"]
    lp_results$se_hot[h+1]    <- sqrt(vv["shock","shock"] +
                                       vv["shock_hot","shock_hot"] +
                                       2 * vv["shock","shock_hot"])
    cat(sprintf("  h=%2d: β(cold)=%.4f (t=%.2f%s), β(hot)=%.4f (t=%.2f%s)\n",
                h,
                lp_results$beta_cold[h+1],
                lp_results$beta_cold[h+1] / lp_results$se_cold[h+1],
                stars(lp_results$beta_cold[h+1] / lp_results$se_cold[h+1]),
                lp_results$beta_hot[h+1],
                lp_results$beta_hot[h+1] / lp_results$se_hot[h+1],
                stars(lp_results$beta_hot[h+1] / lp_results$se_hot[h+1])))
  }, error = function(e) {
    cat(sprintf("  h=%2d: ERROR — %s\n", h, conditionMessage(e)))
  })
}

# ============================================================
# 5. MODEL C: LEVELS REGRESSION
#    Use ln(GDPpc) as the dependent variable instead of growth rate.
#    This directly tests whether temperature levels affect income levels.
#    Controls: country FE (absorbs country baseline), year FE, country trends.
# ============================================================
cat("\n=== Model C: Levels regression (log GDP on temperature levels) ===\n")

fit_C <- feols(
  lnGDP ~ temp + temp2 + prec + prec2 | iso_id[time, time2] + year,
  cluster = ~iso_id, data = dta
)
C_b    <- coef(fit_C)
C_vcov <- vcov(fit_C)

cat(sprintf("Levels regression: temp=%.5f (t=%.2f%s), temp²=%.6f (t=%.2f%s), N=%d\n",
    C_b["temp"],  C_b["temp"] / sqrt(C_vcov["temp","temp"]),
    stars(C_b["temp"] / sqrt(C_vcov["temp","temp"])),
    C_b["temp2"], C_b["temp2"] / sqrt(C_vcov["temp2","temp2"]),
    stars(C_b["temp2"] / sqrt(C_vcov["temp2","temp2"])),
    nobs(fit_C)))
if (C_b["temp2"] < 0 && C_b["temp"] > 0) {
  cat(sprintf("  OptT (levels): %.1f°C\n", -C_b["temp"] / (2 * C_b["temp2"])))
}

# ============================================================
# 6. MODEL D: STATE-DEPENDENT NRK-STYLE (with lags of growth)
#    Δy_it = δ₁τ_it + δ₂(τ_it × T̄_i) + lags(Δy) + lags(τ) | FE
# ============================================================
cat("\n=== Model D: NRK state-dependent + lags ===\n")

# Lags of shock
for (j in 1:3) {
  dta[, paste0("shock_L", j) := shift(shock, j), by = iso]
}

fit_D <- feols(
  growthWDI ~ shock + shock_x_Tbar +
              g_L1 + g_L2 + g_L3 +
              shock_L1 + shock_L2 + shock_L3 |
              iso_id[time, time2] + year,
  cluster = ~iso_id, data = dta
)
D_b    <- coef(fit_D)
D_vcov <- vcov(fit_D)

cat(sprintf("δ₁ (shock):      %.5f (t=%.2f%s)\n",
    D_b["shock"], D_b["shock"] / sqrt(D_vcov["shock","shock"]),
    stars(D_b["shock"] / sqrt(D_vcov["shock","shock"]))))
cat(sprintf("δ₂ (shock×T̄):   %.6f (t=%.2f%s)\n",
    D_b["shock_x_Tbar"], D_b["shock_x_Tbar"] / sqrt(D_vcov["shock_x_Tbar","shock_x_Tbar"]),
    stars(D_b["shock_x_Tbar"] / sqrt(D_vcov["shock_x_Tbar","shock_x_Tbar"]))))

bliss_D <- -D_b["shock"] / D_b["shock_x_Tbar"]
cat(sprintf("Bliss T̄: %.1f°C  |  ME @ T̄=20: %.4f  |  ME @ T̄=28: %.4f\n",
    bliss_D,
    D_b["shock"] + D_b["shock_x_Tbar"] * 20,
    D_b["shock"] + D_b["shock_x_Tbar"] * 28))
cat(sprintf("N = %d\n", nobs(fit_D)))

# ============================================================
# 7. SAVE LP RESULTS
# ============================================================
write.csv(lp_results, "data/output/nrk_local_projections.csv", row.names = FALSE)
cat("\nSaved data/output/nrk_local_projections.csv\n")

# Summary results table
summ <- data.frame(
  model = c("BHM baseline", "A: BHM+lags", "C: Levels OLS", "D: NRK shock+lags"),
  beta1 = c(bhm_b["temp"], A_b["temp"], C_b["temp"], D_b["shock"]),
  se1   = c(sqrt(vcov(fit_bhm)["temp","temp"]), sqrt(A_vcov["temp","temp"]),
            sqrt(C_vcov["temp","temp"]),         sqrt(D_vcov["shock","shock"])),
  beta2 = c(bhm_b["temp2"], A_b["temp2"], C_b["temp2"], D_b["shock_x_Tbar"]),
  se2   = c(sqrt(vcov(fit_bhm)["temp2","temp2"]), sqrt(A_vcov["temp2","temp2"]),
            sqrt(C_vcov["temp2","temp2"]),          sqrt(D_vcov["shock_x_Tbar","shock_x_Tbar"])),
  nobs  = c(nobs(fit_bhm), nobs(fit_A), nobs(fit_C), nobs(fit_D))
)
write.csv(summ, "data/output/nrk_levels_growth_results.csv", row.names = FALSE)
cat("Saved data/output/nrk_levels_growth_results.csv\n")

# ============================================================
# 8. PLOT: Local Projection IRF (cumulative log-GDP response to 1°C shock)
# ============================================================
pdf("figures/ExtendedDataFigs_Input/NRK_IRF.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

# Panel A: unconditional cumulative IRF
h_seq  <- lp_results$h
b_cold <- lp_results$beta_cold
b_hot  <- lp_results$beta_hot
se_cold <- lp_results$se_cold
se_hot  <- lp_results$se_hot

ylim_all <- range(c(b_cold - 1.96*se_cold, b_cold + 1.96*se_cold,
                    b_hot  - 1.96*se_hot,  b_hot  + 1.96*se_hot), na.rm = TRUE)

plot(h_seq, b_cold, type = "l", lwd = 2, col = "steelblue",
     ylim = ylim_all, xlab = "Horizon h (years)",
     ylab = "Cumulative log-GDP response",
     main = "Local Projection IRF\n(Cold vs Hot countries, 1°C shock)")
polygon(c(h_seq, rev(h_seq)),
        c(b_cold - 1.96*se_cold, rev(b_cold + 1.96*se_cold)),
        col = adjustcolor("steelblue", 0.15), border = NA)
lines(h_seq, b_hot,  lwd = 2, col = "tomato")
polygon(c(h_seq, rev(h_seq)),
        c(b_hot - 1.96*se_hot, rev(b_hot + 1.96*se_hot)),
        col = adjustcolor("tomato", 0.15), border = NA)
abline(h = 0, lty = 2, col = "gray50")
legend("topright",
       legend = c(sprintf("Cold (T̄ ≤ %.0f°C)", tbar_med),
                  sprintf("Hot  (T̄ > %.0f°C)", tbar_med)),
       lty = 1, lwd = 2, col = c("steelblue","tomato"), cex = 0.85)

# Panel B: hot minus cold difference
diff_b  <- b_hot - b_cold
diff_se <- sqrt(se_hot^2 + se_cold^2)  # rough (independent assumption)
plot(h_seq, diff_b, type = "l", lwd = 2, col = "purple",
     ylim = range(c(diff_b - 1.96*diff_se, diff_b + 1.96*diff_se), na.rm = TRUE),
     xlab = "Horizon h (years)",
     ylab = "Hot minus Cold IRF",
     main = "Differential Response\n(Hot - Cold countries)")
polygon(c(h_seq, rev(h_seq)),
        c(diff_b - 1.96*diff_se, rev(diff_b + 1.96*diff_se)),
        col = adjustcolor("purple", 0.15), border = NA)
abline(h = 0, lty = 2, col = "gray50")

dev.off()
cat("Saved figures/ExtendedDataFigs_Input/NRK_IRF.pdf\n")

# ============================================================
# 9. FINAL SUMMARY PRINT
# ============================================================
cat("\n")
cat("========================================================\n")
cat("  SUMMARY: NRK-Inspired Model Estimates\n")
cat("========================================================\n\n")
cat("Dependent variable: growthWDI (GDP growth rate)\n\n")

fmt <- function(b, se) {
  t <- b / se
  sprintf("%.5f (t=%.2f%s)", b, t, stars(t))
}

cat("(0) BHM baseline (no lags):\n")
cat(sprintf("    temp  = %s\n", fmt(bhm_b["temp"],  sqrt(vcov(fit_bhm)["temp","temp"]))))
cat(sprintf("    temp² = %s\n", fmt(bhm_b["temp2"], sqrt(vcov(fit_bhm)["temp2","temp2"]))))
cat(sprintf("    OptT  = %.1f°C, ME@20°C = %.4f\n\n",
    -bhm_b["temp"]/(2*bhm_b["temp2"]),
    bhm_b["temp"] + 2*bhm_b["temp2"]*20))

cat("(A) BHM + 3 lags of growth + 3 lags of temp:\n")
cat(sprintf("    temp  = %s\n", fmt(A_b["temp"],  sqrt(A_vcov["temp","temp"]))))
cat(sprintf("    temp² = %s\n", fmt(A_b["temp2"], sqrt(A_vcov["temp2","temp2"]))))
cat(sprintf("    OptT  = %.1f°C, ME@20°C = %.4f\n", -A_b["temp"]/(2*A_b["temp2"]),
    A_b["temp"] + 2*A_b["temp2"]*20))
cat(sprintf("    Cum. ME@20°C (lags 0-3) = %.4f\n\n",
    cum_temp_A + 2*cum_temp2_A*20))

cat("(C) Levels regression (lnGDP ~ temp + temp²):\n")
cat(sprintf("    temp  = %s\n", fmt(C_b["temp"],  sqrt(C_vcov["temp","temp"]))))
cat(sprintf("    temp² = %s\n", fmt(C_b["temp2"], sqrt(C_vcov["temp2","temp2"]))))
if (C_b["temp2"] < 0 && C_b["temp"] > 0) {
  cat(sprintf("    OptT  = %.1f°C\n", -C_b["temp"]/(2*C_b["temp2"])))
}
cat("\n")

cat("(D) NRK shock + lags (δ₁τ + δ₂τ×T̄):\n")
cat(sprintf("    δ₁ (shock)    = %s\n", fmt(D_b["shock"], sqrt(D_vcov["shock","shock"]))))
cat(sprintf("    δ₂ (shock×T̄) = %s\n", fmt(D_b["shock_x_Tbar"], sqrt(D_vcov["shock_x_Tbar","shock_x_Tbar"]))))
cat(sprintf("    Bliss T̄ = %.1f°C\n", bliss_D))

cat("\n(B) Local Projections — cumulative IRF at selected horizons:\n")
cat("    (h=0 omitted: outcome is identically zero by construction)\n")
cat(sprintf("    %-4s  %-12s  %-12s\n", "h", "beta(cold)", "beta(hot)"))
for (h in c(1,2,3,5,10)) {
  idx <- h + 1
  bc_h <- lp_results$beta_cold[idx]
  bh_h <- lp_results$beta_hot[idx]
  sc_h <- lp_results$se_cold[idx]
  sh_h <- lp_results$se_hot[idx]
  if (!is.na(bc_h) && !is.na(sc_h) && sc_h > 0 &&
      !is.na(bh_h) && !is.na(sh_h) && sh_h > 0) {
    cat(sprintf("    h=%-2d  %+.4f (%s)  %+.4f (%s)\n",
        h, bc_h, stars(bc_h / sc_h), bh_h, stars(bh_h / sh_h)))
  } else {
    cat(sprintf("    h=%-2d  %+.4f (NA)   %+.4f (NA)\n", h,
        ifelse(is.na(bc_h), 0, bc_h), ifelse(is.na(bh_h), 0, bh_h)))
  }
}

cat("\nNRK_LevelsGrowth.R complete.\n")
