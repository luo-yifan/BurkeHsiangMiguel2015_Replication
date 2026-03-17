
# NRK_LocalProjections.R
# ============================================================
# NRK-INSPIRED ANALYSIS OF TEMPERATURE EFFECTS ON GDP
# Using BHM panel data (1960-2010)
#
# Addresses three NRK critiques of the BHM quadratic approach:
#
#  (1) NON-LINEARITIES
#      BHM: temp + temp² in FE model — partly between-country identification
#      NRK: shock × T̄_i (state-dependent); effect sign depends on country mean temp
#      Here: replace T² with innovation to temp interacted with T̄_i
#
#  (2) AUTOCORRELATION / TREATMENT VALIDITY
#      BHM: temperature level as treatment — serially correlated, so current
#           temp confounds current + future shocks
#      NRK: use innovation τ_it from nonlinear AR(3) model as exogenous treatment
#      Here: τ_it = residual from T_it = Σ_j γ_j T_{t-j} + Σ_j θ_j T_{t-j}·T̄_i
#            + country FE + year FE  [NRK eq. 4, p = 3 lags]
#
#  (3) LEVELS VS GROWTH
#      BHM: uses annual GDP growth rate as LHS → if significant, implies
#           temperature permanently shifts growth rates (growth effect)
#      NRK: uses LOCAL PROJECTIONS (Jordà 2005) at horizons h = 0,...,H
#           LHS = y_{i,t+h} − y_{i,t−1} (cumulative level change)
#           If β^h ∝ h → permanent growth effect
#           If β^h constant → persistent level shift
#           If β^h → 0    → transitory level effect
#      Here: run LP at h = 0,...,5 with NRK shock as treatment
#
# ============================================================

rm(list = ls())

library(fixest)
library(data.table)

setwd("/Users/yl3699/Downloads/BurkeHsiangMiguel2015_Replication")
dir.create("figures/ExtendedDataFigs_Input", recursive = TRUE, showWarnings = FALSE)
dir.create("data/output", recursive = TRUE, showWarnings = FALSE)

H <- 5   # local projection horizons: h = 0, 1, ..., H

# ============================================================
# 1. LOAD AND PREPARE DATA
# ============================================================
dta <- as.data.table(read.csv("data/input/GrowthClimateDataset.csv"))
setorder(dta, iso_id, year)

dta[, temp  := UDel_temp_popweight]
dta[, temp2 := temp^2]
dta[, prec  := UDel_precip_popweight]
dta[, prec2 := prec^2]
dta[, time  := year - 1985]
dta[, time2 := time^2]

# Country mean temperature: first 30 years of sample (NRK use sample-average; ~1960-1989)
dta[, T_bar := mean(temp[year <= 1989], na.rm = TRUE), by = iso_id]
dta[is.na(T_bar), T_bar := mean(temp, na.rm = TRUE), by = iso_id]   # fallback

cat(sprintf("Panel: %d obs, %d countries, years %d-%d\n",
            nrow(dta), dta[, uniqueN(iso_id)],
            dta[, min(year)], dta[, max(year)]))

# ============================================================
# 2. CREATE LAG VARIABLES
# ============================================================
# Temperature lags (for AR model and LP controls)
for (j in 1:3) {
  dta[, paste0("tL", j)    := shift(temp, j),       by = iso_id]
  dta[, paste0("tLxT", j)  := shift(temp, j) * T_bar, by = iso_id]
}

# GDP growth lags (LP controls — NRK include 3 lags of GDP growth)
for (j in 1:3) {
  dta[, paste0("gL", j) := shift(growthWDI, j), by = iso_id]
}

# ============================================================
# 3. NRK AR MODEL: TEMPERATURE SHOCK IDENTIFICATION
#    T_it = Σ_j γ_j T_{t-j} + Σ_j θ_j T_{t-j}·T̄_i + μ_i + μ_t + τ_it
#    Innovation τ_it is the identified temperature shock.
#
#    Purpose: by projecting out the predictable component of temperature
#    (its own AR dynamics, which vary by country mean), we isolate the
#    unpredictable component. This removes serial-correlation-induced
#    confounding: current temp levels reflect past shocks persisting
#    forward, whereas τ_it is, by construction, uncorrelated with past
#    temperature history (given FEs).
# ============================================================
ar_fml <- temp ~ tL1 + tL2 + tL3 + tLxT1 + tLxT2 + tLxT3 | iso_id + year

ar_fit <- feols(ar_fml, data = dta)
cat(sprintf("\nAR model: N=%d, R²=%.4f  (removes predictable temp component)\n",
            nobs(ar_fit), r2(ar_fit, "r2")))

# Extract within-sample AR(1) persistence ρ̄ (unweighted mean of γ_1 + θ_1·T̄)
# gives intuition for how quickly a 1°C shock fades
b_ar  <- coef(ar_fit)
T_hot  <- 25; T_cold <- 5; T_bliss <- 13
rho_hot   <- b_ar["tL1"] + b_ar["tLxT1"] * T_hot
rho_cold  <- b_ar["tL1"] + b_ar["tLxT1"] * T_cold
rho_bliss <- b_ar["tL1"] + b_ar["tLxT1"] * T_bliss
cat(sprintf("  Lag-1 temp persistence  ρ̂ (T̄=25°C hot):   %.3f\n", rho_hot))
cat(sprintf("  Lag-1 temp persistence  ρ̂ (T̄=13°C bliss): %.3f\n", rho_bliss))
cat(sprintf("  Lag-1 temp persistence  ρ̂ (T̄= 5°C cold):  %.3f\n", rho_cold))

# Add shock to dataset
dta[, tau := NA_real_]
dta[obs(ar_fit), tau := resid(ar_fit)]
dta[, tau_x_Tbar := tau * T_bar]
cat(sprintf("  Shock obs available: %d\n", sum(!is.na(dta$tau))))

# ============================================================
# 4. CONTEMPORANEOUS COMPARISON: BHM vs NRK specification
#    Replicates NRK Table 3 logic on BHM data.
#    Controls: 3 temp lags + 3 GDP growth lags + country FE + year FE
#    (NRK use country + year FE, not country time trends)
# ============================================================
ctrl <- "tL1 + tL2 + tL3 + tLxT1 + tLxT2 + tLxT3 + gL1 + gL2 + gL3"

# Col 1: BHM quadratic (with NRK controls, for clean comparison)
fit_bhm <- feols(
  as.formula(paste("growthWDI ~ temp + temp2 + prec + prec2 +", ctrl, "| iso_id + year")),
  cluster = ~iso_id, data = dta)

# Col 2: NRK state-dependent only (no quadratic)
fit_nrk <- feols(
  as.formula(paste("growthWDI ~ tau + tau_x_Tbar + prec + prec2 +", ctrl, "| iso_id + year")),
  cluster = ~iso_id, data = dta)

# Col 3: Hybrid (both quadratic and state-dependent)
fit_hyb <- feols(
  as.formula(paste("growthWDI ~ temp + temp2 + tau + tau_x_Tbar + prec + prec2 +", ctrl,
                   "| iso_id + year")),
  cluster = ~iso_id, data = dta)

cat("\n=== Table: Contemporaneous GDP growth models ===\n")
cat(sprintf("%-22s  %-14s  %-14s  %-14s\n",
            "Variable", "BHM quadratic", "NRK state-dep", "Hybrid"))
cat(strrep("-", 70), "\n")

print_row <- function(name, var, fits) {
  vals <- sapply(fits, function(f) {
    b  <- tryCatch(coef(f)[var], error = function(e) NA)
    se <- tryCatch(sqrt(vcov(f)[var, var]), error = function(e) NA)
    if (is.na(b)) return("      -       ")
    t  <- b / se
    st <- if (abs(t) > 2.576) "***" else if (abs(t) > 1.960) "**" else
          if (abs(t) > 1.645) "*" else ""
    sprintf("%9.5f%3s  ", b, st)
  })
  cat(sprintf("%-22s  %s\n", name, paste(vals, collapse = "")))
  # print SE row
  se_vals <- sapply(fits, function(f) {
    se <- tryCatch(sqrt(vcov(f)[var, var]), error = function(e) NA)
    if (is.na(se)) return("               ")
    sprintf("  (%7.5f)   ", se)
  })
  cat(sprintf("%-22s  %s\n", "", paste(se_vals, collapse = "")))
}

print_row("temp",    "temp",        list(fit_bhm, fit_nrk, fit_hyb))
print_row("temp²",   "temp2",       list(fit_bhm, fit_nrk, fit_hyb))
print_row("τ (shock)", "tau",       list(fit_bhm, fit_nrk, fit_hyb))
print_row("τ × T̄",  "tau_x_Tbar", list(fit_bhm, fit_nrk, fit_hyb))

cat(sprintf("%-22s  %-14s  %-14s  %-14s\n",
            "N", nobs(fit_bhm), nobs(fit_nrk), nobs(fit_hyb)))

# Implied bliss / optimal temperatures
b_bhm  <- coef(fit_bhm)
b_nrk  <- coef(fit_nrk)
b_hyb  <- coef(fit_hyb)
opt_bhm  <- -b_bhm["temp"]  / (2 * b_bhm["temp2"])
bliss_nrk <- -b_nrk["tau"]   / b_nrk["tau_x_Tbar"]
bliss_hyb <- -b_hyb["tau"]   / b_hyb["tau_x_Tbar"]
cat(sprintf("\nOptimal/bliss temperature: BHM=%.1f°C  NRK=%.1f°C  Hybrid=%.1f°C\n",
            opt_bhm, bliss_nrk, bliss_hyb))

# Marginal effects at 20°C
me_bhm  <- b_bhm["temp"]  + 2 * b_bhm["temp2"]  * 20
me_nrk  <- b_nrk["tau"]   + b_nrk["tau_x_Tbar"] * 20
me_hyb  <- b_hyb["tau"]   + b_hyb["tau_x_Tbar"] * 20
cat(sprintf("ME @ T̄=20°C (pp/year):    BHM=%.4f   NRK=%.4f   Hybrid=%.4f\n",
            me_bhm, me_nrk, me_hyb))

# ============================================================
# 5. LOCAL PROJECTIONS — Level vs Growth Test
#    For each horizon h, regress y_{t+h} - y_{t-1} on shock τ_t
#    (cumulative GDP change from period before shock to h periods after)
#
#    LHS construction: y_{t+h} − y_{t−1} ≈ Σ_{j=0}^{h} growthWDI_{i,t+j}
#    This is the cumulative log-GDP change over the [t, t+h] window.
#
#    Controls are all dated at or before t (pre-determined): 3 lags of temp
#    (×1 and ×T̄), 3 lags of GDP growth.  Country + year FEs absorb
#    country-specific levels and common year effects.
#
#    LEVEL EFFECT:  β^h roughly constant across h → GDP shifts once, stays there
#    GROWTH EFFECT: β^h grows ∝ h+1 → GDP keeps declining each year
#    TRANSITORY:    β^h → 0 as h increases → GDP recovers
# ============================================================

# Construct cumulative GDP change at each horizon h
for (h in 0:H) {
  leads <- lapply(0:h, function(j) shift(dta$growthWDI, j, type = "lead"))
  dta[[paste0("cumgdp_", h)]] <- Reduce("+", leads)
}

# Construct temperature level at h (for temp IRF estimation)
for (h in 1:H) {
  dta[[paste0("temp_Fh", h)]] <- shift(dta$temp, h, type = "lead")
}

# Storage for LP results
lp_res <- data.frame(
  h         = 0:H,
  # GDP IRF coefficients
  b0_hot    = NA_real_,   # β₀^h + β₁^h × 25  (hot country)
  b0_bliss  = NA_real_,   # β₀^h + β₁^h × 13
  b0_cold   = NA_real_,   # β₀^h + β₁^h × 5
  se_hot    = NA_real_,
  se_bliss  = NA_real_,
  se_cold   = NA_real_,
  # Annual-equivalent effect (tests growth vs level)
  ann_hot   = NA_real_,   # b0_hot  / (h+1) — if constant, pure growth effect
  ann_bliss = NA_real_,
  ann_cold  = NA_real_,
  # Temp IRF (to compute CRR)
  alpha_hot  = NA_real_,   # α₀^h + α₁^h × T̄
  alpha_bliss= NA_real_,
  alpha_cold = NA_real_,
  nobs      = NA_integer_
)

lp_fml_base <- paste("~ tau + tau_x_Tbar + prec + prec2 +", ctrl, "| iso_id + year")

cat("\n=== Local Projections: Cumulative GDP response to 1°C shock ===\n")
cat(sprintf("%-4s  %-12s  %-12s  %-12s  %-12s  %-12s  %-12s  %-8s\n",
            "h", "b_hot(SE)", "b_bliss(SE)", "b_cold(SE)",
            "ann_hot", "ann_bliss", "ann_cold", "N"))
cat(strrep("-", 90), "\n")

for (i in seq_len(nrow(lp_res))) {
  h   <- lp_res$h[i]
  lhs <- paste0("cumgdp_", h)

  # GDP local projection
  fml_gdp  <- as.formula(paste(lhs, lp_fml_base))
  fit_gdp  <- tryCatch(
    feols(fml_gdp, cluster = ~iso_id, data = dta),
    error = function(e) NULL)

  if (!is.null(fit_gdp)) {
    b0 <- coef(fit_gdp)["tau"]
    b1 <- coef(fit_gdp)["tau_x_Tbar"]
    V  <- vcov(fit_gdp)[c("tau","tau_x_Tbar"), c("tau","tau_x_Tbar")]

    irf <- function(Tbar) b0 + b1 * Tbar
    se_irf <- function(Tbar) sqrt(max(0, V[1,1] + Tbar^2 * V[2,2] + 2 * Tbar * V[1,2]))

    lp_res$b0_hot[i]   <- irf(T_hot);   lp_res$se_hot[i]   <- se_irf(T_hot)
    lp_res$b0_bliss[i] <- irf(T_bliss); lp_res$se_bliss[i] <- se_irf(T_bliss)
    lp_res$b0_cold[i]  <- irf(T_cold);  lp_res$se_cold[i]  <- se_irf(T_cold)
    lp_res$ann_hot[i]  <- irf(T_hot)   / (h + 1)
    lp_res$ann_bliss[i]<- irf(T_bliss) / (h + 1)
    lp_res$ann_cold[i] <- irf(T_cold)  / (h + 1)
    lp_res$nobs[i]     <- nobs(fit_gdp)
  }

  # Temperature IRF (h=0 is normalized to 1 by construction)
  if (h == 0) {
    lp_res$alpha_hot[i]   <- 1
    lp_res$alpha_bliss[i] <- 1
    lp_res$alpha_cold[i]  <- 1
  } else {
    lhs_t <- paste0("temp_Fh", h)
    fml_t <- as.formula(paste(lhs_t, lp_fml_base))
    fit_t <- tryCatch(
      feols(fml_t, cluster = ~iso_id, data = dta),
      error = function(e) NULL)
    if (!is.null(fit_t)) {
      a0 <- coef(fit_t)["tau"]
      a1 <- coef(fit_t)["tau_x_Tbar"]
      lp_res$alpha_hot[i]   <- a0 + a1 * T_hot
      lp_res$alpha_bliss[i] <- a0 + a1 * T_bliss
      lp_res$alpha_cold[i]  <- a0 + a1 * T_cold
    }
  }

  cat(sprintf("h=%d  %+.4f(%.4f)  %+.4f(%.4f)  %+.4f(%.4f)  %+.4f  %+.4f  %+.4f  %d\n",
              h,
              lp_res$b0_hot[i],   lp_res$se_hot[i],
              lp_res$b0_bliss[i], lp_res$se_bliss[i],
              lp_res$b0_cold[i],  lp_res$se_cold[i],
              lp_res$ann_hot[i],  lp_res$ann_bliss[i], lp_res$ann_cold[i],
              lp_res$nobs[i]))
}

# ============================================================
# 6. LEVEL VS GROWTH SUMMARY
#
#    Cumulative Response Ratios (CRR) at H=5
#    CRR_i = Σ_{h=0}^H β^h_i / Σ_{h=0}^H α^h_i
#    = cumulative GDP effect / cumulative temperature change
#    Interprets β as: "GDP effect per 1°C sustained temperature increase"
#    BHM's long-run damage function assumes CRR is permanent (infinite horizon).
# ============================================================
cat("\n=== Level vs Growth Tests ===\n")

for (Tbar_name in c("hot (T̄=25)", "bliss (T̄=13)", "cold (T̄=5)")) {
  col_b <- switch(Tbar_name,
                  "hot (T̄=25)"    = "b0_hot",
                  "bliss (T̄=13)"  = "b0_bliss",
                  "cold (T̄=5)"   = "b0_cold")
  col_a <- switch(Tbar_name,
                  "hot (T̄=25)"    = "alpha_hot",
                  "bliss (T̄=13)"  = "alpha_bliss",
                  "cold (T̄=5)"   = "alpha_cold")
  col_ann <- switch(Tbar_name,
                    "hot (T̄=25)"   = "ann_hot",
                    "bliss (T̄=13)" = "ann_bliss",
                    "cold (T̄=5)"  = "ann_cold")

  betas  <- lp_res[[col_b]]
  alphas <- lp_res[[col_a]]
  anns   <- lp_res[[col_ann]]

  cum_gdp  <- sum(betas,  na.rm = TRUE)
  cum_temp <- sum(alphas, na.rm = TRUE)
  crr      <- if (!is.na(cum_temp) && cum_temp != 0) cum_gdp / cum_temp else NA

  # Test for growth vs level: compare β^H / β^0
  ratio_H0 <- if (!is.na(betas[1]) && betas[1] != 0) betas[H + 1] / betas[1] else NA

  cat(sprintf("\n  Country type: %s\n", Tbar_name))
  cat(sprintf("    h=0 (contemporaneous GDP effect):  %+.4f pp\n", betas[1]))
  cat(sprintf("    h=5 (5-yr cumulative GDP effect):  %+.4f pp\n", betas[H+1]))
  cat(sprintf("    Ratio β^5 / β^0:  %.2f  (pure growth=%.1f, level=1.0)\n",
              ratio_H0, H + 1))
  cat(sprintf("    Cumul. temp IRF Σα^h:  %.3f  (1=fully transitory)\n", cum_temp))
  cat(sprintf("    CRR (GDP/temp):        %+.4f pp per 1°C sustained\n", crr))
  cat(sprintf("    Annual-equiv trend: h=0:%+.4f  h=3:%+.4f  h=5:%+.4f (decline→level)\n",
              anns[1], anns[4], anns[6]))
}

# ============================================================
# 7. BHM COMPARISON
#    BHM long-run projection for hot country (T̄=25°C):
#    Assumes temperature LEVEL permanently changes, so annual growth
#    effect accumulates forever. This implies GDP effect after T years = T × ME.
# ============================================================
cat("\n=== BHM vs NRK Projected 5-Year GDP Impact (hot country, T̄=25°C) ===\n")
me_bhm_hot <- b_bhm["temp"] + 2 * b_bhm["temp2"] * T_hot
me_nrk_h0  <- lp_res$b0_hot[1]
crr_hot    <- sum(lp_res$b0_hot, na.rm=T) / sum(lp_res$alpha_hot, na.rm=T)

cat(sprintf("  BHM annual ME at T=25°C:    %+.4f pp/year\n", me_bhm_hot))
cat(sprintf("  BHM 5-yr cumulative effect: %+.4f pp (assumes level shift each year)\n",
            5 * me_bhm_hot))
cat(sprintf("  NRK LP h=0 effect:          %+.4f pp\n", me_nrk_h0))
cat(sprintf("  NRK LP h=5 cumulative:      %+.4f pp\n", lp_res$b0_hot[H+1]))
cat(sprintf("  NRK CRR (5yr):              %+.4f pp per °C sustained\n\n", crr_hot))

# ============================================================
# 8. PLOT: IMPULSE RESPONSE FUNCTIONS
# ============================================================
pdf("figures/ExtendedDataFigs_Input/NRK_LocalProjections.pdf", width = 11, height = 8)
par(mfrow = c(2, 3), mar = c(4, 4.5, 3.5, 1.5), oma = c(0, 0, 2, 0))

z90 <- qnorm(0.95)   # 90% CI

# --- Top row: Cumulative GDP IRFs for hot / bliss / cold ---
plot_irf <- function(col_b, col_se, title, col, ylims = NULL) {
  b   <- lp_res[[col_b]]
  se  <- lp_res[[col_se]]
  hs  <- lp_res$h
  yl  <- if (is.null(ylims)) range(c(b - z90*se, b + z90*se), na.rm=TRUE)*1.15 else ylims
  plot(hs, b, type = "l", lwd = 2.5, col = col,
       xlab = "Horizon h (years)", ylab = "Cumulative GDP effect (pp)",
       main = title, ylim = yl, las = 1, xaxt = "n")
  axis(1, at = 0:H)
  polygon(c(hs, rev(hs)), c(b - z90*se, rev(b + z90*se)),
          col = adjustcolor(col, 0.2), border = NA)
  abline(h = 0, lty = 2, col = "grey50")
  # Reference lines: what constant-ME growth effect looks like
  abline(a = b[1], b = b[1], lty = 3, col = "grey30", lwd = 1)
  legend("topleft", c("LP estimate", "90% CI", "Growth effect ref."),
         lty = c(1, NA, 3), pch = c(NA, 15, NA),
         col = c(col, adjustcolor(col, 0.2), "grey30"),
         lwd = c(2.5, 10, 1), cex = 0.75, bg = "white")
}

ylim_hot   <- c(min(lp_res$b0_hot   - z90*lp_res$se_hot,   na.rm=T) * 1.2,
                max(lp_res$b0_hot   + z90*lp_res$se_hot,   na.rm=T) * 1.2)
ylim_bliss <- c(min(lp_res$b0_bliss - z90*lp_res$se_bliss, na.rm=T) * 1.2,
                max(lp_res$b0_bliss + z90*lp_res$se_bliss, na.rm=T) * 1.2)
ylim_cold  <- c(min(lp_res$b0_cold  - z90*lp_res$se_cold,  na.rm=T) * 1.2,
                max(lp_res$b0_cold  + z90*lp_res$se_cold,  na.rm=T) * 1.2)

plot_irf("b0_hot",   "se_hot",   "Hot country (Tbar=25 C)", "#C0392B", ylim_hot)
plot_irf("b0_bliss", "se_bliss", "Bliss-point country (Tbar=13 C)", "#2980B9", ylim_bliss)
plot_irf("b0_cold",  "se_cold",  "Cold country (Tbar=5 C)", "#27AE60", ylim_cold)

# --- Bottom row: Level vs Growth Tests ---

# Panel 4: Annual-equivalent effect β^h / (h+1)
# If constant → pure growth effect; declining → fading (level-like)
hs <- lp_res$h
plot(hs, lp_res$ann_hot, type = "b", pch = 16, col = "#C0392B", lwd = 2,
     xlab = "Horizon h", ylab = "b^h / (h+1)  (pp per year)",
     main = "Annual-Equiv. Effect b^h/(h+1)\n(flat = growth, declining = level)",
     ylim = range(c(lp_res$ann_hot, lp_res$ann_bliss, lp_res$ann_cold), na.rm=T)*c(1.3,1.1),
     las = 1, xaxt = "n")
axis(1, at = 0:H)
lines(hs, lp_res$ann_bliss, type = "b", pch = 17, col = "#2980B9", lwd = 2)
lines(hs, lp_res$ann_cold,  type = "b", pch = 15, col = "#27AE60", lwd = 2)
abline(h = 0, lty = 2, col = "grey50")
legend("topright", c("Hot (Tbar=25 C)", "Bliss (Tbar=13 C)", "Cold (Tbar=5 C)"),
       col = c("#C0392B", "#2980B9", "#27AE60"), lty = 1, pch = c(16,17,15), lwd = 2, cex = 0.8)

# Panel 5: Temperature IRF (persistence of shock itself)
alpha_df <- data.frame(h = 0:H,
                       hot   = lp_res$alpha_hot,
                       bliss = lp_res$alpha_bliss,
                       cold  = lp_res$alpha_cold)
yran <- range(c(alpha_df$hot, alpha_df$bliss, alpha_df$cold), na.rm = TRUE)
plot(alpha_df$h, alpha_df$hot, type = "b", pch = 16, col = "#C0392B", lwd = 2,
     xlab = "Horizon h", ylab = "Temperature level at t+h given shock at t",
     main = "Temperature Persistence IRF\n(from AR model, by country type)",
     ylim = c(0, max(yran, na.rm=T) * 1.15), las = 1, xaxt = "n")
axis(1, at = 0:H)
lines(alpha_df$h, alpha_df$bliss, type = "b", pch = 17, col = "#2980B9", lwd = 2)
lines(alpha_df$h, alpha_df$cold,  type = "b", pch = 15, col = "#27AE60", lwd = 2)
abline(h = 0, lty = 2, col = "grey50")
legend("topright", c("Hot (Tbar=25 C)", "Bliss (Tbar=13 C)", "Cold (Tbar=5 C)"),
       col = c("#C0392B", "#2980B9", "#27AE60"), lty = 1, pch = c(16,17,15), lwd = 2, cex = 0.8)

# Panel 6: Response function at h=0 (NRK contemporaneous)
Tbar_seq <- seq(0, 30, by = 0.5)
b0_h0  <- coef(fit_nrk)["tau"]
b1_h0  <- coef(fit_nrk)["tau_x_Tbar"]
V_h0   <- vcov(fit_nrk)[c("tau","tau_x_Tbar"), c("tau","tau_x_Tbar")]
me_seq  <- b0_h0 + b1_h0 * Tbar_seq
se_seq  <- sqrt(pmax(0, V_h0[1,1] + Tbar_seq^2 * V_h0[2,2] + 2 * Tbar_seq * V_h0[1,2]))

# BHM response for comparison (at T = T̄ for a representative "permanent" year)
bhm_seq <- b_bhm["temp"] + 2 * b_bhm["temp2"] * Tbar_seq

plot(Tbar_seq, me_seq * 100, type = "l", lwd = 2.5, col = "#8E44AD",
     xlab = "Country mean temperature Tbar (deg C)",
     ylab = "Marginal effect of 1°C shock (pp/year)",
     main = "Contemporaneous Response Function\n(NRK shock vs BHM)",
     ylim = range(c((me_seq - z90*se_seq)*100, (me_seq + z90*se_seq)*100,
                     bhm_seq*100), na.rm=T) * c(1.15, 1.05),
     las = 1)
polygon(c(Tbar_seq, rev(Tbar_seq)),
        c((me_seq - z90*se_seq)*100, rev((me_seq + z90*se_seq)*100)),
        col = adjustcolor("#8E44AD", 0.15), border = NA)
lines(Tbar_seq, bhm_seq * 100, lwd = 2, lty = 2, col = "steelblue")
abline(h = 0, lty = 3, col = "grey50")
legend("topright", c("NRK shock x Tbar", "BHM temp + temp^2"),
       col = c("#8E44AD","steelblue"), lty = c(1,2), lwd = 2, cex = 0.85)

mtext("NRK-Inspired Local Projections: BHM Data 1960-2010 (Nath, Ramey, Klenow approach)",
      outer = TRUE, cex = 1.1, font = 2)
dev.off()
cat("Saved figures/ExtendedDataFigs_Input/NRK_LocalProjections.pdf\n")

# ============================================================
# 9. SAVE LP RESULTS
# ============================================================
write.csv(lp_res, "data/output/NRK_LP_results.csv", row.names = FALSE)
cat("Saved data/output/NRK_LP_results.csv\n")

cat("\nNRK_LocalProjections.R complete.\n")
