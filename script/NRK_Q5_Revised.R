
# NRK_Q5_Revised.R
# Question 5 (Revised): NRK-Inspired Critique of BHM
#
# Three upgrades over the original Q5:
#   (i)  Temperature innovations/shocks  — AR(3) innovations tau_it rather than levels
#   (ii) Distributed lags 0-5           — lag profile of shock effects on GDP growth
#   (iii) Growth vs. levels comparison  — both DeltalnGDP (growth) and lnGDP (levels)
#          tested via local projections (Jorda 2005)
#
# Models estimated:
#   Model 1: Distributed lag shock regression (growth outcome)
#            Dy_it = sum_{k=0}^{5} [delta_k * tau_{i,t-k}
#                                 + delta'_k * (tau_{i,t-k} x T_bar_i)]
#                   + controls + FE + eps_it
#
#   Model 2: Local projections (levels-vs-growth test)
#            y_{i,t+h} - y_{i,t-1} = beta_h * tau_it
#                                   + beta'_h * (tau_it x T_bar_i)
#                                   + sum_{k=1}^{3} g_{i,t-k} + FE + eps
#            h = 0 .. 10
#            At h=0: outcome = current growth rate (not constant)
#
#   Model 3: Levels regression
#            lnGDP_it = delta_1 * tau_it + delta_2 * (tau_it x T_bar_i)
#                      + prec + prec2 + FE + eps_it
#
# Outputs:
#   data/output/q5_distributed_lag_profile.csv
#   data/output/q5_local_projections.csv
#   data/output/q5_levels_regression.csv
#   figures/ExtendedDataFigs_Input/Q5_Revised.pdf

rm(list = ls())

library(fixest)
library(data.table)

setwd("/Users/yl3699/Downloads/BurkeHsiangMiguel2015_Replication")
dir.create("data/output", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/ExtendedDataFigs_Input", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. LOAD DATA
# ============================================================
dta <- as.data.table(read.csv("data/input/GrowthClimateDataset.csv"))
setorder(dta, iso, year)

dta[, temp  := UDel_temp_popweight]
dta[, prec  := UDel_precip_popweight]
dta[, prec2 := prec^2]
dta[, time  := year - 1985]
dta[, time2 := time^2]

cat("Loaded:", nrow(dta), "obs,", dta[, uniqueN(iso)], "countries\n")

# ============================================================
# 2. COUNTRY MEAN TEMPERATURE (30-year baseline, as in NRK)
# ============================================================
dta[, T_bar := mean(temp[year <= 1989], na.rm = TRUE), by = iso]
dta[is.na(T_bar), T_bar := mean(temp, na.rm = TRUE), by = iso]  # fallback

# ============================================================
# 3. AR(3) TEMPERATURE INNOVATIONS (NRK eq. 4)
#    T_it = sum_j gamma_j T_{i,t-j}
#          + sum_j theta_j T_{i,t-j} * T_bar_i
#          + mu_i + mu_t + tau_it
# ============================================================
for (j in 1:3) {
  dta[, paste0("temp_L", j) := shift(temp, j), by = iso]
}

ar_fit <- feols(
  temp ~ temp_L1 + temp_L2 + temp_L3 +
         I(temp_L1 * T_bar) + I(temp_L2 * T_bar) + I(temp_L3 * T_bar) |
         iso + year,
  data = dta
)
cat(sprintf("AR model: N=%d, R2=%.4f\n", nobs(ar_fit), r2(ar_fit, "r2")))

dta[, tau := NA_real_]
dta[obs(ar_fit), tau := resid(ar_fit)]
cat("AR innovations (tau) available:", sum(!is.na(dta$tau)), "\n")

# State-dependence interaction
dta[, tau_x_Tbar := tau * T_bar]

# ============================================================
# 4. LOG GDP LEVEL and LAGGED GROWTH RATES
#    lnGDP constructed for Model 3 (levels regression).
#    LP (Model 2) uses rolling sum of actual growthWDI (not cumsum of filled NAs),
#    so NAs in intermediate years correctly propagate as NA in LP outcome.
# ============================================================
dta[, g := growthWDI]  # shorthand (keep NAs for LP — don't fill with zeros)

# lnGDP for Model 3: cumulative sum with NAs treated as 0 (unavoidable approximation)
dta[, g_fill := fifelse(is.na(g), 0, g)]
dta[, lnGDP := cumsum(g_fill) - first(cumsum(g_fill)), by = iso]

# Lagged growth rates for LP controls (k=1,2,3)
for (k in 1:3) {
  dta[, paste0("g_L", k) := shift(g, k), by = iso]
}

# ============================================================
# 5. DISTRIBUTED LAG SHOCK VARIABLES (k = 0 .. 5)
# ============================================================
for (k in 0:5) {
  dta[, paste0("tau_L", k)       := shift(tau, k), by = iso]
  dta[, paste0("tau_xTbar_L", k) := shift(tau_x_Tbar, k), by = iso]
}

# ============================================================
# 6. MODEL 1: DISTRIBUTED LAG SHOCK REGRESSION (growth outcome)
#    Dy_it = sum_{k=0}^{5} [delta_k tau_{t-k} + delta'_k (tau_{t-k} x T_bar)]
#            + prec + prec2 + FE + eps
# ============================================================
lag_vars <- paste0("tau_L", 0:5)
lag_int_vars <- paste0("tau_xTbar_L", 0:5)
all_lag_vars <- c(lag_vars, lag_int_vars)
fe_str <- "iso_id[time, time2] + year"

fml_dl <- as.formula(
  paste("growthWDI ~",
        paste(c(all_lag_vars, "prec", "prec2"), collapse = " + "),
        "|", fe_str)
)

fit_dl <- feols(fml_dl, cluster = ~iso_id, data = dta)
cat("\n=== Model 1: Distributed Lag (growth outcome) ===\n")
print(summary(fit_dl, keep = c("tau_L", "tau_xTbar_L")))

# Extract lag profile
dl_coefs <- coef(fit_dl)
dl_se    <- se(fit_dl)
dl_tstat <- tstat(fit_dl)

lag_profile <- data.table(
  lag = 0:5,
  delta_k      = as.numeric(dl_coefs[paste0("tau_L", 0:5)]),
  se_delta_k   = as.numeric(dl_se[paste0("tau_L", 0:5)]),
  delta_prime_k    = as.numeric(dl_coefs[paste0("tau_xTbar_L", 0:5)]),
  se_delta_prime_k = as.numeric(dl_se[paste0("tau_xTbar_L", 0:5)])
)

# Cumulative effects (sum of delta_k for k=0..lag)
lag_profile[, cumulative_delta := cumsum(delta_k)]

# Reference temperatures: cold=5, moderate=15, hot=25
T_hot  <- 25
T_mod  <- 15
T_cold <- 5
lag_profile[, ME_hot  := delta_k + delta_prime_k * T_hot]
lag_profile[, ME_mod  := delta_k + delta_prime_k * T_mod]
lag_profile[, ME_cold := delta_k + delta_prime_k * T_cold]

cat("\n--- Distributed lag profile ---\n")
print(lag_profile, digits = 4)

fwrite(lag_profile, "data/output/q5_distributed_lag_profile.csv")
cat("Saved: data/output/q5_distributed_lag_profile.csv\n")

# ============================================================
# 7. MODEL 2: LOCAL PROJECTIONS (Jorda 2005)
#    Outcome: forward rolling sum of growth rates
#      lp_y_{h} = sum_{s=0}^{h} growthWDI_{i,t+s}
#    This is h+1 years of cumulative growth starting at t.
#    At h=0: just contemporaneous growth (not constant).
#    NAs in any intermediate growthWDI make lp_y_{h} = NA.
#
#    Regression:
#      lp_y_{h} = beta_h * tau_it + beta'_h * (tau_it x T_bar_i)
#                + g_L1 + g_L2 + g_L3 + prec + prec2 + FE + eps
#
#    Level vs growth test:
#      Growth effect: IRF_h grows proportionally with h
#      Level effect:  IRF_h plateaus (each new year adds ~0)
#      Reverting:     IRF_h returns toward zero
#
#    IRF at T_bar = T_val: beta_h + beta'_h * T_val
#    SE via delta method:
#      sqrt(Var(b) + T^2 * Var(b') + 2*T*Cov(b,b'))
# ============================================================
cat("\n=== Model 2: Local Projections ===\n")

# Build forward rolling sums of growthWDI (using actual NAs, not filled)
# lp_fwd_h = growthWDI_t + growthWDI_{t+1} + ... + growthWDI_{t+h}
lp_results <- list()

for (h in 0:10) {
  # Sum of leads 0..h of growthWDI
  lead_mat <- sapply(0:h, function(s) {
    dta[, shift(growthWDI, s, type = "lead"), by = iso]$V1
  })
  if (h == 0) {
    lp_y <- lead_mat
  } else {
    lp_y <- rowSums(lead_mat, na.rm = FALSE)
  }

  tmp_dt <- copy(dta)
  tmp_dt[, lp_outcome := lp_y]
  tmp_dt <- tmp_dt[
    !is.na(lp_outcome) & !is.na(tau) & !is.na(tau_x_Tbar) &
    !is.na(g_L1) & !is.na(g_L2) & !is.na(g_L3)
  ]

  fml_lp <- as.formula(paste(
    "lp_outcome ~ tau + tau_x_Tbar + g_L1 + g_L2 + g_L3 + prec + prec2 |",
    fe_str
  ))

  fit_lp <- tryCatch(
    feols(fml_lp, cluster = ~iso_id, data = tmp_dt),
    error = function(e) {
      cat("h=", h, "error:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(fit_lp)) next

  cf  <- coef(fit_lp)
  vc  <- vcov(fit_lp)
  b   <- cf["tau"]
  bp  <- cf["tau_x_Tbar"]

  irf_cold <- b + bp * T_cold
  irf_mod  <- b + bp * T_mod
  irf_hot  <- b + bp * T_hot

  var_b   <- vc["tau", "tau"]
  var_bp  <- vc["tau_x_Tbar", "tau_x_Tbar"]
  cov_bbp <- vc["tau", "tau_x_Tbar"]

  se_cold <- sqrt(var_b + T_cold^2 * var_bp + 2 * T_cold * cov_bbp)
  se_mod  <- sqrt(var_b + T_mod^2  * var_bp + 2 * T_mod  * cov_bbp)
  se_hot  <- sqrt(var_b + T_hot^2  * var_bp + 2 * T_hot  * cov_bbp)

  lp_results[[h + 1]] <- data.table(
    h            = h,
    N            = nobs(fit_lp),
    beta_h       = b,
    beta_prime_h = bp,
    irf_cold = irf_cold, se_cold = se_cold, t_cold = irf_cold / se_cold,
    irf_mod  = irf_mod,  se_mod  = se_mod,  t_mod  = irf_mod  / se_mod,
    irf_hot  = irf_hot,  se_hot  = se_hot,  t_hot  = irf_hot  / se_hot
  )

  cat(sprintf(
    "  h=%2d: N=%5d  cold=%+.5f (t=%.2f)  mod=%+.5f (t=%.2f)  hot=%+.5f (t=%.2f)\n",
    h, nobs(fit_lp),
    irf_cold, irf_cold / se_cold,
    irf_mod,  irf_mod  / se_mod,
    irf_hot,  irf_hot  / se_hot
  ))
}

lp_dt <- rbindlist(lp_results)
print(lp_dt, digits = 4)
fwrite(lp_dt, "data/output/q5_local_projections.csv")
cat("Saved: data/output/q5_local_projections.csv\n")

# ============================================================
# 8. MODEL 3: LEVELS REGRESSION
#    lnGDP_it = delta_1 * tau_it + delta_2 * (tau_it x T_bar_i)
#              + prec + prec2 + iso_id + year + eps
#
#    Uses simpler FEs (country + year, no time trends) because
#    lnGDP is itself trend-driven; quadratic trends would absorb
#    nearly all within-country variation, causing near-collinearity.
# ============================================================
cat("\n=== Model 3: Levels Regression ===\n")

fit_lev <- feols(
  lnGDP ~ tau + tau_x_Tbar + prec + prec2 | iso_id + year,
  cluster = ~iso_id, data = dta[!is.na(tau)]
)
print(summary(fit_lev))

lev_cf  <- coef(fit_lev)
lev_se  <- se(fit_lev)
lev_t   <- tstat(fit_lev)

levels_res <- data.table(
  term      = names(lev_cf),
  estimate  = as.numeric(lev_cf),
  std_error = as.numeric(lev_se),
  t_stat    = as.numeric(lev_t)
)
fwrite(levels_res, "data/output/q5_levels_regression.csv")
cat("Saved: data/output/q5_levels_regression.csv\n")

# Bliss point (T_bar at which contemporaneous effect = 0)
d1 <- lev_cf["tau"]
d2 <- lev_cf["tau_x_Tbar"]
bliss_lev <- -d1 / d2
cat(sprintf("\nLevels bliss T_bar = %.1f C\n", bliss_lev))
cat(sprintf("ME at T_bar=20C:  %.5f\n", d1 + d2 * 20))
cat(sprintf("ME at T_bar=28C:  %.5f\n", d1 + d2 * 28))

# ============================================================
# 9. ALSO RUN BASELINE SHOCK REGRESSION (contemporaneous, h=0)
#    for reference in results table
# ============================================================
fit_base_shock <- feols(
  as.formula(paste("growthWDI ~ tau + tau_x_Tbar + prec + prec2 |", fe_str)),
  cluster = ~iso_id, data = dta
)
b0_base <- coef(fit_base_shock)
cat(sprintf("\nBaseline shock (h=0): delta1=%.5f (t=%.2f)  delta2=%.5f (t=%.2f)\n",
            b0_base["tau"],       tstat(fit_base_shock)["tau"],
            b0_base["tau_x_Tbar"], tstat(fit_base_shock)["tau_x_Tbar"]))
cat(sprintf("  Bliss T_bar = %.1f C\n", -b0_base["tau"] / b0_base["tau_x_Tbar"]))

# ============================================================
# 10. PLOT: Four-panel figure
#    Panel A: Lag profile delta_k (growth outcome) at T_hot and T_cold
#    Panel B: Cumulative lag profile
#    Panel C: Local projection IRF — cold vs hot countries
#    Panel D: IRF ratio hot/cold
# ============================================================
pdf("figures/ExtendedDataFigs_Input/Q5_Revised.pdf", width = 11, height = 9)
par(mfrow = c(2, 2), mar = c(4, 4.5, 3, 1), mgp = c(2.5, 0.8, 0))

# --- Panel A: Lag profile ME at three reference temperatures ---
ylim_a <- range(c(lag_profile$ME_hot, lag_profile$ME_mod,
                  lag_profile$ME_cold), na.rm = TRUE) * c(1.5, 1.5)
plot(lag_profile$lag, lag_profile$ME_cold, type = "b", pch = 16, col = "steelblue",
     lwd = 2, xlab = "Lag k (years)", ylab = "Marginal effect on growth (%)",
     main = "(A) Lag profile: ME = delta_k + delta'_k * T_bar",
     ylim = ylim_a, xaxt = "n")
axis(1, at = 0:5)
lines(lag_profile$lag, lag_profile$ME_mod,  type = "b", pch = 15, col = "darkgreen", lwd = 2)
lines(lag_profile$lag, lag_profile$ME_hot,  type = "b", pch = 17, col = "firebrick", lwd = 2)
abline(h = 0, lty = 2, col = "gray50")
legend("topright",
       legend = c(paste0("Cold     (", T_cold, "C)"),
                  paste0("Moderate (", T_mod,  "C)"),
                  paste0("Hot      (", T_hot,  "C)")),
       col = c("steelblue", "darkgreen", "firebrick"),
       pch = c(16, 15, 17), lwd = 2, bty = "n", cex = 0.80)

# --- Panel B: Cumulative lag profile ---
cum_cold <- cumsum(lag_profile$ME_cold)
cum_mod  <- cumsum(lag_profile$ME_mod)
cum_hot  <- cumsum(lag_profile$ME_hot)
ylim_b <- range(c(cum_cold, cum_mod, cum_hot), na.rm = TRUE) * c(1.5, 1.5)
plot(lag_profile$lag, cum_cold, type = "b", pch = 16, col = "steelblue",
     lwd = 2, xlab = "Lag k (years)", ylab = "Cumulative ME on growth (%)",
     main = "(B) Cumulative lag profile",
     ylim = ylim_b, xaxt = "n")
axis(1, at = 0:5)
lines(lag_profile$lag, cum_mod, type = "b", pch = 15, col = "darkgreen", lwd = 2)
lines(lag_profile$lag, cum_hot, type = "b", pch = 17, col = "firebrick", lwd = 2)
abline(h = 0, lty = 2, col = "gray50")
legend("topright",
       legend = c(paste0("Cold     (", T_cold, "C)"),
                  paste0("Moderate (", T_mod,  "C)"),
                  paste0("Hot      (", T_hot,  "C)")),
       col = c("steelblue", "darkgreen", "firebrick"),
       pch = c(16, 15, 17), lwd = 2, bty = "n", cex = 0.80)

# --- Panel C: LP IRF for three reference temperatures ---
if (nrow(lp_dt) > 0) {
  ci_lo_cold <- lp_dt$irf_cold - 1.96 * lp_dt$se_cold
  ci_hi_cold <- lp_dt$irf_cold + 1.96 * lp_dt$se_cold
  ci_lo_mod  <- lp_dt$irf_mod  - 1.96 * lp_dt$se_mod
  ci_hi_mod  <- lp_dt$irf_mod  + 1.96 * lp_dt$se_mod
  ci_lo_hot  <- lp_dt$irf_hot  - 1.96 * lp_dt$se_hot
  ci_hi_hot  <- lp_dt$irf_hot  + 1.96 * lp_dt$se_hot

  ylim_c <- range(c(ci_lo_cold, ci_hi_cold,
                    ci_lo_mod,  ci_hi_mod,
                    ci_lo_hot,  ci_hi_hot), na.rm = TRUE)
  ylim_c <- ylim_c + c(-1, 1) * 0.1 * diff(ylim_c)

  plot(lp_dt$h, lp_dt$irf_cold, type = "b", pch = 16, col = "steelblue", lwd = 2,
       xlab = "Horizon h (years)", ylab = "Cumulative log-GDP response",
       main = "(C) Local Projection IRF",
       ylim = ylim_c, xaxt = "n")
  axis(1, at = 0:10)

  polygon(c(lp_dt$h, rev(lp_dt$h)), c(ci_lo_cold, rev(ci_hi_cold)),
          col = adjustcolor("steelblue",  alpha.f = 0.15), border = NA)
  polygon(c(lp_dt$h, rev(lp_dt$h)), c(ci_lo_mod,  rev(ci_hi_mod)),
          col = adjustcolor("darkgreen",  alpha.f = 0.15), border = NA)
  polygon(c(lp_dt$h, rev(lp_dt$h)), c(ci_lo_hot,  rev(ci_hi_hot)),
          col = adjustcolor("firebrick",  alpha.f = 0.15), border = NA)

  lines(lp_dt$h, lp_dt$irf_mod,  type = "b", pch = 15, col = "darkgreen", lwd = 2)
  lines(lp_dt$h, lp_dt$irf_hot,  type = "b", pch = 17, col = "firebrick", lwd = 2)
  abline(h = 0, lty = 2, col = "gray50")

  legend("topright",
         legend = c(paste0("Cold     (", T_cold, "C)"),
                    paste0("Moderate (", T_mod,  "C)"),
                    paste0("Hot      (", T_hot,  "C)")),
         col = c("steelblue", "darkgreen", "firebrick"),
         pch = c(16, 15, 17), lwd = 2, bty = "n", cex = 0.80)
} else {
  plot.new(); text(0.5, 0.5, "No LP results", cex = 1.5)
}

# --- Panel D: Normalised IRF (level vs. growth effect) ---
if (nrow(lp_dt) > 0) {
  irf0_cold <- lp_dt$irf_cold[lp_dt$h == 0]
  irf0_mod  <- lp_dt$irf_mod[lp_dt$h  == 0]
  irf0_hot  <- lp_dt$irf_hot[lp_dt$h  == 0]

  norm_cold <- lp_dt$irf_cold / abs(irf0_cold)
  norm_mod  <- lp_dt$irf_mod  / abs(irf0_mod)
  norm_hot  <- lp_dt$irf_hot  / abs(irf0_hot)

  ylim_d <- range(c(norm_cold, norm_mod, norm_hot), na.rm = TRUE) * c(1.4, 1.4)
  ylim_d[is.infinite(ylim_d) | is.na(ylim_d)] <- c(-5, 5)

  plot(lp_dt$h, norm_cold, type = "b", pch = 16, col = "steelblue", lwd = 2,
       xlab = "Horizon h (years)", ylab = "Normalised IRF (impact = 1)",
       main = "(D) Normalised IRF: level vs. growth effect",
       ylim = ylim_d, xaxt = "n")
  axis(1, at = 0:10)
  lines(lp_dt$h, norm_mod,  type = "b", pch = 15, col = "darkgreen",  lwd = 2)
  lines(lp_dt$h, norm_hot,  type = "b", pch = 17, col = "firebrick",  lwd = 2)
  abline(h = 0, lty = 2, col = "gray50")
  abline(h = 1, lty = 3, col = "gray40")
  legend("topright", legend = c(paste0("Cold     (", T_cold, "C)"),
                                 paste0("Moderate (", T_mod,  "C)"),
                                 paste0("Hot      (", T_hot,  "C)"),
                                 "Growth effect (IRF grows)",
                                 "Level effect  (IRF = const)"),
         col = c("steelblue", "darkgreen", "firebrick", "gray40", "gray40"),
         pch = c(16, 15, 17, NA, NA), lty = c(1, 1, 1, 3, 2),
         lwd = 2, bty = "n", cex = 0.75)
} else {
  plot.new(); text(0.5, 0.5, "No LP results", cex = 1.5)
}

dev.off()
cat("Saved: figures/ExtendedDataFigs_Input/Q5_Revised.pdf\n")

# ============================================================
# 11. SUMMARY TABLE
# ============================================================
cat("\n============================================================\n")
cat("SUMMARY\n")
cat("============================================================\n")

cat("\nModel 1 — Distributed Lag (growth, k=0..5):\n")
cat(sprintf("  Contemporaneous (k=0): delta_0 = %.5f (t=%.2f), delta'_0 = %.5f (t=%.2f)\n",
            lag_profile$delta_k[1],       lag_profile$delta_k[1] / lag_profile$se_delta_k[1],
            lag_profile$delta_prime_k[1], lag_profile$delta_prime_k[1] / lag_profile$se_delta_prime_k[1]))
cat(sprintf(
  "  Cumulative (k=0..5): sum(delta_k)=%.5f, ME_cold=%.5f, ME_mod=%.5f, ME_hot=%.5f\n",
  sum(lag_profile$delta_k),
  sum(lag_profile$ME_cold), sum(lag_profile$ME_mod), sum(lag_profile$ME_hot)
))

cat("\nModel 2 — Local Projections (growth vs. level test):\n")
if (nrow(lp_dt) > 0) {
  cat(sprintf(
    "  h= 0: cold=%+.5f (t=%.2f), mod=%+.5f (t=%.2f), hot=%+.5f (t=%.2f)\n",
    lp_dt[h==0, irf_cold], lp_dt[h==0, t_cold],
    lp_dt[h==0, irf_mod],  lp_dt[h==0, t_mod],
    lp_dt[h==0, irf_hot],  lp_dt[h==0, t_hot]
  ))
  for (hh in c(1, 2, 3, 5, 10)) {
    row <- lp_dt[h == hh]
    if (nrow(row) == 0) next
    cat(sprintf(
      "  h=%2d: cold=%+.5f (t=%.2f), mod=%+.5f (t=%.2f), hot=%+.5f (t=%.2f)\n",
      hh,
      row$irf_cold, row$t_cold,
      row$irf_mod,  row$t_mod,
      row$irf_hot,  row$t_hot
    ))
  }
  # Level vs growth conclusion
  max_cold <- max(abs(lp_dt$irf_cold), na.rm = TRUE)
  end_cold <- abs(lp_dt$irf_cold[nrow(lp_dt)])
  cat(sprintf("\n  Cold IRF at h=10 vs peak: %.3f of peak => %s\n",
              end_cold / max_cold,
              ifelse(end_cold / max_cold < 0.5, "LEVEL effect (dissipates)", "GROWTH effect (persists)")))
  max_hot <- max(abs(lp_dt$irf_hot), na.rm = TRUE)
  end_hot  <- abs(lp_dt$irf_hot[nrow(lp_dt)])
  cat(sprintf("  Hot  IRF at h=10 vs peak: %.3f of peak => %s\n",
              end_hot / max_hot,
              ifelse(end_hot / max_hot < 0.5, "LEVEL effect (dissipates)", "GROWTH effect (persists)")))
}

cat("\nModel 3 — Levels Regression:\n")
cat(sprintf("  delta_1 = %.5f (t=%.2f), delta_2 = %.5f (t=%.2f)\n",
            lev_cf["tau"],       lev_t["tau"],
            lev_cf["tau_x_Tbar"], lev_t["tau_x_Tbar"]))
cat(sprintf("  Bliss T_bar = %.1f C\n", bliss_lev))

cat("\n============================================================\n")
cat("Done.\n")
