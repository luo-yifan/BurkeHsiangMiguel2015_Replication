
# ConleyAnalysis.R
# Re-estimate BHM baseline with Conley (2008) Spatial HAC standard errors.
# Uses the Frisch-Waugh-Lovell theorem: partial out FEs via fixest::demean(),
# then apply conleyreg to the demeaned (within-FE) data.
#
# Spatial cutoff: 200–2000 km in 100 km increments
# Time cutoff:    Inf (allow correlation across all time periods within a country)
# Kernel:         Bartlett (linear decay to cutoff)
# Coordinates:    Country centroids from ESRI shapefile

rm(list = ls())
library(fixest)
library(conleyreg)
library(sf)
library(data.table)

setwd("/Users/yl3699/Downloads/BurkeHsiangMiguel2015_Replication")
dir.create("figures/ExtendedDataFigs_Input", recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. LOAD AND PREPARE PANEL DATA
# ============================================================
dta <- as.data.table(read.csv("data/input/GrowthClimateDataset.csv"))
dta[, temp  := UDel_temp_popweight]
dta[, temp2 := temp^2]
dta[, prec  := UDel_precip_popweight]
dta[, prec2 := prec^2]
dta[, time  := year - 1985]
dta[, time2 := time^2]

# ============================================================
# 2. COUNTRY CENTROIDS FROM SHAPEFILE
# ============================================================
shp      <- st_read("data/input/shape/country.shp", quiet = TRUE)
ctrs     <- st_centroid(st_geometry(shp))
coords   <- data.frame(
  iso = as.character(shp$GMI_CNTRY),
  lat = st_coordinates(ctrs)[, 2],
  lon = st_coordinates(ctrs)[, 1],
  stringsAsFactors = FALSE
)

# Manual centroids for territories not in shapefile (small islands/territories)
extra <- data.frame(
  iso = c("ADO", "CHI", "HKG", "IMY", "KSV", "MNE", "TMP", "WBG"),
  lat = c(42.55,  49.37, 22.30, 54.24, 42.60, 42.71,  -8.87, 31.95),
  lon = c( 1.57,  -2.20,114.18, -4.48, 21.17, 19.37, 125.73, 35.30),
  stringsAsFactors = FALSE
)
coords <- rbind(coords, extra)
coords <- coords[!duplicated(coords$iso), ]   # keep first centroid per ISO

# Use match() to add lat/lon WITHOUT reordering rows (merge() reorders)
m <- match(dta$iso, coords$iso)
dta$lat <- coords$lat[m]
dta$lon  <- coords$lon[m]
cat("Countries with missing centroid:", sum(is.na(dta$lat) & !is.na(dta$growthWDI)), "\n")

# ============================================================
# 3. BASELINE FEOLS — COEFFICIENTS ARE FIXED ACROSS CUTOFFS
# ============================================================
fit_base <- feols(
  growthWDI ~ temp + temp2 + prec + prec2 | iso_id[time, time2] + year,
  cluster = ~iso_id, data = dta
)
cat("\n=== Baseline (clustered by country) ===\n")
cat(sprintf("beta1 (temp):  %.6f  SE = %.6f  t = %.3f\n",
            coef(fit_base)["temp"],  sqrt(vcov(fit_base)["temp",  "temp"]),
            coef(fit_base)["temp"]  / sqrt(vcov(fit_base)["temp",  "temp"])))
cat(sprintf("beta2 (temp2): %.7f  SE = %.7f  t = %.3f\n",
            coef(fit_base)["temp2"], sqrt(vcov(fit_base)["temp2", "temp2"]),
            coef(fit_base)["temp2"] / sqrt(vcov(fit_base)["temp2", "temp2"])))

# ============================================================
# 4. PARTIAL OUT FIXED EFFECTS VIA FRISCH-WAUGH-LOVELL
#    demean(feols_obj) demeanes ALL model variables (LHS + RHS)
#    by the FE structure, returning within-FE residuals
# ============================================================
cat("\nDemeaning within fixed effects...\n")
dm <- demean(fit_base)         # returns a matrix of demeaned variables
# Column order matches the feols model variables
cat("Demeaned variables:", colnames(dm), "\n")
cat("Dimensions:", dim(dm), "\n")

y_dm    <- dm[, "growthWDI"]
X_dm    <- dm[, c("temp", "temp2", "prec", "prec2"), drop = FALSE]

# Verify FWL: OLS on demeaned data recovers feols coefficients
b_fwl   <- solve(t(X_dm) %*% X_dm) %*% t(X_dm) %*% y_dm
e_fwl   <- as.numeric(y_dm - X_dm %*% b_fwl)
cat(sprintf("\nFWL check — beta1: %.6f  beta2: %.7f (should match feols)\n",
            b_fwl["temp", ], b_fwl["temp2", ]))

# ============================================================
# 5. ASSEMBLE ANALYSIS DATA WITH DEMEANED VARIABLES + COORDS
# ============================================================
# demean() drops NAs. Rebuild using the observations feols actually used.
# obs(fit_base) gives integer row indices into the original dta
obs_used  <- obs(fit_base)
dta_est   <- as.data.frame(dta)[obs_used, ]
dta_est$y_dm     <- y_dm
dta_est$temp_dm  <- X_dm[, "temp"]
dta_est$temp2_dm <- X_dm[, "temp2"]
dta_est$prec_dm  <- X_dm[, "prec"]
dta_est$prec2_dm <- X_dm[, "prec2"]
dta_est$e        <- e_fwl
dta_est$iso_num  <- as.integer(factor(dta_est$iso_id))  # numeric country ID for conleyreg

cat("Estimation obs:", nrow(dta_est), " | Missing lat:", sum(is.na(dta_est$lat)), "\n")
stopifnot(nrow(dta_est) == length(y_dm))  # sanity check

# ============================================================
# 6. CONLEY SEs VIA conleyreg ON FWL-DEMEANED DATA
#    Since FEs are already removed, this is pooled OLS on demeaned vars.
#    unit + time + lag_cutoff=Inf gives full within-unit serial correlation.
#    Bartlett spatial kernel decays linearly to zero at dist_cutoff.
#    intercept=FALSE because demeaning already removed the constant.
# ============================================================
b1 <- b_fwl["temp",  ]
b2 <- b_fwl["temp2", ]

cutoffs <- seq(200, 2000, by = 100)
results <- data.frame(
  cutoff   = cutoffs,
  se_temp  = NA_real_,
  se_temp2 = NA_real_,
  t_temp   = NA_real_,
  t_temp2  = NA_real_
)

cat("\nComputing Conley SEs for", length(cutoffs), "cutoffs...\n")

for (i in seq_along(cutoffs)) {
  d <- cutoffs[i]
  cat(sprintf("  cutoff = %4d km ...", d))
  # conleyreg requires numeric unit and time IDs
  fit_c <- conleyreg(
    formula     = y_dm ~ temp_dm + temp2_dm + prec_dm + prec2_dm,
    data        = dta_est,
    lat         = "lat",
    lon         = "lon",
    dist_cutoff = d,
    unit        = "iso_num",   # numeric country id (factor-converted)
    time        = "year",     # numeric year
    lag_cutoff  = Inf,        # allow correlation across all time periods
    kernel      = "bartlett",
    ncores      = 1,         # avoid NA core detection in constrained environments
    intercept   = FALSE,      # no constant in demeaned regression
    verbose     = FALSE
  )
  # conleyreg returns a coeftest matrix: rows=coefs, cols=Estimate/SE/t/p
  se_temp_c  <- fit_c["temp_dm",  "Std. Error"]
  se_temp2_c <- fit_c["temp2_dm", "Std. Error"]
  results$se_temp[i]  <- se_temp_c
  results$se_temp2[i] <- se_temp2_c
  results$t_temp[i]   <- b1 / se_temp_c
  results$t_temp2[i]  <- b2 / se_temp2_c
  cat(sprintf(" se1=%.5f  t1=%.2f  se2=%.6f  t2=%.2f\n",
              se_temp_c, b1/se_temp_c, se_temp2_c, b2/se_temp2_c))
}

# ============================================================
# 7. PLOT t-STATISTICS vs SPATIAL CUTOFF
# ============================================================
pdf("figures/ExtendedDataFigs_Input/ConleyTstats.pdf", width = 9, height = 5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

# ---- Panel A: beta1 (temp linear) ----
plot(results$cutoff, abs(results$t_temp),
     type = "l", lwd = 2, col = "steelblue",
     xlab = "Spatial cutoff (km)",
     ylab = "|t-statistic|",
     main = expression(paste("Linear temperature term (", beta[1], ")")),
     ylim = c(0, max(abs(results$t_temp), na.rm = TRUE) * 1.15),
     las  = 1)
abline(h = 1.96,  lty = 2, col = "red",    lwd = 1.5)
abline(h = 1.645, lty = 3, col = "orange", lwd = 1.5)
legend("topright", legend = c("|t|", "5% (1.96)", "10% (1.645)"),
       lty = c(1,2,3), col = c("steelblue","red","orange"), lwd = 1.5, cex = 0.85)

# ---- Panel B: beta2 (temp²) ----
plot(results$cutoff, abs(results$t_temp2),
     type = "l", lwd = 2, col = "darkred",
     xlab = "Spatial cutoff (km)",
     ylab = "|t-statistic|",
     main = expression(paste("Quadratic temperature term (", beta[2], ")")),
     ylim = c(0, max(abs(results$t_temp2), na.rm = TRUE) * 1.15),
     las  = 1)
abline(h = 1.96,  lty = 2, col = "red",    lwd = 1.5)
abline(h = 1.645, lty = 3, col = "orange", lwd = 1.5)
legend("topright", legend = c("|t|", "5% (1.96)", "10% (1.645)"),
       lty = c(1,2,3), col = c("darkred","red","orange"), lwd = 1.5, cex = 0.85)

dev.off()
cat("\nSaved figures/ExtendedDataFigs_Input/ConleyTstats.pdf\n")

# ============================================================
# 8. PRINT SUMMARY TABLE
# ============================================================
cat("\n=== Conley SEs: t-statistics by spatial cutoff ===\n")
cat(sprintf("%-12s  %-10s  %-8s  %-10s  %-8s\n",
            "Cutoff (km)", "SE(beta1)", "t(beta1)", "SE(beta2)", "t(beta2)"))
cat(strrep("-", 58), "\n")
for (i in seq_len(nrow(results))) {
  cat(sprintf("%-12d  %-10.5f  %-8.3f  %-10.6f  %-8.3f\n",
              results$cutoff[i], results$se_temp[i], results$t_temp[i],
              results$se_temp2[i], results$t_temp2[i]))
}

# Clustered-SE baseline for comparison
se_clust <- sqrt(diag(vcov(fit_base)))
cat(sprintf("\nClustered (baseline): SE(b1)=%.5f t=%.3f | SE(b2)=%.6f t=%.3f\n",
            se_clust["temp"],  coef(fit_base)["temp"]  / se_clust["temp"],
            se_clust["temp2"], coef(fit_base)["temp2"] / se_clust["temp2"]))

cat("\nConleyAnalysis.R complete.\n")
