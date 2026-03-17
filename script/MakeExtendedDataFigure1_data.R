
# MakeExtendedDataFigure1_data.R
# R equivalent of MakeExtendedDataFigure1.do (data generation portion)
# Generates:
#   data/output/ExtendedDataFig1g.csv  -- country-level marginal effects
#   data/output/ExtendedDataFig1h.csv  -- interacted model marginal effects
#   data/output/ExtendedDataFig1i.csv  -- polynomial + spline predictions
#
# Run this before MakeExtendedDataFigure1.R

rm(list = ls())

library(data.table)
library(fixest)
library(splines)

setwd("/Users/yl3699/Downloads/BurkeHsiangMiguel2015_Replication")
dir.create("data/output", recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------
dta <- as.data.table(read.csv("data/input/GrowthClimateDataset.csv"))
dta[, temp  := UDel_temp_popweight]
dta[, temp2 := temp^2]
dta[, time  := year - 1985]
dta[, time2 := time^2]

# -----------------------------------------------------------------------
# FIG 1G: Country-level marginal effects vs global quadratic model
# Stata: areg growthWDI temp temp2 prec prec2 i.year _yi_* _y2_*, a(iso_id) cl(iso_id)
#        then for each country: reg growthWDI temp precip year year2
# -----------------------------------------------------------------------
cat("Computing Fig1g data (country-level marginals)...\n")

fit_global <- feols(
  growthWDI ~ temp + temp2 + UDel_precip_popweight + UDel_precip_popweight_2 |
    iso_id[time, time2] + year,
  cluster = ~iso_id, data = dta
)
b1 <- coef(fit_global)["temp"]
b2 <- coef(fit_global)["temp2"]

countries <- sort(unique(dta$iso[!is.na(dta$growthWDI) & !is.na(dta$temp)]))
marg_rows <- lapply(countries, function(ct) {
  d_ct <- dta[iso == ct & !is.na(growthWDI) & !is.na(temp)]
  if (nrow(d_ct) < 4) return(NULL)
  fit_ct <- tryCatch(
    lm(growthWDI ~ temp + UDel_precip_popweight + year + I(year^2), data = d_ct),
    error = function(e) NULL
  )
  if (is.null(fit_ct) || !("temp" %in% names(coef(fit_ct)))) return(NULL)
  data.frame(
    iso      = ct,
    b        = coef(fit_ct)["temp"],
    se       = sqrt(diag(vcov(fit_ct)))["temp"],
    meantemp = mean(d_ct$temp, na.rm = TRUE)
  )
})
marg_df <- rbindlist(Filter(Negate(is.null), marg_rows))
marg_df[, fit := b1 + 2 * b2 * meantemp]  # predicted marginal from global model

write.csv(marg_df, "data/output/ExtendedDataFig1g.csv", row.names = FALSE)
cat("Wrote data/output/ExtendedDataFig1g.csv  (", nrow(marg_df), "countries)\n")

# -----------------------------------------------------------------------
# FIG 1H: Interacted model marginal effects
# Model 1: growthWDI ~ temp + temp*mean_temp + prec + prec*mean_prec + FE
# Model 2: + temp*demeaned_mean_income
# Marginal effects at tavg = 1, 3, 5, ..., 29 via delta method
# -----------------------------------------------------------------------
cat("Computing Fig1h data (interacted model marginals)...\n")

dta[, tavg    := mean(temp, na.rm = TRUE), by = iso_id]
dta[, pavg    := mean(UDel_precip_popweight, na.rm = TRUE), by = iso_id]
dta[, tempavg := temp * tavg]
dta[, precavg := UDel_precip_popweight * pavg]

fit_int1 <- feols(
  growthWDI ~ temp + tempavg + UDel_precip_popweight + precavg |
    iso_id[time, time2] + year,
  cluster = ~iso_id, data = dta
)

# Income interaction (demeaned, divided by 1000 for legibility)
dta[, mGDP  := mean(gdpCAP_wdi, na.rm = TRUE), by = iso_id]
dta[, dmGDP := (mGDP - mean(mGDP, na.rm = TRUE)) / 1000]
dta[, tempInc := temp * dmGDP]
dta[, precInc := UDel_precip_popweight * dmGDP]

fit_int2 <- feols(
  growthWDI ~ temp + tempavg + tempInc + UDel_precip_popweight + precavg + precInc |
    iso_id[time, time2] + year,
  cluster = ~iso_id, data = dta
)

b1i <- coef(fit_int1); V1i <- vcov(fit_int1)
b2i <- coef(fit_int2); V2i <- vcov(fit_int2)

fig1h <- data.frame(matrix(NA_real_, nrow = 30, ncol = 8))
names(fig1h) <- c("n1", "n2", "marg1", "cilo1", "cihi1", "marg2", "cilo2", "cihi2")

for (i in seq(1, 29, by = 2)) {
  e1  <- b1i["temp"] + i * b1i["tempavg"]
  v1  <- V1i["temp","temp"] + i^2*V1i["tempavg","tempavg"] + 2*i*V1i["temp","tempavg"]
  e2  <- b2i["temp"] + i * b2i["tempavg"]
  v2  <- V2i["temp","temp"] + i^2*V2i["tempavg","tempavg"] + 2*i*V2i["temp","tempavg"]
  fig1h[i, ] <- c(i-0.2, i+0.2,
                  e1, e1-1.96*sqrt(v1), e1+1.96*sqrt(v1),
                  e2, e2-1.96*sqrt(v2), e2+1.96*sqrt(v2))
}

write.csv(fig1h, "data/output/ExtendedDataFig1h.csv", row.names = FALSE)
cat("Wrote data/output/ExtendedDataFig1h.csv\n")

# -----------------------------------------------------------------------
# FIG 1I: Polynomial (degree 2-7) and spline (3-7 knots) predictions
# Evaluated on grid X = 1/3, 2/3, ..., 99/3  (matches Stata: gen X = _n/3 if _n<100)
# -----------------------------------------------------------------------
cat("Computing Fig1i data (polynomials + splines)...\n")

# Create polynomial columns in data
for (j in 3:7) {
  col_nm <- paste0("temp", j)
  if (!col_nm %in% names(dta)) dta[, (col_nm) := temp^j]
}

X_grid <- (1:99) / 3  # X = 0.333, 0.667, ..., 33.0
poly_preds <- data.frame(X = X_grid)

# Polynomial models: Y2 from quadratic model, Y3 from cubic, ..., Y7 from degree-7
for (k in 2:7) {
  rhs_vars <- c("temp", paste0("temp", 2:k))
  fml <- as.formula(paste(
    "growthWDI ~", paste(rhs_vars, collapse = " + "),
    "+ UDel_precip_popweight + UDel_precip_popweight_2",
    "| iso_id[time, time2] + year"
  ))
  fit_poly <- feols(fml, cluster = ~iso_id, data = dta)
  b_p <- coef(fit_poly)

  # Cumulative polynomial evaluation at X grid (same as Stata's loc add accumulation)
  y_pred <- b_p["temp"] * X_grid
  for (j in 2:k) {
    nm <- paste0("temp", j)
    if (nm %in% names(b_p)) y_pred <- y_pred + b_p[nm] * X_grid^j
  }
  poly_preds[[paste0("Y", k)]] <- as.numeric(y_pred)
  cat("  Polynomial degree", k, "done\n")
}

# Spline models: natural cubic splines with 3-7 knots
# Knot placement matches Stata's mkspline, nknots(k) cubic percentile defaults:
knot_q <- list(
  "3" = c(0.10, 0.50, 0.90),
  "4" = c(0.05, 0.35, 0.65, 0.95),
  "5" = c(0.05, 0.275, 0.50, 0.725, 0.95),
  "6" = c(0.05, 0.23, 0.41, 0.59, 0.77, 0.95),
  "7" = c(0.05, 0.20, 0.35, 0.50, 0.65, 0.80, 0.95)
)

temp_est <- dta$temp[!is.na(dta$growthWDI) & !is.na(dta$temp)]

for (k in 3:7) {
  all_k   <- unname(quantile(temp_est, probs = knot_q[[as.character(k)]], na.rm = TRUE))
  inner_k <- all_k[-c(1L, length(all_k))]
  bound_k <- all_k[c(1L, length(all_k))]

  # Build spline basis on training data
  basis_mat <- ns(dta$temp, knots = inner_k, Boundary.knots = bound_k)
  nc     <- ncol(basis_mat)
  sp_nms <- paste0("sp", k, "_", seq_len(nc))

  # Attach basis columns to a temporary data.table for feols
  dta_sp <- dta[, .(growthWDI, UDel_precip_popweight, UDel_precip_popweight_2,
                     iso_id, year, time, time2)]
  for (j in seq_len(nc)) dta_sp[, (sp_nms[j]) := basis_mat[, j]]

  fml_sp <- as.formula(paste(
    "growthWDI ~", paste(sp_nms, collapse = " + "),
    "+ UDel_precip_popweight + UDel_precip_popweight_2",
    "| iso_id[time, time2] + year"
  ))
  fit_sp <- feols(fml_sp, cluster = ~iso_id, data = dta_sp)
  b_sp   <- coef(fit_sp)

  # Evaluate on X grid using same knots
  basis_grid <- ns(X_grid, knots = inner_k, Boundary.knots = bound_k)
  y_sp <- as.numeric(basis_grid %*% b_sp[sp_nms])
  poly_preds[[paste0("spline_est_", k)]] <- y_sp
  cat("  Spline knots =", k, "done\n")
}

write.csv(poly_preds, "data/output/ExtendedDataFig1i.csv", row.names = FALSE)
cat("Wrote data/output/ExtendedDataFig1i.csv\n")
cat("MakeExtendedDataFigure1_data.R complete.\n")
