
# MakeExtendedDataFigure2_data.R
# R equivalent of MakeExtendedDataFigure2.do (data generation portion)
# Generates:
#   data/output/Effect_Marginals_RichPoor.csv  -- rich/poor marginal effects at T=1..30
#
# Run this before MakeExtendedDataFigure2.R

rm(list = ls())

library(data.table)
library(fixest)

setwd("/Users/yl3699/Downloads/BurkeHsiangMiguel2015_Replication")
dir.create("data/output", recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------
# Load and prepare data
# -----------------------------------------------------------------------
dta <- as.data.table(read.csv("data/input/GrowthClimateDataset.csv"))

dta[, temp  := UDel_temp_popweight]
dta[, temp2 := temp^2]
dta[, prec  := UDel_precip_popweight]
dta[, prec2 := prec^2]
dta[, time  := year - 1985]
dta[, time2 := time^2]

# Rich/poor indicator (bottom 50th percentile of GDP = poor)
dta[, poorWDIppp := as.integer(GDPpctile_WDIppp < 50)]
dta[GDPpctile_WDIppp == 0 | is.na(GDPpctile_WDIppp), poorWDIppp := NA]

# Interaction terms
dta[, temppoor  := temp  * (poorWDIppp == 1)]
dta[, temp2poor := temp2 * (poorWDIppp == 1)]
dta[, precpoor  := prec  * (poorWDIppp == 1)]
dta[, prec2poor := prec2 * (poorWDIppp == 1)]

# -----------------------------------------------------------------------
# For each outcome, run rich/poor interacted regression, then compute
# marginal effects at T = 1..30 via delta method
# -----------------------------------------------------------------------
# Column: c' * V * c  where c is gradient vector
quad_var <- function(V, nms, c_vec) {
  # c_vec is a named numeric vector; compute c' V c
  idx <- match(names(c_vec), nms)
  ok  <- !is.na(idx)
  c_sub <- c_vec[ok]
  V_sub <- V[idx[ok], idx[ok]]
  as.numeric(t(c_sub) %*% V_sub %*% c_sub)
}

vars <- c("growthWDI", "AgrGDPgrowthCap", "NonAgrGDPgrowthCap")

out <- data.frame(x = 1:30)

for (v in vars) {
  cat("Running regression for", v, "...\n")
  fml <- as.formula(paste(
    v, "~ temp + temp2 + temppoor + temp2poor +",
    "prec + prec2 + precpoor + prec2poor |",
    "iso_id[time, time2] + year"
  ))
  fit <- feols(fml, cluster = ~iso_id, data = dta[!is.na(poorWDIppp)])
  b   <- coef(fit)
  V   <- vcov(fit)
  nms <- names(b)

  br  <- numeric(30); bep <- numeric(30); bc  <- numeric(30)
  ser <- numeric(30); sep <- numeric(30); sec <- numeric(30)
  pr  <- numeric(30); pp  <- numeric(30); pc  <- numeric(30)

  for (t in 1:30) {
    # Rich marginal: b_temp + 2*t*b_temp2
    # Gradient: c = (1, 2t, 0, 0) for (temp, temp2, temppoor, temp2poor)
    e_r   <- b["temp"] + 2 * t * b["temp2"]
    c_r   <- c(temp = 1, temp2 = 2 * t)
    v_r   <- quad_var(V, nms, c_r)
    se_r  <- sqrt(max(v_r, 0))

    # Poor marginal: b_temp + b_temppoor + 2*t*(b_temp2 + b_temp2poor)
    # Gradient: (1, 2t, 1, 2t)
    e_p   <- b["temp"] + b["temppoor"] + 2 * t * (b["temp2"] + b["temp2poor"])
    c_p   <- c(temp = 1, temp2 = 2 * t, temppoor = 1, temp2poor = 2 * t)
    v_p   <- quad_var(V, nms, c_p)
    se_p  <- sqrt(max(v_p, 0))

    # Contrast (rich - poor) = -(b_temppoor + 2*t*b_temp2poor)
    # Gradient: (0, 0, -1, -2t)
    e_c   <- -(b["temppoor"] + 2 * t * b["temp2poor"])
    c_c   <- c(temppoor = -1, temp2poor = -2 * t)
    v_c   <- quad_var(V, nms, c_c)
    se_c  <- sqrt(max(v_c, 0))

    br[t]  <- e_r;   ser[t] <- se_r;  pr[t] <- 2 * pnorm(-abs(e_r / se_r))
    bep[t] <- e_p;   sep[t] <- se_p;  pp[t] <- 2 * pnorm(-abs(e_p / se_p))
    bc[t]  <- e_c;   sec[t] <- se_c;  pc[t] <- 2 * pnorm(-abs(e_c / se_c))
  }

  out[[paste0(v, "_br")]]  <- br
  out[[paste0(v, "_bp")]]  <- bep
  out[[paste0(v, "_bc")]]  <- bc
  out[[paste0(v, "_ser")]] <- ser
  out[[paste0(v, "_sep")]] <- sep
  out[[paste0(v, "_sec")]] <- sec
  out[[paste0(v, "_pr")]]  <- pr
  out[[paste0(v, "_pp")]]  <- pp
  out[[paste0(v, "_pc")]]  <- pc

  cat("  Done:", v, "\n")
}

write.csv(out, "data/output/Effect_Marginals_RichPoor.csv", row.names = FALSE)
cat("Wrote data/output/Effect_Marginals_RichPoor.csv (", nrow(out), "rows x", ncol(out), "cols)\n")
cat("MakeExtendedDataFigure2_data.R complete.\n")
