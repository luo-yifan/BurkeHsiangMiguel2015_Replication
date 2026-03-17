##############################################################################
# ExtendPanel2022.R
# Extends the BHM (Burke, Hsiang & Miguel 2015) panel from 1960-2013 to
# 1960-2022 using:
#   - Climate data: Weighted Climate Dataset (Gortan et al. 2024) from Figshare
#       CRU TS 4.07 and ERA5 pop-weighted temperature + precipitation
#   - GDP data: World Bank WDI (NY.GDP.PCAP.KD.ZG), downloaded via wbstats
#
# Three regression panels:
#   (1) BHM original   1960-2013, UDel climate, BHM GDP
#   (2) Extended CRU   1960-2022, CRU TS climate, BHM+WB GDP
#   (3) Extended ERA5  1960-2022, ERA5 climate,   BHM+WB GDP
#   (4) Hybrid CRU     1960-2022, UDel 1960-2013 + CRU 2014-2022, BHM+WB GDP
##############################################################################

library(fixest)
library(data.table)

setwd("/Users/yl3699/Downloads/BurkeHsiangMiguel2015_Replication")

##############################################################################
# 1. LOAD BHM ORIGINAL DATA
##############################################################################
cat("=== Loading BHM original data ===\n")
bhm <- as.data.table(
  read.csv("data/input/GrowthClimateDataset.csv")
)

# Keep only rows with non-missing growthWDI and climate data
bhm_base <- bhm[!is.na(growthWDI) & !is.na(UDel_temp_popweight) &
                  !is.na(UDel_precip_popweight)]
cat("BHM observations:", nrow(bhm_base), "\n")
cat("BHM years:", min(bhm_base$year), "-", max(bhm_base$year), "\n")
cat("BHM countries:", length(unique(bhm_base$iso)), "\n")

##############################################################################
# 2. PARSE WCD MONTHLY CLIMATE FILES -> ANNUAL PANEL
##############################################################################
cat("\n=== Parsing WCD climate files ===\n")

parse_wcd <- function(filepath, var_name, agg_fn = mean) {
  # Read wide format: rows = XYYYYMM months, cols = ISO3 country codes
  dt <- fread(filepath)
  setnames(dt, "Date", "month_code")
  # Parse year from XYYYYMM: X190101 -> 1901-01
  dt[, year := as.integer(substr(month_code, 2, 5))]
  # Drop non-year columns
  country_cols <- setdiff(names(dt), c("month_code", "year"))
  # Remove synthetic region codes (Z-prefixed and ATA=Antarctica)
  country_cols <- country_cols[!grepl("^Z", country_cols) & country_cols != "ATA"]
  dt_sub <- dt[, c("year", country_cols), with = FALSE]
  # Aggregate to annual
  dt_ann <- dt_sub[, lapply(.SD, agg_fn, na.rm = TRUE), by = year]
  # Melt to long
  dt_long <- data.table::melt(dt_ann, id.vars = "year",
                              variable.name = "iso3c", value.name = var_name)
  # Replace NaN with NA (NaN can appear from mean/sum of all-NA months)
  vals <- dt_long[[var_name]]
  vals[is.nan(vals)] <- NA
  set(dt_long, j = var_name, value = vals)
  return(dt_long)
}

# Temperature: annual mean (average of 12 monthly values)
cat("  CRU temperature...\n")
cru_tmp <- parse_wcd("data/input/WCD/gadm0_cru_tmp_pop_2015_monthly.csv",
                     "temp_cru", mean)
cat("  ERA5 temperature...\n")
era_tmp <- parse_wcd("data/input/WCD/gadm0_era_tmp_pop_2015_monthly.csv",
                     "temp_era5", mean)

# Precipitation: annual sum (sum of 12 monthly values in mm/month -> mm/year)
cat("  CRU precipitation...\n")
cru_pre <- parse_wcd("data/input/WCD/gadm0_cru_pre_pop_2015_monthly.csv",
                     "prec_cru", sum)
cat("  ERA5 precipitation...\n")
era_pre <- parse_wcd("data/input/WCD/gadm0_era_pre_pop_2015_monthly.csv",
                     "prec_era5", sum)

# UDel precipitation (for hybrid: covers 1960-2017)
cat("  UDel temperature...\n")
dela_tmp <- parse_wcd("data/input/WCD/gadm0_dela_tmp_pop_2015_monthly.csv",
                      "temp_dela", mean)

# Merge climate sources
cat("  Merging climate sources...\n")
clim <- merge(cru_tmp, cru_pre, by = c("year", "iso3c"), all = TRUE)
clim <- merge(clim, era_tmp, by = c("year", "iso3c"), all = TRUE)
clim <- merge(clim, era_pre, by = c("year", "iso3c"), all = TRUE)
clim <- merge(clim, dela_tmp, by = c("year", "iso3c"), all = TRUE)
cat("Climate panel:", nrow(clim), "rows,", length(unique(clim$iso3c)),
    "countries,", min(clim$year), "-", max(clim$year), "\n")

# Filter to 1960-2022
clim <- clim[year >= 1960 & year <= 2022]

##############################################################################
# 3. ISO CODE MAPPING (BHM iso -> WCD iso3c)
##############################################################################
# Most codes match. Map legacy/alternative BHM codes to ISO3
iso_map <- c(
  ADO = "AND",  # Andorra
  CHI = "CHN",  # Channel Islands -> use NA (no match)
  HKG = "HKG",  # Hong Kong (in WCD)
  IMY = NA,     # Isle of Man
  KSV = "XKX",  # Kosovo
  MAC = "MAC",  # Macao (in WCD)
  ROM = "ROU",  # Romania (old code)
  TMP = "TLS",  # East Timor
  WBG = "PSE",  # West Bank & Gaza
  ZAR = "COD"   # DR Congo (old Zaire)
)

# Verify which of these mappings exist in WCD
wcd_isos <- unique(clim$iso3c)
cat("\nISO mapping check:\n")
for (old in names(iso_map)) {
  new <- iso_map[old]
  if (!is.na(new)) {
    cat(sprintf("  %s -> %s: %s\n", old, new,
                ifelse(new %in% wcd_isos, "found", "NOT FOUND")))
  }
}

##############################################################################
# 4. VALIDATE: COMPARE CRU vs UDel FOR OVERLAP PERIOD
##############################################################################
cat("\n=== Validating CRU vs UDel temperature (1960-2013 overlap) ===\n")

# Create a mapped iso column in clim
clim[, iso_bhm := as.character(iso3c)]
# Apply reverse mapping (WCD iso3c -> BHM iso)
rev_map <- c(AND = "ADO", ROU = "ROM", TLS = "TMP", PSE = "WBG", COD = "ZAR")
for (wcd_code in names(rev_map)) {
  clim[iso3c == wcd_code, iso_bhm := rev_map[wcd_code]]
}

# Merge BHM UDel with CRU for overlap check
val_bhm <- bhm[!is.na(UDel_temp_popweight),
               .(iso, year, temp_udel = UDel_temp_popweight,
                 prec_udel = UDel_precip_popweight)]
val_clim <- clim[year >= 1960 & year <= 2013,
                 .(iso_bhm, year, temp_cru, prec_cru, temp_era5, temp_dela)]
setnames(val_clim, "iso_bhm", "iso")

val <- merge(val_bhm, val_clim, by = c("iso", "year"))
val <- val[!is.na(temp_udel) & !is.na(temp_cru)]
cat("Validation obs (BHM + CRU overlap):", nrow(val), "\n")
cat(sprintf("UDel vs CRU temp: r=%.4f, RMSE=%.3f C\n",
            cor(val$temp_udel, val$temp_cru, use = "complete.obs"),
            sqrt(mean((val$temp_udel - val$temp_cru)^2, na.rm = TRUE))))
cat(sprintf("UDel vs ERA5 temp: r=%.4f, RMSE=%.3f C\n",
            cor(val$temp_udel, val$temp_era5, use = "complete.obs"),
            sqrt(mean((val$temp_udel - val$temp_era5)^2, na.rm = TRUE))))
cat(sprintf("UDel vs DELA temp: r=%.4f, RMSE=%.3f C\n",
            cor(val$temp_udel, val$temp_dela, use = "complete.obs"),
            sqrt(mean((val$temp_udel - val$temp_dela)^2, na.rm = TRUE))))
cat(sprintf("UDel vs CRU prec: r=%.4f, RMSE=%.1f mm/yr\n",
            cor(val$prec_udel, val$prec_cru, use = "complete.obs"),
            sqrt(mean((val$prec_udel - val$prec_cru)^2, na.rm = TRUE))))

##############################################################################
# 5. LOAD WORLD BANK GDP DATA
##############################################################################
cat("\n=== Loading World Bank GDP data ===\n")
wb <- fread("data/input/WCD/wb_gdp_1960_2022.csv")
setnames(wb, "NY.GDP.PCAP.KD.ZG", "growthWDI_wb")
wb[, growthWDI_wb := growthWDI_wb / 100]  # Convert % to fraction (match BHM)
wb <- wb[!is.na(growthWDI_wb) & date >= 1960 & date <= 2022]
setnames(wb, "date", "year")
cat("WB GDP rows:", nrow(wb), "\n")
cat("WB countries with 2022 data:",
    nrow(unique(wb[year == 2022, .(iso3c)])), "\n")

##############################################################################
# 6. BUILD EXTENDED PANELS
##############################################################################
cat("\n=== Building extended panels ===\n")

# ---- Helper: attach iso_bhm to WB data ----
wb[, iso_bhm := as.character(iso3c)]
for (wcd_code in names(rev_map)) {
  wb[iso3c == wcd_code, iso_bhm := rev_map[wcd_code]]
}

# ---- Panel A: BHM original (baseline) ----
# Use BHM data as-is; iso_id numeric FE already present
panel_bhm <- bhm_base[!is.na(growthWDI) & !is.na(UDel_temp_popweight) &
                        !is.na(UDel_precip_popweight),
                      .(iso, year, growthWDI,
                        temp = UDel_temp_popweight,
                        prec = UDel_precip_popweight,
                        temp2 = UDel_temp_popweight^2,
                        prec2 = UDel_precip_popweight^2,
                        iso_id, GDPpctile_WDIppp, time, time2)]
cat("Panel BHM original:", nrow(panel_bhm), "obs\n")

# ---- Function: build an extended panel with given climate source ----
build_extended <- function(tmp_col, prec_col, label) {
  # Step 1: BHM rows 1960-2013 (use BHM GDP, can optionally override climate)
  bhm_rows <- merge(
    bhm_base[!is.na(growthWDI), .(iso, year, growthWDI, iso_id, GDPpctile_WDIppp)],
    clim[year >= 1960 & year <= 2013,
         .(iso_bhm, year, temp = get(tmp_col), prec = get(prec_col))],
    by.x = c("iso", "year"), by.y = c("iso_bhm", "year"), all.x = TRUE
  )
  # Keep original BHM climate where WCD is missing
  bhm_rows_full <- merge(
    bhm_rows,
    bhm_base[, .(iso, year, temp_udel = UDel_temp_popweight,
                 prec_udel = UDel_precip_popweight)],
    by = c("iso", "year"), all.x = TRUE
  )
  bhm_rows_full[is.na(temp), temp := temp_udel]
  bhm_rows_full[is.na(prec), prec := prec_udel]
  bhm_rows_full[, c("temp_udel", "prec_udel") := NULL]
  bhm_rows_full <- bhm_rows_full[!is.na(temp) & !is.na(prec) & !is.na(growthWDI)]

  # Step 2: New rows 2014-2022 from WB GDP + WCD climate
  new_clim <- clim[year >= 2014 & year <= 2022,
                   .(iso = iso_bhm, year,
                     temp = get(tmp_col), prec = get(prec_col))]
  new_gdp  <- wb[year >= 2014 & year <= 2022,
                 .(iso = iso_bhm, year, growthWDI = growthWDI_wb)]
  new_rows <- merge(new_gdp, new_clim, by = c("iso", "year"))
  new_rows <- new_rows[!is.na(growthWDI) & !is.na(temp) & !is.na(prec)]
  # Restrict to countries in BHM panel
  bhm_ctrs <- unique(bhm_base$iso)
  new_rows <- new_rows[iso %in% bhm_ctrs]
  # Add iso_id from BHM mapping
  iso_id_map <- unique(bhm_base[, .(iso, iso_id)])
  new_rows   <- merge(new_rows, iso_id_map, by = "iso")
  # Add GDPpctile_WDIppp (use last BHM value for each country)
  gdp_pct <- bhm_base[!is.na(GDPpctile_WDIppp),
                      .(GDPpctile_WDIppp = GDPpctile_WDIppp[which.max(year)]),
                      by = iso]
  new_rows <- merge(new_rows, gdp_pct, by = "iso", all.x = TRUE)

  # Step 3: Stack and compute time/time2 (years since country's first observation)
  combined <- rbindlist(list(bhm_rows_full, new_rows), use.names = TRUE, fill = TRUE)
  combined <- combined[order(iso, year)]
  combined[, time := year - min(year), by = iso]
  combined[, time2 := time^2]
  combined[, temp2 := temp^2]
  combined[, prec2 := prec^2]
  combined <- combined[!is.na(growthWDI) & !is.na(temp) & !is.na(prec)]
  cat(sprintf("Panel %s: %d obs, %d countries, %d-%d\n",
              label, nrow(combined), length(unique(combined$iso)),
              min(combined$year), max(combined$year)))
  return(combined)
}

# Build panels using CRU and ERA5
panel_cru  <- build_extended("temp_cru",  "prec_cru",  "Extended-CRU")
panel_era5 <- build_extended("temp_era5", "prec_era5", "Extended-ERA5")

# Hybrid: use UDel climate 1960-2013 (BHM original) + CRU 2014-2022
build_hybrid <- function() {
  # BHM 1960-2013 rows with original UDel climate
  bhm_rows <- bhm_base[!is.na(growthWDI) & !is.na(UDel_temp_popweight) &
                          !is.na(UDel_precip_popweight),
                        .(iso, year, growthWDI, iso_id, GDPpctile_WDIppp,
                          temp = UDel_temp_popweight, prec = UDel_precip_popweight)]
  # 2014-2022 with CRU climate + WB GDP
  new_clim <- clim[year >= 2014 & year <= 2022,
                   .(iso = iso_bhm, year, temp = temp_cru, prec = prec_cru)]
  new_gdp  <- wb[year >= 2014 & year <= 2022,
                 .(iso = iso_bhm, year, growthWDI = growthWDI_wb)]
  new_rows <- merge(new_gdp, new_clim, by = c("iso", "year"))
  new_rows <- new_rows[!is.na(growthWDI) & !is.na(temp) & !is.na(prec)]
  new_rows <- new_rows[iso %in% unique(bhm_base$iso)]
  iso_id_map <- unique(bhm_base[, .(iso, iso_id)])
  new_rows   <- merge(new_rows, iso_id_map, by = "iso")
  gdp_pct <- bhm_base[!is.na(GDPpctile_WDIppp),
                      .(GDPpctile_WDIppp = GDPpctile_WDIppp[which.max(year)]),
                      by = iso]
  new_rows <- merge(new_rows, gdp_pct, by = "iso", all.x = TRUE)
  combined <- rbindlist(list(bhm_rows, new_rows), use.names = TRUE, fill = TRUE)
  combined <- combined[order(iso, year)]
  combined[, time := year - min(year), by = iso]
  combined[, time2 := time^2]
  combined[, temp2 := temp^2]
  combined[, prec2 := prec^2]
  combined <- combined[!is.na(growthWDI) & !is.na(temp) & !is.na(prec)]
  cat(sprintf("Panel Hybrid (UDel+CRU): %d obs, %d countries, %d-%d\n",
              nrow(combined), length(unique(combined$iso)),
              min(combined$year), max(combined$year)))
  return(combined)
}

panel_hybrid <- build_hybrid()

##############################################################################
# 7. RUN BHM BASELINE REGRESSION ON ALL FOUR PANELS
##############################################################################
cat("\n=== Running BHM baseline regressions ===\n")

run_bhm <- function(panel, label) {
  fit <- feols(
    growthWDI ~ temp + temp2 + prec + prec2 | iso_id[time, time2] + year,
    cluster = ~iso_id, data = panel,
    warn = FALSE
  )
  b1 <- coef(fit)["temp"]
  b2 <- coef(fit)["temp2"]
  se1 <- se(fit)["temp"]
  se2 <- se(fit)["temp2"]
  opt_t <- -b1 / (2 * b2)
  me_20 <- b1 + 2 * b2 * 20
  n_obs <- nobs(fit)
  n_cty <- length(unique(panel$iso))
  yr_range <- paste(min(panel$year), max(panel$year), sep = "-")
  list(
    label  = label, yr_range = yr_range, n_obs = n_obs, n_cty = n_cty,
    b1 = b1, se1 = se1, t1 = b1 / se1,
    b2 = b2, se2 = se2, t2 = b2 / se2,
    opt_t = opt_t, me_20 = me_20
  )
}

res1 <- run_bhm(panel_bhm,   "BHM original (UDel, 1960-2013)")
res2 <- run_bhm(panel_cru,   "Extended CRU (1960-2022)")
res3 <- run_bhm(panel_era5,  "Extended ERA5 (1960-2022)")
res4 <- run_bhm(panel_hybrid,"Hybrid UDel+CRU (1960-2022)")

##############################################################################
# 8. PRINT COMPARISON TABLE
##############################################################################
cat("\n")
cat("==========================================================================\n")
cat("BHM PANEL EXTENSION TO 2022: COMPARISON TABLE\n")
cat("==========================================================================\n")
cat(sprintf("%-35s %6s %5s %7s %7s %7s %7s %6s %7s\n",
            "Specification", "N", "Cty", "b1(temp)", "se1",
            "b2(temp2)", "se2", "Opt.T", "ME@20"))
cat(paste(rep("-", 90), collapse = ""), "\n")
for (r in list(res1, res2, res3, res4)) {
  cat(sprintf("%-35s %6d %5d  %7.5f %7.5f  %8.6f %8.6f %6.1f %8.5f\n",
              r$label, r$n_obs, r$n_cty,
              r$b1, r$se1, r$b2, r$se2,
              r$opt_t, r$me_20))
}
cat(paste(rep("-", 90), collapse = ""), "\n")
cat("Notes: b1=temp coeff, b2=temp^2 coeff, Opt.T=-b1/(2*b2) optimal temp (C),\n")
cat("       ME@20=marginal effect at 20C. All panels use country quad time trends + year FE,\n")
cat("       cluster SE by country.\n\n")

# Significance stars
sig_stars <- function(t) {
  if (abs(t) > 3.29) return("***")
  if (abs(t) > 2.58) return("**")
  if (abs(t) > 1.96) return("*")
  return("")
}

cat("Significance:\n")
for (r in list(res1, res2, res3, res4)) {
  cat(sprintf("  %-35s: temp %s%s, temp2 %s%s\n",
              r$label,
              formatC(r$b1, format = "f", digits = 5),
              sig_stars(r$t1),
              formatC(r$b2, format = "f", digits = 6),
              sig_stars(r$t2)))
}

##############################################################################
# 9. INCOME HETEROGENEITY ON EXTENDED PANEL
##############################################################################
cat("\n=== Income heterogeneity (rich vs poor) on extended panels ===\n")

run_heterog <- function(panel, label) {
  sub <- panel[!is.na(GDPpctile_WDIppp)]
  sub[, poor := as.integer(GDPpctile_WDIppp < 50)]
  fit <- feols(
    growthWDI ~ temp + temp2 + prec + prec2 +
      I(temp * poor) + I(temp2 * poor) + I(prec * poor) + I(prec2 * poor) |
      iso_id[time, time2] + year,
    cluster = ~iso_id, data = sub, warn = FALSE
  )
  cf <- coef(fit)
  b1r <- cf["temp"];  b2r <- cf["temp2"]
  b1p <- cf["temp"] + cf["I(temp * poor)"]
  b2p <- cf["temp2"] + cf["I(temp2 * poor)"]
  cat(sprintf("  %-35s: Rich opt.T=%5.1f C, Poor opt.T=%5.1f C\n",
              label, -b1r/(2*b2r), -b1p/(2*b2p)))
}

run_heterog(panel_bhm,   "BHM original")
run_heterog(panel_cru,   "Extended CRU")
run_heterog(panel_era5,  "Extended ERA5")
run_heterog(panel_hybrid,"Hybrid UDel+CRU")

##############################################################################
# 10. SAVE EXTENDED PANELS AND RESULTS
##############################################################################
cat("\n=== Saving outputs ===\n")

# Save the extended CRU panel (most natural extension)
panel_cru_save <- panel_cru[, .(iso, year, growthWDI, temp, temp2, prec, prec2,
                                 iso_id, GDPpctile_WDIppp, time, time2)]
write.csv(panel_cru_save, "data/output/extended_panel_cru_1960_2022.csv",
          row.names = FALSE)
cat("Saved extended CRU panel:", nrow(panel_cru_save), "rows\n")

panel_era5_save <- panel_era5[, .(iso, year, growthWDI, temp, temp2, prec, prec2,
                                   iso_id, GDPpctile_WDIppp, time, time2)]
write.csv(panel_era5_save, "data/output/extended_panel_era5_1960_2022.csv",
          row.names = FALSE)
cat("Saved extended ERA5 panel:", nrow(panel_era5_save), "rows\n")

# Summary results to CSV
results_df <- rbind(
  data.frame(spec = res1$label, n_obs = res1$n_obs, n_cty = res1$n_cty,
             b_temp = res1$b1, se_temp = res1$se1, b_temp2 = res1$b2,
             se_temp2 = res1$se2, opt_T = res1$opt_t, ME_20 = res1$me_20),
  data.frame(spec = res2$label, n_obs = res2$n_obs, n_cty = res2$n_cty,
             b_temp = res2$b1, se_temp = res2$se1, b_temp2 = res2$b2,
             se_temp2 = res2$se2, opt_T = res2$opt_t, ME_20 = res2$me_20),
  data.frame(spec = res3$label, n_obs = res3$n_obs, n_cty = res3$n_cty,
             b_temp = res3$b1, se_temp = res3$se1, b_temp2 = res3$b2,
             se_temp2 = res3$se2, opt_T = res3$opt_t, ME_20 = res3$me_20),
  data.frame(spec = res4$label, n_obs = res4$n_obs, n_cty = res4$n_cty,
             b_temp = res4$b1, se_temp = res4$se1, b_temp2 = res4$b2,
             se_temp2 = res4$se2, opt_T = res4$opt_t, ME_20 = res4$me_20)
)
write.csv(results_df, "data/output/panel_extension_results.csv", row.names = FALSE)
cat("Saved regression comparison to data/output/panel_extension_results.csv\n")

##############################################################################
# 11. PLOT: RESPONSE CURVES FOR ALL FOUR SPECS
##############################################################################
cat("\n=== Generating comparison figure ===\n")

pdf("figures/ExtendedDataFigs_Input/PanelExtension2022.pdf",
    width = 10, height = 6)

# Refit to get full vcov for CI
make_response <- function(panel) {
  fit <- feols(
    growthWDI ~ temp + temp2 + prec + prec2 | iso_id[time, time2] + year,
    cluster = ~iso_id, data = panel, warn = FALSE
  )
  temp_grid <- seq(-5, 35, by = 0.5)
  b1 <- coef(fit)["temp"]; b2 <- coef(fit)["temp2"]
  V  <- vcov(fit)[c("temp", "temp2"), c("temp", "temp2")]
  y  <- b1 * temp_grid + b2 * temp_grid^2
  # SE via delta method: se = sqrt(t^2*V11 + t^4*V22 + 2*t^3*V12)
  se_y <- sapply(temp_grid, function(t) {
    g <- c(t, t^2)
    sqrt(max(0, t(g) %*% V %*% g))
  })
  data.frame(temp = temp_grid, effect = y, se = se_y)
}

resp1 <- make_response(panel_bhm)
resp2 <- make_response(panel_cru)
resp3 <- make_response(panel_era5)
resp4 <- make_response(panel_hybrid)

# Normalize so effect = 0 at optimal temperature of BHM original
t_star <- -res1$b1 / (2 * res1$b2)
norm_at <- function(df) {
  y0 <- with(df, approx(temp, effect, xout = t_star)$y)
  df$effect <- df$effect - y0
  df
}
resp1 <- norm_at(resp1)
resp2 <- norm_at(resp2)
resp3 <- norm_at(resp3)
resp4 <- norm_at(resp4)

ylim_range <- range(c(resp1$effect - 1.645*resp1$se,
                       resp1$effect + 1.645*resp1$se,
                       resp2$effect, resp3$effect, resp4$effect), na.rm = TRUE)
ylim_range <- c(max(ylim_range[1], -0.20), min(ylim_range[2], 0.10))

cols <- c("black", "#0072B2", "#D55E00", "#009E73")
ltys  <- c(1, 2, 3, 4)

par(mar = c(4.5, 4.5, 3, 8))
plot(resp1$temp, resp1$effect, type = "n",
     xlim = c(-5, 35), ylim = ylim_range,
     xlab = "Annual mean temperature (C)",
     ylab = "Change in GDP per capita growth (pp)",
     main = "BHM Temperature-Growth Response: Original vs Extended to 2022",
     las = 1, cex.axis = 0.9, cex.lab = 1.0)

# Shaded CI for BHM original
polygon(c(resp1$temp, rev(resp1$temp)),
        c(resp1$effect - 1.645*resp1$se, rev(resp1$effect + 1.645*resp1$se)),
        col = adjustcolor("grey70", alpha.f = 0.5), border = NA)

abline(h = 0, lty = 2, col = "grey50")
abline(v = t_star, lty = 3, col = "grey60")

lines(resp1$temp, resp1$effect, col = cols[1], lwd = 2, lty = ltys[1])
lines(resp2$temp, resp2$effect, col = cols[2], lwd = 2, lty = ltys[2])
lines(resp3$temp, resp3$effect, col = cols[3], lwd = 2, lty = ltys[3])
lines(resp4$temp, resp4$effect, col = cols[4], lwd = 2, lty = ltys[4])

legend("topright", inset = c(-0.32, 0),
       legend = c(
         sprintf("BHM original\n(UDel, 1960-2013, N=%d)", res1$n_obs),
         sprintf("Ext. CRU TS\n(1960-2022, N=%d)", res2$n_obs),
         sprintf("Ext. ERA5\n(1960-2022, N=%d)", res3$n_obs),
         sprintf("Hybrid UDel+CRU\n(1960-2022, N=%d)", res4$n_obs)
       ),
       col = cols, lty = ltys, lwd = 2,
       cex = 0.78, bty = "n", xpd = TRUE)

mtext("Grey band = 90% CI for BHM original. Curves normalized to 0 at BHM optimal T.",
      side = 1, line = 3.5, cex = 0.75)

dev.off()
cat("Saved figure: figures/ExtendedDataFigs_Input/PanelExtension2022.pdf\n")

cat("\n=== Done ===\n")
