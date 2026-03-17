
# MakeExtendedDataFigure3_data.R
# R equivalent of MakeExtendedDataFigure3.do (data generation portion)
# Runs 6 regressions comparing DJO and BHM specifications.
#
# Output:
#   data/output/DJO/Coefficients_DJOcompare.csv  (mod, b1, b2 for 6 models)
#
# Requires: install.packages(c("haven", "fixest", "data.table"))
# Run this before MakeExtendedDataFigure3.R

rm(list = ls())

library(haven)
library(data.table)
library(fixest)

setwd("/Users/yl3699/Downloads/BurkeHsiangMiguel2015_Replication")
dir.create("data/output/DJO", recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------
# 1. READ AND PREPARE DJO DATA
# -----------------------------------------------------------------------
raw <- as.data.table(zap_labels(read_dta("data/input/DellJonesOlken/climate_panel.dta")))

# Restrict to <= 2003 (as in DJO)
djo <- raw[year <= 2003]

# Numeric group IDs for fixed effects and clustering
djo[, cc_num     := .GRP, by = fips60_06]
djo[, parent_num := .GRP, by = parent]

# Compute log GDP (WDI) and growth rate
djo[, lngdpwdi := log(gdpLCU)]
setorderv(djo, c("fips60_06", "year"))
djo[, lngdpwdi_lag := data.table::shift(lngdpwdi, 1), by = fips60_06]
djo[, g := lngdpwdi - lngdpwdi_lag]   # log growth, same units as growthWDI

# Drop countries with < 20 years of GDP data
djo[, nobs := sum(!is.na(g)), by = fips60_06]
djo <- djo[nobs >= 20]

# Create poor = bottom half of initial log GDP per capita (one value per country)
initgdp <- djo[!is.na(lnrgdpl_t0), .SD[1], by = fips60_06][, .(fips60_06, lnrgdpl_t0)]
med_lgdp <- median(log(initgdp$lnrgdpl_t0), na.rm = TRUE)
initgdp[, poor_djo := as.integer(log(lnrgdpl_t0) <= med_lgdp)]
djo <- merge(djo, initgdp[, .(fips60_06, poor_djo)], by = "fips60_06", all.x = TRUE)

# Region variable and region-year FE ID
region_cols <- c("_MENA", "_SSAF", "_LAC", "_WEOFF", "_EECA", "_SEAS")
djo[, region := "other"]
for (rc in intersect(region_cols, names(djo))) {
  djo[get(rc) == 1, region := rc]
}
djo[, regionyear := paste0(region, year)]
djo[, rynum := .GRP, by = regionyear]

# Climate variables
setnames(djo, c("wtem", "wpre"), c("temp", "prec"))
djo[, prec  := prec / 10]    # rescale to meters (as in Stata: replace prec = prec/10)
djo[, temp2 := temp^2]
djo[, prec2 := prec^2]

# Country time trends (anchored at 1949 like DJO)
djo[, time  := year - 1949]
djo[, time2 := time^2]

# Rename country_code to iso
setnames(djo, "country_code", "iso_djo")
djo[, iso := as.character(iso_djo)]

# -----------------------------------------------------------------------
# 2. DJO REGRESSIONS (models 1-3): on DJO data
# -----------------------------------------------------------------------
cat("Running DJO regressions (models 1-3)...\n")

# mod 1: DJO original spec
#   cgmreg g temp temp2 RY* i.cc_num, cluster(parent_num rynum)
#   region×year FEs (rynum) + country FEs (cc_num), no trends, no precip
m1 <- feols(g ~ temp + temp2 | cc_num + rynum,
            cluster = ~parent_num + rynum, data = djo)

# mod 2: DJO + quadratic country time trends + precip
#   cgmreg g temp temp2 prec prec2 RY* i.cc_num _yi_* _y2_*, cluster(parent_num rynum)
m2 <- feols(g ~ temp + temp2 + prec + prec2 | cc_num[time, time2] + rynum,
            cluster = ~parent_num + rynum, data = djo)

# mod 3: year FE + country trends + precip (closest to BHM baseline)
#   cgmreg g temp temp2 prec prec2 i.year i.cc_num _yi_* _y2_*, cluster(parent_num rynum)
m3 <- feols(g ~ temp + temp2 + prec + prec2 | cc_num[time, time2] + year,
            cluster = ~parent_num + rynum, data = djo)

# -----------------------------------------------------------------------
# 3. BHM REGRESSIONS (models 4-6): on BHM data
# -----------------------------------------------------------------------
cat("Running BHM regressions (models 4-6)...\n")

bhm <- as.data.table(read.csv("data/input/GrowthClimateDataset.csv"))
bhm[, temp  := UDel_temp_popweight]
bhm[, temp2 := temp^2]
bhm[, prec  := UDel_precip_popweight / 1000]
bhm[, prec2 := prec^2]
bhm[, time  := year - 1985]
bhm[, time2 := time^2]
bhm[, poor  := as.integer(GDPpctile_WDIppp < 50)]

# mod 4: our baseline on full WDI sample with >= 20 observations
#   cgmreg growthWDI temp temp2 prec prec2 i.year _yi_* _y2_* i.iso_id if wdinomiss>=20, cluster(iso_id)
m4 <- feols(growthWDI ~ temp + temp2 + prec + prec2 | iso_id[time, time2] + year,
            cluster = ~iso_id, data = bhm[wdinomiss >= 20])

# DJO estimation sample: countries present in both DJO data and BHM data
djo_isos <- unique(djo$iso)
bhm_djo  <- bhm[iso %in% djo_isos]

# mod 5: our baseline restricted to DJO sample
#   cgmreg growthWDI temp temp2 prec prec2 i.year _yi_* _y2_* i.iso_id if DJOmerge==3, cluster(iso_id)
m5 <- feols(growthWDI ~ temp + temp2 + prec + prec2 | iso_id[time, time2] + year,
            cluster = ~iso_id, data = bhm_djo)

# mod 6: DJO FE structure (continent×year + poor×year), DJO sample, our growth + our temp
#   cgmreg growthWDI temp temp2 prec prec2 _cy_* _py_* i.iso_id if DJOmerge==3, cluster(iso_id)
#   _cy_* = continent×year FE  =>  year^continent
#   _py_* = poor×year FE       =>  year^poor
m6 <- feols(growthWDI ~ temp + temp2 + prec + prec2 | iso_id + year^continent + year^poor,
            cluster = ~iso_id, data = bhm_djo[!is.na(poor)])

# -----------------------------------------------------------------------
# 4. COLLECT COEFFICIENTS AND SAVE CSV
# -----------------------------------------------------------------------
extract_b <- function(fit) {
  b <- coef(fit)
  c(b1 = unname(b["temp"]), b2 = unname(b["temp2"]))
}

models <- list(m1, m2, m3, m4, m5, m6)
coef_df <- data.frame(
  mod = 1:6,
  b1  = sapply(models, function(f) extract_b(f)["b1"]),
  b2  = sapply(models, function(f) extract_b(f)["b2"])
)

write.csv(coef_df, "data/output/DJO/Coefficients_DJOcompare.csv", row.names = FALSE)
cat("Wrote data/output/DJO/Coefficients_DJOcompare.csv\n")
print(coef_df)
cat("MakeExtendedDataFigure3_data.R complete.\n")
