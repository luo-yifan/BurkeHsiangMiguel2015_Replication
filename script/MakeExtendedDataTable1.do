

* Script to make Extended Data Table 1

clear all
set mem 1G
set matsize 10000
set maxvar 10000


postfile robust mod b1 b2 using data/output/Coefficients_robustness, replace
use data/input/GrowthClimateDataset, clear
xtset iso_id year
gen temp = UDel_temp_popweight
gen temp2 = temp*temp
gen precip = UDel_precip_popweight/1000 //rescaling so precip coeffs woudl be legible
gen precip2 = precip*precip
cap est drop reg*
loc j = 1
reg growthWDI temp temp2 precip precip2 i.year _yi_* _y2_* i.iso_id , cluster(iso_id)  //full sample
est sto reg`j'
estadd sca opt = _b[temp]/-2/_b[temp2]
post robust (`j') (_b[temp]) (_b[temp2])
loc j = `j' + 1
reg growthWDI temp temp2 precip precip2 i.year _yi_* _y2_* i.iso_id if wdinomiss>=20, cluster(iso_id)  //drop countries with <20 obs
est sto reg`j'
estadd sca opt = _b[temp]/-2/_b[temp2]
post robust (`j') (_b[temp]) (_b[temp2])
loc j = `j' + 1
reg growthWDI temp temp2 precip precip2 i.year _yi_* _y2_* i.iso_id if  oil==0, cluster(iso_id)  //Drop oil countries
est sto reg`j'
estadd sca opt = _b[temp]/-2/_b[temp2]
post robust (`j') (_b[temp]) (_b[temp2])
loc j = `j' + 1
reg growthWDI temp temp2 precip precip2 i.year _yi_* _y2_* i.iso_id if  iso~="USA" & iso~="CHN", cluster(iso_id)  //drop US and China
est sto reg`j'
estadd sca opt = _b[temp]/-2/_b[temp2]
post robust (`j') (_b[temp]) (_b[temp2])
loc j = `j' + 1
xi i.year*i.continent, pref(_cy_)  //make continent-yr FE
reg growthWDI temp temp2 precip precip2 _cy_* _yi_* _y2_* i.iso_id , cluster(iso_id)  //add continent year FE
est sto reg`j'
estadd sca opt = _b[temp]/-2/_b[temp2]
post robust (`j') (_b[temp]) (_b[temp2])
loc j = `j' + 1
reg growthWDI temp temp2 precip precip2 _cy_* i.iso_id , cluster(iso_id)  //continent-year FE and no time trends
est sto reg`j'
estadd sca opt = _b[temp]/-2/_b[temp2]
post robust (`j') (_b[temp]) (_b[temp2])
loc j = `j' + 1
reg growthWDI temp temp2 precip precip2 _yi_* _y2_* i.iso_id, cluster(iso_id)  // no year FE
est sto reg`j'
estadd sca opt = _b[temp]/-2/_b[temp2]
post robust (`j') (_b[temp]) (_b[temp2])
loc j = `j' + 1
reg growthWDI temp temp2 precip precip2 i.year _yi_* i.iso_id, cluster(iso_id)  //linear time trends
est sto reg`j'
estadd sca opt = _b[temp]/-2/_b[temp2]
post robust (`j') (_b[temp]) (_b[temp2])
loc j = `j' + 1
reg growthWDI L.growthWDI temp temp2 precip precip2 _yi_* _y2_* i.iso_id , cluster(iso_id)  // one lag of growth
est sto reg`j'
estadd sca opt = _b[temp]/-2/_b[temp2]
post robust (`j') (_b[temp]) (_b[temp2])
loc j = `j' + 1
reg growthWDI L(1/3).growthWDI temp temp2 precip precip2 _yi_* _y2_* i.iso_id , cluster(iso_id) //3 lags of growth
est sto reg`j'
estadd sca opt = _b[temp]/-2/_b[temp2]
post robust (`j') (_b[temp]) (_b[temp2])
loc j = `j' + 1
reg rgdpCAPgr temp temp2 precip precip2 _yi_* _y2_* i.iso_id , cluster(iso_id)  //Penn World Tables data
est sto reg`j'
estadd sca opt = _b[temp]/-2/_b[temp2]
post robust (`j') (_b[temp]) (_b[temp2])

lab var temp "Temp."
lab var temp2 "Temp. sq."
lab var precip "Precip."
lab var precip2 "Precip. sq."
esttab reg* using figures/ExtendedDataFigs_Input/Table1.tex, drop(_yi_* _y2_* *.year *.iso_id _cy_* L*)  ///
	replace b(%10.4f %10.4f %10.4f %10.4f %10.4f) se l sfmt(%10.0f %10.3f %10.2f) ///
	scalars("N Observations" "r2 R squared" "opt Optimum") ///
	star(* 0.10 ** 0.05 *** 0.01) compress nonotes width(\hsize) ///
	mtitles("Base" "$>$20yrs" "No oil" "No US/China" "ContYr FE" "ContYr + noTrend" "No YrFE" "LinearTime" "LDV 1lag" "LDV 3lags" "PWT") 
postclose robust
use data/output/Coefficients_robustness, clear
outsheet using data/output/Coefficients_robustness.csv, comma replace
