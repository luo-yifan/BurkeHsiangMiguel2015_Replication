
* RESULTS TO CONSTRUCT RESPONSE FUNCTIONS IN FIGURE 2
//		Note: users will need to install parmest command (e.g. "ssc install parmest")

clear all
set mem 1G
set matsize 10000
set maxvar 10000
set more off

// ---------- ESTIMATE GLOBAL RESPONSE WITH BASELINE REGRESSION SPECIFICATION   
// ---------- THEN WRITE OUT FOR CONSTUCTION OF FIGURE 2, panel A

use data/input/GrowthClimateDataset, clear
gen temp = UDel_temp_popweight
reg growthWDI c.temp##c.temp UDel_precip_popweight UDel_precip_popweight_2 i.year _yi_* _y2_* i.iso_id, cluster(iso_id)
	mat b = e(b)
	mat b = b[1,1..2] //save coefficients
	di _b[temp]/-2/_b[c.temp#c.temp]
loc min -5
margins, at(temp=(`min'(1)35)) post noestimcheck level(90)
parmest, norestore level(90)
split parm, p("." "#")
ren parm1 x
destring x, replace
replace x = x + `min' - 1  
drop parm* 
outsheet using data/output/estimatedGlobalResponse.csv, comma replace  //writing out results for R
use data/input/GrowthClimateDataset, clear
keep UDel_temp_popweight Pop TotGDP growthWDI GDPpctile_WDIppp continent iso countryname year
outsheet using data/output/mainDataset.csv, comma replace
clear
svmat b
outsheet using data/output/estimatedCoefficients.csv, comma replace


// -----------------  DATA FOR FIGURE 2, PANELS B, D, E  -----------------------

loc vars growthWDI AgrGDPgrowthCap NonAgrGDPgrowthCap 
foreach var of loc vars  {
use data/input/GrowthClimateDataset, clear
drop _yi_* _y2_* time time2
gen time = year - 1985
gen time2 = time^2
qui xi i.iso_id*time, pref(_yi_)  //linear country time trends
qui xi i.iso_id*time2, pref(_y2_) //quadratic country time trend
qui drop _yi_iso_id* 
qui drop _y2_iso_id* 
gen temp = UDel_temp_popweight 
gen poorWDIppp = (GDPpctile_WDIppp<50)
replace poorWDIppp=. if GDPpctile_WDIppp==.
gen interact = poorWDIppp
qui reg `var' interact#c.(c.temp##c.temp UDel_precip_popweight UDel_precip_popweight_2)  _yi_* _y2_* i.year i.iso_id, cl(iso_id) //PPP WDI baseline
loc min 0
margins, over(interact) at(temp=(`min'(1)30)) post noestimcheck force level(90)
parmest, norestore level(90)
split parm, p("." "#")
ren parm1 x
ren parm3 interact
destring x interact, replace
replace x = x - `min' - 1  
drop parm* 
gen model = "`var'"
if "`var'"=="growthWDI" {
	save data/output/EffectHeterogeneity, replace
	}
	else {
	append using data/output/EffectHeterogeneity
	save data/output/EffectHeterogeneity, replace
	}
}
use data/output/EffectHeterogeneity, clear
outsheet using data/output/EffectHeterogeneity.csv, comma replace


// -----------------  DATA FOR FIGURE 2, PANEL C  -----------------------

use data/input/GrowthClimateDataset, clear
drop _yi_* _y2_* time time2
gen time = year - 1985
gen time2 = time^2
qui xi i.iso_id*time, pref(_yi_)  //linear country time trends
qui xi i.iso_id*time2, pref(_y2_) //quadratic country time trend
qui drop _yi_iso_id* 
qui drop _y2_iso_id* 
gen temp = UDel_temp_popweight 
gen early = year<1990
gen interact = early
qui reg growthWDI interact#c.(c.temp##c.temp UDel_precip_popweight UDel_precip_popweight_2)  _yi_* _y2_* i.year i.iso_id, cl(iso_id) 
loc min 0
margins, over(interact) at(temp=(`min'(1)30)) post noestimcheck level(90)
parmest, norestore level(90)
split parm, p("." "#")
ren parm1 x
ren parm3 interact
destring x interact, replace
replace x = x - `min' - 1  
drop parm* 
outsheet using data/output/EffectHeterogeneityOverTime.csv, comma replace
