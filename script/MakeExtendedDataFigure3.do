
* SCRIPT TO GENERATE DATA FOR EXTENDED DATA FIGURE 3
*  USER NEEDS TO INSTALL USER-WRITTEN PROGRAM "ESTOUT"  (type "findit estout")
*  USER NEEDS TO ADD /SCRIPT DIRECTORY TO ADOPATH  (type "adopath + parent_directories/script") so Stata can find cgmreg.ado

clear all
set mem 1G
set matsize 10000
set maxvar 10000

cap mkdir data/output/DJO

adopath + "script/"


// CREATE DJO DATA FILE (FOLLOWING THEIR DO FILE) SO WE CAN COMPARE DATA

* This sample is the sample of all countries 
* Init GDP is defined based on your first year in the data
* Must have at least 20 years of GDP data
use data/input/DellJonesOlken/climate_panel, clear

* restrict to 2003
keep if year <= 2003

encode parent, g(parent_num)
gen lgdp=ln(rgdpl)
encode fips60_06, g(cc_num)
sort country_code year
tsset cc_num year
g lngdpwdi = ln(gdpLCU)

*calculate GDP growth
gen temp1 = l.lngdpwdi
gen g=lngdpwdi-temp1
replace g = g * 100 
drop temp1
summarize g

* Drop if less than 20 yrs of GDP data
g tempnonmis = 1 if g != .
replace tempnonmis = 0 if g == .
bys fips60_06: egen tempsumnonmis = sum(tempnonmis)
drop if tempsumnonmis  < 20
	
preserve
keep if lnrgdpl_t0 < . 
bys fips60_06: keep if _n == 1 
xtile initgdpbin = ln(lnrgdpl_t0), nq(2)
keep fips60_06 initgdpbin
tempfile tempxtile
save `tempxtile',replace
restore

merge n:1 fips60_06 using `tempxtile', nogen
tab initgdpbin, g(initxtilegdp)

preserve
keep if wtem50 < . 
bys fips60_06: keep if _n == 1 
xtile initwtem50bin = wtem50 , nq(2)
keep fips60_06 initwtem50bin
save `tempxtile',replace
restore

merge n:1 fips60_06 using `tempxtile', nogen
tab initwtem50bin, g(initxtilewtem)

tsset
foreach Y in wtem wpre  {
	gen `Y'Xlnrgdpl_t0 =`Y'*lnrgdpl_t0 
	for var initxtile*: gen `Y'_X =`Y'*X	
	label var `Y'Xlnrgdpl_t0 "`Y'.*inital GDP pc"
	for var initxtile*: label var `Y'_X "`Y'* X"
}
	tab year, gen (yr)
local numyears = r(r) - 1
//their region x year FE and poor x year FE
	foreach X of num 1/`numyears' {
			foreach Y in MENA SSAF LAC WEOFF EECA SEAS {
				quietly gen RY`X'X`Y'=yr`X'*_`Y'
				quietly tab RY`X'X`Y'
			}
			quietly gen RYPX`X'=yr`X'*initxtilegdp1
		}


g region=""
foreach X in _MENA   _SSAF   _LAC    _WEOFF  _EECA _SEAS {
	replace region="`X'" if `X'==1
}

g regionyear=region+string(year)
encode regionyear, g(rynum)	

rename initxtilegdp1 poor
keep poor g cc_num _MENA _SSAF _LAC _WEOFF _EECA _SEAS wtem wpre year fips RY* country_code country parent_num rynum
ren fips fips_cntry
ren country_code iso
ren country countryDJO
ren wtem temp //rename temp and precip variables so we can get them in same table
ren wpre prec
replace prec = prec/10 //re-scaling to meters
replace g = g/100 //to match our units
save data/output/DJO/DJOdata, replace

merge 1:1 iso year using data/input/GrowthClimateDataset

// NOW RUN REGRESSIONS ON DJO DATA, VARYING THE CONTROLS AND SAVING THE COEFFICIENTS AND THE FITS
cap est drop reg*
gen X = _n/3 if _n<90
cap postutil clear
postfile djo mod b1 b2 using data/output/DJO/Coefficients_DJOcompare, replace

//their original specification, global sample, adding a quadratic
gen temp2 = temp^2
cgmreg g temp temp2 RY* i.cc_num, cluster(parent_num rynum)
gen Ydjo1 = X*_b[temp] + X^2*_b[temp2]
est sto reg1
estadd sca opt = _b[temp]/-2/_b[temp2]
post djo (1) (_b[temp]) (_b[temp2])

//adding a quadratic and precip. no difference from without precip
gen prec2 = prec^2
cgmreg g temp temp2 prec prec2 RY* i.cc_num, cluster(parent_num rynum)
gen Ydjo2 = X*_b[temp] + X^2*_b[temp2]

// keeping continent-yr FE but adding a quadratic country trend
cap drop time time2
gen time = year - 1949
xi i.cc_num*time, pref(_yi_)  //generate country trends
gen time2 = time*time
qui xi i.cc_num*time2, pref(_y2_) //quadratic country time trend
drop _yi_cc_num* _y2_cc_num*
cgmreg g temp temp2 prec prec2 RY* i.cc_num _yi_* _y2_*, cluster(parent_num rynum)
gen Ydjo3 = X*_b[temp] + X^2*_b[temp2]
est sto reg3
estadd sca opt = _b[temp]/-2/_b[temp2]
post djo (2) (_b[temp]) (_b[temp2])

// year FE and trends, closest to ours
cgmreg g temp temp2 prec prec2 i.year i.cc_num _yi_* _y2_*, cluster(parent_num rynum)
gen Ydjo4 = X*_b[temp] + X^2*_b[temp2]
est sto reg4
estadd sca opt = _b[temp]/-2/_b[temp2]
post djo (3) (_b[temp]) (_b[temp2])


//calculate effects relative to a year at 20C and save
forvalues i = 1/4 {
	qui summ Ydjo`i' if X==20
	replace Ydjo`i' = Ydjo`i' - r(mean)
	}
keep Ydjo* X
drop if _n>90
save data/output/DJO/DJOresults, replace

//get the DJO estimation sample
use data/output/DJO/DJOdata, clear
qui reg g temp RY* i.cc_num
gen DJOsample = e(sample)
ren temp tempDJO
keep g temp year DJOsample iso countryDJO poor
save data/output/DJO/DJOsample, replace

//NOW SIMILAR EXERCISE WITH OUR DATA, MAKING OUR SAMPLE AND ESTIMATION LOOK MORE AND MORE LIKE THEIRS
use data/input/GrowthClimateDataset, clear
gen X = _n/3 if _n<90
gen temp = UDel_temp_popweight
gen temp2 = temp*temp
gen prec = UDel_precip_popweight/1000  //rescaling to meters to make coefficients legible
gen prec2 = prec^2
//first our base result
cgmreg growthWDI temp temp2 prec prec2 i.year _yi_* _y2_* i.iso_id if wdinomiss>=20, cluster(iso_id)
est sto reg5
estadd sca opt = _b[temp]/-2/_b[temp2]
gen Ybase = X*_b[temp] + X^2*_b[temp2]
post djo (4) (_b[temp]) (_b[temp2])

//restrict sample to DJO sample as best we can
merge 1:1 iso year using data/output/DJO/DJOsample, gen(DJOmerge)
br year countryname isocode iso countryDJO if DJOmerge==2
cgmreg growthWDI temp temp2 prec prec2 i.year _yi_* _y2_* i.iso_id if DJOmerge==3, cluster(iso_id)
est sto reg6
estadd sca opt = _b[temp]/-2/_b[temp2]
gen Yr = X*_b[temp] + X^2*_b[temp2]
post djo (5) (_b[temp]) (_b[temp2])
// restrict sample to DJO sample as best we can and use same FE
qui xi i.year*i.continent, pref(_cy_)  //continent x year FE
qui drop _cy_cont* _cy_year*
//gen poor = (GDPpctile_WDIppp<50)
//replace poor=. if GDPpctile_WDIppp==.
qui xi i.year*i.poor, pref(_py_)  //generating poor x year FE using their definition of poor countries
drop _py_poor _py_year_*
cgmreg growthWDI temp temp2 prec prec2 _cy_* _py_* i.iso_id if DJOmerge==3, cluster(iso_id)  //DJO sample and FE, our data
est sto reg7
estadd sca opt = _b[temp]/-2/_b[temp2]
gen Yrr = X*_b[temp] + X^2*_b[temp2]
post djo (6) (_b[temp]) (_b[temp2])
cgmreg g temp temp2 prec prec2 _cy_* _py_* i.iso_id if DJOmerge==3, cluster(iso_id)  //DJO sample and FE, their growth data our temperature data
est sto reg8
estadd sca opt = _b[temp]/-2/_b[temp2]
drop temp temp2
ren tempDJO temp
gen temp2 = temp*temp
cgmreg growthWDI temp temp2 prec prec2 _cy_* _py_* i.iso_id if DJOmerge==3, cluster(iso_id)  //DJO sample and FE, our growth data their temperature data
est sto reg9
estadd sca opt = _b[temp]/-2/_b[temp2]

lab var temp "Temperature"
lab var temp2 "Temperature squared"
lab var prec "Precip."
lab var prec2 "Precip. squared"
esttab reg* using figures/SI/TableS3.tex, drop(RY* *.cc_* _yi_* _y2_* *.year *.iso_id _cy_* _py_*)  ///
	replace b(%10.4f %10.4f %10.4f %10.4f %10.4f) se l sfmt(%10.0f %10.3f %10.2f) ///
	scalars("N Observations" "r2 R squared" "opt Optimum") ///
	star(* 0.10 ** 0.05 *** 0.01) compress nonotes width(\hsize) ///
	mtitles("base" "+trend" "yearFE+trend" "base" "DJOsample" "DJO FE" "DJO growth" "DJO temp") 
postclose djo


use data/output/DJO/Coefficients_DJOcompare, clear
outsheet using data/output/DJO/Coefficients_DJOcompare.csv, comma replace

