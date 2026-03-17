

*  Script to make Extended Data Table 2

clear all
set mem 1G
set matsize 10000
set maxvar 10000

use data/input/GrowthClimateDataset, clear
gen temp = UDel_temp_popweight 
gen prec = UDel_precip_popweight
gen poorWDIppp = (GDPpctile_WDIppp<50)
drop if GDPpctile_WDIppp==. // drop countries w/o PPP data 
gen temp2 = temp*temp
gen temppoor = temp*(poorWDIppp==1)
gen temp2poor = temp2*(poorWDIppp==1)
gen prec2 = prec*prec
gen precpoor = prec*(poorWDIppp==1)
gen prec2poor = prec2*(poorWDIppp==1)
xi i.year*i.continent, pref(_cy_)  //make continent-yr FE
loc c1 "_yi_* _y2_* i.year"
loc c2 "_yi_* _y2_* i.poorWDIppp#i.year"
loc c3 "_yi_* _y2_* i.year if wdinomiss>=20 "
loc c4 "_yi_* _y2_* i.poorWDIppp#i.year if wdinomiss>=20"
loc c5 "_yi_* _y2_* _cy_*"
loc c6 "_cy_*"
cap est drop reg*
forvalues z = 1/6 {
	areg growthWDI temp temp2 temppoor temp2poor prec prec2 precpoor prec2poor `c`z'', a(iso_id) cl(iso_id)
	est sto reg`z'
	qui lincom _b[temp] + _b[temppoor]
	estadd sca linb = r(estimate)
	estadd sca linse = r(se)
	qui lincom _b[temp2] + _b[temp2poor]
	estadd sca qb = r(estimate)
	estadd sca qse = r(se)
	estadd sca optr = _b[temp]/-2/_b[temp2]
	estadd sca optp = (_b[temp] + _b[temppoor])/(-2*(_b[temp2] + _b[temp2poor]))
	}	
lab var temp "Temperature ($\beta_1$)"
lab var temp2 "Temperature sq. ($\beta_2$)"
lab var temppoor "Temperature * poor ($\beta_3$)"
lab var temp2poor "Temperature sq * poor ($\beta_4$)"
esttab reg* using figures/ExtendedDataFigs_Input/Table2.tex, drop(prec* _yi_* _y2_* *.year _cy_* _cons)  ///
	replace b(%10.4f %10.4f %10.4f %10.4f %10.4f) se l sfmt(%10.0f %10.3f %10.4f %10.4f %10.4f %10.4f %10.1f %10.1f) ///
	scalars("N Observations" "r2 R squared" "linb $\beta_1 + \beta_3$" "linse se($\beta_1 + \beta_3$)" "qb $\beta_2 + \beta_4$" "qse se($\beta_2 + \beta_4$)" "optr Rich-country optimum (C)" "optp Poor-country optimum (C)") ///
	star(* 0.10 ** 0.05 *** 0.01) compress nonotes width(\hsize) substitute(\_ _) nocons ///
	mtitles("Base" "poor-yr FE" "$>$20Yrs" "$>$20Yrs+poor-yr FE" "ContYr FE" "ContYr + noTrend") 
