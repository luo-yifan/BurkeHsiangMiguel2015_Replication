
* Calculations for panels in Extended Data Figure 1
* 	users need to install user-written commands vallist and estout
*		ssc install vallist
*		ssc install estout


clear all
set mem 1G
set matsize 10000
set maxvar 10000

cap mkdir figures/SI
global plots "figures/ExtendedDataFigs_Input"


// Make sub-panels for Extended Data Fig 1
use data/input/GrowthClimateDataset, clear
gen temp = UDel_temp_popweight
gen year2 = year*year
loc ctys ISL FRA USA VNM MLI
foreach ct of loc ctys {
cap drop growthresid tempresid
qui reg growthWDI UDel_precip_popweight UDel_precip_popweight year year2 if iso=="`ct'"
qui predict growthresid if e(sample), resid
qui reg temp UDel_precip_popweight year year2 if iso=="`ct'"
qui predict tempresid if e(sample), resid
tw (lfitci growthresid tempresid, legend(off) xtitle("temperature deviations (C)")) ///
	(scatter growthresid tempresid, ms(circle) mc(black) title("`ct'") saving($plots/g`ct'.gph, replace) )
}
graph combine $plots/gISL.gph $plots/gFRA.gph $plots/gUSA.gph $plots/gVNM.gph $plots/gMLI.gph, rows(2) xsize(12) ysize(10) 
graph export $plots/Fig1_PanelB.pdf, as(pdf) replace
foreach ct of loc ctys {
	rm $plots/g`ct'.gph
	}

*  Now calculate response of each individual country in time series regression, for bottom panel of figure
use data/input/GrowthClimateDataset, clear
gen temp = UDel_temp_popweight
gen temp2 = temp*temp
areg growthWDI temp temp2 UDel_precip_popweight UDel_precip_popweight_2 i.year _yi_* _y2_*, a(iso_id) cl(iso_id)
loc b1 = _b[temp]
loc b2 = _b[temp2]
gen year2 = year*year
keep if e(sample)
postfile marg str3(iso) b se meantemp using data/output/Marg_Effects, replace
vallist iso
loc ctys `r(list)'
xtset iso_id year
foreach ct of loc ctys {
	di "`ct'"
	qui summ temp if iso=="`ct'" & growthWDI<.
	loc mn = r(mean)
	qui reg growthWDI temp UDel_precip_popweight year year2 if iso=="`ct'"
	post marg ("`ct'") (_b[temp]) (_se[temp]) (`mn')
	}	
postclose marg

use data/output/Marg_Effects, clear
gen fit = `b1' + 2*`b2'*meantemp
outsheet using data/output/ExtendedDataFig1g.csv, comma replace


// Now make last panel, interacting 
use data/input/GrowthClimateDataset, clear
gen temp = UDel_temp_popweight
gen temp2 = temp*temp
gen prec = UDel_precip_popweight
gen prec2 = prec*prec
egen tavg = mean(temp), by(iso_id)
egen pavg = mean(prec), by(iso_id)
gen tempavg = temp*tavg
gen precavg = prec*pavg
cap est drop reg*
areg growthWDI temp tempavg prec precavg i.year _yi_* _y2_*, a(iso_id) cluster(iso_id) //so this looks very much like the quadratic model, as expected
est sto reg1
gen marg1=.
gen se1=.
forvalues i = 1(2)30 {
	qui lincom temp + `i'*tempavg
	replace marg1 = r(estimate) in `i'
	replace se1 = r(se) in `i'
	}
// now run adding interaction with mean GDP/cap. demeaning the GDP/cap variable so that temperature terms can be interpreted at average income
egen mGDP = mean(gdpCAP_wdi), by(iso_id)
qui summ mGDP
g dmGDP = (mGDP - r(mean))/1000  //demean so temp coeffs are interpretable, and divide by 1000 to coefficients are visible
gen tempInc = temp*dmGDP
gen precInc = prec*dmGDP

areg growthWDI temp tempavg tempInc prec precavg precInc i.year _yi_* _y2_*, a(iso_id) cluster(iso_id) 
est sto reg2
gen marg2=.
gen se2=.
forvalues i = 1(2)30 {
	qui lincom temp + `i'*tempavg
	replace marg2 = r(estimate) in `i'
	replace se2 = r(se) in `i'
	}
	
//now make comparison plot	
gen cihi1 = marg1 + 1.96*se1
gen cilo1 = marg1 - 1.96*se1
gen cihi2 = marg2 + 1.96*se2
gen cilo2 = marg2 - 1.96*se2
gen n1 = _n-0.2 if _n<31
gen n2 = _n+0.2 if _n<31
preserve
keep n* ci* marg*
keep in 1/30
outsheet using data/output/ExtendedDataFig1h.csv, comma replace
restore

* Estimate a few additional interacted models with alternate FE (poor-by-year FE or continent-by-year FE) or log income, and write out as Table S1

xi i.year*i.continent, pref(_cy_)  //make continent-yr FE
areg growthWDI temp tempavg tempInc prec precavg precInc _cy_* _yi_* _y2_*, a(iso_id) cluster(iso_id) 
est sto reg3
gen lmGDP = log(mGDP)
qui summ lmGDP
g dlmGDP = lmGDP - r(mean)
gen tempIncL = temp* dlmGDP
gen precIncL = prec* dlmGDP
areg growthWDI temp tempavg tempIncL prec precavg precIncL i.year _yi_* _y2_*, a(iso_id) cluster(iso_id) 
est sto reg4

lab var temp "$T_{it}$"
lab var tempavg "$T_{it}*\bar{T_i}$"
lab var tempInc "$T_{it}*\bar{Y_i}$"
lab var tempIncL "$T_{it}*log(\bar{Y_i})$"
esttab reg* using figures/SI/TableS1.tex, drop(_yi_* _y2_* *.year _cy_* prec* _cons)  ///
	replace b(%10.4f ) se l sfmt(%10.0f %10.3f) ///
	scalars("N Observations" "r2 R squared" ) ///
	star(* 0.10 ** 0.05 *** 0.01) compress nonotes width(\hsize) substitute(\_ _) 
	


* Make panels i-j that shows how model performs under more flexible functional forms.  first with higher-order polynomials, then splines
use data/input/GrowthClimateDataset, clear
gen X = _n/3 if _n<100
gen temp = UDel_temp_popweight
loc add _b[temp]*X
forvalues i = 2/7 {
	gen temp`i' = UDel_temp_popweight^`i'
	qui reg growthWDI temp* UDel_precip_popweight UDel_precip_popweight_2 i.year _yi_* _y2_* i.iso_id, cl(iso_id)
	loc add  "`add' + _b[temp`i']*X^`i'"
	di "`add'"
	gen Y`i' = `add'
	}

// add the same thing running the restricted cubic spline with different numbers of knots
cap drop spline*
forvalues i = 3/7 {
mkspline spline`i'_ = temp, nknots(`i') cubic 
mat knot_vector = r(knots) // store vector of knots
loc knot_list ""
forvalues k = 1/`i' { 
	loc a = knot_vector[1,`k']
	loc knot_list "`knot_list' `a'"
}
loc Jsegments = `i'-1 // one less segment than there are knots
qui reg growthWDI spline`i'_* UDel_precip_popweight UDel_precip_popweight_2 i.year _yi_* _y2_* i.iso_id, cl(iso_id)
drop spline`i'*
mkspline spline`i'_ = X, cubic knots(`knot_list') // generate the same set of splines as those used in the regression
// compute the effect of X on Y by summing predictions for the splines
loc knot_inner_product "0"
forvalues j = 1/`Jsegments' {
	loc knot_inner_product "`knot_inner_product' + _b[spline`i'_`j']*spline`i'_`j'"
}
gen spline_est_`i' = `knot_inner_product'
}

keep X spline_est* Y?*
keep in 1/100
outsheet using data/output/ExtendedDataFig1i.csv, comma replace
