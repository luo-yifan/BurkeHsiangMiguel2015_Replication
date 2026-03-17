
* Create bottom panels of Extended Data Figure 2. 

clear all
set mem 1G
set matsize 10000
set maxvar 10000

global plots "figures/ExtendedDataFigs_Input"


// DATA FOR PANEL C
use data/input/GrowthClimateDataset, clear
xtset iso_id year
gen temp = UDel_temp_popweight
gen temp2 = temp^2
cap postutil clear
postfile lags nlag b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 b11 b12 using data/output/Coefficients_lags, replace
gen x = _n if _n<31
loc vals 0 1 3 5
foreach j of loc vals {
cap drop cilo cihi
qui reg growthWDI L(0/`j').(temp temp2 UDel_precip_popweight UDel_precip_popweight_2) i.year _yi_* _y2_* i.iso_id, cl(iso_id)
mat b = e(b)
gen Y`j'est=.
gen Y`j'se=.
forvalues i = 1/30 {
	if `j'==0 {
		qui lincom temp +  2*`i'*temp2
		}
	if `j'==1 {
		qui lincom temp + L.temp +  2*`i'*(temp2 + L.temp2)
		}
	if `j'==3 {
		qui lincom temp + L.temp + L2.temp + L3.temp + 2*`i'*(temp2 + L.temp2 + L2.temp2 + L3.temp2)
		}
	if `j'==5 {
		qui lincom temp + L.temp + L2.temp + L3.temp + L4.temp + L5.temp + 2*`i'*(temp2 + L.temp2 + L2.temp2 + L3.temp2 + L4.temp2 + L5.temp2)
		}
	qui replace Y`j'est= r(estimate) in `i'
	qui replace Y`j'se = r(se) in `i'
	}
gen cilo = Y`j'est - 1.96*Y`j'se
gen cihi = Y`j'est + 1.96*Y`j'se
	tw (rarea cilo cihi x, color(bluishgray)  legend(off) ylabel(,angle(horizontal))) ///
	(line Y`j'est x , lc(black) yline(0, lp(dash)) saving($plots/g`j'.gph, replace)) 
	post lags (`j') (b[1,1]) (b[1,2]) (b[1,3]) (b[1,4]) (b[1,5]) (b[1,6]) (b[1,7]) (b[1,8]) (b[1,9]) (b[1,10]) (b[1,11]) (b[1,12])
	// going to be saving stuff that we don't need for lag<5, but easiest way to do it
	}
graph combine $plots/g0.gph $plots/g1.gph $plots/g3.gph $plots/g5.gph, xsize(11) ysize(10) iscale(0.7) ycomm
graph export $plots/Fig2_PanelC.pdf, as(pdf) replace 
postclose lags
foreach j of loc vals {
	rm $plots/g`j'.gph
	}
* make matrix of estimated marginal effects at different points in the temperature distribution
keep x Y0est Y0se Y1est Y1se Y3est Y3se Y5est Y5se
keep if x==5 | x==10 | x==15 | x==20 | x == 25 | x == 30
//xpose, clear
//drop in 1
mkmat Y*, mat(res)
frmttable using $plots/TableS2.tex, statmat(res) substat(1) sdec(4) tex replace fra ///
	ctitles("", "Lags = 0", "Lags = 1", "Lags = 3", "Lags = 5") ///
	rtitles("5\C" \ "" \ "10\C" \ "" \ "15\C" \ "" \ "20\C" \"" \ "25\C" \ "" \ "30\C") ///
	hlines(11010101010101) spacebef(1101010101010)
	
	
* DATA FOR REST OF PANELS
use data/input/GrowthClimateDataset, clear
gen temp = UDel_temp_popweight 
gen prec = UDel_precip_popweight
gen poorWDIppp = (GDPpctile_WDIppp<50)
replace poorWDIppp=. if GDPpctile_WDIppp==.
gen temp2 = temp*temp
gen temppoor = temp*(poorWDIppp==1)
gen temp2poor = temp2*(poorWDIppp==1)
gen prec2 = prec*prec
gen precpoor = prec*(poorWDIppp==1)
gen prec2poor = prec2*(poorWDIppp==1)
gen T = _n if _n<=30
loc vars growthWDI AgrGDPgrowthCap NonAgrGDPgrowthCap
foreach var of loc vars  {
qui reg `var' temp temp2 temppoor temp2poor prec prec2 precpoor prec2poor  _yi_* _y2_* i.year i.iso_id if poorWDIppp<., cl(iso_id) 
gen `var'_br=. 
gen `var'_bp=.
gen `var'_bc=.
gen `var'_ser=.
gen `var'_sep=.
gen `var'_sec=.
gen `var'_pr=.
gen `var'_pp=.
gen `var'_pc=.
forvalues t = 1/30 {
	qui lincom (temp + 2*temp2*`t')  // test rich country marginal effect
	qui replace `var'_br = r(estimate) in `t'
	qui replace `var'_ser = r(se) in `t'
	qui replace `var'_pr =  2*ttail(e(df_r),abs(r(estimate)/r(se))) in `t'  //pvalue
	qui lincom (temp + temppoor + 2*temp2*`t' + 2*temp2poor*`t')  // test poor country marginal effect
	qui replace `var'_bp = r(estimate) in `t'
	qui replace `var'_sep = r(se) in `t'
	qui replace `var'_pp =  2*ttail(e(df_r),abs(r(estimate)/r(se)))  in `t' //pvalue
	qui lincom -(temppoor + 2*temp2poor*`t')  //difference in slopes: rich minus poor
	qui replace `var'_bc = r(estimate) in `t'
	qui replace `var'_sec = r(se) in `t'
	qui replace `var'_pc =  2*ttail(e(df_r),abs(r(estimate)/r(se)))  in `t' //pvalue
	}
di "`var'"
}
keep *_b? *_se? *_p?
gen x = _n
drop if x>30
outsheet using data/output/Effect_Marginals_RichPoor.csv, comma replace
