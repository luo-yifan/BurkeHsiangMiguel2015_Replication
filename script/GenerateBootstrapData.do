

* BOOTSTRAP THE REGRESSION, WITH AND WITHOUT LAGS, AND ASSUMING EITHER A POOLED OR DIFFERENTIATED RICH/POOR RESPONSE FUNCTION 
*	Sampling countries, with replacement

clear all
set mem 1G
set matsize 10000
set maxvar 10000

cap mkdir data/output/bootstrap

//  BOOTSTRAP MAIN SPECIFICATION 1000 times.  First run is baseline specification (i.e. not resampled); next thousand are resampled
cap postutil clear
postfile boot run temp temp2 prec prec2 using data/output/bootstrap/bootstrap_noLag, replace
set seed 8675309
use data/input/GrowthClimateDataset, clear
qui gen UDel_temp_popweight_2 = UDel_temp_popweight^2
qui reg growthWDI UDel_temp_popweight UDel_temp_popweight_2 UDel_precip_popweight UDel_precip_popweight_2 i.year _yi_* _y2_* i.iso_id 
post boot (0) (_b[UDel_temp_popweight]) (_b[UDel_temp_popweight_2]) (_b[UDel_precip_popweight]) (_b[UDel_precip_popweight_2])
// now bootstrap model, sampling countries
forvalues nn = 1/1000 {
use data/input/GrowthClimateDataset, clear
bsample, cl(iso_id)  //draw a sample of countries with replacement
qui gen UDel_temp_popweight_2 = UDel_temp_popweight^2
	qui reg growthWDI UDel_temp_popweight UDel_temp_popweight_2 UDel_precip_popweight UDel_precip_popweight_2 i.year _yi_* _y2_* i.iso_id 
	post boot (`nn') (_b[UDel_temp_popweight]) (_b[UDel_temp_popweight_2]) (_b[UDel_precip_popweight]) (_b[UDel_precip_popweight_2])
di "`nn'"  //print a counter to watch progress
}
postclose boot
use data/output/bootstrap/bootstrap_noLag, clear
outsheet using data/output/bootstrap/bootstrap_noLag.csv, comma replace


//ZERO LAG MODEL ALLOWING RICH/POOR TO RESPOND DIFFERENTLY
cap postutil clear
postfile boot run temp temppoor temp2 temp2poor prec precpoor prec2 prec2poor using data/output/bootstrap/bootstrap_richpoor, replace
set seed 8675309
use data/input/GrowthClimateDataset, clear
qui gen UDel_temp_popweight_2 = UDel_temp_popweight^2
qui gen poor = (GDPpctile_WDIppp<50)
qui replace poor=. if GDPpctile_WDIppp==.
qui reg growthWDI poor#c.(UDel_temp_popweight UDel_temp_popweight_2 UDel_precip_popweight UDel_precip_popweight_2) i.year _yi_* _y2_* i.iso_id 
mat b = e(b)
post boot (0) (b[1,1]) (b[1,2]) (b[1,3]) (b[1,4]) (b[1,5]) (b[1,6]) (b[1,7]) (b[1,8])
// now bootstrap model, sampling countries
forvalues nn = 1/1000 {
	use data/input/GrowthClimateDataset, clear
	bsample, cl(iso_id)  //draw a sample of countries with replacement
	qui gen UDel_temp_popweight_2 = UDel_temp_popweight^2
	qui gen poor = (GDPpctile_WDIppp<50)
	qui replace poor=. if GDPpctile_WDIppp==.
	qui reg growthWDI poor#c.(UDel_temp_popweight UDel_temp_popweight_2 UDel_precip_popweight UDel_precip_popweight_2) i.year _yi_* _y2_* i.iso_id 
	mat b = e(b)
	post boot (`nn') (b[1,1]) (b[1,2]) (b[1,3]) (b[1,4]) (b[1,5]) (b[1,6]) (b[1,7]) (b[1,8])
	di "`nn'"
}
postclose boot
use data/output/bootstrap/bootstrap_richpoor, clear
outsheet using data/output/bootstrap/bootstrap_richpoor.csv, comma replace


//POOLED MODEL WITH LAGS
postfile boot run temp L1temp L2temp L3temp L4temp L5temp temp2 L1temp2 L2temp2 L3temp2 L4temp2 L5temp2  using data/output/bootstrap/bootstrap_5Lag, replace
set seed 8675309
use data/input/GrowthClimateDataset, clear
xtset iso_id year
qui gen UDel_temp_popweight_2 = UDel_temp_popweight^2
qui reg growthWDI L(0/5).(UDel_temp_popweight UDel_temp_popweight_2 UDel_precip_popweight UDel_precip_popweight_2) i.year _yi_* _y2_* i.iso_id
mat b = e(b)
post boot (0) (b[1,1]) (b[1,2]) (b[1,3]) (b[1,4]) (b[1,5]) (b[1,6]) (b[1,7]) (b[1,8]) (b[1,9]) (b[1,10]) (b[1,11]) (b[1,12])
// now bootstrap model, sampling countries. 
forvalues nn = 1/1000 {
	use data/input/GrowthClimateDataset, clear
	bsample, cl(iso_id) idcluster(id) //draw a sample of countries with replacement
	xtset id year  //need to use the new cluster variable it creates. 
	qui gen UDel_temp_popweight_2 = UDel_temp_popweight^2	
	qui reg growthWDI L(0/5).(UDel_temp_popweight UDel_temp_popweight_2 UDel_precip_popweight UDel_precip_popweight_2) i.year _yi_* _y2_* i.iso_id
	mat b = e(b)
	post boot (`nn') (b[1,1]) (b[1,2]) (b[1,3]) (b[1,4]) (b[1,5]) (b[1,6]) (b[1,7]) (b[1,8]) (b[1,9]) (b[1,10]) (b[1,11]) (b[1,12])
	di "`nn'"
}
postclose boot
use data/output/bootstrap/bootstrap_5Lag, clear
gen tlin = temp + L1temp + L2temp + L3temp + L4temp + L5temp
gen tsq = temp2 + L1temp2 + L2temp2 + L3temp2 + L4temp2 + L5temp2
outsheet using data/output/bootstrap/bootstrap_5Lag.csv, comma replace



//MODEL WITH RICH POOR AND 5 LAGS
postfile boot run temp temppoor L1temp L1temppoor L2temp L2temppoor L3temp L3temppoor L4temp L4temppoor L5temp L5temppoor ///
	temp2 temp2poor L1temp2 L1temp2poor L2temp2 L2temp2poor L3temp2 L3temp2poor L4temp2 L4temp2poor L5temp2 L5temp2poor ///
	using data/output/bootstrap/bootstrap_richpoor_5lag, replace
set seed 8675309
use data/input/GrowthClimateDataset, clear
xtset iso_id year
qui gen UDel_temp_popweight_2 = UDel_temp_popweight^2
qui gen poor = (GDPpctile_WDIppp<50)
qui replace poor=. if GDPpctile_WDIppp==.
qui reg growthWDI poor#c.(L(0/5).(UDel_temp_popweight UDel_temp_popweight_2 UDel_precip_popweight UDel_precip_popweight_2)) i.year _yi_* _y2_* i.iso_id 
mat b = e(b)
post boot (0) (b[1,1]) (b[1,2]) (b[1,3]) (b[1,4]) (b[1,5]) (b[1,6]) (b[1,7]) (b[1,8]) (b[1,9]) (b[1,10]) (b[1,11]) (b[1,12]) (b[1,13]) (b[1,14]) (b[1,15]) (b[1,16]) (b[1,17]) (b[1,18]) (b[1,19]) (b[1,20]) (b[1,21]) (b[1,22]) (b[1,23]) (b[1,24])
// now bootstrap model, sampling countries. 
forvalues nn = 1/1000 {
	use data/input/GrowthClimateDataset, clear
	bsample, cl(iso_id) idcluster(id) //draw a sample of countries with replacement
	qui xtset id year  //need to use the new cluster variable it creates. 
	qui gen UDel_temp_popweight_2 = UDel_temp_popweight^2	
	qui gen poor = (GDPpctile_WDIppp<50)
	qui replace poor=. if GDPpctile_WDIppp==.
	qui reg growthWDI poor#c.(L(0/5).(UDel_temp_popweight UDel_temp_popweight_2 UDel_precip_popweight UDel_precip_popweight_2)) i.year _yi_* _y2_* i.iso_id 
	mat b = e(b)
	post boot (`nn') (b[1,1]) (b[1,2]) (b[1,3]) (b[1,4]) (b[1,5]) (b[1,6]) (b[1,7]) (b[1,8]) (b[1,9]) (b[1,10]) (b[1,11]) (b[1,12]) (b[1,13]) (b[1,14]) (b[1,15]) (b[1,16]) (b[1,17]) (b[1,18]) (b[1,19]) (b[1,20]) (b[1,21]) (b[1,22]) (b[1,23]) (b[1,24])
	di "`nn'"
}
postclose boot
use data/output/bootstrap/bootstrap_richpoor_5lag, clear
gen tlin = temp + L1temp + L2temp + L3temp + L4temp + L5temp
gen tlinpoor = temppoor + L1temppoor + L2temppoor + L3temppoor + L4temppoor + L5temppoor
gen tsq = temp2 + L1temp2 + L2temp2 + L3temp2 + L4temp2 + L5temp2
gen tsqpoor = temp2poor + L1temp2poor + L2temp2poor + L3temp2poor + L4temp2poor + L5temp2poor
outsheet using data/output/bootstrap/bootstrap_richpoor_5lag.csv, comma replace


