
# Script to construct our main impact projections under future climate change
# MB, March 2015

# Main steps in the script:
#   1. Assembles the necessary scenarios on future population and income growth.  These come from the
#       Shared Socioeconomic Pathways (SSPs), as well as from a "baseline" scenario that fixes future growth rates
#       in absence of climate change at historical rates as.  These scenarios are saved as lists.
#   2. Reads in projected country-level temperature change.  These are calculated in the getTemperatureChange.R script
#   3. Projects changes in per capita GDP for different historical regression models, and for each of these, 
#         different population/income projections.  
# The main output is country- and year-specific per capita GDP projections to 2100, with and without climate change, 
#     for multiple regression models and multiple future pop/income scenarios.  These are written out. 


rm(list = ls())

require(fields)
require(classInt)
require(plotrix)
require(data.table)
require(dplyr)
require(raster)
library(maps)
library(ggplot2)
library(xtable)
library(reshape2)
"%&%"<-function(x,y)paste(x,y,sep="")  #define a function for easy string pasting

dir.create("data/output/projectionOutput")

# FIRST DEFINE FUNCTION TO INTERPOLATE SSP POPULATION AND GROWTH DATASETS.  THESE COME AS 5-YR PROJECTIONS. we want opbs for each year for each projection, so linearly interpolate between 5-yr estimates
#   writing this as a function that spits out a data frame, where the first three columns gives the scenario and country, and the rest give the projections by year
#   growth projections only go through 2095, so repeating the 2095 value for 2096-2099
ipolate <- function(mat) {
  mat1 <- array(dim=c(dim(mat)[1],length(yrs)))
  ys <- as.numeric(unlist(strsplit(names(mat),"X")))
  est <- seq(2010,2100,5)  #the 5yr estimates in the SSP dataset
  for (i in 1:length(yrs)) {
    y = yrs[i]
    if ("X"%&%y %in% names(pop) == T) {  #if the year falls on the 5-yr interval, use their point estimate. otherwise interpolate between nearest endpoints
      mat1[,i]  <- as.numeric(mat[,which(names(mat)=="X"%&%y)])
    } else {
      z <- y-est
      yl = est[which(z==min(z[z>0]))]  #the 5-year endpoint lower than the year
      y5 = yl+5  #the next endpoint
      el <- as.numeric(mat[,which(names(mat)=="X"%&%yl)])  #values at lower endpoint
      eu <- as.numeric(mat[,which(names(mat)=="X"%&%y5)]) #values at upper endpoint
      if (y > max(ys,na.rm=T)) {  mat1[,i] <- el   #this is to account for growth projections ending in 2095  
      }  else { mat1[,i] <- el + (eu-el)*(y-yl)/5 }
    }
  } 
  mat1 <- data.frame(mat[,1:3],mat1)
  names(mat1)[4:dim(mat1)[2]] <- yrs
  levels(mat1$Region)[levels(mat1$Region)=="COD"] <- "ZAR"  #our code for the DRC
  levels(mat1$Region)[levels(mat1$Region)=="ROU"] <- "ROM"  #our code for Romania  
  return(mat1)
}

# Get baseline mean growth rate and temperature, using 1980-2010. 
dta <- read.csv("data/output/mainDataset.csv")
gdpCap = dta$TotGDP/dta$Pop
dta <- data.frame(dta,gdpCap)
mt <- dta %>%   #the following few lines gets the average temperature in each country for the years we want, using dplyr
  filter(year>=1980 & is.na(UDel_temp_popweight)==F & is.na(growthWDI)==F) %>% 
  group_by(iso) %>% 
  summarize(meantemp = mean(UDel_temp_popweight,na.rm=T), basegrowth = mean(growthWDI, na.rm=T), gdpCap = mean(gdpCap,na.rm=T))
mt <- as.data.frame(mt)
yrs <- 2010:2099

# We will also run our own scenario, where future growth rates are just set at their historical average for each country, and pop projections come from UN
#   Need to get these in same shape as SSP above so they can easily be fed into the calculations below
# First process the UN population data
pop <- read.csv("data/input/populationData/WPP2012_DB02_POPULATIONS_ANNUAL.csv")
# reshape to wide, just keeping yrs 2010 and after
pop <- pop[pop$Time%in%2010:2100,]  
pop <- pop[,c(1,5,9)]
pop <- reshape(pop,v.names="PopTotal",timevar="Time",idvar="LocID",direction="wide")
pop <- pop[pop$LocID<900,]
cod <- read.csv("data/input/populationData/WPP2012_F01_LOCATIONS.csv")  #this file matches location ID in UN database to our ISO codes, with excpetion of DRC and ROM, as in the SSP
cod <- cod[cod$LocID<900,]
iso <- as.character(cod$ISO3_Code)
iso[iso=="COD"] <- "ZAR"
iso[iso=="ROU"] <- "ROM"
LocID <- cod$LocID
cod <- data.frame(LocID, iso)
cod <- cod[iso!="GRL",]  #dropping greenland because not in SSP
poploc <- merge(pop,cod,by="LocID")
popprojb <- merge(mt,poploc,by="iso")  #165 countries which we have both baseline data and pop projections for
popprojb <- popprojb[,names(popprojb)!="LocID"]
names(popprojb)[5:dim(popprojb)[2]] <- yrs
popprojb[5:dim(popprojb)[2]] <- popprojb[5:dim(popprojb)[2]]/1000 #units will be in millions to match SSP
basegrowth <- matrix(rep(popprojb$basegrowth,length(yrs)),ncol=length(yrs),byrow=F) #same growth rate for every year, which is the baseline growth rate
colnames(basegrowth) <- yrs
basegrowth <- cbind(popprojb[,1:4],basegrowth)

popProjections <- NULL  #initialize list that we will will with population projections for each scenario
growthProjections <- NULL   #same for growth
popProjections[[1]] <- popprojb[]  
growthProjections[[1]] <- basegrowth

# ADD IN PROJECTIONS FROM SSP
# read in data and interpolate
pop <- read.csv("data/input/SSP/SSP_PopulationProjections.csv")
levels(pop$Scenario)[levels(pop$Scenario)=="SSP4d_v9_130115"] <- "SSP4_v9_130115"  #renaming one of the scenarios slightly so the loop works
growth <- read.csv("data/input/SSP/SSP_GrowthProjections.csv")
pop1 <- ipolate(pop)  #warning here is just from stringsplit function
growth1 <- ipolate(growth)
growth1[,names(growth1)%in%yrs] = growth1[,names(growth1)%in%yrs]/100

#   First we merge countries in historical database with the growth and pop projections from SSP, restricted to the scenario we want
# we are using growth projections from OECD, which are the only ones with data for every country; population projections are from IIASA
popSSP <- merge(mt,pop1,by.x="iso",by.y="Region")  #merge our data and SSP for population
# length(unique(popproj$iso))  #165 countries that we can match in our data to SSP data
growthSSP <- merge(mt,growth1,by.x="iso",by.y="Region")

for (scen in 1:5) {  #now add each scenario to the list
  # projections for economic growth - using OECD, because they have projections for every country
  pgrow <- growthSSP$Model=="OECD Env-Growth" & growthSSP$Scenario=="SSP"%&%scen%&%"_v9_130325"
  growthProjections[[scen+1]] <- growthSSP[pgrow,]
  
  # population projections from IIASA
  ppop <- popSSP$Scenario=="SSP"%&%scen%&%"_v9_130115"
  popProjections[[scen+1]] <- popSSP[ppop,]
}

# SAVE THIS SCENARIO DATA TO BE USED IN CONSTRUCTION OF DAMAGE FUNCTION
save(popProjections,file="data/output/projectionOutput/popProjections.Rdata")
save(growthProjections,file="data/output/projectionOutput/growthProjections.Rdata")

# Finally, read in projections of future temperature change, generated by the getTemperatureChange.R script
#   These are population-weighted country level projections averaged across all CMIP5 models
Tchg <- read.csv("data/input/CCprojections/CountryTempChange_RCP85.csv")
Tchg <- merge(popProjections[[1]][,1:3],Tchg,by.x="iso",by.y="GMI_CNTRY")
Tchg <- Tchg[Tchg$CNTRY_NAME%in%c("West Bank","Gaza Strip","Bouvet Island")==F,]
Tchg_named <- setNames(Tchg$Tchg, Tchg$iso)  # named vector for per-scenario country lookup
Tchg <- Tchg$Tchg # just keep the vector of temperature changes, since sort order is now correct


#####################################################################################################################
#  NOW WE ARE GOING TO RUN FUTURE PROJECTIONS OF IMPACTS BASED ON DIFFERENT HISTORICAL MODELS,
#     1. pooled model, zero lag
#     2. pooled model, 5 lags
#     3. rich/poor model, zero lag
#     4. rich/poor model, 5 lags
#  	Calculating impacts relative to base period of 1980-2010; calculating each growth/pop scenario for each model
#   Writing results for each model out, so we can easily make plots of multiple scenarios if desired.  
#   For each we write out:  
#       (i) country specific GDP/cap for climate change and no climate change scenarios
#       (ii) global average GDP/cap, and total global GDP, for climate change and no climate change scenarios 
#####################################################################################################################

# set temperature change for all runs
dtm <- Tchg   #the country-specific changes
scens <- c("base","SSP"%&%1:5)
ccd <- dtm/length(yrs)  #rate of increase in temperature per year.  

####################################################################################################
# POOLED MODEL WITH NO LAGS
####################################################################################################

prj <- read.csv("data/output/bootstrap/bootstrap_noLag.csv")  #bootstrapped projections, all obs
np = dim(prj)[1]  #number of bootstrap replicates we ran. the first one is the baseline model

# now loop over SSP scenarios and our own scenario. just doing base, SSP3, SSP5 for now, since those appear to be the relevant ones
for (scen in c(1,4,6)) {

  growthproj <- growthProjections[[scen]]  #select growth projections
  popproj <- popProjections[[scen]]  # select population projections
  # Align countries across growth, population, and temperature-change data
  iso_g <- as.character(growthproj[, 1]); iso_p <- as.character(popproj$iso)
  common_iso <- Reduce(intersect, list(iso_g, iso_p, names(Tchg_named)))
  growthproj <- growthproj[iso_g %in% common_iso, ]
  popproj    <- popproj[match(common_iso, iso_p), ]
  ccd        <- Tchg_named[common_iso] / length(yrs)

  basegdp = popproj$gdpCap  #baseline GDP/cap
  temp <- popproj$meantemp  #baseline temperature.
  
  GDPcapCC = GDPcapNoCC = array(dim=c(dim(growthproj)[1],length(yrs),np))  #array to fill with GDP/cap for each country
  dimnames(GDPcapCC) <- dimnames(GDPcapNoCC) <- list(growthproj[,1],yrs,1:np)
  GDPcapCC[,1,] = GDPcapNoCC[,1,] = basegdp  #initialize with baseline per cap GDP
  tots = array(dim=c(np,length(yrs),4))  #array to hold average global per cap GDP and total global GDP across scenarios, with and without climate change
  dimnames(tots) <- list(1:np,yrs,c("avgGDPcapCC","avgGDPcapNoCC","TotGDPCC","TotGDPNoCC"))
  
  for (tt in 1:np) {  #looping over bootstrap estimates
    bg = prj$temp[tt]*temp + prj$temp2[tt]*temp*temp  #this finds the predicted growth level for each country's temperature for the particular bootstrap run
    for (i in 2:length(yrs)) {
      j = i - 1
      y = yrs[i]
      basegrowth <- growthproj[,which(names(growthproj)==y)]  #growth rate without climate change
      GDPcapNoCC[,i,tt] = GDPcapNoCC[,j,tt]*(1+basegrowth)  #last year's per cap GDP times this years growth rate, as projected by scenario
      newtemp = temp+j*ccd
      dg = prj$temp[tt]*newtemp + prj$temp2[tt]*newtemp*newtemp  #predicted growth under new temperature
      dg[newtemp>30] = prj$temp[tt]*30 + prj$temp2[tt]*30*30  #constrain response to response at 30C if temp goes above that.  this is so we are not projecting out of sample
      
      diff = dg - bg  #difference between predicted baseline growth and predicted growth under new temp
      GDPcapCC[,i,tt] = GDPcapCC[,j,tt]*(1+basegrowth + diff)  #last year's GDPcap (w/ climate change) times climate-adjusted growth rate for this year
      
      #now calculate global average per cap GDP, weighting by population
      wt = popproj[,which(names(popproj)==y)]  #population weights 
      tots[tt,i,1] <- round(weighted.mean(GDPcapCC[,i,tt],wt),3)  # per cap GDP, with climate change
      tots[tt,i,2] <- round(weighted.mean(GDPcapNoCC[,i,tt],wt),3)  # per cap GDP, no climate change
      #total GDP with and without climate change. multiplying by 1e6 because population is in millions
      tots[tt,i,3] <- sum(GDPcapCC[,i,tt]*wt*1e6)  #with climate change
      tots[tt,i,4] <- sum(GDPcapNoCC[,i,tt]*wt*1e6) #without climate change
    }
  }
  #write out scenario specific results
  GDPcapCC <- round(GDPcapCC,3) #round to nearest dollar
  GDPcapNoCC <- round(GDPcapNoCC,3)  #ditto
  save(GDPcapCC,file="data/output/projectionOutput/GDPcapCC_pooled_"%&%scens[scen]%&%".Rdata")
  save(GDPcapNoCC,file="data/output/projectionOutput/GDPcapNoCC_pooled_"%&%scens[scen]%&%".Rdata")
  save(tots,file="data/output/projectionOutput/GlobalChanges_pooled_"%&%scens[scen]%&%".Rdata")
  GDPcapCC <- GDPcapNoCC <- NULL
  print(scen)
  
}


####################################################################################################
# RICH/POOR MODEL WITH NO LAGS
####################################################################################################

#   We have different response functions for countries that were above or below median income in the historical sample.  Countries move onto
#     the rich-country response function in the future if their income rises above the historical median income (and vice versa if it falls).
#     Depending where they are in the temperature distribution when this happens, the marginal effect of temperature on growth could change sign

prj <- read.csv("data/output/bootstrap/bootstrap_richpoor.csv")  #interacted model, all obs
np = dim(prj)[1]  #number of bootstrap replicates we ran. the first one is the one from the baseline model

for (scen in c(1,4,6)) {

  growthproj <- growthProjections[[scen]]
  popproj <- popProjections[[scen]]
  # Align countries across growth, population, and temperature-change data
  iso_g <- as.character(growthproj[, 1]); iso_p <- as.character(popproj$iso)
  common_iso <- Reduce(intersect, list(iso_g, iso_p, names(Tchg_named)))
  growthproj <- growthproj[iso_g %in% common_iso, ]
  popproj    <- popproj[match(common_iso, iso_p), ]
  ccd        <- Tchg_named[common_iso] / length(yrs)

  basegdp = popproj$gdpCap  #baseline GDP/cap
  medgdp <- median(basegdp)
  temp <- popproj$meantemp  #baseline temperature.
  
  GDPcapCC = GDPcapNoCC = array(dim=c(dim(growthproj)[1],length(yrs),np))  #array to fill with GDP/cap for each country
  dimnames(GDPcapCC) <- dimnames(GDPcapNoCC) <- list(growthproj[,1],yrs,1:np)
  GDPcapCC[,1,] = GDPcapNoCC[,1,] = basegdp  #initialize with baseline per cap GDP
  tots = array(dim=c(np,length(yrs),4))  #array to hold average global per cap GDP and total global GDP across scenarios, with and without climate change
  dimnames(tots) <- list(1:np,yrs,c("avgGDPcapCC","avgGDPcapNoCC","TotGDPCC","TotGDPNoCC"))
  
  for (tt in 1:np) {  #looping over bootstrap estimates
    for (i in 2:length(yrs)) {
      j = i - 1
      y <- yrs[i]
      #baseline growth rate from SSP
      basegrowth <- growthproj[,which(names(growthproj)==y)]
      
      #was country poor last year in climate change world?
      poor <- GDPcapCC[,j,tt]<=medgdp
      
      #calculate appropriate predicted growth rates in world without climate change
      bg = prj$temp[tt]*temp + prj$temp2[tt]*temp*temp  #predicted growth level for rich countries absent climate change
      bg[poor] = prj$temppoor[tt]*temp[poor] + prj$temp2poor[tt]*temp[poor]*temp[poor]  #predicted growth level for poor countries absent climate change
      
      GDPcapNoCC[,i,tt] = GDPcapNoCC[,j,tt]*(1+basegrowth)  #last year's per cap GDP times this years growth rate
      newtemp = temp+j*ccd
      dg = prj$temp[tt]*newtemp + prj$temp2[tt]*newtemp*newtemp  #predicted growth under new temperature for rich countries
      dg[poor] = prj$temppoor[tt]*newtemp[poor] + prj$temp2poor[tt]*newtemp[poor]*newtemp[poor]  #predicted growth for poor countries
      
      dg[newtemp>30 & poor==0] = prj$temp[tt]*30 + prj$temp2[tt]*30*30  #constrain response to response at 30C if temp goes above that for rich countries.  this is so we are not projecting out of sample
      dg[newtemp>30 & poor==1] = prj$temppoor[tt]*30 + prj$temp2poor[tt]*30*30  #constrain response to response at 30C for poor countries  
      diff = dg - bg  #difference between predicted baseline growth and predicted growth under new temp
      GDPcapCC[,i,tt] = GDPcapCC[,j,tt]*(1+basegrowth + diff)  #last year's GDPcap (w/ climate change) times climate adjusted growth rate for this year
      
      #now calculate global average per cap GDP, weighting by population
      wt = popproj[,which(names(popproj)==y)]  #population weights 
      tots[tt,i,1] <- round(weighted.mean(GDPcapCC[,i,tt],wt),3)
      tots[tt,i,2] <- round(weighted.mean(GDPcapNoCC[,i,tt],wt),3)
      #total GDP with and without climate change. multiplying by 1e6 because population is in millions
      tots[tt,i,3] <- sum(GDPcapCC[,i,tt]*wt*1e6)  #with climate change
      tots[tt,i,4] <- sum(GDPcapNoCC[,i,tt]*wt*1e6) #without climate change
    }
  }
  #write out scenario specific results
  GDPcapCC <- round(GDPcapCC,3) #round to nearest dollar
  GDPcapNoCC <- round(GDPcapNoCC,3)  #ditto
  save(GDPcapCC,file="data/output/projectionOutput/GDPcapCC_richpoor_"%&%scens[scen]%&%".Rdata")
  save(GDPcapNoCC,file="data/output/projectionOutput/GDPcapNoCC_richpoor_"%&%scens[scen]%&%".Rdata")
  save(tots,file="data/output/projectionOutput/GlobalChanges_richpoor_"%&%scens[scen]%&%".Rdata")
  GDPcapCC <- GDPcapNoCC <- NULL
  print(scen)
  
}


####################################################################################################
# POOLED MODEL WITH 5 LAGS
####################################################################################################

prj <- read.csv("data/output/bootstrap/bootstrap_5lag.csv")  #bootstrapped projections of model with 5 lags
np = dim(prj)[1]  #number of bootstrap replicates we ran. the first one is the one from the baseline model

# now loop over SSP scenarios and our own scenario. just doing base, SSP3, SSP5 for now, since those appear to be the relevant ones
for (scen in c(1,4,6)) {

  growthproj <- growthProjections[[scen]]
  popproj <- popProjections[[scen]]
  # Align countries across growth, population, and temperature-change data
  iso_g <- as.character(growthproj[, 1]); iso_p <- as.character(popproj$iso)
  common_iso <- Reduce(intersect, list(iso_g, iso_p, names(Tchg_named)))
  growthproj <- growthproj[iso_g %in% common_iso, ]
  popproj    <- popproj[match(common_iso, iso_p), ]
  ccd        <- Tchg_named[common_iso] / length(yrs)

  basegdp = popproj$gdpCap  #baseline GDP/cap
  temp <- popproj$meantemp  #baseline temperature.
  
  GDPcapCC = GDPcapNoCC = array(dim=c(dim(growthproj)[1],length(yrs),np))  #array to fill with GDP/cap for each country
  dimnames(GDPcapCC) <- dimnames(GDPcapNoCC) <- list(growthproj[,1],yrs,1:np)
  GDPcapCC[,1,] = GDPcapNoCC[,1,] = basegdp  #initialize with baseline per cap GDP
  tots = array(dim=c(np,length(yrs),4))  #array to hold average global per cap GDP and total global GDP across scenarios, with and without climate change
  dimnames(tots) <- list(1:np,yrs,c("avgGDPcapCC","avgGDPcapNoCC","TotGDPCC","TotGDPNoCC"))
  
  for (tt in 1:np) {  #looping over bootstrap estimates
    bg = prj$tlin[tt]*temp + prj$tsq[tt]*temp*temp  #this finds the predicted growth level for each country's temperature for the particular bootstrap run
    for (i in 2:length(yrs)) {
      j = i - 1
      y = yrs[i]
      basegrowth <- growthproj[,which(names(growthproj)==y)]
      GDPcapNoCC[,i,tt] = GDPcapNoCC[,j,tt]*(1+basegrowth)  #last year's per cap GDP times this years growth rate, as projected by scenario
      newtemp = temp+j*ccd
      dg = prj$tlin[tt]*newtemp + prj$tsq[tt]*newtemp*newtemp  #predicted growth under new temperature
      dg[newtemp>30] = prj$tlin[tt]*30 + prj$tsq[tt]*30*30  #constrain response to response at 30C if temp goes above that.  this is so we are not projecting out of sample
      
      diff = dg - bg  #difference between predicted baseline growth and predicted growth under new temp
      GDPcapCC[,i,tt] = GDPcapCC[,j,tt]*(1+basegrowth + diff)  #last year's GDPcap (w/ climate change) times climate-adjusted growth rate for this year
      
      #now calculate global average per cap GDP, weighting by population
      wt = popproj[,which(names(popproj)==y)]  #population weights 
      tots[tt,i,1] <- round(weighted.mean(GDPcapCC[,i,tt],wt),3)
      tots[tt,i,2] <- round(weighted.mean(GDPcapNoCC[,i,tt],wt),3)
      #total GDP with and without climate change. multiplying by 1e6 because population is in millions
      tots[tt,i,3] <- sum(GDPcapCC[,i,tt]*wt*1e6)  #with climate change
      tots[tt,i,4] <- sum(GDPcapNoCC[,i,tt]*wt*1e6) #without climate change
    }
  }
  #write out scenario specific results
  GDPcapCC <- round(GDPcapCC,3) #round to nearest dollar
  GDPcapNoCC <- round(GDPcapNoCC,3)  #ditto
  save(GDPcapCC,file="data/output/projectionOutput/GDPcapCC_pooled5lag_"%&%scens[scen]%&%".Rdata")
  save(GDPcapNoCC,file="data/output/projectionOutput/GDPcapNoCC_pooled5lag_"%&%scens[scen]%&%".Rdata")
  save(tots,file="data/output/projectionOutput/GlobalChanges_pooled5lag_"%&%scens[scen]%&%".Rdata")
  GDPcapCC <- GDPcapNoCC <- NULL
  print(scen)
  
}


####################################################################################################
# RICH/POOR MODEL WITH 5 LAGS
####################################################################################################

prj <- read.csv("data/output/bootstrap/bootstrap_richpoor_5lag.csv")  #interacted model, all obs
np = dim(prj)[1]  #number of bootstrap replicates we ran. the first one is the one from the baseline model

for (scen in c(1,4,6)) {

  growthproj <- growthProjections[[scen]]
  popproj <- popProjections[[scen]]
  # Align countries across growth, population, and temperature-change data
  iso_g <- as.character(growthproj[, 1]); iso_p <- as.character(popproj$iso)
  common_iso <- Reduce(intersect, list(iso_g, iso_p, names(Tchg_named)))
  growthproj <- growthproj[iso_g %in% common_iso, ]
  popproj    <- popproj[match(common_iso, iso_p), ]
  ccd        <- Tchg_named[common_iso] / length(yrs)
  
  basegdp = popproj$gdpCap  #baseline GDP/cap
  medgdp <- median(basegdp)
  temp <- popproj$meantemp  #baseline temperature.  
  
  GDPcapCC = GDPcapNoCC = array(dim=c(dim(growthproj)[1],length(yrs),np))  #array to fill with GDP/cap for each country
  dimnames(GDPcapCC) <- dimnames(GDPcapNoCC) <- list(growthproj[,1],yrs,1:np)
  GDPcapCC[,1,] = GDPcapNoCC[,1,] = basegdp  #initialize with baseline per cap GDP
  tots = array(dim=c(np,length(yrs),4))  #array to hold average global per cap GDP and total global GDP across scenarios, with and without climate change
  dimnames(tots) <- list(1:np,yrs,c("avgGDPcapCC","avgGDPcapNoCC","TotGDPCC","TotGDPNoCC"))
  
  for (tt in 1:np) {  #looping over bootstrap estimates
    for (i in 2:length(yrs)) {
      j = i - 1
      y <- yrs[i]
      #baseline growth rate from SSP
      basegrowth <- growthproj[,which(names(growthproj)==y)]
      
      #was country poor last year in climate change world?
      poor <- GDPcapCC[,j,tt]<=medgdp
      
      #calculate appropriate predicted growth rates in world without climate change
      bg = prj$tlin[tt]*temp + prj$tsq[tt]*temp*temp  #predicted growth level for rich countries absent climate change
      bg[poor] = prj$tlinpoor[tt]*temp[poor] + prj$tsqpoor[tt]*temp[poor]*temp[poor]  #predicted growth level for poor countries absent climate change
      
      GDPcapNoCC[,i,tt] = GDPcapNoCC[,j,tt]*(1+basegrowth)  #last year's per cap GDP times this years growth rate
      newtemp = temp+j*ccd
      dg = prj$tlin[tt]*newtemp + prj$tsq[tt]*newtemp*newtemp  #predicted growth under new temperature for rich countries
      dg[poor] = prj$tlinpoor[tt]*newtemp[poor] + prj$tsqpoor[tt]*newtemp[poor]*newtemp[poor]  #predicted growth for poor countries
      
      dg[newtemp>30 & poor==0] = prj$tlin[tt]*30 + prj$tsq[tt]*30*30  #constrain response to response at 30C if temp goes above that for rich countries.  this is so we are not projecting out of sample
      dg[newtemp>30 & poor==1] = prj$tlinpoor[tt]*30 + prj$tsqpoor[tt]*30*30  #constrain response to response at 30C for poor countries  
      
      diff = dg - bg  #difference between predicted baseline growth and predicted growth under new temp
      GDPcapCC[,i,tt] = GDPcapCC[,j,tt]*(1+basegrowth + diff)  #last year's GDPcap (w/ climate change) times climate adjusted growth rate for this year
      
      #now calculate global average per cap GDP, weighting by population
      wt = popproj[,which(names(popproj)==y)]  #population weights 
      tots[tt,i,1] <- round(weighted.mean(GDPcapCC[,i,tt],wt),3)
      tots[tt,i,2] <- round(weighted.mean(GDPcapNoCC[,i,tt],wt),3)
      #total GDP with and without climate change. multiplying by 1e6 because population is in millions
      tots[tt,i,3] <- sum(GDPcapCC[,i,tt]*wt*1e6)  #with climate change
      tots[tt,i,4] <- sum(GDPcapNoCC[,i,tt]*wt*1e6) #without climate change
    }
  }
  #write out scenario specific results
  GDPcapCC <- round(GDPcapCC,3) #round to nearest dollar
  GDPcapNoCC <- round(GDPcapNoCC,3)  #ditto
  save(GDPcapCC,file="data/output/projectionOutput/GDPcapCC_richpoor5lag_"%&%scens[scen]%&%".Rdata")
  save(GDPcapNoCC,file="data/output/projectionOutput/GDPcapNoCC_richpoor5lag_"%&%scens[scen]%&%".Rdata")
  save(tots,file="data/output/projectionOutput/GlobalChanges_richpoor5lag_"%&%scens[scen]%&%".Rdata")
  GDPcapCC <- GDPcapNoCC <- NULL
  print(scen)
  
}

