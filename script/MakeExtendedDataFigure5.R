
# Make Extended data FIGURE 5 

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

#load functions we need
weightProj <- dget("script/weightProj.R")
layOut <- dget("script/layOut.R")

##########################################################################################
# FIRST GENERATE GLOBAL UNCERTAINTY PLOTS
##########################################################################################

yrs <- 2010:2099
scens <- c("base","SSP"%&%1:5)
runs <- c("pooled","pooled5lag","richpoor","richpoor5lag")
rn <- c("Pooled SR", "Pooled LR", "Differentiated SR", "Differentiated LR")
pl <- NULL  #will hold plots
yz = seq(-100,100,0.5)  #vector of impacts over which densities should be evaluated
for (j in runs)  {
  for (scen in c(1,4,6)) {
    load("data/output/projectionOutput/GlobalChanges_"%&%j%&%"_"%&%scens[scen]%&%".Rdata")
    chgs <- tots[,,1]/tots[,,2]-1
    chgs <- chgs*100
    chgs[,1] <- 0  # first year is end of baseline, so no change
    pp <- weightProj(chgs,scens[scen]%&%", "%&%rn[which(runs==j)],textsize=16,namesize=6)
    pl[[length(pl)+1]] <- pp
    print(j)
  }}

##########################################################################################
# GENERATE INEQUALITY PLOTS
##########################################################################################

load(file="data/output/projectionOutput/popProjections.Rdata")
load(file="data/output/projectionOutput/growthProjections.Rdata")
yrs <- 2010:2099

popproj <- popProjections[[6]]
basegdp <- popproj$gdpCap
ee <- ecdf(basegdp)
pct <- ee(basegdp)  #so this gives percentile of each country in terms of base gdp
qnt <- rep(1,length(pct))  #this will be a vector indicating baseline income quintile
for (i in 2:5) {
  n <- seq(0,0.8,0.2)[i]
  qnt[pct>=n & pct< (n+0.2)] <- i
}

scen=6 # SSP5
for (j in runs) {
  load("data/output/projectionOutput/GDPcapCC_"%&%j%&%"_"%&%scens[scen]%&%".Rdata")
  load("data/output/projectionOutput/GDPcapNoCC_"%&%j%&%"_"%&%scens[scen]%&%".Rdata")
  
  chgs <- array(dim=c(5,length(yrs)))  #now loop over quintiles and calculate population-weighted impacts by year for each quintile
  for (b in 1:5) {  
    zz <- (qnt==b)
    wt = popproj[zz,which(names(popproj)%in%yrs)]  #population weights 
    chgs[b,] <- diag(t(wt)%*%GDPcapCC[zz,,1])/diag(t(wt)%*%GDPcapNoCC[zz,,1]) - 1  #quicker way to get the weighted averages for each year in the region without looping over years
  }
  chg <- data.frame(t(chgs)*100,yrs)
  df <- melt(chg,id="yrs")
  pl[[length(pl)+1]] <- ggplot() + geom_line(data=df,aes(x=yrs,y=value,group=variable,colour=variable),size=1) + 
    geom_abline(intercept=0,slope=0,size=0.1) + 
    scale_colour_manual(values=rev(c("#fef0d9","#fdcc8a","#fc8d59","#e34a33","#b30000"))) + 
    theme_classic() +
    xlab("years") + ylab("% change in GDP/cap") +
    theme(legend.position="none", panel.grid.minor = element_blank()) +
    scale_y_continuous(limits=c(-90,30),breaks=seq(-100,75,25)) +
    scale_x_continuous(limits=c(2010,2100),breaks=seq(2020,2100,20)) +
    theme(text = element_text(size=16)) +
  annotate("text",label=scens[scen]%&%", "%&%rn[which(runs==j)],x=2010,y=25,size=6,hjust=0)
}

# COMBINE PANELS AND WRITE OUT THE FIGURE
pdf(file="figures/ExtendedDataFigs_Input/Figure5.pdf",width=16,height=16)
layOut(list(pl[[1]],1,1),
       list(pl[[2]],1,2),
       list(pl[[3]],1,3),
       list(pl[[13]],1,4),
       list(pl[[4]],2,1),
       list(pl[[5]],2,2),
       list(pl[[6]],2,3),
       list(pl[[14]],2,4),
       list(pl[[7]],3,1),
       list(pl[[8]],3,2),
       list(pl[[9]],3,3),
       list(pl[[15]],3,4),
       list(pl[[10]],4,1),
       list(pl[[11]],4,2),
       list(pl[[12]],4,3),
       list(pl[[16]],4,4)
)
dev.off()

