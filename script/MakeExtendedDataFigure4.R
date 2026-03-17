
#  Make Extended Data Figure 4

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

yrs <- 2010:2099
scens <- c("base","SSP"%&%1:5)
scen=6 # SSP5

iso <- read.csv("data/input/iso/ISOcountryCodes.csv")
region <- read.csv("data/input/iso/ISOregionCodes.csv")
iso <- merge(iso,region,by="sub.region.code")
load("data/output/projectionOutput/popProjections.Rdata")

runs <- c("pooled","pooled5lag","richpoor","richpoor5lag")
scen=6 # SSP5
popproj <- popProjections[[scen]]
nms <- as.character(iso$alpha.3)
nms[nms=="COD"] <- "ZAR"  #match to our weird codes
nms[nms=="ROU"] <- "ROM"  
mm <- match(popproj$iso,nms)
rg <- iso[mm,]
rr <- unique(rg$metaRegion)
rrs <- c("SAsia","SSA", "Europe", "MENA", "LAmer", "Ocea.", "SEAsia", "NAmer", "CEAsia")

# get point estimate of impact for each region for pooled SSP5, and sort order of plots positive -> negative as in Fig4
dmg <- c()
j="pooled"
load("data/output/projectionOutput/GDPcapCC_"%&%j%&%"_"%&%scens[scen]%&%".Rdata")
load("data/output/projectionOutput/GDPcapNoCC_"%&%j%&%"_"%&%scens[scen]%&%".Rdata")
for (i in 1:length(rr)) {
  zz <- which(rg$metaRegion==rr[i])
  wt = popproj[zz,which(names(popproj)%in%yrs)]  #population weights 
  chgs <- array(dim=c(dim(GDPcapCC)[3],dim(GDPcapCC)[2]))
  b=1
  dd <- diag(t(wt)%*%GDPcapCC[zz,,b])/diag(t(wt)%*%GDPcapNoCC[zz,,b]) - 1
  dmg = c(dmg,dd[length(dd)]) 
}
sr <- sort(dmg,index.return=T,decreasing=T)


pl <- NULL  #list to fill with plots
for (j in runs)  {  #loop over historical models
  load("data/output/projectionOutput/GDPcapCC_"%&%j%&%"_"%&%scens[scen]%&%".Rdata")
  load("data/output/projectionOutput/GDPcapNoCC_"%&%j%&%"_"%&%scens[scen]%&%".Rdata")

  # to build regional impact projections, we need population projections for the relevant SSP/scenario to use as weights. so need to have run the top lines that generate these
  for (i in sr$ix) {
    zz <- which(rg$metaRegion==rr[i])
    wt = popproj[zz,which(names(popproj)%in%yrs)]  #population weights 
    chgs <- array(dim=c(dim(GDPcapCC)[3],dim(GDPcapCC)[2]))
    for (b in 1:dim(chgs)[1]) {
      chgs[b,] <- diag(t(wt)%*%GDPcapCC[zz,,b])/diag(t(wt)%*%GDPcapNoCC[zz,,b]) - 1  #quicker way to get the weighted averages for each year in the region without looping over years
    }
    chgs <- chgs*100
    yz = seq(-100,100,1)  #vector of impacts over which densities should be evaluated
    pl[[length(pl)+1]] <- weightProj(chgs,as.character(rrs[i]),textsize=7,namesize=5)
    print(as.character(rr[i]))
  } #end region loop
} #end response function loop
  

   png(file="figures/ExtendedDataFigs_Input/Fig4.png",width=10,height=18,bg="transparent",pointsize=3,res=200,units="in")
  layOut(list(pl[[1]],1,1),
         list(pl[[2]],2,1),
         list(pl[[3]],3,1),
         list(pl[[4]],4,1),
         list(pl[[5]],5,1),
         list(pl[[6]],6,1),
         list(pl[[7]],7,1),
         list(pl[[8]],8,1),
         list(pl[[9]],9,1),
         list(pl[[10]],1,2),
         list(pl[[11]],2,2),
         list(pl[[12]],3,2),
         list(pl[[13]],4,2),
         list(pl[[14]],5,2),
         list(pl[[15]],6,2),
         list(pl[[16]],7,2),
         list(pl[[17]],8,2),
         list(pl[[18]],9,2),
         list(pl[[19]],1,3),
         list(pl[[20]],2,3),
         list(pl[[21]],3,3),
         list(pl[[22]],4,3),
         list(pl[[23]],5,3),
         list(pl[[24]],6,3),
         list(pl[[25]],7,3),
         list(pl[[26]],8,3),
         list(pl[[27]],9,3),
         list(pl[[28]],1,4),
         list(pl[[29]],2,4),
         list(pl[[30]],3,4),
         list(pl[[31]],4,4),
         list(pl[[32]],5,4),
         list(pl[[33]],6,4),
         list(pl[[34]],7,4),
         list(pl[[35]],8,4),
         list(pl[[36]],9,4)
  )
  dev.off()

