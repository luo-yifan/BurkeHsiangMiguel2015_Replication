
#  SCRIPT TO MAKE FIGURE 4 IN BURKE, HSIANG, MIGUEL

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


####################################################################
# TOP PANEL: MAP OF COUNTRY-LEVEL DAMAGES IN 2100 for figure 4
####################################################################

j="pooled"
yrs <- 2010:2099
scens <- c("base","SSP"%&%1:5)
scen=6

load("data/output/projectionOutput/GDPcapCC_"%&%j%&%"_"%&%scens[scen]%&%".Rdata")
load("data/output/projectionOutput/GDPcapNoCC_"%&%j%&%"_"%&%scens[scen]%&%".Rdata")
proj <- (GDPcapCC[,90,1]/GDPcapNoCC[,90,1] - 1)*100
proj[proj>100] = 100
rg <- range(proj,na.rm=T)

library(sf)
world <- st_read('data/input/shape/country.shp', quiet = TRUE)
data <- st_drop_geometry(world)

# now merge country codes
mm <- match(data$GMI_CNTRY,dimnames(GDPcapCC)[[1]])  #gives the row number in data[] for each row in the country vector of projections
proj <- proj[mm]

brks = seq(-100,100,1)  #where to set the color breaks
quant = classIntervals(proj,style="fixed",fixedBreaks = brks)
setwhite = 0  #which value you want to be white
col.pal = designer.colors(length(brks),col=c("red","white","lightblue"),x=c(0,which(brks==setwhite)/length(brks),1))
col.plot = findColours(quant,col.pal)

toplot = data[,3]!="Antarctica"  # make plot without Antarctica

pdf(file="figures/MainFigs_Input/Figure4_TopPanel.pdf",width=12,height=6)
par(mar=c(0,0,0,0))
plot(st_geometry(world[toplot,]),col=col.plot[toplot],ylim=c(-80,80))
# now add legend
rb = c(-30,-75,30,-65)
gradient.rect(rb[1],rb[2],rb[3],rb[4],col=smoothColors(col.pal),gradient="x",border="black")
tx <- seq(-100,100,50)
txw = which(brks%in%tx)
zz = seq(rb[1],rb[3],length=length(tx))
text(zz,rb[2]-2,as.character(tx),pos=1,cex=1)
text((rb[1]+ rb[3])/2,rb[2]-13,"% change in GDP/cap",cex=1)
segments(zz,rb[2]-1,zz,rb[2])
dev.off()


####################################################################
#  REGIONAL PROJECTIONS FOR FIGURE 4, SSP5 AND POOLED MODEL. WE GENERATE THESE FOR THE FULL SUITE OF PROJECTIONS BELOW
####################################################################

#load functions we need
weightProj <- dget("script/weightProj.R")
layOut <- dget("script/layOut.R")

j="pooled"
yrs <- 2010:2099
scens <- c("base","SSP"%&%1:5)
scen=6 # SSP5

iso <- read.csv("data/input/iso/ISOcountryCodes.csv")
region <- read.csv("data/input/iso/ISOregionCodes.csv")
iso <- merge(iso,region,by="sub.region.code")
load("data/output/projectionOutput/GDPcapCC_"%&%j%&%"_"%&%scens[scen]%&%".Rdata")
load("data/output/projectionOutput/GDPcapNoCC_"%&%j%&%"_"%&%scens[scen]%&%".Rdata")
load("data/output/projectionOutput/popProjections.Rdata")
popproj <- popProjections[[scen]]
nms <- as.character(iso$alpha.3)
nms[nms=="COD"] <- "ZAR"  #match to our weird codes
nms[nms=="ROU"] <- "ROM"  
sum(nms%in%popproj[,1])
mm <- match(popproj[,1],nms)
rg <- iso[mm,]
rr <- unique(rg$metaRegion)

# get point estimate of impact for each region, and sort order of plots positive -> negative
dmg <- c()
for (i in 1:length(rr)) {
  zz <- which(rg$metaRegion==rr[i])
  wt = popproj[zz,which(names(popproj)%in%yrs)]  #population weights 
  chgs <- array(dim=c(dim(GDPcapCC)[3],dim(GDPcapCC)[2]))
  b=1
  dd <- diag(t(wt)%*%GDPcapCC[zz,,b])/diag(t(wt)%*%GDPcapNoCC[zz,,b]) - 1
  dmg = c(dmg,dd[length(dd)])  #quicker way to get the weighted averages for each year in the region without looping over years
}
sr <- sort(dmg,index.return=T,decreasing=T)


# to build regional impact projections, we need population projections for the relevant SSP/scenario to use as weights. so need to have run the top lines that generate these
pl <- NULL  #list to fill with plots
#for (i in 1:length(rr)) {
for (i in sr$ix) {
  zz <- which(rg$metaRegion==rr[i])
  wt = popproj[zz,which(names(popproj)%in%yrs)]  #population weights 
  chgs <- array(dim=c(dim(GDPcapCC)[3],dim(GDPcapCC)[2]))
  for (b in 1:dim(chgs)[1]) {
    chgs[b,] <- diag(t(wt)%*%GDPcapCC[zz,,b])/diag(t(wt)%*%GDPcapNoCC[zz,,b]) - 1  #quicker way to get the weighted averages for each year in the region without looping over years
  }
  chgs <- chgs*100
  yz = seq(-100,100,0.1)  #vector of impacts over which densities should be evaluated
  pll <- weightProj(chgs,as.character(rr[i]),textsize=16,namesize=6)
  pl[[length(pl)+1]] <- pll + theme(axis.text.x=element_blank(),
                                    axis.text.y=element_blank(),
                                    axis.title.x=element_blank(),
                                    axis.title.y=element_blank())
  print(as.character(rr[i]))
}

# now make panels
png(file="figures/MainFigs_Input/Figure4_BottomPanels_v2.png",width=12,height=11,units="in",res=200)
layOut(list(pl[[1]],1,1),
       list(pl[[2]],1,2),
       list(pl[[3]],1,3),
       list(pl[[4]],2,1),
       list(pl[[5]],2,2),
       list(pl[[6]],2,3),
       list(pl[[7]],3,1),
       list(pl[[8]],3,2),
       list(pl[[9]],3,3)
)
dev.off()



## Statistics on assumed average growth rates in different years in SSP3 and SSP5, for SI
load("data/output/projectionOutput/popProjections.Rdata")
load("data/output/projectionOutput/growthProjections.Rdata")

yy <- which(names(growthProjections[[6]])=="2050")
weighted.mean(growthProjections[[6]][,yy],popProjections[[6]][,yy])  #pop-weighted growth rate in 2050, ssp5
weighted.mean(growthProjections[[4]][,yy],popProjections[[4]][,yy])  # same for ssp3

yy <- which(names(growthProjections[[6]])=="2090")
weighted.mean(growthProjections[[6]][,yy],popProjections[[6]][,yy])  #pop-weighted growth rate in 2090, ssp5
weighted.mean(growthProjections[[4]][,yy],popProjections[[4]][,yy])  # same for ssp3