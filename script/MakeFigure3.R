
#  Make Figure 3

#  Script also computes reported values on the percentage of countries poorer in absolute and relative terms

rm(list = ls())

require(fields)
require(classInt)
require(plotrix)
require(data.table)
require(dplyr)
require(raster)
library(maps)
library(xtable)
library(reshape2)
"%&%"<-function(x,y)paste(x,y,sep="")  #define a function for easy string pasting


scens <- c("base","SSP"%&%1:5)
yrs <- 2010:2099

pdf(file="figures/MainFigs_Input/Figure3.pdf",height=15,width=15)
par(mar=c(5,4,4,4),xaxs="i",yaxs="i",mfrow=c(1,2))
scen=6
load("data/output/projectionOutput/GDPcapCC_pooled_"%&%scens[scen]%&%".Rdata")
load("data/output/projectionOutput/GDPcapNoCC_pooled_"%&%scens[scen]%&%".Rdata")
cx = 1.4 #scaling factor for labels
plot(1,type="n",xlim=c(1920,2100),ylim=c(4.5,14),axes=F,ylab="GDP per capita (1000 $)",xlab="",las=1,cex.lab=cx)
ll <- c(100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000)
axis(2,at=log(ll),labels=ll/1000,las=1,cex.axis=cx)
axis(4,at=log(ll),labels=ll/1000,las=1,cex.axis=cx)
axis(1,at=c(1921,1970,2010,2050,2099),labels=c("2100","2050","today","2050","2100"),cex.axis=cx)
for (i in 1:dim(GDPcapCC)[1])  {
  lines(yrs,log(GDPcapNoCC[i,,1]),col="grey")
  lines(rev(yrs)-89,log(GDPcapCC[i,,1]),col="red")
}
abline(v=2010)

scen=4
load("data/output/projectionOutput/GDPcapCC_pooled_"%&%scens[scen]%&%".Rdata")
load("data/output/projectionOutput/GDPcapNoCC_pooled_"%&%scens[scen]%&%".Rdata")
cx = 1.4 #scaling factor for labels
plot(1,type="n",xlim=c(1920,2100),ylim=c(4.5,14),axes=F,ylab="GDP per capita (1000 $)",xlab="",las=1,cex.lab=cx)
ll <- c(100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000)
axis(2,at=log(ll),labels=ll/1000,las=1,cex.axis=cx)
axis(4,at=log(ll),labels=ll/1000,las=1,cex.axis=cx)
axis(1,at=c(1921,1970,2010,2050,2099),labels=c("2100","2050","today","2050","2100"),cex.axis=cx)
for (i in 1:dim(GDPcapCC)[1])  {
  lines(yrs,log(GDPcapNoCC[i,,1]),col="grey")
  lines(rev(yrs)-89,log(GDPcapCC[i,,1]),col="red")
}
abline(v=2010)
dev.off()

#statistics in paper on number of countries above/below thresholds
for (scen in c(4,6)) {
  load("data/output/projectionOutput/GDPcapCC_pooled_"%&%scens[scen]%&%".Rdata")
  load("data/output/projectionOutput/GDPcapNoCC_pooled_"%&%scens[scen]%&%".Rdata")
  z1 <- round(sum(GDPcapCC[,90,1]<GDPcapCC[,1,1])/dim(GDPcapCC)[1]*100,2)
  print("%countries absolutely poorer, "%&%scens[scen]%&%" = "%&%z1)  #percent of countries that are poorer in 2100 under CC relative to today
  z2 <- round(sum(GDPcapCC[,90,1]<GDPcapNoCC[,90,1])/dim(GDPcapCC)[1]*100,2)
  print("%countries relatively poorer, "%&%scens[scen]%&%" = "%&%z2)  #percent of countries that are poorer in 2100 under CC relative to today
}

