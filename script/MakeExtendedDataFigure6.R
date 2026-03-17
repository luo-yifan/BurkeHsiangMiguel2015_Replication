
# Make Extended Data Figure 6

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


runs <- c("pooled","pooled5lag","richpoor","richpoor5lag") #historical models
scens <- c("base","SSP"%&%1:5)  #socioeconomic scenarios
yrs <- 2010:2099

iam <- read.csv("data/input/IAMdata/ProcessedKoppData.csv")
mods <- c("DICE","FUND","PAGE")
iam <- iam[iam$IAM%in%mods,]
inc <- sort(unique(iam$T[iam$T<=6 & iam$T>1]))  #temperatures under which some IAM was run
incs <- c(0.8,inc)  #adding in a zero increase relative to today (=0.8 increase relative to preindustrial)

pdf(file="figures/ExtendedDataFigs_Input/Fig6.pdf",height=6,width=8)
layout(matrix(c(1,1,1,2,1,1,1,3,1,1,1,4),byrow=T,nrow=3))
par(mar=c(5,4,1,1))


#####################################################
#  Panel A:  
#####################################################
colz = c("black","red","orange","blue")
ltyz = c("dashed","dotted","solid")
plot(1,type="n",xlim=c(0.8,5),ylim=c(-90,50),xlab="temperature increase by 2100 (C)", ylab="% change in global GDP", las=1,cex.lab=2)
abline(h=0)
for (j in 1:length(runs))  {
  load("data/output/projectionOutput/DamageFunction_"%&%runs[j]%&%".Rdata")
  cc <- (tots[,,1]/tots[,,2] - 1)*100
  for (i in 1:3) 
    lines(incs,cc[,i],col=colz[j],lty=ltyz[i],lwd=2)
  }
  

#####################################################
#  Ratio panels  
#####################################################

par(mar=c(3,3,1,1))

for (mod in mods) {
plot(1,type="n",xlim=c(0.8,5),ylim=c(0,1),xlab="temp. increase by 2100 (C)", ylab="IAM estimate/our estimate",las=1)
abline(h=seq(0.2,0.8,0.2),lty=2,lwd=0.8,col="grey")
abline(h=c(0,1))
for (j in 1:length(runs))  {
  load("data/output/projectionOutput/DamageFunction_"%&%runs[j]%&%".Rdata")
  cc <- (tots[,,1]/tots[,,2] - 1)*100
  cc <- data.frame(cc,incs)
  toplot <- merge(cc,iam[iam$IAM==mod,],by.x="incs",by.y="T")
  for (i in 1:3) {
    rr <- toplot$damages*(-1)/toplot[,i+1]
    lines(toplot$incs,rr,col=colz[j],lty=ltyz[i],lwd=1)
}}
}

dev.off()