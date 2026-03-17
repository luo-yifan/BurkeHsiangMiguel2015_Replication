
# Code to construct Figure 5 in Burke, Hsiang, Miguel 2015

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

pl <- NULL  #list that we will fill with plots
runs <- c("pooled","pooled5lag","richpoor","richpoor5lag") #historical models
scens <- c("base","SSP"%&%1:5)  #socioeconomic scenarios
yrs <- 2010:2099

#load functions we need
weightProj <- dget("script/weightProj.R")
layOut <- dget("script/layOut.R")

#initialize list to hold plot panels
pl <- NULL  


############################################################################
# Panel A: global projections for SSP5, pooled model
############################################################################

scen=6 # SSP5
j <- "pooled"
load("data/output/projectionOutput/GlobalChanges_"%&%j%&%"_"%&%scens[scen]%&%".Rdata")
chgs <- tots[,,1]/tots[,,2]-1
chgs <- chgs*100
chgs[,1] <- 0  # first year is end of baseline, so no change
yz = seq(-100,100,0.1)  #vector of impacts over which densities should be evaluated
pl[[1]] <- weightProj(chgs,"",ylo=-75,yhi=60,textsize=18,annotate_y = 50)


############################################################################
# Panel B: global projections for SSP5, all historical models 
############################################################################

scen=6
imp <- NULL
for (j in runs)  {
    load("data/output/projectionOutput/GlobalChanges_"%&%j%&%"_"%&%scens[scen]%&%".Rdata")
    chgs <- tots[,,1]/tots[,,2]-1
    chgs <- chgs*100
    chgs[,1] <- 0  # first year is end of baseline, so no change
    imp <- cbind(imp,chgs[1,])
  }
colnames(imp) <- runs
imp1 <- data.frame(yrs,imp)
df <- melt(imp1,id="yrs")
p <- ggplot() + geom_line(data=df,aes(x=yrs,y=value,group=variable,colour=variable),size=1) + 
  geom_abline(intercept=0,slope=0,size=0.1) + 
  scale_colour_manual(values=c("black","red","orange","blue")) + 
  theme_classic()  +  
  xlab("years") + ylab("% change in GDP/cap") +
  theme(panel.background=element_rect(fill="white",color="black"), legend.position="none", panel.grid.minor = element_blank()) +
  scale_y_continuous(limits=c(-75,60),breaks=seq(-100,75,25)) +
  scale_x_continuous(limits=c(2010,2100),breaks=seq(2020,2100,20)) +
  theme(text = element_text(size=18))
pl[[2]] <- p


  
  
############################################################################
# Panel C: damages by quintile for SSP5, pooled model
############################################################################

scen=6 # SSP5
j="pooled"
load("data/output/projectionOutput/GDPcapCC_"%&%j%&%"_"%&%scens[scen]%&%".Rdata")
load("data/output/projectionOutput/GDPcapNoCC_"%&%j%&%"_"%&%scens[scen]%&%".Rdata")
load("data/output/projectionOutput/popProjections.Rdata")
load("data/output/projectionOutput/growthProjections.Rdata")
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
chgs <- array(dim=c(5,length(yrs)))  #now loop over quintiles and calculate population-weighted impacts by year for each quintile
for (b in 1:5) {  
  zz <- (qnt==b)
  wt = popproj[zz,which(names(popproj)%in%yrs)]  #population weights 
  chgs[b,] <- diag(t(wt)%*%GDPcapCC[zz,,1])/diag(t(wt)%*%GDPcapNoCC[zz,,1]) - 1  #quicker way to get the weighted averages for each year in the region without looping over years
}
chg <- data.frame(t(chgs)*100,yrs)
df <- melt(chg,id="yrs")
pl[[3]] <- ggplot() + geom_line(data=df,aes(x=yrs,y=value,group=variable,colour=variable),size=1) + 
  geom_abline(intercept=0,slope=0,size=0.1) + 
  scale_colour_manual(values=rev(c("#fef0d9","#fdcc8a","#fc8d59","#e34a33","#b30000"))) + 
  theme_classic()  +  
  xlab("years") + ylab("% change in GDP/cap") +
  theme(panel.background=element_rect(fill="white",color="black"), legend.position="none", panel.grid.minor = element_blank()) +
  scale_y_continuous(limits=c(-90,30),breaks=seq(-100,75,25)) +
  scale_x_continuous(limits=c(2010,2100),breaks=seq(2020,2100,20)) +
  theme(text = element_text(size=18))


############################################################################
# Panel D: global damage function for SSP5, all historical models + uncertainty on pooled model
############################################################################

iam <- read.csv("data/input/IAMdata/ProcessedKoppData.csv")
iam$damages <- -1*iam$damages  #flipping signs to match everything else we have
mods <- c("DICE","FUND","PAGE")
iam1 <- iam[iam$IAM%in%mods,]

inc <- sort(unique(iam1$T[iam1$T<=6 & iam1$T>1]))  #temperatures under which some IAM was run
ll <- c(0.8,inc)  #adding in a zero increase relative to today (=0.8 increase relative to preindustrial)
scen = 6
chgs <- array(dim=c(1001,length(ll)))
for (dt in 1:length(ll)) {
  load("data/output/projectionOutput/DamageFunction_pooled_"%&%scens[scen]%&%"_"%&%ll[dt]%&%".Rdata")
  chgs[,dt] <- (tots[,90,1]/tots[1,90,2] - 1)*100  
}
cihi=est=cilo=ci1=ci2=c()
for (j in 1:length(ll)) {
  zz <- quantile(chgs[,j],probs=c(0.05,0.25,0.75,0.95))
  cilo = c(cilo,zz[1])
  cihi = c(cihi,zz[4])
  ci1  = c(ci1,zz[2]) #bottom of IQR
  ci2 = c(ci2,zz[3]) #top of IQR
  est = c(est,chgs[1,j])
}
dta1 <- data.frame(id=rep(1,length(ll)*2),x=c(ll,rev(ll)),y=c(cilo,rev(cihi)))
dta2 <- data.frame(id=rep(2,length(ll)*2),x=c(ll,rev(ll)),y=c(ci1,rev(ci2)))
dta <- rbind(dta1,dta2)  #our CI
est <- data.frame(est,ll)

#bring in estimates using all historical models
runs <- c("pooled","pooled5lag","richpoor","richpoor5lag") #historical models
imp <- NULL
for (j in runs)  {
  load("data/output/projectionOutput/DamageFunction_"%&%j%&%".Rdata")
  chg <- (tots[,,1]/tots[,,2] - 1)*100
  chg <- data.frame(ll,chg[,which(colnames(chg)=="SSP5")],rep(j,length(ll)))
  imp <- rbind(imp,chg)
}

names(imp) <- c("T","damages","model")
# #write this out because people are going to want it
# names(imp1)[1] <- "degC_abovePreIndustrial"
# write.csv(imp1,file="data/output/DamageFunctionPointEstimates.csv",row.names=F)

# now make plot
pl[[4]] <- ggplot(dta, aes(x=x,y=y)) + geom_polygon(aes(fill=id,group=id)) +
  geom_abline(intercept=0,slope=0,size=0.1) + 
  xlab("temperature change (C)") + ylab("% change in global GDP") +
  theme_classic()  +  
  scale_fill_gradient(low="lightblue",high="skyblue") + 
  theme(panel.background=element_rect(fill="white",color="black"), legend.position="none", panel.grid.minor = element_blank()) +
  theme(text = element_text(size=18)) +
  coord_cartesian(ylim=c(-90,50),xlim=c(0.8,5)) +
  scale_y_continuous(breaks=seq(-75,50,25)) +  
  geom_line(data=imp,aes(x=T,y=damages,group=model,color=model),size=1) +
  geom_line(data=iam1,aes(x=T,y=damages,group=IAM,colour=IAM),size=1,linetype="dashed") +
  scale_colour_manual(values=c("violetred1","purple","brown","black","red","orange","blue")) 


# NOW WRITE COMBINED FIGURE OUT
pdf(file="figures/MainFigs_Input/Figure5.pdf",height=10,width=10)
layOut(list(pl[[1]],1,1),
       list(pl[[2]],1,2),
       list(pl[[3]],2,1),
       list(pl[[4]],2,2)
)
dev.off()


