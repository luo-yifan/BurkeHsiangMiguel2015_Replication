
# SCRIPT TO MAKE FIGURE 2
# This file calls data created in GenerateFigure2Data.do

rm(list = ls())

require(fields)
require(classInt)
require(plotrix)
require(dplyr)
"%&%"<-function(x,y)paste(x,y,sep="")  #define a function for easy string pasting

pdf(file="figures/MainFigs_Input/Figure2.pdf",width=10,height=5.5,useDingbats=F)

mat <- matrix(c(1,1,2,3,1,1,4,5),nrow=2,byrow=T)
layout(mat)
ax = 1.5  #scaling for axes
par(mar=c(4,4,2,1))

########################################################
#  Panel A
########################################################

resp <- read.csv("data/output/estimatedGlobalResponse.csv")
dta <- read.csv("data/output/mainDataset.csv")
smpl <- is.na(dta$growthWDI)==F & is.na(dta$UDel_temp_popweight)==F   #main estimation sample
coef <- read.csv("data/output/estimatedCoefficients.csv")

# center response at optimum
x = resp$x
mx = max(resp$estimate)
est = resp$estimate - mx
min90 = resp$min90 - mx
max90 = resp$max90 - mx

ctys = c('USA','CHN','DEU','JPN','IND','NGA','IDN','BRA','FRA','GBR')
ctt = c('US','CHN',"GER","JPN",'IND','NGA','INDO','BRA','FRA','UK')

#initialize plot
plot(1,xlim=c(-2,30),ylim=c(-0.4,0.1),type="n",las=1,cex.axis=1.3)

# add vertical average temperature lines for selected countries
for (j in 1:length(ctys)) {
  tt = mean(dta$UDel_temp_popweight[dta$iso==ctys[j]],na.rm=T)
  segments(tt,-0.23,tt,0.15,lwd=0.5)
}

# plot CI and main effect
polygon(c(x,rev(x)),c(min90,rev(max90)),col="lightblue",border=NA)
lines(x,est,lwd=2)

# now add histograms at bottom
# first calculate percent of population and global gdp produced at each temperature bin, for our estimation sample
bins = seq(-7,30,0.5)
histtemp = dta$UDel_temp_popweight[smpl]
histpop = dta$Pop[smpl]
histgdp = dta$TotGDP[smpl]
pop = gdp = c()
for (j in 1:(length(bins)-1)) {
  lo = bins[j]
  hi = bins[j+1]
  pop = c(pop,sum(histpop[histtemp>=lo & histtemp<hi]))
  gdp = c(gdp,sum(histgdp[histtemp>=lo & histtemp<hi]))
}
pop = pop/sum(pop)
gdp = gdp/sum(gdp)

#parameters that set where histograms go
dis = 0.055
base = -0.3

# now make histograms
#temperature
zz <- hist(histtemp,plot=F,breaks=bins)
cts = zz$counts/max(zz$counts)*0.05  #sets the height of the tallest bar to 0.05
rect(bins,base,bins+0.5,base+cts,col="red")
# pop
cts = pop/max(pop)*0.05
rect(bins,base-dis*(1),bins+0.5,base-dis*(1)+cts,col="grey")
# gdp
cts = gdp/max(gdp)*0.05
rect(bins,base-dis*(2),bins+0.5,base-dis*(2)+cts,col="black")



########################################################
#  Panels b
########################################################
resp <- read.csv("data/output/EffectHeterogeneity.csv")
poor <- dta$GDPpctile_WDIppp<50
rich <- dta$GDPpctile_WDIppp>=50

resp <- resp[resp$x>=5,]  #dropping estimates below 5C, since so little poor country exposure down there
mods = unique(as.character(resp$model))

m <- "growthWDI"
plot(1,xlim=c(5,30),ylim=c(-0.35,0.1),type="n",las=1,cex.axis=1.3,cex.lab=1.3,ylab="",xlab="")
smp = resp$model==m & resp$interact==1  #poor countries
xx = resp$x[smp]
mx = max(resp$estimate[smp])
est = resp$estimate[smp] - mx
min90 = resp$min90[smp] - mx
max90 = resp$max90[smp] - mx

polygon(c(xx,rev(xx)),c(min90,rev(max90)),col="lightblue",border=NA)
lines(xx,est,lwd=2,col="steelblue3")

# now add rich countries
smp = resp$model==m & resp$interact==0  #rich countries
xx = resp$x[smp]
mx = max(resp$estimate[smp])
est = resp$estimate[smp] - mx
lines(xx,est,lwd=2,col="red")

# now add histograms of temperature exposures at the base
bins = seq(-7,30,0.5)
poortemp = dta$UDel_temp_popweight[smpl==T & poor==T]
richtemp = dta$UDel_temp_popweight[smpl==T & rich==T]
base = -0.3
zz <- hist(richtemp,plot=F,breaks=bins)
cts = zz$counts/max(zz$counts)*0.05  #sets the height of the tallest bar to 0.05
rect(bins,base,bins+0.5,base+cts,border="red",col="NA")
base = -0.35
zz1 <- hist(poortemp,plot=F,breaks=bins)
cts = zz1$counts/max(zz1$counts)*0.05
rect(bins,base,bins+0.5,base+cts,col="lightblue")


########################################################
#  Panel c
########################################################
resp <- read.csv("data/output/EffectHeterogeneityOverTime.csv")
early <- dta$year<1990

smp = resp$interact==1  #early period
xx = resp$x[smp]
mx = max(resp$estimate[smp])
est = resp$estimate[smp] - mx	
min90 = resp$min90[smp] - mx
max90 = resp$max90[smp] - mx

plot(1,xlim=c(5,30),ylim=c(-0.35,0.1),type="n",las=1,cex.axis=1.3,cex.lab=1.3,ylab="",xlab="")
polygon(c(xx,rev(xx)),c(min90,rev(max90)),col="lightblue",border=NA)
lines(xx,est,lwd=2,col="steelblue3")
# now add point estimate for later period
smp = resp$interact==0  #poor countries
xx = resp$x[smp]
mx = max(resp$estimate[smp])
est = resp$estimate[smp] - mx  
lines(xx,est,lwd=2,col="red")

# now add histograms of temperature exposures at the base
bins = seq(-7,30,0.5)
earlytemp = dta$UDel_temp_popweight[smpl==T & early==T]
latetemp = dta$UDel_temp_popweight[smpl==T & early==F]
base = -0.3
zz <- hist(earlytemp,plot=F,breaks=bins)
cts = zz$counts/max(zz$counts)*0.05  #sets the height of the tallest bar to 0.05
rect(bins,base,bins+0.5,base+cts,border="red",col=NA)
base = -0.35
zz1 <- hist(latetemp,plot=F,breaks=bins)
cts = zz1$counts/max(zz1$counts)*0.05
rect(bins,base,bins+0.5,base+cts,col="lightblue")


########################################################
#  Panels d, e
########################################################

resp <- read.csv("data/output/EffectHeterogeneity.csv")
poor <- dta$GDPpctile_WDIppp<50
rich <- dta$GDPpctile_WDIppp>=50
resp <- resp[resp$x>=5,]  #dropping estimates below 5C, because so little poor country exposure there
mods = unique(as.character(resp$model))
toplot=c("AgrGDPgrowthCap","NonAgrGDPgrowthCap")

for (m in toplot) {
  plot(1,xlim=c(5,30),ylim=c(-0.35,0.1),type="n",las=1,cex.axis=1.3,cex.lab=1.3,ylab="",xlab="")
  smp = resp$model==m & resp$interact==1  #poor countries
  xx = resp$x[smp]
  mx = max(resp$estimate[smp])
  est = resp$estimate[smp] - mx  
  min90 = resp$min90[smp] - mx
  max90 = resp$max90[smp] - mx
  
  polygon(c(xx,rev(xx)),c(min90,rev(max90)),col="lightblue",border=NA)
  lines(xx,est,lwd=2,col="steelblue3")
  # now add rich countries
  smp = resp$model==m & resp$interact==0  #poor countries
  xx = resp$x[smp]
  mx = max(resp$estimate[smp])
  est = resp$estimate[smp] - mx  
  lines(xx,est,lwd=2,col="red")
  
}

dev.off()

