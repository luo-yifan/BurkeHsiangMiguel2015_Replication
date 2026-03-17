
#  Script to create panels of Extended Data Figure 1

rm(list = ls())
"%&%"<-function(x,y)paste(x,y,sep="")  #define a function for easy string pasting

##################################################################
#   FIGURE TO SHOW HOW INDIVIDUAL COUNTRIES MAP INTO GLOBAL RESPONSE FUNCTION
##################################################################

resp <- read.csv("data/output/estimatedGlobalResponse.csv")
dta <- read.csv("data/output/mainDataset.csv")
smpl <- is.na(dta$growthWDI)==F & is.na(dta$UDel_temp_popweight)==F   #main estimation sample
coef <- read.csv("data/output/estimatedCoefficients.csv")

# center response at 20C
x = resp$x
at20 = resp$estimate[x==20]
est = resp$estimate - at20
min90 = resp$min90 - at20
max90 = resp$max90 - at20

ctys = c('RUS','FRA','USA','VNM','MLI')

#initialize plot
pdf(file="figures/ExtendedDataFigs_Input/Fig1_PanelA.pdf",width=5,height=4)

plot(1,xlim=c(-2,30),ylim=c(-0.25,0.1),type="n",las=1,xlab="temperature",ylab="growth rate")
abline(h=0,lwd=1)

# add vertical average temperature lines and ranges for selected countries
bot = -0.25
top = 0.1
for (j in 1:length(ctys)) {
  tt = mean(dta$UDel_temp_popweight[dta$iso==ctys[j]],na.rm=T)
  rg = range(dta$UDel_temp_popweight[dta$iso==ctys[j]],na.rm=T)
  rect(rg[1],bot,rg[2],top,col="grey",border=NA)
  #  segments(tt,bot,tt,top,lwd=0.5)
}

# plot CI and main effect
polygon(c(x,rev(x)),c(min90,rev(max90)),col="lightblue",border=NA)
lines(x,est,lwd=2)
dev.off()



##################################################################
# 	MAP OF MARGINAL EFFECTS OF WARMING
##################################################################

resp <- read.csv("data/output/estimatedGlobalResponse.csv")
dta <- read.csv("data/output/mainDataset.csv")
coef <- read.csv("data/output/estimatedCoefficients.csv")

library(dplyr)
library(sf)
dta1 <- dta  # dta is already a data.frame; tbl_dt() is defunct
mt <- dta1 %>%   #the following few lines gets the average temperature in each country for the proscribed years
  filter(year>=1980) %>%
  group_by(iso) %>%
  summarize(meantemp = mean(UDel_temp_popweight,na.rm=T))
mt <- as.data.frame(mt)

world <- st_read('data/input/shape/country.shp', quiet = TRUE)
data <- st_drop_geometry(world)

# now merge country codes
mm <- match(data$GMI_CNTRY,mt$iso)  #gives the row number in data[] for each row in mt

# calculate projected change in growth rate for 1C increase in temperature
tt = mt$meantemp[mm]
# base = coef$b1*tt + coef$b2*tt^2
# t1 = coef$b1*(tt+1) + coef$b2*(tt+1)^2
# proj = (t1 - base)*100
proj = (coef$b1 + 2*coef$b2*tt)*100
rg <- range(proj,na.rm=T)

brks = as.numeric(round(seq(round(rg[1]-0.1,1),round(rg[2]+0.1,1),0.1),1))  #where to set the color breaks
quant = classIntervals(proj,style="fixed",fixedBreaks = brks)
setwhite = 0  #which value you want to be white
col.pal = designer.colors(length(brks),col=c("red","white","lightblue"),x=c(0,which(brks==setwhite)/length(brks),1))
col.plot = findColours(quant,col.pal)

toplot = data[,3]!="Antarctica"
pdf(file="figures/ExtendedDataFigs_Input/Fig1_PanelG.pdf",width=12,height=6)
par(mar=c(0,0,0,0))
plot(st_geometry(world[toplot,]),col=col.plot[toplot],ylim=c(-80,80))  #subset sf object to exclude Antarctica
# now add legend
rb = c(-30,-75,30,-65)
gradient.rect(rb[1],rb[2],rb[3],rb[4],col=smoothColors(col.pal),gradient="x",border="black")
tx = unique(round(brks))
tx = tx[tx>rg[1] & tx<rg[2]]
txw = which(brks%in%tx)
zz = seq(rb[1],rb[3],length=length(brks))
text(zz[txw],rb[2]-2,as.character(tx),pos=1,cex=1)
text((rb[1]+ rb[3])/2,rb[2]-13,"ppt effect on growth rate",cex=1)
segments(zz[txw],rb[2]-1,zz[txw],rb[2])
dev.off()


####################################################################################################################################
## PLOT PANELS H & I & J & K:  MARGINAL EFFECTS FROM TIME SERIES REGRESSIONS + Interacted model + more flexible functional forms
####################################################################################################################################

pdf(file="figures/ExtendedDataFigs_Input/Fig1_PanelHIJK.pdf",width=10,height=8,useDingbats = F)
par(mfrow=c(2,2))

marg <- read.csv("data/output/ExtendedDataFig1g.csv")
cilo = marg$b - 2*marg$se
cihi = marg$b + 2*marg$se
mt <- marg$meantemp
out <- marg$fit<cilo | marg$fit>cihi  #countries where CI does not overlap predicted marginal effect from quadratic model
sum(out)  #count these countries

cty <- c("USA","VNM","MLI","ISL","FRA")
ctys <- which(marg$iso%in%cty)

plot(marg$meantemp,marg$b,pch=21,bg="grey",las=1,xlab="mean temp", ylab = "marginal effect (dy/dT)",cex=0.8)
colz=rep("grey",length(cilo))
segments(mt,cilo,mt,cihi,lwd=0.5,col=colz)
pcol <- rep("grey",length(cilo))
pcol[ctys] <- "orange"
points(marg$meantemp,marg$b,pch=21,bg=pcol,cex=1)
abline(a=0.013456, b=2*-0.0005026,col="black",lwd=2)  #first derivative of global estimated response function
text(mt[ctys],marg$b[ctys],marg$iso[ctys],pos=4,cex=0.6,offset=0.3)


#  PLOT INTERACTED MODEL, PANEL H
eff <- read.csv("data/output/ExtendedDataFig1h.csv")
#pdf(file="figures/ExtendedDataFigs_Input/Fig1_PanelH.pdf",width=6,height=5,useDingbats = F)
plot(1,type="n",xlab="mean temp", ylab = "marginal effect (dy/dT)", xlim= c(1,30), ylim=c(-0.025,0.025), las=1)
abline(h=0)
segments(eff$n1,eff$cilo1,eff$n1,eff$cihi1,lwd=2,col="blue")
segments(eff$n2,eff$cilo2,eff$n2,eff$cihi2,lwd=2,col="orange")
points(eff$n1,eff$marg1,pch=21,bg="blue",cex=1)
points(eff$n2,eff$marg2,pch=21,bg="orange",cex=1)


# PANEL I, J LOOKING AT MORE FLEXIBLE FUNCITONAL FORMS
fnc <- read.csv("data/output/ExtendedDataFig1i.csv")
fnc <- fnc[fnc$X<=30 & is.na(fnc$X)==F,]
colz = c("black","red","orange","black","red","orange")
ltyz = c(1,1,1,2,2,2)
#polynomials
plot(1,type="n",xlab="mean temp", ylab = "growth rate", xlim= c(1,30), ylim=c(-0.15,0.07), las=1)
for (i in 2:7) {
  ww <- which(names(fnc)=="Y"%&%i)
  yy <- fnc[,ww]
  y <- yy - yy[fnc$X==20]  #plot relative to year at 20C
  lines(fnc$X,y,lwd=2,col=colz[i-1],lty=ltyz[i-1])  
}

#splines
plot(1,type="n",xlab="mean temp", ylab = "growth rate", xlim= c(1,30), ylim=c(-0.15,0.07), las=1)
for (i in 3:7) {
  ww <- which(names(fnc)=="spline_est_"%&%i)
  yy <- fnc[,ww]
  y <- yy - yy[fnc$X==20]  #plot relative to year at 20C
  lines(fnc$X,y,lwd=2,col=colz[i-1],lty=ltyz[i-1])  
}
# add quadratic for comparison
ww <- which(names(fnc)=="Y"%&%2)
yy <- fnc[,ww]
y <- yy - yy[fnc$X==20]  #plot relative to year at 20C
lines(fnc$X,y,lwd=2,col="black")  

dev.off()
