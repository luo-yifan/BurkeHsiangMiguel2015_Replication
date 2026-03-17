
#  Script to make Extended Data Table 3
#   Table give point estimates of damages and percentiles of distribution

rm(list = ls())

library(xtable)
"%&%"<-function(x,y)paste(x,y,sep="")  #define a function for easy string pasting

runs <- c("pooled","pooled5lag","richpoor","richpoor5lag")
scens <- c("base","SSP"%&%1:5)
yrs <- 2010:2099

nr = 9  #number of rows per scenario
runs <- c("pooled","pooled5lag","richpoor","richpoor5lag")
fill <- matrix(nrow=4*nr,ncol=3)
sz = c(1,4,6)  #scenarios to loop over
for (rr in 1:4) {
  j <- runs[rr]
  for (sc in 1:3) {
    load("data/output/projectionOutput/GlobalChanges_"%&%j%&%"_"%&%scens[sz[sc]]%&%".Rdata")
    chgs <- tots[,,1]/tots[,,2]-1
    chgs <- chgs*100
    fill[(rr-1)*nr+1,sc] <- chgs[1,90]
    fill[(rr-1)*nr+2:6,sc] <- quantile(chgs[,90],probs=c(0.05,0.25,0.5,0.75,0.95))
    fill[(rr-1)*nr+7,sc] <- sum(chgs[,90]<0)/length(chgs[,90])*100  #number of obs <0
    fill[(rr-1)*nr+8,sc] <- sum(chgs[,90]<(-10))/length(chgs[,90])*100  #number of obs < -10%
    fill[(rr-1)*nr+9,sc] <- sum(chgs[,90]<(-20))/length(chgs[,90])*100  #number of obs < -20%
  }
}
nms1 <- c("Pooled model",rep("",nr-1),"Pooled with 5 lags",rep("",nr-1), "Rich/poor model",rep("",nr-1),"Rich/poor with 5 lags",rep("",nr-1))
nms2 <- rep(c("point est.","5th","25th","50th","75th","95th","% runs < 0", "% runs < -10", "% runs < -20"),4)
out<-data.frame(nms1,nms2,fill)
names(out) <- c("model","estimate","baseline growth","SSP3","SSP5")
write.csv(out,file="figures/ExtendedDataFigs_Input/Table3.csv")
tbl <- xtable(out)
digits(tbl) <- 0
print(tbl,include.rownames=FALSE,hline.after=nr*0:4)
