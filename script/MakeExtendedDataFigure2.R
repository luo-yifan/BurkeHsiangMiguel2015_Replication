

#   SCRIPT TO MAKE EXTENDED DATA FIGURE 2

"%&%" <- function(x, y) paste(x, y, sep = "")  # string concatenation helper

pdf(file="figures/ExtendedDataFigs_Input/Fig2_Panel.pdf",width=12,height=8,useDingbats=F)
matr <- cbind(matrix(1:9,byrow=T,nrow=3),matrix(10:18,byrow=T,nrow=3))
layout(matr)
par(mar=c(4,4,3,1))

rp <- read.csv("data/output/Effect_Marginals_RichPoor.csv")
x=1:30
zz = c("rich","poor","rich - poor")
for (v in c("growthWDI","AgrGDPgrowthCap","NonAgrGDPgrowthCap")) {
  for (j in c("r","p","c")) {
    b = rp[,which(names(rp)==v%&%"_b"%&%j)]
    se = rp[,which(names(rp)==v%&%"_se"%&%j)]
    cilo = b - 1.96*se
    cihi = b + 1.96*se
    plot(1,ylim=c(-0.04,0.04),xlim=c(1,30),cex.axis=1.2,las=1,xlab="temperature",ylab=v,main=zz[which(c("r","p","c")==j)]) 
    abline(h=0)
    abline(v=c(10,20,30),lty=2,lwd=0.5)
    polygon(c(x,rev(x)),c(cilo,rev(cihi)),col="lightblue",border=NA)
    lines(x,b,lwd=2)
  }
}


#plot of p-values
zz = c("rich","poor","rich - poor")
for (v in c("growthWDI","AgrGDPgrowthCap","NonAgrGDPgrowthCap")) {
  for (j in c("r","p","c")) {
    z = which(names(rp)==v%&%"_p"%&%j)
    plot(1:30,rp[,z],ylim=c(0,1),cex.axis=1.2,las=1,xlab="temperature",ylab="p-value",main=zz[which(c("r","p","c")==j)],pch=19,yaxs='i') 
    abline(h=seq(0,1,0.2),lty=2,lwd=0.5)
    abline(v=c(10,20,30),lty=2,lwd=0.5)
    abline(h=c(0.05,0.1),col="red")
  }
}
dev.off()

