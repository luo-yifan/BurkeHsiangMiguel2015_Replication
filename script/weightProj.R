
# Define function that makes visually weighted projections.

# This takes the output from the historical bootstrap, calculates the density of impacts for each year across the bootstrap runs, and makes an image plot that shows
#     how these densities evolve over time. 
# Function requires two arguments:  the N x K matrix of bootstrap projections (N=#bootstraps and K = years), and a name for the country/region for plotting purposes that get's plotted in upper left, 
#     User can also specify a low and high CI that the plotting can be restricted to (default is 5th-95th percentiles), and some sizes for the text labels

weightProj <- function(mat,nm,lo=0.05,hi=0.95,textsize=10,namesize=5,ylo=-100,yhi=75, annotate_y=70) {  #two required arguments
  # first construct dataframe of weights. putting it in dataframe for ggplot()
  gd <- c()
  for (i in 1:length(yrs)) {
    qt<-quantile(mat[,i],probs=c(lo,hi))  
    keep <- mat[,i]>=qt[1] & mat[,i]<=qt[2]  
    z <- density(mat[keep,i],from=-100,to=100,n=length(yz),adjust=0.6)
    zz <- (z$y - min(z$y))/max(z$y)
    yr <- rep(yrs[i],length(yz))
    add = data.frame(yr,yz,zz)
    gd <- rbind(gd,add)
  }
  yy <- mat[1,]  #use first row as point estimate
  mn <- data.frame(yrs,yy)
  
  # now make plot 
  p <- ggplot(gd, aes(x=yr,y=yz)) + geom_tile(aes(fill=zz)) + scale_fill_gradientn(colours=c("white","red")) + 
    geom_line(data=mn,aes(x=yrs,y=yy),size=1) +
    geom_abline(intercept=0,slope=0,size=0.1) + 
    xlab("year") + ylab("% change in GDP/cap") +
    theme_classic()  +  
    theme(legend.position="none") +
    scale_y_continuous(limits=c(ylo,yhi),breaks=seq(-100,75,25)) +
    scale_x_continuous(limits=c(2010,2100),breaks=seq(2020,2100,20)) +
    theme(text = element_text(size=textsize)) + 
    annotate("text",x=2010,y=annotate_y,label=nm,size=namesize,hjust=0)
  return(p)
}