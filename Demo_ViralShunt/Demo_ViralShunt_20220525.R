### SUPPLEMENTARYINFORMATION_RCODE - Full Code
### Last Updated May 25, 2022
rm(list = ls())
library(lmodel2)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(splines)
library(fields)

# Set working directory
Path_root='D:/data/SEATS'
setwd(paste(Path_root,"/Demo_ViralShunt",sep=''))
source('ViralShunt_Xcor.R') # Functions for conducting cross-correlation analysis and associated statistical fittings
seed.t=21591
set.seed(seed.t)

# Loading high-frequency microbial data collected from South-East Asia Time Series station (i.e. SEATS; 18 oN, 116 oE)
# The data users should contact the database administrator (fkshiah@gate.sinica.edu.tw) and strictly follow the terms of use below
# Terms of Use
##### I. You must acknowledge the use of content.
##### II. Monitoring data is made available for use in activities of a non-profit nature only.
##### III. Users must contact the SEATS database administrator (fkshiah@gate.sinica.edu.tw) before using monitoring data for any publications, 
#####      including conference presentations as well as handouts and presentation materials for meetings such as committees and councils. 
#####      Co-authorship may be required for some publications depending on the manner in which the monitoring data is to be used.

plot.dat <- read.csv("SEATS_MicrobialData.csv") # read in plot-level data with raw vars
date.name=as.character(plot.dat[,'Date_first'])
site.name=unique(date.name)
n.cruise=length(site.name)
site.name2=paste(1:n.cruise,
                        month.abb[as.numeric(unlist(lapply(site.name,function(x){strsplit(x,'/')[[1]][2]})))],
                        as.numeric(unlist(lapply(site.name,function(x){strsplit(x,'/')[[1]][1]}))),sep='. ')
color.all=c('blue3','green2','red3',"olivedrab","cyan1","gray50","orange3","turquoise4","violet")
date.pty=color.all[plot.dat[,'cruise_number']]


win.graph(20, 10); par(mfrow=c(1,2))
# Figure 1A
x='Virus';y='BB';
plot(log(plot.dat[,y]) ~ log(plot.dat[,x]),ylim=range(log(plot.dat[,y])[!is.na(plot.dat[,x])],na.rm=T),
     cex.main=1.3, font.main=1, pch=21, col="black", main="", bg=date.pty, lwd=1.5, cex.lab=1.4, cex=1.2, cex.axis=1.2,
     xlab="log Virus abundance", ylab="log Bacterial biomass")

lm.t=apply(lmodel2(log(plot.dat[,y]) ~ log(plot.dat[,x]),  nperm=1000)$regression.results[2,-1],2,as.numeric)
lty1=1;if(lm.t['P-perm (1-tailed)']>0.05){lty1=2}
abline(as.numeric(lm.t[c('Intercept','Slope')]),lwd=2,lty=lty1)
r2t=cor(log(plot.dat[,x]),log(plot.dat[,y]),use='complete.obs')^2
pv=round(lm.t['P-perm (1-tailed)'],2); if(pv<0.01){pv='<0.01'}
nf=sum(!apply(is.na(cbind(plot.dat[,c(x,y)])),1,any))
  
se.slope=sd(log(plot.dat[,y]),na.rm=T)/sd(log(plot.dat[,x]),na.rm=T)*((1-r2t)/nf)^0.5

legend(-1.6,1.5,
       paste('Slope =',round(lm.t[c('Slope')],2),'+/-',round(se.slope,2)
             ,'\nR2 =',round(r2t,2)
             ,'\np =',pv),
       bty='n',cex=1.2)
legend('topleft',site.name2,col=color.all,pch=rep(16),bty='n')

# Figure 1B
x='Virus';y='SGR'
plot(log(plot.dat[,y]) ~ log(plot.dat[,x]),ylim=range(log(plot.dat[,y])[!is.na(plot.dat[,x])],na.rm=T),
     cex.main=1.3, font.main=1, pch=21, col="black", main="", bg=date.pty, lwd=1.5, cex.lab=1.4, cex=1.2, cex.axis=1.2,
     xlab="log Virus abundance", ylab="log Bacterial specific growth rate")

lm.t=apply(lmodel2(log(plot.dat[,y]) ~ log(plot.dat[,x]),  nperm=1000)$regression.results[2,-1],2,as.numeric)
lty1=1;if(lm.t['P-perm (1-tailed)']>0.05){lty1=2}
abline(as.numeric(lm.t[c('Intercept','Slope')]),lwd=2,lty=lty1)
r2t=cor(log(plot.dat[,x]),log(plot.dat[,y]),use='complete.obs')^2
pv=round(lm.t['P-perm (1-tailed)'],2); if(pv<0.01){pv='<0.01'}
nf=sum(!apply(is.na(cbind(plot.dat[,c(x,y)])),1,any))

se.slope=sd(log(plot.dat[,y]),na.rm=T)/sd(log(plot.dat[,x]),na.rm=T)*((1-r2t)/nf)^0.5

legend(-1.5,-2.7,
       paste('Slope =',round(lm.t[c('Slope')],2),'+/-',round(se.slope,2)
             ,'\nR2 =',round(r2t,2)
             ,'\np =',pv),
       bty='n',cex=1.2)




###############################################################################
### Cross-correlation analysis for depth-integrated data
cvbm.d=NULL
da.na.d=NULL
dat.intd=NULL
dcri=4 # Two adjacent measures separated by more than 4 hours needs to insert NA 
lgmx=4 # performing cross-correlation with lags at most 4*3=12 hours 
for(i in 1:n.cruise){
  dat.t=filter(plot.dat,cruise_number==i)
  utime=unique(dat.t[,'time'])
  # depth-integrated values observed at each time point
  dat.intd.t=data.frame(cruise_number=i,time=utime,
                        apply(dat.t[,-c(1:7)],2,function(y){z=tap(dp=-dat.t[,'depth'],z=dat.t[,'time'],y,dpL=0,dpd=100,aves=T);return(z)}),# depth-integrated values
                        aggregate(dat.t[,c('BB','SGR','Virus')],list(dat.t[,'time']),function(x){sd(x,na.rm=T)/length(x)})[,-1]) # standard errors of depth-integrated mean at every time points 
  colnames(dat.intd.t)=c('cruise_number','time','BB','SGR','Virus','BB_se','SGR_se','Virus_se')
  dat.intd=rbind(dat.intd,dat.intd.t)
  
  dtime=diff(dat.intd.t[,'time'])
  x.t=as.numeric(dat.intd.t[,'Virus'])
  y.t=as.numeric(dat.intd.t[,'SGR'])
  out=xcor.ts(dtime,x.t,y.t,dcri=4,lgmx=4)
  da.na.d=rbind(da.na.d,data.frame(cruise_number=i,out[[1]]))
  cvbm.d=rbind(cvbm.d,c(i,out[[2]]))
  x.t=y.t=dtime=NULL
}# end of i
colnames(cvbm.d)=c('cruise_number',paste('lag',c(0:lgmx)*3,sep='-'))
colnames(da.na.d)[-c(1:2)]=c('Virus','SGR','SGR_lag6','Virus_TR','SGR_TR','SGR_TR_lag6')

###############################################################################
### Cross-correlation analysis (discrete depth)
cvbm=NULL
da.na=NULL
dcri=4 # Two adjacent measures separated by more than 4 hours needs to insert NA 
lgmx=4 # performing cross-correlation with lags at most 4*3=12 hours 
for(i in 1:n.cruise){
  dat.t=filter(plot.dat,cruise_number==i)
  dep=sort(unique(dat.t[,'depth']),decreasing = T)
  for(j in 1:length(dep)){
    dat.p=filter(dat.t,depth==dep[j])
    # Time-averaged value at each depth
    if(nrow(dat.p)>=3){
      dtime=diff(dat.p[,'time'])
      x.t=as.numeric(dat.p[,'Virus'])
      y.t=as.numeric(dat.p[,'SGR'])
      out=xcor.ts(dtime,x.t,y.t,dcri=4,lgmx=4)
      da.na=rbind(da.na,data.frame(cruise_number=i,depth=rep(dep[j]),out[[1]]))
      cvbm=rbind(cvbm,c(i,dep[j],out[[2]]))
    }
    x.t=y.t=dtime=NULL
  }# end of j
}# end of i

colnames(cvbm)=c('cruise_number','depth',paste('lag',c(0:lgmx)*3,sep='-'))
colnames(da.na)[-c(1:3)]=c('Virus','SGR','SGR_lag6','Virus_TR','SGR_TR','SGR_TR_lag6')

######################################################################
# Fig 2. The diel changes of depth-averaged values
win.graph(90,60)
par(mfcol=c(3,3),mar=c(4,4,2,4))
for(i in 1:n.cruise){
  dat.DepInt=filter(dat.intd,cruise_number==i) # Short time series of depth-integrated values
  v.up=dat.DepInt[,'Virus']+dat.DepInt[,'Virus_se'];v.lo=dat.DepInt[,'Virus']-dat.DepInt[,'Virus_se']; # depth-integrated values +/- S.E.M. 
  bm.up=dat.DepInt[,'SGR']+dat.DepInt[,'SGR_se'];bm.lo=dat.DepInt[,'SGR']-dat.DepInt[,'SGR_se']; # depth-integrated values +/- S.E.M.
  plot(Virus~time,dat.DepInt,ylim=range(c(v.lo,v.up),na.rm=T),type='n',xlab='Elapsed time (hours)',ylab='Virus abundance')
  polygon(c(dat.DepInt[,'time'],dat.DepInt[nrow(dat.DepInt):1,'time']),c(v.up,v.lo[length(v.lo):1]),col=t_col('grey',perc = 50),border=NA)
  lines(Virus~time,dat.DepInt,type='l',lwd=2)
  par(new=T)
  plot(SGR~time,dat.DepInt,axes=F,ylab='',xlab='',ylim=range(c(bm.lo,bm.up),na.rm=T),type='n')
  polygon(c(dat.DepInt[,'time'],dat.DepInt[nrow(dat.DepInt):1,'time']),c(bm.up,bm.lo[length(bm.lo):1]),col=t_col('pink',perc = 50) ,border=NA)
  lines(SGR~time,dat.DepInt,type='l',lwd=2,col='red')
  legend('topright',site.name2[i],bty='n')
  axis(4, col="red", col.axis="red")
  mtext("Bacterial specific growth rate", side=4, line=2,cex=0.7, col='red')
}


#########################################################################
### Fig. 3 Statistical relationship between virus abundance and 6-hour lagged SGR 
win.graph(20, 10); par(mfrow=c(1,2))
x='Virus';y='SGR_lag6'
da.t=da.na
plot(log(da.t[,y]) ~ log(da.t[,x]),ylim=range(log(da.t[,y]),na.rm=T),
     cex.main=1.3, font.main=1, pch=21, col="black", main="", bg=color.all[da.t[,'cruise_number']], lwd=1.5, cex.lab=1.4, cex=1.2, cex.axis=1.2,
     xlab="log Virus abundance", ylab="log Bacterial specific growth")


lm.t=apply(lmodel2(log(da.t[,y]) ~ log(da.t[,x]),  nperm=1000)$regression.results[2,-1],2,as.numeric)
lty1=1;if(lm.t['P-perm (1-tailed)']>0.05){lty1=2}
abline(as.numeric(lm.t[c('Intercept','Slope')]),lwd=2,lty=lty1)
r2t=cor(log(da.t[,x]),log(da.t[,y]),use='complete.obs')^2
pv=round(lm.t['P-perm (1-tailed)'],2); if(pv<0.01){pv='< 0.01'}
nf=sum(!apply(is.na(cbind(da.t[,c(x,y)])),1,any))

se.slope=sd(log(da.t[,y]),na.rm=T)/sd(log(da.t[,x]),na.rm=T)*((1-r2t)/nf)^0.5

legend('topleft',paste('R2=',round(r2t,2)
                        ,'\np ',pv),bty='n',cex=1.2)

x='Virus_TR';y='SGR_TR_lag6'
plot(da.t[,y] ~ da.t[,x],ylim=range(da.t[,y],na.rm=T),
     cex.main=1.3, font.main=1, pch=21, col="black", main="", bg=color.all[da.t[,'cruise_number']], lwd=1.5, cex.lab=1.4, cex=1.2, cex.axis=1.2,
     xlab="Empirical viral turnover rate (h-1)", ylab="Empirical bacterial turnover rate (h-1)")


lm.t=apply(lmodel2(da.t[,y] ~ da.t[,x],  nperm=1000)$regression.results[2,-1],2,as.numeric)
lty1=1;if(lm.t['P-perm (1-tailed)']>0.05){lty1=2}
abline(as.numeric(lm.t[c('Intercept','Slope')]),lwd=2,lty=lty1)
r2t=cor(da.t[,x],da.t[,y],use='complete.obs')^2
pv=round(lm.t['P-perm (1-tailed)'],2); if(pv<0.01){pv='< 0.01'}
nf=sum(!apply(is.na(cbind(da.t[,c(x,y)])),1,any))

se.slope=sd(da.t[,y],na.rm=T)/sd(da.t[,x],na.rm=T)*((1-r2t)/nf)^0.5

legend('topleft',paste('R2=',round(r2t,2)
                       ,'\np ',pv),bty='n',cex=1.2)
legend('bottomright',site.name2,col=color.all,pch=rep(16),bty='n')


##################################################################################################
### Fig. 4 Violin plots for cross-correlation coefficients with fitted statistical relationships between time-lags and cross-correlation coefficients
lagh=seq(0,12,3)
# Depth-integrated results
cvbm.viod=data.frame(lag=as.numeric(rep(lagh,each=n.cruise)),value=c(cvbm.d[,-1]))
# Discrete depths
cvbm.viot=data.frame(lag=as.numeric(rep(lagh,each=nrow(cvbm))),depth=as.numeric(rep(cvbm[,'depth'],length(lagh))),value=c(cvbm[,-c(1:2)]))
cvbm.viod=cvbm.viod[!is.na(cvbm.viod[,'value']),];cvbm.viot=cvbm.viot[!is.na(cvbm.viot[,'value']),];

g1=plotF(da.t=cvbm.viod) 
g2=plotF(da.t=cvbm.viot)
# Mixed layer
gs=plotF(da.t=filter(cvbm.viot,abs(depth)<50))
# Subsurface layer
gd=plotF(da.t=filter(cvbm.viot,abs(depth)>=50))

win.graph(60,60)
grid.arrange(g1, g2, gs, gd, nrow = 2)
