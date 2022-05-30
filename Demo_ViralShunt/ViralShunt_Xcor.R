##########################################################################
### Function for cross-correlation analysis
xcor.ts=function(dtime,x.t,y.t,dcri=4,lgmx=4,y.lag=2){
  dtime.na=c(0,which(dtime>dcri),length(x.t))
  e.n=sum(apply(!is.na(cbind(x.t,y.t)),1,all))
  # Filling up NA values in the short time series when two adjacent measures are separated more than 4 hours
  if(length(dtime.na)>2){  
    x.na=NULL
    y.na=NULL
    for(k in 2:length(dtime.na)){
      x.na=c(x.na,x.t[(dtime.na[k-1]+1):dtime.na[k]],NA)
      y.na=c(y.na,y.t[(dtime.na[k-1]+1):dtime.na[k]],NA)
    }
    x.na=x.na[-length(x.na)]
    y.na=y.na[-length(y.na)]
    x.t=x.na
    y.t=y.na
  }
  n.t=length(x.t)
  
  
  x.i=as.numeric(cinter(x.t));y.i=as.numeric(cinter(y.t))
  n.t=length(x.i)
  
  x.r=c(log(x.i[-1]/x.i[-n.t]),NA)/3
  y.r=c(log(y.i[-1]/y.i[-n.t]),NA)/3
  
  
  da.na.t=data.frame(time = seq(0,3*(length(x.t)-1),3),
                     x = x.t, y = y.t, y_lag = c(y.t[-c(1:(y.lag))],rep(NA,y.lag)),
                     xr = x.r, yr = y.r, yr_lag=c(y.r[-c(1:(y.lag))],rep(NA,y.lag)))
  
  
  dn=n.t-(lgmx+3)
  
  # Cross-correlation analysis
  if(e.n>=lgmx){ # Select short time series without too many NAs
    ccf1=ccf(x=x.i,y=y.i,plot=F,type = c("correlation"),
             lag.max=lgmx,na.action = na.pass, ylab = "cross-correlation")
    # Extract cross-correlation coefficients
    ccv.t=as.numeric(ccf1[c(-lgmx:0)]$acf)
    if(dn<0){ccv.na=rep(NA,1+lgmx);ccv.na[(1+abs(dn)):(1+lgmx)]=ccv.t[(1+abs(dn)):(1+lgmx)];ccv.t=ccv.na}
    cvbm.t=ccv.t[(lgmx+1):1]
  }else{cvbm.t=rep(NA,lgmx+1)}
  
  return(list(da.na.t,cvbm.t))  
}


####################################################################################
## Plot violin plot and conduct statistical fitting
plotF=function(da.t){
  newx=seq(0,12,0.01)
  newx2=seq(0,6,0.01)
  da.t2=filter(da.t,lag<=6)
  
  # Statistical fitting based on simple linear regression and regression splines 
  lm.t=lm(value~lag,da.t2)
  fit0 <- lm(value~1, da.t)
  fit2 <- lm(value~bs(lag,degree = 3), da.t)
  
  (stat.t=c(R2_linear=summary(lm.t)$r.square,
            p_linear=summary(lm.t)$coefficients[2,4],
            R2_cubic=cor(predict(fit2),da.t[!is.na(da.t[,'value']),'value'])^2,
            p_cubic=anova(fit0,fit2)[['Pr(>F)']][2],n=nrow(filter(da.t,lag==0))))
  
  pv=c(stat.t[2],stat.t[4]);pv2=round(pv,2);pv2[pv2<0.01]='<0.01';pv2[pv2!='<0.01']=paste('=',pv2[pv2!='<0.01'],sep='')
  
  # Model predictions 
  dat1= data.frame(x=newx2, y=predict(lm.t,data.frame(lag=newx2)))
  dat2= data.frame(x=newx, y=predict(fit2,data.frame(lag=newx)))
  
  # Violin plot + fitted results
  gi = ggplot() + 
    geom_violin(data = da.t, mapping = aes(x = lag, y = value, group = lag), trim=T, fill="gray")+
    scale_x_continuous(breaks = seq(0,12,3)) +
    scale_y_continuous(limits = c(-0.85, 0.85))+
    geom_boxplot(data = da.t, mapping = aes(x = lag, y = value, group = lag), width=0.3) +
    geom_line(data = dat1, mapping = aes(x = x, y = y),color='red', size=1, linetype = 1+(pv[1]>0.05)*1+(pv[1]>0.1)*1)+
    geom_line(data = dat2, mapping = aes(x = x, y = y),color='blue', size=1, linetype = 1+(pv[2]>0.05)*1+(pv[2]>0.1)*1)+
    geom_hline(yintercept=0,col='grey')+
    geom_text(aes(x=10, y=-0.7, label = paste('R2=',round(stat.t[1],2),', p',pv2[1],sep='')),color='red')+
    geom_text(aes(x=10, y=-0.8, label = paste('R2=',round(stat.t[3],2),', p',pv2[2],sep='')),color='blue')+
    labs(title="",x="Time lag", y = "Cross-correlation coefficient")+
    theme_classic()
  return(gi)
  
}


#######Function for cubic-spline interpolation
cinter=function(x){
  x=as.matrix(x)
  n=nrow(x)
  y=x
  for(i in 1:ncol(x)){
    x.t=x[,i]
    id=which(is.na(x.t))
    if(length(id)!=0){
      xtt=cbind(1:n,x.t)    
      inp=fields::splint(xtt[-id,1],xtt[-id,2],1:n)
      ####exclude extrapolation 
      nak=NULL
      for(j in 1:length(id)){      
        ind=which((id-j)==0)
        if(length(ind)==0){break}else{nak=c(nak,id[ind])}      
      }
      for(j in 1:length(id)){      
        ind2=which((id-n-j+1)==0)
        if(length(ind2)==0){break}else{nak=c(nak,id[ind2])}      
      }
      id=setdiff(id,nak)
      if(length(id)!=0){x.t[id]=inp[id];y[,i]=x.t}      
    }
  }
  return(y)
}

# Function for depth-integration; dp = depth of integration y :data; z: time; dpL: minimal depth; dpd: maximal depth
tap=function(dp,y,z,dpL,dpd,aves=T){          
  dp=as.numeric(as.matrix(dp))
  uz=unique(z)
  ddy=cbind(abs(dp),y)
  tl=NULL
  for(i in 1:length(uz)){
    if(all(is.na(ddy[,2]))){
      tl=c(tl,NA)
    }else{
      dy=matrix(ddy[which(z==uz[i]),],ncol=2)
      dy=matrix(dy[dy[,1]<=dpd&dy[,1]>=dpL,],ncol=2)
      dy=matrix(dy[order(dy[,1]),],ncol=2)
      a=which(is.na(dy[,2]))
      if(length(a)!=0){dy=dy[-a,]}
      n=length(dy)
      if(n!=0){
        if(n>2){
          x.n1 = dy[,1][-1]
          x.n0 = dy[,1][-n/2]
          y.n1 = dy[,2][-1]
          y.n0 = dy[,2][-n/2]
          if(aves){ccc=0.5*t(x.n1-x.n0)%*%(y.n1+y.n0)/(dy[,1][n/2]-dy[,1][1])}else{
            ccc=0.5*t(x.n1-x.n0)%*%(y.n1+y.n0)}
          if(all(dy[,1]==unique(dy[,1])[1])){ccc=mean(dy[,2],na.rm=T)}
        }else{ccc = dy[2]}
      }else(ccc=NA)
      if(all(is.na(dy))){ccc=NA}
      tl=c(tl,ccc)
      dy=NULL
    }
  }
  return(tl)
}

## Make transparent colors; Mark Gardener 2015; www.dataanalytics.org.uk
t_col <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  invisible(t.col)
}

