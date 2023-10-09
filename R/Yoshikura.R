Yoshikura <- function(data,t){
  pointestimate <- function(data,t)
  {
    # the number of confirmed cases when there is first death
    X0=data[data$day <=t & data$Death!=0,]$Conf

    # the number of deaths when there is first death
    Y0=data[data$day <=t & data$Death!=0,]$Death

    if (length(Y0)<=1) {CFR=NA}
    else{
      C=numeric(0); D=numeric(0)
      rec=numeric(0)
      h=1; i=1
      repeat{if (Y0[i+1]!=Y0[i]) {D[h]=Y0[i]; rec[h]=i; h=h+1}
        i=i+1
        if (i+1>length(Y0)) break}
      for(i in 1: length(rec)){C[i]=X0[rec[i]]}
      X=log(C)
      Y=log(D)
      if (length(Y)<=1) {CFR=NA}
      else {
        lm.sol<-lm(Y~X)
        Intercep=summary(lm.sol)[[4]][1]
        slope=summary(lm.sol)[[4]][2]
        N0=exp(-Intercep/slope)
        CFR=((max(C)/N0)^slope)/max(C)
      }
    }
    return(CFR)
  }
  # boostrap method to calculate variance
  # The time when there is first death
  if (is.na(pointestimate(data,t))) {upper=NA;lower=NA;var=NA}
  else {
    firsttime=data[data$Death!=0,]$day[1]
    time=firsttime
    # Sampling
    samp=numeric(0)
    i=1
    repeat{
      samp[i]=pointestimate(data,t=time)
      time=time+1;i=i+1
      if (time>t) break
    }
    Length=length(samp)
    # resampling K times
    s=numeric(0)
    k=500
    for (i in 1:k)
    {
      resamples=sample(samp,length(samp), replace = T)
      s[i]=mean(resamples,na.rm = T)
    }
    CFR=pointestimate(data,t)
    var=var(s,na.rm = T)
    upper=CFR+1.96*sd(s)
    lower=CFR-1.96*sd(s)
  }
  list("CFR"=CFR,"lower"=lower,"upper"=upper,"var"=var)
}
