stepcov.p<-function(cov,p,data)
{

  chen.cov<-function(p,cov,data)
  {
    beta<-p
    sf<-data
    sut<-sort(unique(sf$time))
    sf2<-sf
    sd<-data.frame(time=sut,p1=0,p0=0,event=0,risk=0)

    l<-rep(0,length(sut))
    ln<-length(cov)
    beta0<-beta[1]
    beta1<-beta[2:(ln+1)]
    beta2<-beta[(ln+2):(2*ln+1)]

    for(i in 1:length(sut))
    {
      sf1<-sf[sf$time==sut[i],]
      num1<-0;num2<-0
      for (j in 1:length(sf1$time))
      {
        if (sf1$dc[j]==1)   num1<-num1+1
        if (sf1$dc[j]==2)   num2<-num2+1
      }
      sd$p1[i]<-num1;sd$p0[i]<-num2;sd$event[i]<-num1+num2

      sf2<-sf2[sf2$time!=sut[i],]
      sf3<-rbind(sf2,sf1[sf1$dc!=0,],sf1[sf1$dc==0,])

      sf4<-sf1[sf1$dc==1,]

      sf5<-sf1[sf1$dc==2,]

      sf31<-subset(sf3,select=cov);sf31[is.na(sf31)]<-0
      sf41<-subset(sf4,select=cov);sf41[is.na(sf41)]<-0
      sf51<-subset(sf5,select=cov);sf51[is.na(sf51)]<-0

      l[i]<--(num1+num2)*log(sum(exp(beta1%*%t(sf31)))+sum(exp(beta0+beta2%*%t(sf31))))+sum(beta1%*%t(sf41))+num2*beta0+sum(beta2%*%t(sf51)) #轮廓似然函数
    }

    return(-sum(l))

  }


  beta<-p
  data<-data

  T<-length(cov)

  beta0<-beta[1]
  beta1<-beta[2:(T+1)]
  beta2<-beta[(T+2):(2*T+1)]

  dn<-attributes(data)$deaths

  if(T==1) {cov2=cov;p1=beta}

  if(T>1)
  {

    aic1<-function(cov2,beta10,beta20)
    {
      j<-1

      Aic1<-numeric(0)


      if(all(cov2=="0"))
      {
        repeat
        {cov1=cov[j]
        p1=append(beta0,c(beta1[j],beta2[j]))
        result=nlm(chen.cov,p=p1,cov=cov1,data=data)
        Aic1[j]<-(2*result$minimum)+length(cov1)*log(dn)
        j=j+1;if(j>T) break
        }
      }

      if(all(cov2!="0"))
      {
        repeat
        {
          if(all(cov2!=cov[j]))
          {
            cov1=append(cov2,cov[j])
            p1=append(beta0,c(beta10,beta1[j],beta20,beta2[j]))
            result=nlm(chen.cov,p=p1,cov=cov1,data=data)
            Aic1[j]<-(2*result$minimum)+length(cov1)*log(dn)
          }
          j=j+1;if(j>T) break
        }
      }

      length(Aic1)=T

      return(Aic1)
    }


    order1<-numeric(0);cov2<-character(0);beta10<-numeric(0);beta20<-numeric(0);Aic<-numeric(0)
    Aic1<-aic1(cov2="0")
    for(j in 1:T)
    {
      Aic[j]<-min(Aic1[!is.na(Aic1)])

      if(j==1)
      {
        order1[j]<-which.min(Aic1)
        beta10[j]<-beta1[order1[j]]
        beta20[j]<-beta2[order1[j]]
        cov2[j]<-cov[order1[j]]
      }

      if(j!=1)
      {
        if(Aic[j]>Aic[j-1]) 	{length(Aic)=j-1;cov2=cov2;break}
        if(Aic[j]<=Aic[j-1])
        {
          order1[j]<-which.min(Aic1)
          beta10[j]<-beta1[order1[j]]
          beta20[j]<-beta2[order1[j]]
          cov2[j]<-cov[order1[j]]
        }

        if(j==T)	{Aic=Aic;cov2=cov2;break}
      }

      Aic1<-aic1(cov2,beta10,beta20)

    }
  }
  pc<-c(beta0,beta10,beta20)

  out=nlm(chen.cov,p=pc,cov=cov2,data=data,hessian=T)

  estimate=out$estimate

  sd=solve(out$hessian)

  return(list(BICcr=Aic,cov=cov2,estimate=estimate,var=sd,time=attributes(data)$stime,deaths=attributes(data)$deaths))
}
