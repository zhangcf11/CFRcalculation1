#' Estimate Coefficient values In Individual Survival Data
#' @description Estimate coefficient values of covariate.
#' @usage chencov.p(data, cov, p)
#' @param p starting parameter values;p[1] for beta0,p[2:(n+1)] for beta1 with respect to death,p[(n+2):(2*n+1)] for beta2 with respect to cure.
#' @param cov covariate
#' @param data individual survival data transformed by surdata.
#'
#' @return A list containing the following components:
#' \item{method}{the source.}
#' \item{cov}{ covariate.}
#' \item{estimate}{ coefficient values.}
#' \item{var}{ variance matrix}
#' \item{time}{ study time}
#' \item{deaths}{ number of death}
#' @references Zheng Chen. Estimating the case fatality raye using a constant cure-death hazard ration. Lifetime Data Anal. 2009,15,316-329.
#' @examples
#' data(individual)
#' #get individual data
#' sudata<-surdata(data=individual,t=60,type=1)
#' #data for calculation
#' ## Not run:
#' chencov.p(sudata,cov=c("sex","age"),p=c(0.2,0.4,0.5,0.4,0.6))
#' ## End(Not run)
#' #cov:covariate; p:starting parameter values.
chencov.p<-function(p,cov,data)
{
  pc<-p
  data1<-data
  covc<-cov
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
  out=nlm(chen.cov,p=pc,cov=covc,data=data1,hessian=T)
  estimate=out$estimate
  sd=solve(out$hessian)
  return(list(cov=covc,estimate=estimate,var=sd,time=attributes(data1)$stime,deaths=attributes(data1)$deaths))
}
