#' Calculating CFR Using Parametric Method In Individual Data
#' @description Calculating CFR using gamma parametric method in individual survival data. See references for details.
#' @usage out.donga(data, p)
#' @param data individual survival data transformed by surdata.
#' @param p starting parameter values. p[1:2] for the mean and variance of death,p[3:4] for the mean and variance of cure,p[5] for CFR.
#'
#' @return A list containing the following components:
#' \item{method}{the source.}
#' \item{estimation}{the estimate at the study time,which contains Time(study time),CFR,Var(variance) and CI.}
#' \item{p}{parameter values}
#' \item{Var}{variance matrix}
#' @references Christl A Donnelly. Epidemiological determinants of spread of causal agent of severe acute respiratory syndrome in Hong Kong. The Lancet,2003,361:1761-1832.
#'
#' @examples
#' data(individual)
#' #get data
#' data1<-surdata(data=individual,t=60,type=1)
#' ##survival data for calculation
out.donga<-function(data,p)
{
  donga<-function(data,p)
  {
    sf<-data
    b1<-data.frame(time=sf$time,dc=sf$dc)
    n<-length(b1$time)
    l<-rep(1,n)
    for (i in 1:n)
    {
      if (b1$dc[i]==1) l[i]<-p[5]*(pgamma(b1$time[i]+1,p[1],p[2])-pgamma(b1$time[i],p[1],p [2]))
      else if (b1$dc[i]==2) l[i]<-(1-p[5])*(pgamma(b1$time[i]+1,p[3],p[4])-pgamma(b1$time[i],p[3],p[4]))
      else if (b1$dc[i]==0) l[i]<-(p[5]*(1-pgamma(b1$time[i],p[1],p[2]))+(1-p[5])*(1-pgamma(b1$time[i],p[3],p[4])))
    }
    ll<-sum(log(l))
    return(-ll)
  }


  out<-nlm(donga,p=p,data=data,hessian=T)
  var<-solve(out$hessian)
  p<-out$estimate[5]
  sd1<-var[5,5]
  time<-attributes(data)$stime
  up<-p+1.96*sqrt(sd1)
  low<-p-1.96*sqrt(sd1)
  df<-data.frame(Time=time,CFR=p,Var=sd1,Up95=up,Low95=low)
  return(list(method="Donnelly(LANCET,2003)",estimation=df,p=out$estimate,Var=var))
}
