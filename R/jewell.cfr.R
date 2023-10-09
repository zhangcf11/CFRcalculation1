#' Modified KM Method For CFR In Individual Survival Data
#' @description Calculate CFR using modified Kaplan-Meier(KM) method, including method-A(death cumulative incidence) and method-B(death cumulative incidence divided by death and cure cumulative incidence). See references for details
#' @usage jewell.cfr(data)
#' @param data individual survival data transformed by surdata.

#' @return A list containing the following components:
#' \describe{
#' \item{method}{ the source.}
#' \item{CFR_a}{the interpretation of method A.}
#' \item{CFR_b}{the interpretation of method B.}
#' \item{estimation}{the estimate at the study time,which contains method,Time(study time),CFR,Greedvar(Greedwood variance) and CI.}
#' }
#' @references  Nicholas P.Jewell. Non-parametric estimation of the case fatality ratio with competing risks data:An application to Severe Acute Respiratory Syndrome(SARS), Statist.Med. 2007,26,1982-1998.
#'
#' @examples
#' data(individual)
#' data1<-surdata(data=individual,t=60,type=1)
#' jewell.cfr(data1)
#' #calculate CFR
jewell.cfr<-function(data){
  sf<-data
  sf3<-chdata(data)
  n3<-length(sf3$risk)
  for(i in 1:n3)
  {
    if(i==1) {sf3$kmd[i]<-sf3$h1[i];sf3$kmc[i]<-sf3$h2[i]}
    else
    {
      sf3$kmd[i]<-(sf3$h1[i])*(sf3$survival[i-1])
      sf3$kmc[i]<-(sf3$h2[i])*(sf3$survival[i-1])

    }
  }
  f1<-sum(sf3$kmd)
  f2<-sum(sf3$kmc)
  CFR_a<-f1
  CFR_b<-f1/(f1+f2)

  l<-length(sf3$time)
  F01<-rep(0,l);F02<-rep(0,l)
  M<-matrix(0,(l-1),l)
  for(i in 1:(l-1)){
    sf4<-sf3[(i+1):l,]
    for(j in (i+1):l){
      M[i,j]<-prod(1-sf3$event[(i+1):j]/sf3$risk[(i+1):j])
    }
    F01[i]<-sum((sf4$h1)*M[i,i:(l-1)])
    F02[i]<-sum((sf4$h2)*M[i,i:(l-1)])
  }
  sur<-c(1,sf3$survival[1:(l-1)])
  cov_gw<--sum(sur^2*(1-F01)*F02*(sf3$risk-1)*sf3$p1/(sf3$risk^3))-sum(sur^2*F01*(1-F02)*(sf3$risk-1)*sf3$p2/(sf3$risk^3))
  var1_gw<-sum((sur*F01)^2*(sf3$risk-1)*sf3$event/(sf3$risk^3))+sum(sur^2*(1-2*F01)*(sf3$risk-1)*sf3$p1/(sf3$risk^3))
  var0_gw<-sum((sur*F02)^2*(sf3$risk-1)*sf3$event/(sf3$risk^3))+sum(sur^2*(1-2*F02)*(sf3$risk-1)*sf3$p2/(sf3$risk^3))

  Var_a<-var1_gw
  Var_b<-(f2^2*var1_gw+f1^2*var0_gw-2*f1*f2*cov_gw)/(f1+f2)^4  #CFR_b方差

  time<-attributes(sf)$stime
  upa<-CFR_a+1.96*sqrt(Var_a)
  upb<-CFR_b+1.96*sqrt(Var_b)
  lowa<-CFR_a-1.96*sqrt(Var_a)
  lowb<-CFR_b-1.96*sqrt(Var_b)
  out<-data.frame(Method=c("CFR_a","CFR_b"),time=c(time,time),CFR=c(CFR_a,CFR_b),
                  Greedvar=c(Var_a,Var_b),Up95=c(upa,upb),Low95=c(lowa,lowb))
  return(list(method="Jewell(Statist.Med.,2007)",CFR_a="cumulative death incidence",CFR_b="cumulative death incidence divided by total cumulative incidence",estimation=out))
}
