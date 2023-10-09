#' Calculating CFR Using KM method In Individual Data
#' @description Calculating CFR using KM method in individual survival data. See references for details.
#' @usage out.donkm(data, p)
#' @param data individual survival data transformed by surdata.
#' @param p starting parameter values. p[1:2] for the mean and variance of death,p[3:4] for the mean and variance of cure,p[5] for CFR.
#'
#' @return A list containing the following components:
#' \item{method}{the source.}
#' \item{estimation}{the estimate at the study time,which contains Time or study time,CFR,Var orvariance and CI.}
#' \item{p}{parameter values}
#' \item{Var}{variance matrix}
#' @references Christl A Donnelly. Epidemiological determinants of spread of causal agent of severe acute respiratory syndrome in Hong Kong. The Lancet,2003,361:1761-1832.

out.donkm<-function(data,p){


  donkm<-function(data,p)
  {
    sf<-data
    sf2<-chdata(data)
    sf2$kmdeath<-1-cumprod(1-sf2$h1)

    sf2$kmlive<-1-cumprod(1-sf2$h2)

    n<-length(sf2$time)


    lg<-function(x) {ifelse(x==0,0,log(x))}

    l<-numeric(0)
    attach(sf2)
    for(i in 1:n)
    {
      if(i==1) l[i]=p1[i]*lg(p*kmdeath[i])+p2[i]*lg((1-p)*kmlive[i])+censor[i]*lg(1)
      else     l[i]=p1[i]*lg(p*(kmdeath[i]-kmdeath[i-1]))+p2[i]*lg((1-p)*(kmlive[i]-kmlive[i-1]))+censor[i]*lg(p*(1-kmdeath[i-1])+(1-p)*(1-kmlive[i-1]))

    }
    detach(sf2)
    ll<-sum(l)
    return(-ll)
  }

  out<-nlm(donkm,p=p,data=data,hessian=T)
  var<-solve(out$hessian)
  p<-out$estimate
  time<-attributes(data)$stime
  df<-data.frame(time=time,cfr=p,var=var) #矩阵表示
  return(list(method="Donnelly(LANCET,2003) using KM method",estimation=df,p=out$estimate,Var=var))
}

