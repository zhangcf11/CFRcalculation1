#' Calculate CFR In Individual Survival Data
#' @description Calculate CFR using a constant cure-death hazard ratio,  See references for details
#' @usage chen.cfr(data)
#' @param data individual survival data transformed by surdata
#'
#' @return A list containing the following components:
#' \item{method}{the source}
#' \item{estimation}{ estimate of CFR and variance}
#' @references Zheng Chen. Estimating the case fatality raye using a constant cure-death hazard ration. Lifetime Data Anal. 2009,15,316-329.
chen.cfr<-function(data)
{
  sf<-data
  sf2<-chdata(data)
  a<-sum(sf2$p2)/sum(sf2$p1)
  sda<-sum(sf2$p2)*sum(sf2$p2+sf2$p1)/(sum(sf2$p1))^3
  cfr<-1/(1+a)
  sd<-sda/(1+a)^4
  time<-attributes(sf)$stime
  out<-data.frame(time=time,cfr=cfr,var=sd)
  return(list(method="constant cure-death hazard ratio",estimation=out))
}
