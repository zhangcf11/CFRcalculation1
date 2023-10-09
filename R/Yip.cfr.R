#' Calculating CFR Using Chain Multinomial In Summarized Data
#' @description Calculating CFR chain multinomial in summarized data. See references for details.
#' @usage Yip.cfr(data,t,B)
#' @param data summarized survival data
#' @param t the study time.
#' @param B the window width of a kernel function.
#' @return A list containing the following components:
#' \item{CFR}{CFR.}
#' \item{lower}{the lowe of CFR's  CI}
#' \item{upper}{the upper of CFR's CI.}
#' \item{var}{ the variance of CFR}

Yip.cfr<- function(data,t,B){
  Kb=function(x,b) # construct kernel function
  {
    if (abs(x/b)<1) y=0.75*(1-(x/b)^2)/b
    else y=0
    return(y)
  }
  ai=c(NA,data[1:(nrow(data)-1),]$a)# the number of cases time at (t-1)

  data$ai=ai
  data$p1=data$d/data$ai
  data$p2=data$c/data$ai
  sum=0;sum1=0;sum2=0;sumv1=0;sumv2=0;sumcov=0
  for(s in data$day[2]:t )
  {
    sum=Kb(t-s,B)+sum
    sum1=Kb(t-s,B)*data[data$day==s,]$p1+sum1
    sum2=Kb(t-s,B)*data[data$day==s,]$p2++sum2
    sumv1=(Kb(t-s,B)^2)*data[data$day==s,]$p1*(1-data[data$day==s,]$p1)/data[data$day==s,]$ai+sumv1
    sumv2=(Kb(t-s,B)^2)*data[data$day==s,]$p2*(1-data[data$day==s,]$p2)/data[data$day==s,]$ai+sumv2
    sumcov=(Kb(t-s,B)^2)*data[data$day==s,]$p1*data[data$day==s,]$p2/data[data$day==s,]$ai+sumcov
  }
  P1=sum1/sum
  P2=sum2/sum
  varP1=sumv1/(sum^2) # variance of P1
  varP2=sumv2/(sum^2) # variance of P2
  covP1P2=-sumcov/(sum^2) # covariance of P1 and P2
  theta=P1/P2
  vartheta=(varP1/(P2^2))+(P1^2)*varP2/(P2^4)-2*P1*covP1P2/(P2^3) # variance of theta
  CFR=theta/(1+theta)
  var=vartheta/((1+theta)^4) # variance of CFR
  sd=sqrt(var)
  upper=CFR+1.96*sd
  lower=CFR-1.96*sd
  return(list("CFR"=CFR,"lower"=lower,"upper"=upper,"Var"=var))
}
