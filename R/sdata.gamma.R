#' Simulate Individual Data Using Weibull Distribution
#' @description Simulate individual data using Weibull distribution
#' @usage sdata.gamma(num, cfr, para, seed = Sys.Date())
#' @param num number of size.
#' @param cfr the set cfr.
#' @param para the parameter values. para[1:2] for shape and scale of death,para[3:4] for shape and scale of cure.
#' @param seed the default seed is current time.
#' @usage sdata.wei(num, cfr, para, seed = Sys.Date())
#'
#' @return A list containing the following components:
#' \item{id}{patient id number.}
#' \item{time1}{entry time.}
#' \item{time2}{the time of cure or death.}
#' \item{time}{survival time.}
#' \item{dc}{indicator of outcome.1-death,2-cure.}
#' @examples
#' sdata.gamma(120,0.2,para=c(25,300,35,300)) #for example
sdata.gamma<-function(num,cfr,para,seed=Sys.Date()){

  seed<-as.numeric(seed)

  set.seed(seed);rdata<-runif(3*num)

  id<-c(1:num)

  dc<-rdata[(2*num+1):(3*num)]

  num1<-0
  num2<-0
  for(i in 1:num) {
    if (dc[i]<=cfr)
    {dc[i]<-1
    num1<-num1+1}
    else {dc[i]<-2;num2<-num2+1}
  }

  scale1<-para[1]/para[2]
  shape1<-para[1]^2/para[2]
  scale2<-para[3]/para[4]
  shape2<-para[3]^2/para[4]

  set.seed(2*seed+1);dtau<-c(ceiling(rgamma(num1,shape1,scale1)))
  set.seed(3*seed+1);ctau<-c(ceiling(rgamma(num2,shape2,scale2)))

  tau<-c(0);dd<-0;cc<-0
  for(i in 1:num){
    if(dc[i]==1){dd<-dd+1;tau[i]<-dtau[dd]}
    else {cc<-cc+1;tau[i]<-ctau[cc]}
  }

  set.seed(4*seed+1);ti<-c(ceiling(rgamma(num,shape=6.319032,scale=(0.16101)^(-1))))

  ti<-ti-min(ti)+1
  tend<-c(ti+tau-1)

  dat<-data.frame(time1=ti,time2=tend,time=tau,dc=dc)
  dat<-dat[order(dat$time1),]
  data<-data.frame(id=id,time1=dat$time1,time2=dat$time2,time=dat$time,dc=dat$dc)

  return(data)
}
