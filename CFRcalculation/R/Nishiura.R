Nishiura <- function(data,t,m,sd,distr){
  if (is.element("Bhat", installed.packages()[,1])==FALSE){
    install.packages("Bhat")
  }
  require(Bhat) #same as library statement
  library(Bhat)
  # estimate the exponential growth rate, r
  C=function(time)
  {
    Conf=data[data$day==time,]$Conf
    return(Conf)
  }
  CI=0
  for (i in data$day[1]:(t-1))
  {CI=CI+C(i)}

  f=function(r)
  {
    -CI*(1-exp(-r))+(C(t)-C(data$day[1]))*exp(-r)
  }

  r=uniroot(f,lower=0,upper=1)[[1]]
  # factor of underestimation
  # 1-exp, 2-gamma
  if (distr==1)
  {u=1/(1+r*m)}
  if (distr==2)
  {# v  coefficient of variation
    v=sd/m
    u=(1+r*m*v^2)^(-1/(v^2))
  }
  Conf=data[data$day==t,]$Conf
  Death=data[data$day==t,]$Death
  nlogf=function(x) # negative loglikelihood
  {
    lnL=-(Death*log(x)+(u*Conf-Death)*log(1-x))
    return(lnL)
  }
  x <- list(label=c("p"),est=c(0.1),low=c(0),upp=c(0.9))
  q  <-  dfp(x,f=nlogf)
  x$est<- q$est
  if (is.na(try( plkhci(x,nlogf,'p',prob=0.95))[2]))
  {CFR=NA;lower=NA;upper=NA;var=NA}
  else {
    interval <- plkhci(x,nlogf,'p',prob=0.95)
    lower=interval[1]
    upper=interval[2]
    CFR=x$est
  }
  list("CFR"=CFR,"lower"=lower,"upper"=upper,"var"=NA)
}
