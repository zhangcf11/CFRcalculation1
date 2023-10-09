Chen <- function(data,t){
  sumc=sum(data[data$day<=t,]$c)
  sumd=sum(data[data$day<=t,]$d)
  theta=sumc/sumd
  CFR=1/(1+theta)
  avartheta=sumc*(sumc+sumd)/((sumd)^3);
  avar=avartheta/((1+theta)^4)
  sd=avar^(1/2)
  upper=CFR+1.96*sd
  lower=CFR-1.96*sd
  list("CFR"=CFR,"lower"=lower,"upper"=upper,"var"=avar)
}
