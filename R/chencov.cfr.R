#' Obtain CFR In Individual Survival Data
#' @description  Obtain CFR with covariate value for chencov.p. See references for details.
#' @usage chencov.cfr(object, covv)
#' @param object a model object from chencov.p.
#' @param covv covariate value.
#'
#' @return A list containing the following components:
#' \item{method}{the source.}
#' \item{estimation}{the estimate at the study time,which contains Time(study time),CFR and Var(variance)}
#' @references Zheng Chen. Estimating the case fatality raye using a constant cure-death hazard ration. Lifetime Data Anal. 2009,15,316-329.
#' @examples
#' data(individual)
#' #get individual data
#' sudata<-surdata(data=individual,t=60,type=1)
#' #data for calculation,study time:60
#' ## Not run:
#' obj<-chencov.p(sudata,cov=c("sex","age"),p=c(0.2,0.4,0.5,0.4,0.6))
#' ## End(Not run)
#' #cov:covariate; p:starting parameter values.
#' ## Not run:
#' chencov.cfr(object=obj,covv=c(0,1))
#' ## End(Not run)
#' #calculating CFR with covariate value c(0,1);sex=0 for male,age=1 for more than 60 years old.

chencov.cfr<-function(object,covv)
{
  beta<-object$estimate
  T<-length(object$cov)
  pvar<-object$var
  time<-object$time
  beta0<-beta[1]
  beta1<-beta[2:(T+1)]
  beta2<-beta[(T+2):(2*T+1)]
  theta<-exp(beta0-sum(beta1*covv)+sum(beta2*covv))
  cfr<-1/(1+theta)
  logvar<-t(c(1,-covv,covv))%*%pvar%*%c(1,-covv,covv)
  var.theta<-logvar*(theta)^2
  var<-var.theta/(1+theta)^4
  time<-object$time
  result<-data.frame(time=time,cfr=cfr,var=var)
  return(result)
}
