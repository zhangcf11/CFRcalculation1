#' Title Transform Individual Data By Survival Time
#' @description  Transform Individual Data By Survival Time.
#' @usage chdata(data)
#' @param data individual survival data transformed by surdata
#'
#'
#' @examples
#' data(individual)
#' #get data
#' data1<-surdata(data=individual,t=60,type=1)
#' #study time:60;type:1,cout time.
#' chdata(data=data1)
chdata<-function(data)
{
  sf<-data

  sut<-sort(unique(sf$time))

  sd<-data.frame(time=sut,p1=0,p2=0,event=0,censor=0,risk=0)

  for(i in 1:length(sut))
  {

    sf1<-sf[sf$time==sut[i],]

    sd$censor[i]<-length(sf1[sf1$dc==0,]$time)

    sd$p1[i]<-length(sf1[sf1$dc==1,]$time)

    sd$p2[i]<-length(sf1[sf1$dc==2,]$time)

    sd$event[i]<-length(sf1[sf1$dc!=0,]$time)

    sd$risk[i]<-length(sf[sf$time>=sut[i],]$time)
  }

  sd$h1<-sd$p1/sd$risk
  sd$h2<-sd$p2/sd$risk
  sur<-1-sd$event/sd$risk
  sd$survival<-cumprod(sur)

  return(sd)
}
