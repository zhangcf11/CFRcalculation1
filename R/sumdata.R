#' Title
#' @description Conversion of individual data to aggregated data
#' @param data simulation data,for example sdata.gamma and sdata.wei
sumdata<-function(data)
{
  mday<-max(data$time2)
  time<-numeric(0);n<-numeric(0);d<-numeric(0);c<-numeric(0);N<-numeric(0);D<-numeric(0);C<-numeric(0);a<-numeric(0)
  for(i in 1:mday)
  {
    time[i]<-i
    n[i]<-length(subset(data,time1==i)$dc)
    N[i]<-length(subset(data,time1<=i)$dc)
    d[i]<-length(subset(data,time2==i&dc==1)$dc)
    D[i]<-length(subset(data,time2<=i&dc==1)$dc)
    c[i]<-length(subset(data,time2==i&dc==2)$dc)
    C[i]<-length(subset(data,time2<=i&dc==2)$dc)
  }
  a<-N-D-C
  data<-data.frame(day=time,n=n,d=d,c=c,Conf=N,Death=D,Recov=C,a=a)
  data
}
