#' Title Transform Individual Data for Calculation
#' @description  Transform collected individual survival data into which can be used to calculate CFR.
#' @usage surdata(data, t, type)
#' @param data the collected individual data.
#' @param t the study time.
#' @param type the type of collected data.type=1 for cout time(e.g. 1,2....), type=2 for calendar time, type=3 for survival time.
#'
#' @return The returned data.frame contains original information only prior to cut point time(study time).
#' @examples
#' data(individual)
#' #get data
#' surdata(data=individual,t=60,type=1)
#' #study time:60;type:1,cout time.
surdata <-
  function(data,t,type)
  {

    ###############################################
    if(type!=3)
    {
      if(type==1)
      {
        sf<-data
        sf$time2<-sf$time2-min(sf$time1)+1
        sf$time1<-sf$time1-min(sf$time1)+1
        sf<-sf[sf$time1<t,]

        {if(max(sf$time2)<=t) sf<-sf
          else{
            sf[sf$time2>t,]$dc<-0
            sf[sf$time2>t,]$time2<-t
          }

          sf$time<-sf$time2-sf$time1+1
        }
      }

      if(type==2)
      {	t<-as.Date(t)
      sf<-data
      sf$time1<-as.Date(sf$time1)
      sf$time2<-as.Date(sf$time2)
      sf<-sf[sf$time1<t,]
      {if(max(sf$time2)<=t) sf<-sf
        else{
          sf[sf$time2>t,]$dc<-0
          sf[sf$time2>t,]$time2<-t

        }
      }

      sf$time<-as.numeric(difftime(as.Date(sf$time2),as.Date(sf$time1),units="days"))+1
      }
    }

    ############################################################################
    if(type==3)
    {
      sf<-data
      {if(max(sf$time)<=t) sf<-sf
        else	{
          sf[sf$time>t,]$dc<-0
          sf[sf$time>t,]$time<-t
        }
      }
    }

    n<-length(sf$dc)
    death<-length(sf[sf$dc==1,]$dc)
    cure<-length(sf[sf$dc==2,]$dc)
    sf["id"]<-1:n
    attr(sf,"cure")<-cure
    attr(sf,"deaths")<-death
    attr(sf,"Sample size")<-n
    attr(sf,"stime")<-t

    return(sf)
  }

