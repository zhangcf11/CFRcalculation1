lamtest <- function(data,T){
  num <- dim(data)[1]
  if (num>=T){
    deltaD <- NULL
    deltaR <- NULL
    for (i in 1:num)
    { if (i==1)
    {deltaD[i] <- data$Death[i]
    deltaR[i] <- data$Recov[i]}
      else
      {deltaD[i] <- data$Death[i]- data$Death[i-1]
      deltaR[i] <-  data$Recov[i]- data$Recov[i-1]}
    }
    data$deltaD <- deltaD
    data$deltaR <- deltaR
    Z <- 0
    for (i in 1:T)
    {if (data$a[i]!=0)
    {Z <- (data$Death[i]/data$a[i])*data$deltaR[i] - (data$Recov[i]/data$a[i])*data$deltaD[i] + Z}
    }
    F <- function(x,ind)
    {s <- 0
    for (i in 1:x) {
      if (data$a[i]!=0 & ind==1){ # ind=1 for death
        s <- (1/data$a[i])*data$deltaD[i]}
      if (data$a[i]!=0 & ind==2){ # ind=2 for recov
        s <- (1/data$a[i])*data$deltaR[i]}
    }
    return(s)
    }
    cum <- 0
    for (i in 1:T)
    {
      {if (data$a[i]!=0 & i > 1)
        cum <- data$deltaR[i]*((data$Death[i-1]/data$a[i])-F(x=T,ind=1)+F(x=i,ind=1))^2 +
          data$deltaD[i]*((data$Recov[i-1]/data$a[i])-F(x=T,ind=2)+F(x=i,ind=2))^2 + cum
      }
    }
    Sig <- sqrt(cum)
    V <- Z/Sig
    p <- 2*(1-pnorm(abs(V)))
    list("Z"=V,"P_valve"=p)
  }
  else return("T should be smaller than or equal to the last day of the pandemic")
}
