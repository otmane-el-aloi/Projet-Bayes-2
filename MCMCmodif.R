MHWeibull_REffect <- function(data, iter = 5 * 10^4,
                              prop.sd = c(1,1,0.5,0.5,0.5,0.5,0.1,0.3, rep(1,38))){
  
  ## Data is n x 9 matrix with time, age, status,sex,desease, censoring indicator
  time1 <-data[,1]
  time2<-data[,2]
  cens1 <- data[,3]
  cens2 <-data[,4]
  age1<-data[,5]
  age2<-data[,6]
  sex<-data[,7]
  d1<-data[,8]
  d2<-data[,9]
  d3<-data[,10]
  
  ## Define a "suitable" initial state for the chain
  init <- c(0,0,0,0,0,0,1,1, rep(0,38))
  acc.rates <- rep(0, 8+38)
  chain <- matrix(NA, iter + 1, 8+38)
  chain[1,] <- init
  
  
  for (i in 1:iter){
    current <- chain[i,]
    
    #Modifying alpha, beta_age,beta_sex,beta_1,beta_2,beta_3
    for (j in 1:6){
      prop<-current
      #Proposition
      prop[j]<-rnorm(1,current[j],prop.sd[j])

      #The top 
      scale1=exp(prop[1]+prop[2]*age1+prop[3]*sex+prop[4]*d1+prop[5]*d3+prop[6]*d3+prop[9:46])
      scale2=exp(prop[1]+prop[2]*age2+prop[3]*sex+prop[4]*d1+prop[5]*d2+prop[6]*d3+prop[9:46])   
      shape=rep(prop[8],38)

      top <- dnorm(prop[j],0,100, log = TRUE) + sum(dweibull(time1,shape,scale1, log=TRUE)+dweibulll(time2,shape,scale2, log=TRUE)+__cens1__pweibull(time1,shape,scale1, log=TRUE, lower.tail = FALSE)+__cens2__pweibull(time2,shape,scale2, log=TRUE, lower.tail = FALSE))
      
      #The bottom
      scale1=exp(current[1]+current[2]*age1+current[3]*sex+current[4]*d1+current[5]*d2+current[6]*d3+curent[9:46])
      scale2=exp(current[1]+current[2]*age2+current[3]*sex+current[4]*d1+current[5]*d2+current[6]*d3+curent[9:46])
      shape=rep(current[8],38)
      
      bottom <- dnorm(current[j],0,100, log = TRUE) + sum(dweibull(time1,shape,scale1, log=TRUE)+dweibulll(time2,shape,scale2, log=TRUE)+__cens1__pweibull(time1,shape,scale1, log=TRUE, lower.tail = FALSE)+__cens2__pweibull(time2,shape,scale2, log=TRUE, lower.tail = FALSE))
      
      alpha <- exp(top-bottom)
      
      if (runif(1) < alpha){
        current <- prop
        acc.rates[j] <- acc.rates[j] + 1
      }
    }

    #Modfifying sigma2 (conjugate prior)
    current[7] <- 1/rgamma(1,10^(-3)+38/2, rate=10^(-3)+sum(current[9:46]^2)/2)
  
    
    #Modfifying r
    #Proposition
    prop[8]<-rnorm(1,current[8],prop.sd[8])

    #The top 
    scale1=exp(prop[1]+prop[2]*age1+prop[3]*sex+prop[4]*d1+prop[5]*d2+prop[6]*d3+prop[9:46])
    scale2=exp(prop[1]+prop[2]*age2+prop[3]*sex+prop[4]*d1+prop[5]*d2+prop[6]*d3+prop[9:46])
    shape=rep(prop[8],38)
    
    top <- dgamma(prop[8],1,10^(3), log = TRUE) + sum(dweibull(time1,shape,scale1, log=TRUE)+dweibulll(time2,shape,scale2, log=TRUE)+__cens1__pweibull(time1,shape,scale1, log=TRUE, lower.tail = FALSE)+__cens2__pweibull(time2,shape,scale2, log=TRUE, lower.tail = FALSE))
    
    #The bottom
    scale1=exp(current[1]+current[2]*age1+current[3]*sex+current[4]*d1+current[5]*d2+current[6]*d3+current[9:46])
    scale2=exp(current[1]+current[2]*age2+current[3]*sex+current[4]*d1+current[5]*d2+current[6]*d3+current[9:46])
    shape=rep(current[8],38)
    
    bottom <- dgamma(current[8],1,10^3, log = TRUE) + sum(dweibull(time1,shape,scale1, log=TRUE)+dweibulll(time2,shape,scale2, log=TRUE)+__cens1__pweibull(time1,shape,scale1, log=TRUE, lower.tail = FALSE)+__cens2__pweibull(time2,shape,scale2, log=TRUE, lower.tail = FALSE))
  
    alpha <- exp(top-bottom)
  
    if (runif(1) < alpha){
      current <- prop
      acc.rates[8] <- acc.rates[8] + 1
    }
  
    #Modifying bi (on s'intéresse seulement au patient i)
    for (j in 9:46){
      prop<-current
      #Proposition
      prop[j]<-rnorm(1,current[j],current[7])
      
      #The top 
      scale1=exp(prop[1]+prop[2]*age1[j-8]+prop[3]*sex[j-8]+prop[4]*d1[j-8]+prop[5]*d3[j-8]+prop[6]*d3[j-8]+prop[j])
      scale2=exp(prop[1]+prop[2]*age2[j-8]+prop[3]*sex[j-8]+prop[4]*d1[j-8]+prop[5]*d2[j-8]+prop[6]*d3[j-8]+prop[j])   
      shape=rep(prop[8],2)
      
      top <- dnorm(prop[j],0,100, log = TRUE) + sum(dweibull(time1[j-8],shape,scale1, log=TRUE)+dweibulll(time2[j-8],shape,scale2, log=TRUE)+__cens1__pweibull(time1[j-8],shape,scale1, log=TRUE, lower.tail = FALSE)+__cens2__pweibull(time2[j-8],shape,scale2, log=TRUE, lower.tail = FALSE))
      
      #The bottom
      scale1=exp(current[1]+current[2]*age1[j-8]+current[3]*sex[j-8]+current[4]*d1[j-8]+current[5]*d3[j-8]+current[6]*d3[j-8]+current[j])
      scale2=exp(current[1]+current[2]*age2[j-8]+current[3]*sex[j-8]+current[4]*d1[j-8]+current[5]*d2[j-8]+current[6]*d3[j-8]+current[j])
      shape=rep(current[8],2)
      
      bottom <- dnorm(prop[j],0,100, log = TRUE) + sum(dweibull(time1[j-8],shape,scale1, log=TRUE)+dweibulll(time2[j-8],shape,scale2, log=TRUE)+__cens1__pweibull(time1[j-8],shape,scale1, log=TRUE, lower.tail = FALSE)+__cens2__pweibull(time2[j-8],shape,scale2, log=TRUE, lower.tail = FALSE))
      
      alpha <- exp(top-bottom)
      
      if (runif(1) < alpha){
        current <- prop
        acc.rates[j] <- acc.rates[j] + 1
      }
    }
    
    chain[i+1,] <- current
  }
  return(list(chain = chain, acc.rates = acc.rates / iter))
}