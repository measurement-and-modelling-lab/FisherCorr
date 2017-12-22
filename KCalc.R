k <- 6
n <- 10


k2Calc <- function(k,n){
  a <- (n-4)/2
  b <- (n-3)/2
  sum <- 0
  ans <- 0
  if(k%%2 == 0){
    for (i in 1:a){
      sum <- sum + k^(-2)
    }
    sum <- sum*(1/2)
    ans <- ((pi^2)/12) - sum
    return(ans)
  }
  else{
    for (i in 1:b){
      sum <- sum + (2*k-1)^(-2)
    }
    sum <- sum*2
    ans <- ((pi^2)/4) - sum
    return(ans)
    
  }
  
}

k4Calc <- function(k,n){
  a <- (n-4)/2
  b <- (n-3)/2
  sum <- 0
  ans <- 0
  if(k%%2 == 0){
    for (i in 1:a){
      sum <- sum + k^(-4)
    }
    sum <- sum*(3/4)
    ans <- ((pi^4)/120) - sum
    return(ans)
  }
  else{
    for (i in 1:b){
      sum <- sum + (2*k-1)^(-4)
    }
    sum <- sum*12
    ans <- ((pi^4)/8) - sum
    return(ans)
    
  }
  
}

k6Calc <- function(k,n){
  a <- (n-4)/2
  b <- (n-3)/2
  sum <- 0
  ans <- 0
  if(k%%2 == 0){
    for (i in 1:a){
      sum <- sum + k^(-6)
    }
    sum <- sum*(15/4)
    ans <- ((pi^6)/252) - sum
    return(ans)
  }
  else{
    for (i in 1:b){
      sum <- sum + (2*k-1)^(-6)
    }
    sum <- sum*240
    ans <- ((pi^6)/4) - sum
    return(ans)
    
  }
  
}

k8Calc <- function(k,n){
  a <- (n-4)/2
  b <- (n-3)/2
  sum <- 0
  ans <- 0
  if(k%%2 == 0){
    for (i in 1:a){
      sum <- sum + k^(-8)
    }
    sum <- sum*(315/8)
    ans <- ((pi^8)/240) - sum
    return(ans)
  }
  else{
    for (i in 1:b){
      sum <- sum + (2*k-1)^(-8)
    }
    sum <- sum*10080
    ans <- (17*pi^8)/16 - sum
    return(ans)
    
  }
  
}





