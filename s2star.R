#Hotelling series expansion formulae
s2star <- function(rho,n,k,S){
  
  v1Calc <- function(rho, n){
    x1 <- rho/2 + ((5*rho + rho^3) / (8*(n-1)))
    x2 <- (1*rho + 2*rho^3 +3*rho^5) / (16*(n-1)^2)
    x3 <- (83* rho + 13*rho^3 - 27* rho^5+ 75* rho^7) / (128*(n-1)^3)
    x4 <- (143*rho + 20*rho^3 +138*rho^5 - 780*rho^7 + 735*rho^9) / (265*(n-1)^4)
    x5 <- (625*rho + 113*rho^3 - 990*rho^5 + 14250*rho^7 - 33075*rho^9 + 19845*rho^11) / (1024*(n-1)^5)
    
    v1 = 1/(n-1) * (x1+x2+x3+x4+x5)
    return(v1)
  }
  
  v2Calc <- function(rho, n){
    x1 <- (8-rho^2)/(4*(n-1))
    x2 <- (88- 9*rho^2 -9*rho^4)/(24*(n-1)^2)
    x3 <- (384 -19*rho^2 + 2* rho^4 - 75* rho^6)/(64*(n-1)^3)
    x4 <- (16256 - 225* rho^2 - 375* rho^4 + 8025*rho^6 -11025*rho^8)/(1920*(n-1)^4)
    x5 <- (5120 +113*rho^2 +20*rho^4 - 7194*rho^6 + 26460*rho^8- 19845*rho^10)/(512*(n-1)^5)
    
    v2 <- 1/(n-1) * (1 + x1+x2+x3+x4+x5)
    return(v2)
  }
  v3Calc <- function(rho, n){
    x1 <- (3*rho)/(2*(n-1))
    x2 <- (39* rho + 6*rho^3)/(8*(n-1)^2)
    x3 <- (362*rho +45*rho^3 +63*rho^5)/(32*(n-1)^3)
    x4 <- (2809*rho +296*rho^3 - 327*rho^5 +1170*rho^7)/(128*(n-1)^4)
    x5 <- (47461*rho +5720*rho^3 + 11190*rho^5-73350*rho^7 + 77175*rho^9)/(1280*(n-1)^5)
    
    v3 <- 1/(n-1) * (x1+x2+x3+x4+x5)
    return(v3)
  }
  v4Calc <- function(rho, n){
    x1 <- 3/(n-1)
    x2 <- (28 - 3 *rho^3)/(2*(n-1)^2)
    x3 <- (736 - 84* rho^2 - 51* rho^4)/(16*(n-1)^3)
    x4 <- (31744- 3016*rho^2 - 864*rho^4 - 3480*rho^6)/(256*(n-1)^4)
    x5 <- (543616 - 41310* rho^2 - 17595*rho^4 + 87300* rho^6 - 165375*rho^8)/(1920*(n-1)^5)
    
    v4 <- 1/(n-1) * (x1+x2+x3+x4+x5)
    return(v4)
  }
  
  v1 <- v1Calc(rho,n)
  v2 <- v2Calc(rho,n)
  v3 <- v3Calc(rho,n)
  v4 <- v4Calc(rho,n)
  
  aStar <- (1/(n-3)) * sqrt((2/(v4-(v2)^2)))
  bStar <- (k*(k-1)/2)*(1-(aStar*(n-3)*v2))
  
  return(aStar*S+bStar)
}
