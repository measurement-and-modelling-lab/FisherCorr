function (j, k, h, m, R, moments) {
    ## j and k are the row and column number of first correlation from the hypothesis matrix
    ## h and m are the row and column number of second correlation from the hypothesis matrix
    ## R is the correlation matrix if using ADF and the OLS matrix if using TSADF
    ## moments are the fourth order moments, i.e. values on kurtosis
    ## output is the covariance of the two correlations

    FRHO <- dget("./multicorr/FRHO.R")

    term1 <- FRHO(j,k,h,m,moments)
    term2 <- (1/4)*R[j,k]*R[h,m]*(FRHO(j,j,h,h,moments) + FRHO(k,k,h,h,moments) + FRHO(j,j,m,m,moments) + FRHO(k,k,m,m,moments))
    term3 <- (1/2)*R[j,k]*(FRHO(j,j,h,m,moments) + FRHO(k,k,h,m,moments))
    term4 <- (1/2)*R[h,m]*(FRHO(j,k,h,h,moments) + FRHO(j,k,m,m,moments))
    Cov <- (term1 + term2 - term3 - term4)

    return(Cov)
}
