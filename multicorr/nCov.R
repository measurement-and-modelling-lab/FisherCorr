function (j, k, h, m, R) {
    ## j and k are the row and column number of one correlation from the hypothesis matrix
    ## h and m are the row and column number of another (possibly the same) correlation from the hypothesis matrix
    ## R is the correlation matrix if using ADF and the OLS matrix if using TSADF
    ## output is the covariance of the two correlations

    term1 <- ((R[j,h] - R[j,k]*R[k,h])*(R[k,m] - R[k,h]*R[h,m]))
    term2 <- ((R[j,m] - R[j,h]*R[h,m])*(R[k,h] - R[k,j]*R[j,h]))
    term3 <- ((R[j,h] - R[j,m]*R[m,h])*(R[k,m] - R[k,j]*R[j,m]))
    term4 <- ((R[j,m] - R[j,k]*R[k,m])*(R[k,h] - R[k,m]*R[m,h]))
    Cov <- (0.5)*(term1 + term2 + term3 + term4)

    return(Cov)
}
