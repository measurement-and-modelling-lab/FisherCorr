function (j, k, h, m, M) {
    ## Called by adfCov.R
    ## j and k are the row and column number of one correlation from the hypothesis matrix
    ## h and m are the row and column number of another (possibly the same) correlation from the hypothesis matrix
    ## Returns the kurtosis for two the correlations
    
  findpos <- dget("./multicorr/findpos.R")
  temp <- findpos(j, k, h, m)
  fpho <- M[temp]

  return(fpho)
}
