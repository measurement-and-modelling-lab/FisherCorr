function(j, k, h, m) {
    ## Called by FRHO.R
    ## j and k are the row and column number of one correlation from the hypothesis matrix
    ## h and m are the row and column number of another (possibly the same) correlation from the hypothesis matrix
    ## Returns the index number of the kurtosis for the two correlations in the moments list

    ordered.indices <- sort(c(j, k, h, m), decreasing = TRUE)

    a <- ordered.indices[1]
    b <- ordered.indices[2]
    c <- ordered.indices[3]
    d <- ordered.indices[4]

    post <- (a-1)*a*(a+1)*(a+2)/24 + (b-1)*b*(b+1)/6 + c*(c-1)/2 + d

    return(post)
}
