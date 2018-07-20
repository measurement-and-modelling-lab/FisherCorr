function(matrix,n){
    
    fisherz <- dget("./multicorr/fisherz.R")

    sum <- 0
    for (i in 1:nrow(matrix))
        for (j in 1:ncol(matrix)) {
            if (j < i) {
                sum <- sum + (fisherz(matrix[[i,j]]))^2
            }
        }
    S <- (n-3) * sum

    return (S)
}
