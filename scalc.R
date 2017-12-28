fisherz <- dget("fisherz.R")
scalc <- function(matrix,n){
                               
sum <- 0

for (i in 1:nrow(matrix))
  for (j in 1:ncol(matrix)) {
    if (j < i) {
      sum <- sum + (fisherz(matrix[[i,j]]))^2
    }
  }

return ((n-3)*sum)
}
