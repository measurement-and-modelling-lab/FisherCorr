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
matrix <- matrix(c(.1,.2,.3,.4,
                   .4,.5,.6,.5,
                   .1,.2,.4,.6,
                   .7,.8,.8,.8), nrow = 4, ncol = 4)  
n <- 20

print(sCalc(matrix,n))