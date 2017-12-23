z <- dget("fisherz.R")
makesymmetric <- dget("MakeSymmetricMatrix.R")
makedelta <- dget("MakeDeltaFromHypothesis.R")

data <- matrix(c(
1,0.08,0.03,
0.08,1,0.11,
0.03,0.11,1), nrow=3, ncol=3)

hypothesis <- matrix(c(1, 1, 1,
                       2, 3, 3,
					   1, 1, 2,
					   1, 1, 1,
					   0, 0, 0), nrow=3, ncol=5)

N = 50

delta <- makedelta(hypothesis)
k <- nrow(delta)
q <- ncol(delta)

rows <- nrow(hypothesis)
correlations <- c(0)
for (jj in 1:rows) {
	j <- hypothesis[jj,2]
	k <- hypothesis[jj,3]
	correlations[jj] <- data[j,k]
}

rhostar <- hypothesis[,5]
gammaLS <- solve(t(delta)%*%delta)%*%t(delta)%*%(correlations - rhostar)
rhoLS <- delta%*%gammaLS + rhostar
Rlist <- data

rows <- nrow(hypothesis)
for (jj in 1:rows) {
	j <- hypothesis[jj,2]
	k <- hypothesis[jj,3]
	Rlist[j,k] <- rhoLS[jj]
	Rlist[k,j] <- rhoLS[jj]
}

Psi <- matrix(0, nrow=rows, ncol=rows)
for (jj in 1:rows) {
    for (kk in 1:jj) {
        j <- hypothesis[jj,2]
        k <- hypothesis[jj,3]
        h <- hypothesis[kk,2]
        m <- hypothesis[kk,3]

  			term1 <- ((Rlist[j,h] - Rlist[j,k]*Rlist[k,h])*(Rlist[k,m] - Rlist[k,h]*Rlist[h,m]))
  			term2 <- ((Rlist[j,m] - Rlist[j,h]*Rlist[h,m])*(Rlist[k,h] - Rlist[k,j]*Rlist[j,h]))
  			term3 <- ((Rlist[j,h] - Rlist[j,m]*Rlist[m,h])*(Rlist[k,m] - Rlist[k,j]*Rlist[j,m]))
  			term4 <- ((Rlist[j,m] - Rlist[j,k]*Rlist[k,m])*(Rlist[k,h] - Rlist[k,m]*Rlist[m,h]))
  			Psi[jj,kk] <- 0.5*(term1 + term2 + term3 + term4)
  			Psi[kk,jj] <- 0.5*(term1 + term2 + term3 + term4)
	}
}

sigmaLS <- Psi

gammaGLS <- solve(t(delta)%*%solve(sigmaLS)%*%delta)%*%(t(delta)%*%solve(sigmaLS)%*%(correlations - rhostar))

rhoGLS <- delta%*%gammaGLS + rhostar

Rlist2 <- data
rows <- nrow(hypothesis)
for (jj in 1:rows) {
	j <- hypothesis[jj,2]
	k <- hypothesis[jj,3]
	Rlist2[j,k] <- rhoGLS[jj]
}

SLS <- matrix(0, nrow=rows, ncol=rows)
for (jj in 1:rows) {
    for (kk in 1:jj) {
        j <- hypothesis[jj,2]
        k <- hypothesis[jj,3]
        h <- hypothesis[kk,2]
        m <- hypothesis[kk,3]
		SLS[jj,kk] <- sigmaLS[jj,kk]/(1 - Rlist2[j,k]^2)*(1 - Rlist2[h,m]^2)
	}
}

X2 <- (N - 3)*t((z(correlations) - z(rhoGLS)))%*%solve(SLS)%*%(z(correlations) - z(rhoGLS))
cat('X2 =',X2,'\n')

p <- dchisq(X2, df=(k-q))
cat('p-value =',p,'\n')
