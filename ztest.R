z <- dget("fisherz.R")
makedelta <- dget("MakeDeltaFromHypothesis.R")
adfCov <- dget("adfCov.r")
compute4thOrderMoments <- dget("compute4thOrderMoments.r")
findpos <- dget("findpos.r")
FRHO <- dget("FRHO.r")
makecorr <- dget("makecorr.R")

estimationmethod <- 'TSGLS'
datatype <- 'rawdata'
deletion <- 'pairwise'
N <- 50

data <- read.csv(file='data1.csv',head=FALSE,sep=",")
data <- as.matrix(data)

hypothesis <- read.csv(file='hypothesis1.csv',head=FALSE,sep=",")
hypothesis <- as.matrix(hypothesis)

# Produce the correlation matrix using raw data, with pairwise deletion if requested
if (datatype == 'rawdata'){
	N <- nrow(data)
	output <- makecorr(data, deletion)
	data <- output[[1]]
    moments <- output[[2]]
}

delta <- makedelta(hypothesis)

hypothesis_rows <- nrow(hypothesis)
correlations <- c(0)
for (jj in 1:hypothesis_rows) {
	j <- hypothesis[jj,2]
	k <- hypothesis[jj,3]
	correlations[jj] <- data[j,k]
}

rhostar <- hypothesis[,5]

if (is.null(delta)) {
	rhoLS <- rhostar
} else {
	gammaLS <- solve(t(delta)%*%delta)%*%t(delta)%*%(correlations - rhostar)
	rhoLS <- delta%*%gammaLS + rhostar
}

Rlist <- data

if (estimationmethod %in% c('TSGLS', 'TSADF')) { # I think this conditional gives us GLS and TSGLS
	for (jj in 1:hypothesis_rows) {
		j <- hypothesis[jj,2]
		k <- hypothesis[jj,3]
		Rlist[j,k] <- rhoLS[jj]
		Rlist[k,j] <- rhoLS[jj]
	}
}

Psi <- matrix(0, nrow=hypothesis_rows, ncol=hypothesis_rows)
for (jj in 1:hypothesis_rows) {
    for (kk in 1:jj) {
        j <- hypothesis[jj,2]
        k <- hypothesis[jj,3]
        h <- hypothesis[kk,2]
        m <- hypothesis[kk,3]

		if (estimationmethod %in% c('GLS', 'TSGLS')) {
			term1 <- ((Rlist[j,h] - Rlist[j,k]*Rlist[k,h])*(Rlist[k,m] - Rlist[k,h]*Rlist[h,m]))
			term2 <- ((Rlist[j,m] - Rlist[j,h]*Rlist[h,m])*(Rlist[k,h] - Rlist[k,j]*Rlist[j,h]))
			term3 <- ((Rlist[j,h] - Rlist[j,m]*Rlist[m,h])*(Rlist[k,m] - Rlist[k,j]*Rlist[j,m]))
			term4 <- ((Rlist[j,m] - Rlist[j,k]*Rlist[k,m])*(Rlist[k,h] - Rlist[k,m]*Rlist[m,h]))
			Psi[jj,kk] <- 0.5*(term1 + term2 + term3 + term4)
			Psi[kk,jj] <- 0.5*(term1 + term2 + term3 + term4)
		} else {
			term1 <- FRHO(j,k,h,m,moments)
			term2 <- 1/4*Rlist[j,k]*Rlist[h,m]*(FRHO(j,j,h,h,moments) + FRHO(k,k,h,h,moments) + FRHO(j,j,m,m,moments) + FRHO(k,k,m,m,moments))
			term3 <- 1/2*Rlist[j,k]*(FRHO(j,j,h,m,moments) + FRHO(k,k,h,m,moments))
			term4 <- 1/2*Rlist[h,m]*(FRHO(j,k,h,h,moments) + FRHO(j,k,m,m,moments))
        	Psi[jj,kk] <- (term1 + term2 - term3 - term4)
        	Psi[kk,jj] <- (term1 + term2 - term3 - term4)
		}
	}
}

sigmaLS <- Psi

if (is.null(delta)) {
	e <- correlations - rhostar
} else {
	gammaGLS <- solve(t(delta)%*%solve(sigmaLS)%*%delta)%*%(t(delta)%*%solve(sigmaLS)%*%(correlations - rhostar))
	rhoGLS <- delta%*%gammaGLS + rhostar
	e <- z(correlations) - z(rhoGLS)
}

Rlist2 <- data
for (jj in 1:hypothesis_rows) {
	j <- hypothesis[jj,2]
	k <- hypothesis[jj,3]
	Rlist2[j,k] <- rhoGLS[jj]
}

SLS <- matrix(0, nrow=hypothesis_rows, ncol=hypothesis_rows)
for (jj in 1:hypothesis_rows) {
    for (kk in 1:jj) {
        j <- hypothesis[jj,2]
        k <- hypothesis[jj,3]
        h <- hypothesis[kk,2]
        m <- hypothesis[kk,3]
		SLS[jj,kk] <- sigmaLS[jj,kk]/(1 - Rlist2[j,k]^2)*(1 - Rlist2[h,m]^2)
		#SLS[kk,jj] <- sigmaLS[jj,kk]/(1 - Rlist2[j,k]^2)*(1 - Rlist2[h,m]^2) # if this is added then the result is different from Multicor---shouldn't it be symmetric though?
	}
}

X2 <- (N - 3)*t(e)%*%solve(SLS)%*%e

k <- nrow(data)
if (is.null(delta)) {
	q <- 0
} else {
	q <- ncol(delta)
}

p <- pchisq(X2, df=(k-q), lower.tail = FALSE)

cat('X2 =',X2,'\n')
cat('p-value =',p,'\n')
