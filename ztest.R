function (data, N, hypothesis, datatype, estimationmethod, deletion) {

	z <- dget("fisherz.R")
	makedelta <- dget("MakeDeltaFromHypothesis.R")
	adfCov <- dget("adfCov.r")
	compute4thOrderMoments <- dget("compute4thOrderMoments.r")
	findpos <- dget("findpos.r")
	FRHO <- dget("FRHO.r")
	makecorr <- dget("makecorr.R")
	tablegen <- dget("tablegen.r")
	errorcheck <- dget("errorcheck.r")
	
	
	if (deletion == 'listwise') { # apply listwise deletion
	    temp1 <- suppressWarnings(as.numeric(data))
	    temp1 <- matrix(temp1, nrow=nrow(data), ncol=ncol(data))
	    data <- temp1[complete.cases(temp1),]
	}

	error <- errorcheck(data, datatype, hypothesis,deletion)
	if (error == TRUE) {
	  return(invisible())
	}

	# Produce the correlation matrix using raw data, with pairwise deletion if requested
	if (datatype == 'rawdata'){
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
	  rhoGLS <- rhostar
		e <- correlations - rhostar
	} else {
		gammaGLS <- solve(t(delta)%*%solve(sigmaLS)%*%delta)%*%(t(delta)%*%solve(sigmaLS)%*%(correlations - rhostar))
		rhoGLS <- delta%*%gammaGLS + rhostar
		e <- z(correlations) - z(rhoGLS)
		
		# variance of parameter tag estimates---differs from WBCORR but that might be correct
		nMatrix <- diag(rep(N, hypothesis_rows))
		OmegaHatInverse <- sqrt(nMatrix)%*%solve(Psi)%*%sqrt(nMatrix)
		Psi <- t(delta)%*%OmegaHatInverse
		covgamma <- solve(Psi%*%delta)
		covgamma <- round(covgamma, 3)
		
		# confidence interval on parameter tag estimates
		gammaGLS_ci <- lapply(gammaGLS, function(x) {
		  UL <- z(x) + 1.96*sqrt(1/(N-3))
		  UL <- tanh(UL)
		  UL <- round(UL, 3)
		  LL <- z(x) - 1.96*sqrt(1/(N-3))
		  LL <- tanh(LL)
		  LL <- round(LL, 3)
		  return(paste0('[',LL,', ',UL,']'))
		})
		gammaGLS_ci <- unlist(gammaGLS_ci)
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
	X2 <- round(X2, 3)

	k <- nrow(data)
	if (is.null(delta)) {
		q <- 0
	} else {
		q <- ncol(delta)
	}

	p <- pchisq(X2, df=(k-q), lower.tail = FALSE)
	p <- round(p, 3)


	printfunction <- function () {
	  
	  cat('<br><div style="line-height: 175%; margin-left:15px"><b>Hypothesis Matrix</b></div>', sep="")
	  hypothesis <- rbind(c("Group", "Row", "Column", "Parameter Tag", "Fixed Value"), hypothesis)
	  tablegen(hypothesis,TRUE)

	  cat('<br><div style="line-height: 175%; margin-left:15px"><b>Input Correlation Matrix (N=', N, ')</b></div>', sep="")
	  data <- round(data, 3)
	  tablegen(data,FALSE)

	  if (!(is.null(delta))) {
	    cat('<br><div style="line-height: 175%; margin-left:15px"><b>', estimationmethod, 'Parameter Estimates</b></div>')
	    title <- c('Parameter Tag', 'Estimate', 'Std. Error', '95% Confidence Interval')
	    gammaGLS <- round(gammaGLS, 3)
	    numbers <- c(1:length(gammaGLS))
	    gammaGLS <- cbind(numbers, gammaGLS, covgamma, gammaGLS_ci)
	    gammaGLS <- rbind(title, gammaGLS)
	    tablegen(gammaGLS,TRUE)
	  }
	  
	  cat('<br><div style="line-height: 175%; margin-left:15px"><b>Significance Test Results</b></div>', sep="")
	  results <- matrix(c('Chi Square', X2, '&nbsp;&nbsp;df', k-q, '&nbsp;&nbsp;&nbsp;&nbsp;Sig.', p), nrow=2, ncol=3)
	  tablegen(results,TRUE)
	}
	
	x <- capture.output(printfunction())
	cat(x)

}
