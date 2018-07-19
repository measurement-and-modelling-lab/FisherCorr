function (data, N, hypothesis, datatype, estimationmethod, deletion) {

    ## Import functions
    fisherTransform <- dget("./multicorr/fisherz.R")
    adfCov <- dget("./multicorr/adfCov.R")
    compute4thOrderMoments <- dget("./multicorr/compute4thOrderMoments.R")
    findpos <- dget("./multicorr/findpos.R")
    FRHO <- dget("./multicorr/FRHO.R")
    errorcheck <- dget("./multicorr/errorcheck.R")
    assess_range <- dget("./multicorr/assess_range.R")
    assess_mvn <- dget("./multicorr/assess_mvn.R")
    MultivariateSK <- dget("./multicorr/MultivariateSK.R")


    ## If the upper triangle of a correlation matrix is empty, make the matrix symmetric
    ## Otherwise, check whether the matrix is symmetric and if so return an error
    if (datatype == "correlation") {
        upper.triangle <- data[upper.tri(data)]
        lower.triangle <- data[lower.tri(data)]
        if (all(is.na(upper.triangle))) {
            data[upper.tri(data)] <- lower.triangle
        } else if (!all(lower.triangle == upper.triangle)) {
            stop("Correlation matrix is not symmetric.")
        }
    }


    ## Apply listwise deletion
    if (deletion == 'listwise') { ## apply listwise deletion
        temp1 <- suppressWarnings(as.numeric(data))
        temp1 <- matrix(temp1, nrow=nrow(data), ncol=ncol(data))
        data <- temp1[complete.cases(temp1),]
        N <- nrow(data)
    }


    ## Correct N if deletion is pairwise
    if (deletion == 'pairwise') {
        ns <- colSums(!is.na(data))
        hmean <- 1/mean(1/ns)
        hmean <- round(hmean, 1)
        N <- hmean
    }

    ## Error checking
    errorcheck(data, datatype, hypothesis,deletion)


    ## Renumber parameter tags if a number is skipped
    parameter.tags <- hypothesis[hypothesis[,4] != 0, 4]
    if (max(parameter.tags) > length(unique(parameter.tags))) {
        hypothesis[hypothesis[,4] != 0, 4] <- as.numeric(as.factor(parameter.tags))
    }


    ## Assess multivariate normality using Yuan, Lambert & Fouladi (2004) if using pairwise deletion, Mardia (1970) otherwise
    if (datatype == 'rawdata') {
        if (deletion == 'pairwise') {
            MardiaSK <- list(assess_range(list(data)), assess_mvn(list(data)))
        } else {
            MardiaSK <- MultivariateSK(list(data))
        }
    } else {
        MardiaSK <- NA
    }


    ## Produce the correlation matrix using raw data, with pairwise deletion if requested
    if (datatype == 'rawdata'){
        temp <- suppressWarnings(as.numeric(data)) ## Convert anything that can't be interpreted as a number to NA
        if (deletion == "pairwise" & !(NA %in% temp)) {
            stop("You can't use pairwise deletion without missing data.")
        }
        R <- cor(data, use="pairwise")
        moments <- compute4thOrderMoments(data)
    } else {
        R <- data
    }


    ## Check that the correlation matrices are positive definite
    eigen_values <- eigen(R)[[1]]
    if (TRUE %in% (eigen_values <= 0)) {
        stop('Data matrix is not positive definite.')
    }


    ## Check if there are any parameter tags; if there are, create delta
    delta <- sapply(1:max(hypothesis[,4]), function(p) {
        as.numeric(p == hypothesis[,4])
    })
    no.parameters <- max(delta) == 0


    ## getVecR
    hypothesis_rows <- nrow(hypothesis)
    correlations <- c(0)
    for (jj in 1:hypothesis_rows) {
        j <- hypothesis[jj,2]
        k <- hypothesis[jj,3]
        correlations[jj] <- R[j,k]
    }

    
    ## Create a vector of fixed values
    hypothesis[hypothesis[,4] != 0, 5] <- 0
    rhostar <- hypothesis[,5]

    if (is.null(delta)) {
        rhoLS <- rhostar
    } else {
        gammaLS <- solve(t(delta)%*%delta)%*%t(delta)%*%(correlations - rhostar)
        rhoLS <- delta%*%gammaLS + rhostar
    }


    ## Create OLS matrix
    R.OLS <- R
    if (estimationmethod %in% c('TSGLS', 'TSADF')) { ## I think this conditional gives us GLS and TSGLS
        for (jj in 1:hypothesis_rows) {
            j <- hypothesis[jj,2]
            k <- hypothesis[jj,3]
            R.OLS[j,k] <- rhoLS[jj]
            R.OLS[k,j] <- rhoLS[jj]
        }
    }


    ## Calculate correlation covariance
    Psi <- matrix(0, nrow=hypothesis_rows, ncol=hypothesis_rows)
    for (jj in 1:hypothesis_rows) {
        for (kk in 1:jj) {
            j <- hypothesis[jj,2]
            k <- hypothesis[jj,3]
            h <- hypothesis[kk,2]
            m <- hypothesis[kk,3]

            if (estimationmethod %in% c('GLS', 'TSGLS')) {
                term1 <- ((R.OLS[j,h] - R.OLS[j,k]*R.OLS[k,h])*(R.OLS[k,m] - R.OLS[k,h]*R.OLS[h,m]))
                term2 <- ((R.OLS[j,m] - R.OLS[j,h]*R.OLS[h,m])*(R.OLS[k,h] - R.OLS[k,j]*R.OLS[j,h]))
                term3 <- ((R.OLS[j,h] - R.OLS[j,m]*R.OLS[m,h])*(R.OLS[k,m] - R.OLS[k,j]*R.OLS[j,m]))
                term4 <- ((R.OLS[j,m] - R.OLS[j,k]*R.OLS[k,m])*(R.OLS[k,h] - R.OLS[k,m]*R.OLS[m,h]))
                Psi[jj,kk] <- 0.5*(term1 + term2 + term3 + term4)
                Psi[kk,jj] <- 0.5*(term1 + term2 + term3 + term4)
            } else {
                term1 <- FRHO(j,k,h,m,moments)
                term2 <- 1/4*R.OLS[j,k]*R.OLS[h,m]*(FRHO(j,j,h,h,moments) + FRHO(k,k,h,h,moments) + FRHO(j,j,m,m,moments) + FRHO(k,k,m,m,moments))
                term3 <- 1/2*R.OLS[j,k]*(FRHO(j,j,h,m,moments) + FRHO(k,k,h,m,moments))
                term4 <- 1/2*R.OLS[h,m]*(FRHO(j,k,h,h,moments) + FRHO(j,k,m,m,moments))
                Psi[jj,kk] <- (term1 + term2 - term3 - term4)
                Psi[kk,jj] <- (term1 + term2 - term3 - term4)
            }
        }
    }
    
    sigmaLS <- Psi

    parameters <- unique(hypothesis[, 4])
    parameters <- parameters[parameters > 0]
    parameters.length <- length(parameters)

    if (no.parameters) {
        rhoGLS <- rhostar
        e <- correlations - rhostar
        estimates.table <- NA
    } else {
        gammaGLS <- solve(t(delta)%*%solve(sigmaLS)%*%delta)%*%(t(delta)%*%solve(sigmaLS)%*%(correlations - rhostar))
        rhoGLS <- delta%*%gammaGLS + rhostar
        e <- fisherTransform(correlations) - fisherTransform(rhoGLS)
        
        ## variance of parameter tag estimates---differs from WBCORR but that might be correct
        nMatrix <- diag(rep(N, hypothesis_rows))
        OmegaHatInverse <- sqrt(nMatrix)%*%solve(Psi)%*%sqrt(nMatrix)
        Psi <- t(delta)%*%OmegaHatInverse
        covgamma <- solve(Psi%*%delta)
        covgamma <- covgamma[,1]

        ## Confidence interval calculation
        gammaGLS_ci <- c()
        for (i in 1:parameters.length) {
            parameter <- parameters[i]
            corrected_alpha <- 0.05/parameters.length
            critical_value <- qnorm(1-corrected_alpha/2)
            point.estimate <- gammaGLS[i]

            UL <- fisherTransform(point.estimate) + critical_value*sqrt(1/(N-3))
            UL <- tanh(UL)
            LL <- fisherTransform(point.estimate) - critical_value*sqrt(1/(N-3))
            LL <- tanh(LL)


            UL <- round(UL, 3)
            LL <- round(LL, 3)
            gammaGLS_ci[i] <- paste0('[',LL,', ',UL,']')

        }
            covgamma <- round(covgamma, 3)
            gammaGLS <- round(gammaGLS, 3)
            estimates.table <- cbind(parameters, gammaGLS, covgamma, gammaGLS_ci)
            colnames(estimates.table) <- c("Parameter Tags", "Point Estimate", "Std. Error", "Confidence Interval")
    }

    R.GLS <- R
    for (jj in 1:hypothesis_rows) {
        j <- hypothesis[jj,2]
        k <- hypothesis[jj,3]
        R.GLS[j,k] <- rhoGLS[jj]
        R.GLS[k,j] <- rhoGLS[jj]
    }

    SLS <- matrix(0, nrow=hypothesis_rows, ncol=hypothesis_rows)
    for (jj in 1:hypothesis_rows) {
        for (kk in 1:jj) {
            j <- hypothesis[jj,2]
            k <- hypothesis[jj,3]
            h <- hypothesis[kk,2]
            m <- hypothesis[kk,3]
            SLS[jj,kk] <- sigmaLS[jj,kk]/(1 - R.GLS[j,k]^2)*(1 - R.GLS[h,m]^2)
            SLS[kk,jj] <- sigmaLS[jj,kk]/(1 - R.GLS[j,k]^2)*(1 - R.GLS[h,m]^2) ## Adding this makes the output differ from Multicorr; probably a bug in the original
        }
    }


    ## Calculate test statistic and p value
    chisquare <- (N - 3)*t(e)%*%solve(SLS)%*%e
    k <- ncol(R)
    q <- parameters.length
    p <- pchisq(chisquare, df=(k-q), lower.tail = FALSE)
    p <- round(p, 3)
    if (p == 0) {
        p <- '< .001'
    }
    sigtable <- matrix(c(chisquare, k-q, p), nrow=1)
    sigtable <- round(sigtable, 3)
    colnames(sigtable) <- c("Chi Square", "df", "pvalue")
    
    ## Test whether the hypothesis is the identity hypothesis
    ## If so, run superior S test
    test.matrix <- R
    for (i in 1:nrow(hypothesis)) {
        if (hypothesis[i,4] == 0) {
            row <- hypothesis[i,2]
            col <- hypothesis[i,3]
            test.matrix[row,col] <- hypothesis[i,5]
            test.matrix[col,row] <- hypothesis[i,5]
        }
    }
    identity.matrix <- diag(k)
    if (all(test.matrix == identity.matrix)) {
        identity <- TRUE
        scalc <- dget("scalc.R")
        s2star <- dget("s2star.R")
        S <- scalc(data, N)
        S <- s2star(0, N, nrow(data), S)
        S.p <- pchisq(S, df=k, lower.tail=FALSE)
        S <- round(S, 3)
        S.p <- round(Sp, 3)
        S.result <- c(S, k, S.p)
    } else {
        S.result <- NA
    }

    output <- list(hypothesis, N, R, R.OLS, estimates.table, sigtable, R.GLS, S.result, MardiaSK)
    return(output)

}
