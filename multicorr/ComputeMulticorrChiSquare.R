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
    adfCov <- dget("./multicorr/adfCov.R")
    nCov <- dget("./multicorr/nCov.R")


    ## If the upper triangle of a correlation matrix is empty, make the matrix symmetric
    ## Otherwise, check whether the matrix is symmetric and if not return an error
    if (datatype == "correlation") {
        current.upper.triangle <- data[upper.tri(data)]
        symmetric.upper.triangle <- t(data)[upper.tri(t(data))]
        if (all(is.na(current.upper.triangle))) {
            data[upper.tri(data)] <- symmetric.upper.triangle
        } else if (!all(current.upper.triangle == symmetric.upper.triangle)) {
            stop("Correlation matrix is not symmetric.")
        }
    }


    ## Error checking
    errorcheck(data, datatype, hypothesis, deletion, N)


    ## Renumber parameter tags if a number is skipped
    parameter.tags <- hypothesis[hypothesis[,4] != 0, 4]
    if (length(parameter.tags) > 0) {
        if (max(parameter.tags) > length(unique(parameter.tags))) {
            hypothesis[hypothesis[,4] != 0, 4] <- as.numeric(as.factor(parameter.tags))
        }
    }


    ## Apply listwise deletion
    if (deletion == 'listwise') { ## apply listwise deletion
        data <- data[complete.cases(data),]
        N <- nrow(data)
    }


    ## Correct N if deletion is pairwise
    if (deletion == 'pairwise') {
        ns <- colSums(!is.na(data))
        hmean <- 1/mean(1/ns)
        hmean <- round(hmean, 1)
        N <- hmean
    }


    ## Assess multivariate normality using Yuan, Lambert & Fouladi (2004) if using pairwise deletion, Mardia (1970) otherwise
    if (datatype == 'rawdata') {
        if (deletion == 'pairwise') {
            temp <- assess_range(list(data))
            MardiaSK <- list(temp[[1]], assess_mvn(list(data)))
            missing <- temp[[2]]
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
    eigen.values <- eigen(R)[[1]]
    if (TRUE %in% (eigen.values <= 0)) {
        stop('Data matrix is not positive definite.')
    }


    ## Check if there are any parameter tags; if there are, create delta
    delta <- sapply(1:max(hypothesis[,4]), function(p) {
        as.numeric(p == hypothesis[,4])
    })
    no.parameters <- max(hypothesis[,4]) == 0


    ## Create a vector of the correlations referenced in the hypothesis
    hypothesis.length <- nrow(hypothesis)
    correlations <- c(0)
    for (jj in 1:hypothesis.length) {
        j <- hypothesis[jj,2]
        k <- hypothesis[jj,3]
        correlations[jj] <- R[j,k]
    }


    ## Create a column matrix of fixed values
    hypothesis[hypothesis[,4] != 0, 5] <- 0
    rhostar <- hypothesis[,5,drop=FALSE]


    ## Create column matrix of least squares estimates
    if (no.parameters) {
        rhoLS <- rhostar
    } else {
        gammaLS <- solve(t(delta)%*%delta)%*%t(delta)%*%(correlations - rhostar)
        rhoLS <- delta%*%gammaLS + rhostar
    }


    ## Create OLS estimates
    R.OLS <- R
    for (jj in 1:hypothesis.length) {
        j <- hypothesis[jj,2]
        k <- hypothesis[jj,3]
        R.OLS[j,k] <- rhoLS[jj]
        R.OLS[k,j] <- rhoLS[jj]
    }


    ## Calculate correlation covariance for each pair of correlations from the hypothesis matrix
    Psi <- matrix(0, nrow=hypothesis.length, ncol=hypothesis.length)
    for (jj in 1:hypothesis.length) {
        for (kk in 1:jj) {
            j <- hypothesis[jj,2]
            k <- hypothesis[jj,3]
            h <- hypothesis[kk,2]
            m <- hypothesis[kk,3]

            if (estimationmethod %in% c('GLS', 'TSGLS')) {
                Psi[jj,kk] <- nCov(j, k, h, m, R)
                Psi[kk,jj] <- Psi[jj,kk]
            } else {
                Psi[jj,kk] <- adfCov(j, k, h, m, R, moments)
                Psi[kk,jj] <- Psi[jj,kk]
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
        nMatrix <- diag(rep(N, hypothesis.length))
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

            weight <- nrow(hypothesis[hypothesis[,4] == parameter,])

            UL <- fisherTransform(point.estimate) + critical_value*sqrt(1/(weight*N-3))
            UL <- tanh(UL)
            LL <- fisherTransform(point.estimate) - critical_value*sqrt(1/(weight*N-3))
            LL <- tanh(LL)

            UL <- round(UL, 3)
            LL <- round(LL, 3)
            gammaGLS_ci[i] <- paste0('[',LL,', ',UL,']')

        }
            covgamma <- round(covgamma, 3)
            gammaGLS <- round(gammaGLS, 3)
            estimates.table <- cbind(parameters, gammaGLS, covgamma, gammaGLS_ci)
            colnames(estimates.table) <- c("Parameter Tags", "Point Estimate", "Std. Error", paste0(100 - corrected_alpha * 100, "% Confidence Interval"))
    }


    R.GLS <- R
    for (jj in 1:hypothesis.length) {
        j <- hypothesis[jj,2]
        k <- hypothesis[jj,3]
        R.GLS[j,k] <- rhoGLS[jj]
        R.GLS[k,j] <- rhoGLS[jj]
    }

    SLS <- matrix(0, nrow=hypothesis.length, ncol=hypothesis.length)
    for (jj in 1:hypothesis.length) {
        for (kk in 1:jj) {
            j <- hypothesis[jj,2]
            k <- hypothesis[jj,3]
            h <- hypothesis[kk,2]
            m <- hypothesis[kk,3]
            SLS[jj,kk] <- sigmaLS[jj,kk]/(1 - R.GLS[j,k]^2)*(1 - R.GLS[h,m]^2)
            SLS[kk,jj] <- SLS[jj,kk] ## Adding this makes the output differ from Multicorr; probably a bug in the original
        }
    }


    ## Calculate test statistic and p value
    chisquare <- (N - 3)*t(e)%*%solve(SLS)%*%e
    k <- ncol(R)
    q <- parameters.length
    p <- pchisq(chisquare, df=(k-q), lower.tail = FALSE)


    ## Round and assemble output
    source("./multicorr/pRound.R")
    chisquare <- round(chisquare, 3)
    p <- pRound(p)
    sigtable <- matrix(c(chisquare, k-q, p), nrow=1)
    colnames(sigtable) <- c("Chi Square", "df", "pvalue")


    ## Test whether the hypothesis is the identity hypothesis
    ## If so, run superior S test
    identity.matrix <- diag(k)
    if (all(R.GLS == identity.matrix)) {

        ## Run S test
        scalc <- dget("./multicorr/scalc.R")
        s2star <- dget("./multicorr/s2star.R")
        S <- scalc(R, N)
        S <- s2star(0, N, k, S)
        Sp <- pchisq(S, df=k, lower.tail=FALSE)

        ## Round and assemble output
        source("./multicorr/pRound.R")
        S <- round(S, 3)
        k <- round(k, 3)
        Sp <- pRound(Sp)
        S.result <- matrix(c(S, k, Sp), nrow=1, ncol=3)
        colnames(S.result) <- c("S", "df", "pvalue")

    } else {
        S.result <- NA
    }

    output <- list(hypothesis, N, R, R.OLS, estimates.table, sigtable, S.result, MardiaSK)

    return(output)
}
