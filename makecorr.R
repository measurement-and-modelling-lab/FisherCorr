function (data, deletion) {
	RList <- list()
	moments <-list()

	temp <- suppressWarnings(as.numeric(data))
	data <- matrix(temp, nrow=nrow(data), ncol=ncol(data))

	if (deletion == 'pairwise' && NA %in% temp) {
	  RList <- cov(data, use='pairwise')
	} else {
	  RList <- cov(data)
	}

	RList <- cov2cor(RList)
	moments <- compute4thOrderMoments(data)

	output <- list(RList, moments)
	return(output)
}
