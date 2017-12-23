function (data, datatype, hypothesis, deletion) {

	if (ncol(hypothesis) != 5) {
		cat('<br>Error: The hypothesis matrix has the wrong number of columns.')
		return(invisible(TRUE))
	}

    fixed <- hypothesis[,5]
    first4columns <- hypothesis[,-5]

    if (is.numeric(hypothesis) == FALSE) {
        cat('<br>Error: The hypothesis matrix has a non-numeric entry.')
        return(invisible(TRUE))
    } else if ((length(fixed[fixed > 1]) + length(fixed[fixed < -1])) > 0) {
        cat('<br>Error: The hypothesis matrix has a fixed value that is less than -1 or greater than 1.', sep="")
        return(invisible(TRUE))
    } else if (FALSE %in% ((first4columns - floor(first4columns)) == 0)) {
        cat('<br>Error: The hypothesis matrix has a non-integer where it shouldn\'t.', sep="")
        return(invisible(TRUE))
    } else if (TRUE %in% (first4columns < 0)) {
        cat('<br>Error: The hypothesis matrix has a negative number where it shouldn\'t.', sep="")
        return(invisible(TRUE))
    }

    for (i in 1:nrow(hypothesis)) {
        if (!(hypothesis[i,2] %in% 1:ncol(data))) {
            cat('<br>Error: Row ', i, ' of the hypothesis matrix references a non-existent variable.', sep="")
            return(invisible(TRUE))
        } else  if (!(hypothesis[i,3] %in% 1:ncol(data))) {
            cat('<br>Error: Row ', i, ' of the hypothesis matrix references a non-existent variable.', sep="")
            return(invisible(TRUE))
        }
    }

	if (datatype == 'rawdata') {
		if (nrow(data) <= ncol(data)) {
			cat('<br>Error: Either data matrix #', jj, ' is a correlation matrix, or it is a raw data matrix where n <= p.', sep="")
			return(invisible(TRUE))
		} else if (deletion == 'pairwise') {
			m <- data
			data <- m[rowSums(is.na(m))!=ncol(m), ]
			colCheck <- m[,colSums(is.na(m))!=nrow(m)]
			if(ncol(colCheck) != ncol(m)){
			  cat('<br>Error: data matrix #', jj, ' has at least one empty column.', sep="")
			  return(invisible(TRUE))
			}
		}
	}

	if (datatype == 'correlation') {
		if ((length(data[data > 1]) + length(data[data < -1])) > 0) {
		  cat('<br>Error: The correlation matrix has a value that is less than -1 or greater than 1.')
		  return(invisible(TRUE))
		} else if (nrow(data) != ncol(data)) {
			cat('<br>Error: The correlation matrix is not square.')
			return(invisible(TRUE))
		}
	}

	if (deletion == 'no') {
		if (TRUE %in% is.na(data)) {
			cat('<br>Error: Data matrix #', jj, ' has at least one empty entry.', sep="")
			return(invisible(TRUE))
		} else if (is.numeric(data) == FALSE) {
			cat('<br>Error: Data matrix #', jj, ' has at least one non-numeric entry.', sep="")
			return(invisible(TRUE))
		}
	}

    return(invisible(FALSE))

}
