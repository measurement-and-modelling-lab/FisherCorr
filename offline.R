## ## Load functions
ComputeMulticorrChiSquare <- dget("./multicorr/ComputeMulticorrChiSquare.R")
tablegen <- dget("./multicorr/tablegen.R")
RoundPercentile <- dget("./multicorr/RoundPercentile.R")


## Stipulate data type
cat("\nWhat type of data will you be using? 1. Raw data 2. Correlation data\n")
datatype <- readline(prompt="")
if (datatype == 1) {
    methods <- c("GLS", "TSGLS", "ADF", "TSADF")
    cat_methods <- "1. GLS  2. TSGLS  3. ADF  4. TSADF"
    datatype <- "rawdata"
} else if (datatype == 2) {
    methods <- c("GLS", "TSGLS")
    cat_methods <- "1. GLS  2. TSGLS"
    datatype <- "correlation"
} else {
    stop("Invalid data type.")
}


## Specify input files
cat("\nInput your files in the following format: data.csv;hypothesis.csv\n")
files <- readline(prompt="")
files <- strsplit(files, ";")[[1]]
data.length <- length(files) - 1
if (data.length != 1) {
    stop("Your files should consist of one data file and one hypothesis file.")
}


## Read data file(s)
filename <- files[[1]]
file.exists <- file.exists(filename)
if (file.exists) {
    data <- read.csv(file=filename,head=FALSE, sep=",")
    data <- as.matrix(data)
} else {
    stop("Data file does not exist.")
}


## Read hypothesis file
hypothesis.file <- files[[2]]
file.exists <- file.exists(hypothesis.file)
if (file.exists) {
    hypothesis <- read.csv(file=hypothesis.file,head=FALSE)
    hypothesis <- as.matrix(hypothesis)
} else {
    stop("Hypothesis file does not exist")
}


## If the hypothesis matrix doesn't have a group column, add one.
if (ncol(hypothesis) == 4) {
    hypothesis <- cbind(1, hypothesis)
}


## Generate sample size list
if (datatype == "rawdata") {
    N <- nrow(data)
} else {
    cat("\nInput the sample size for your data:\n")
    N <- readline(prompt="")
}


## Select estimation methods to be included
cat("\nWhat estimation method would you like to use?", cat_methods, "\n")
method.index <- readline(prompt="")
if (method.index %in% 1:length(methods)) {
    method.index <- as.numeric(method.index)
    estimation.method <- methods[method.index]
} else {
    stop("Invalid estimation method.")
}


## Choose how to deal with missing data
if (datatype == "rawdata") {
    cat("\nHow would you like to deal with missing data? 1. Listwise  2. Pairwise  3. There is no missing data\n")
    deletion.index <- readline(prompt="")
    if (!(deletion.index %in% c("1", "2", "3"))) {
        stop("Invalid deletion method.")
    } else {
        deletion.index <- as.numeric(deletion.index)
        deletion <- c("listwise", "pairwise", "nodeletion")[deletion.index]
    }
} else {
    deletion <- "nodeletion"
}


## Run the test
output <- ComputeMulticorrChiSquare(data, N, hypothesis, datatype, estimation.method, deletion)


## Print the original hypothesis matrix
cat("\nInput Hypothesis Matrix\n\n")
colnames(hypothesis) <- c("Group", "Row", "Column", "Parameter Tag", "Fixed Value")
tablegen(hypothesis, TRUE)


## Print the amended hypothesis matrix
hypothesis.amended <- output[[1]]
if (!all(hypothesis == hypothesis.amended)) {
    cat("\nAmended Hypothesis Matrix\n\n")
    colnames(hypothesis.amended) <- c("Group", "Row", "Column", "Parameter Tag", "Fixed Value")
    tablegen(hypothesis.amended, TRUE)
}


## Print the correlation matrices
N <- output[[2]]
R <- output[[3]]
cat("\nInput Correlation Matrix (N = ", N, ")\n\n", sep="")
R <- round(R, 3)
tablegen(R, FALSE)


## Print the OLS estimates
R.OLS <- output[[4]]
cat("\nOLS Matrix (N = ", N, ")\n\n", sep="")
R.OLS <- round(R.OLS, 3)
tablegen(R.OLS, FALSE)


## Print the parameter estimates
estimates.table <- output[[5]]
estimates.table <- round(estimates.table, 3)
if (!identical(NA, estimates.table)) {
    cat("\n", estimation.method, "Parameter Estimates\n\n")
    tablegen(estimates.table, TRUE)
}


## Print the significance of the test
sigtable <- output[[6]]
sigtable[,1] <- round(sigtable[,1], 3)
sigtable[,3] <- RoundPercentile(sigtable[,3])

cat("\nSignificance Test Results\n\n")
tablegen(sigtable, TRUE)


## Print the significance of S test
S.result <- output[[6]]
if (is.matrix(S.result)) {
    S.result[,1] <- round(S.result[,1], 3)
    S.result[,3] <- RoundPercentile(S.result[,3])
    cat("\nSignificance Test Results\n\n")
    tablegen(S.result, TRUE)
}


## Print MVN test
if (datatype == "rawdata") {
    MardiaSK <- output[[8]]
    if (deletion == "pairwise") {

        range.table <- MardiaSK[[1]]
        range.table[,4] <- round(range.table[,4], 3)
        range.table[,5] <- RoundPercentile(range.table[,5])

        normality.table <- MardiaSK[[2]]
        normality.table[,2] <- round(normality.table[,2], 3)
        normality.table[,3] <- RoundPercentile(normality.table[,3])
        
        cat("\nAssessment of the Distribution of the Observed Marginals\n\n\n")
        tablegen(range.table, TRUE)

        cat("\nAssessment of Multivariate Normality\n\n")
        tablegen(normality.table, TRUE)

    } else {

        skew.table <- MardiaSK[[1]]
        skew.table[,-5] <- round(skew.table[,-5], 3)
        skew.table[,5] <- RoundPercentile(skew.table[,5])

        kurt.table <- MardiaSK[[2]]
        kurt.table[,-4] <- round(kurt.table[,-4], 3)
        kurt.table[,4] <- RoundPercentile(kurt.table[,4])
        
        cat("\nAssessment of Multivariate Normality\n\n")
        tablegen(skew.table, TRUE)
        cat("\n")
        tablegen(kurt.table, TRUE)
    }
}
