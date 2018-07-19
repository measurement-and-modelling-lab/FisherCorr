## ## Load functions
## compute4thOrderMoments <- dget("./wbcorr/compute4thOrderMoments.R")
ztest <- dget("./multicorr/ztest.R")
## errorcheck <- dget("./wbcorr/errorcheck.R")
tablegen <- dget("./multicorr/tablegen.R")


## ## Stipulate data type
## cat("\nWhat type of data will you be using? 1. Raw data 2. Correlation data\n")
## datatype <- readline(prompt="")
## if (datatype == 1) {
##     methods <- c("GLS", "TSGLS", "ADF", "TSADF")
##     cat_methods <- "1. GLS  2. TSGLS  3. ADF  4. TSADF"
##     datatype <- "rawdata"
## } else if (datatype == 2) {
##     methods <- c("GLS", "TSGLS")
##     cat_methods <- "1. GLS  2. TSGLS"
##     datatype <- "correlation"
## } else {
##     stop("Invalid data type.")
## }


## ## Specify input files
## cat("\nInput your files in the following format: data.csv;hypothesis.csv\n")
## files <- readline(prompt="")
## files <- strsplit(files, ";")[[1]]
## data.length <- length(files) - 1
## if (length(files) != 1) {
##     stop("Your files should consist of one data file and one hypothesis file.")
## }


## ## Read data file(s)
## filename <- files[[1]]
## file.exists <- file.exists(filename)
## if (file.exists) {
##     data <- read.csv(file=filename,head=FALSE, sep=",")
##     data <- as.matrix(data)
## } else {
##     stop("Data file does not exist.")
## }


## ## Read hypothesis file
## hypothesis.file <- files[[2]]
## file.exists <- file.exists(hypothesis.file)
## if (file.exists) {
##     hypothesis <- read.csv(file=hypothesis.file,head=FALSE)
##     hypothesis <- as.matrix(hypothesis)
## } else {
##     stop("Hypothesis file does not exist")
## }


## ## Generate sample size list
## if (datatype == "rawdata") {
##     N <- nrow(data)
## } else {
##     cat("\nInput the sample size for your data:\n")
##     N <- readline(prompt="")
## }


## ## Select estimation methods to be included
## cat("\nWhat estimation method would you like to use?", cat_methods, "\n")
## method.index <- readline(prompt="")
## if (method.index %in% 1:length(methods)) {
##     method.index <- as.numeric(method.index)
##     estimation.method <- methods[method.index]
## } else {
##     stop("Invalid estimation method.")
## }


## ## Choose how to deal with missing data
## if (datatype == "rawdata") {
##     cat("\nHow would you like to deal with missing data? 1. Listwise  2. Pairwise  3. There is no missing data\n")
##     deletion.index <- readline(prompt="")
##     if (!(deletion.index %in% c("1", "2", "3"))) {
##         stop("Invalid deletion method.")
##     } else {
##         deletion.index <- as.numeric(deletion.index)
##         deletion <- c("listwise", "pairwise", "nodeletion")[deletion.index]
##     }
## } else {
##     deletion <- "nodeletion"
## }




data <- read.csv(file="example1_data1.csv",head=FALSE, sep=",")
data <- as.matrix(data)

hypothesis <- read.csv(file="example1_hypothesis.csv",head=FALSE, sep=",")
hypothesis <- as.matrix(hypothesis)

N <- nrow(data)
hypothesis <- matrix(c(1,1,1,2,3,3,1,1,2,1,1,1,0,0,0), nrow=3)
datatype <- "rawdata"
estimation.method <- "TSADF"
deletion <- "nodeletion"

## Run the test
output <- ztest(data, N, hypothesis, datatype, estimation.method, deletion)


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
if (!identical(NA, estimates.table)) {
    cat("\n", estimation.method, "Parameter Estimates\n\n")
    tablegen(estimates.table, TRUE)
}


## Print the significance of the test
sigtable <- output[[6]]
cat("\nSignificance Test Results\n\n")
tablegen(sigtable, TRUE)


## Print MVN test
if (datatype == "rawdata") {
    MardiaSK <- output[[9]]
    if (deletion == "pairwise") {
        cat("\nAssessment of the Distribution of the Observed Marginals\n\n\n")
        tablegen(MardiaSK[[1]], TRUE)
        cat("\nAssessment of Multivariate Normality\n\n")
        tablegen(MardiaSK[[2]], TRUE)
    } else {
        cat("\nAssessment of Multivariate Normality\n\n")
        tablegen(MardiaSK[[1]], TRUE)
        cat("\n")
        tablegen(MardiaSK[[2]], TRUE)
    }
}
