require(htmlTable) || install.packages(htmlTable)
require(shiny) || install.packages(shiny)
require(shinythemes) || install.packages(shinythemes)

shinyServer(function(input, output, session) {
    options(shiny.maxRequestSize = 30 * 1024^2 )

    observe({
        options(shiny.sanitize.errors = FALSE)
    })
    
    output$estimationmethodInput <- renderUI({
        html_ui <- " "
        if (input$datatype == "rawdata") {
            html_ui <- paste0(radioButtons("estimationmethod", "Estimation method:", c("GLS" = "GLS", "TSGLS" = "TSGLS", "ADF" = "ADF", "TSADF" = "TSADF")))
        } else {
            html_ui <- paste0(radioButtons("estimationmethod", "Estimation method:", c("GLS" ="GLS", "TSGLS" = "TSGLS")))
        }
        HTML(html_ui)
    })

    ## Produce one data/file input for each group
    output$datafileInput <- renderUI({ 
        if (input$samples %in% 1:8) {
            html_ui <- " "
            if (input$datatype == 'correlation') {
                for (i in 1:input$samples){
                    html_ui <- paste0(html_ui,
                                      div(style="display: inline-block;vertical-align:top; height: 70px;",
                                          numericInput(paste0("samplesize",i), paste0("N for Data File #",i,":"), "", min = 1, max = NA, step = 1)),
                                      div(style="height: 65px;", fileInput(paste0("datafile",i), label=paste0("Data file #",i,":"))), '<hr>')
                }
            } else {
                for (i in 1:input$samples) {
                    html_ui <- paste0(html_ui, div(style="height: 65px;", fileInput(paste0("datafile",i), label=paste0("Data file #",i,":"))))
                    if (i < input$samples) {
                        html_ui <- paste0(html_ui, '<br>')
                    }
                }
            }
            HTML(html_ui)
        }
    })
    
    wbsctOutput <- eventReactive(input$runButton, {

        ## Read functions
        ComputeMulticorrChiSquare <- dget("./multicorr/ComputeMulticorrChiSquare.R")


        ## Stipulate data type
        datatype <- input$datatype


        ## Stipulate estimation method
        estimation.method <- input$estimationmethod


        ## Import data files
        validate(need("input$datafile", "")) ## Check that data file i exists

        ## Check that R can read the data file as a .csv
        tryCatch({
            read.csv(file=input$datafile[[4]], head=FALSE)
        }, warning = function(w) {
            stop("There was a problem reading one of your .csv files. You may need to add a blank line to the end of the file.")
        }, error = function(e) {
            stop("There was a problem reading one of your .csv files. You may need to add a blank line to the end of the file.")
        })

        data <- as.matrix(read.csv(file=input$datafile[[4]], head=FALSE))

        if (ncol(data) > 16) {
            stop("The web version of MML-Multicorr does not support more than 16 variables.")
        }


        ## Check that hypothesis file is readable as a .csv
        validate(need(input$hypothesisfile, ""))
        tryCatch({
            read.csv(file=input$hypothesisfile[[4]], head=FALSE)
        }, warning = function(w) {
            stop("There was a problem reading your hypothesis file. You may need to add a blank line to the end of the file.")
        }, error = function(e) {
            stop("There was a problem reading your hypothesis file. You may need to add a blank line to the end of the file.")
        })


        ## Read the hypothesis file
        hypothesis <- as.matrix(read.csv(file=input$hypothesisfile[[4]], head=FALSE))


        ## Import N (calculate if raw data) for each group
        validate(need(input$samplesize, ""))
        if (datatype == "correlation") {
            N <- input$samplesize
        } else {
            N <- nrow(data)
        }


        ## Define deletion method
        if (datatype == "rawdata") {
            if (estimation.method %in% c("ADF","TSADF")) {
                deletion <- input$adfdeletion
            } else if (estimation.method %in% c("GLS","TSGLS")) {
                deletion <- input$glsdeletion
            }
        } else {
            deletion <- "nodeletion"
        }


        ## Run the test
        output <- ComputeMulticorrChiSquare(data, N, hypothesis, datatype, estimation.method, deletion)


        html.output <- ""
        #NList <- output[[6]] ## Import amended sample sizes


        ## Print the original hypothesis matrix
        hypothesis.colnames <- c("Group", "Row", "Column", "Parameter Tag", "Fixed Value")
        html.output <- paste0(html.output, htmlTable(hypothesis, align="c", caption="Input Hypothesis Matrix", header=hypothesis.colnames))


        ## Print the amended hypothesis matrix, if changes were made
        hypothesis.amended <- output[[1]]
        if (!all(hypothesis == hypothesis.amended)) {
            html.output <- paste0(html.output, htmlTable(hypothesis.amended, align="c", caption="Amended Hypothesis Matrix", header=hypothesis.colnames))
        }

        N <- output[[2]]

        ## Print the correlation matrices
        R <- output[[3]]
        variables <- nrow(R)
        rownames(R) <- colnames(R) <- lapply(1:variables, function(i) "")
        R <- round(R, 3)
        labels <- paste0("<b>X<sub>", 1:variables, "</sub></b>")
        caption <- paste0("Input Correlation Matrix (N=", N, ")")
        html.output <- paste0(html.output, htmlTable(R, align="r", caption=caption, rnames=labels, header=labels, align.header="r", css.cell = "padding-left: .5em; padding-right: .2em;"))


        ## Print the OLS matrices
        R.OLS <- output[[4]]
        rownames(R.OLS) <- colnames(R.OLS) <- lapply(1:variables, function(i) "")
        R.OLS <- round(R.OLS, 3)
        labels <- paste0("<b>X<sub>", 1:variables, "</sub></b>")
        caption <- paste0("OLS Matrix (N=", N, ")")
        html.output <- paste0(html.output, htmlTable(R.OLS, align="r", caption=caption, rnames=labels, header=labels, align.header="r", css.cell = "padding-left: .5em; padding-right: .2em;"))


        ## Print the parameter estimates
        gammahatDisplay <- output[[5]]
        if (is.matrix(gammahatDisplay)) {
            gammahatDisplay <- gammahatDisplay[order(gammahatDisplay[,1]), , drop=FALSE] ## Order the estimates by parameter tag
            if (!identical(NA, gammahatDisplay)) {
                header <- paste0(estimation.method, " Parameter Estimates")
                html.output <- paste0(html.output, htmlTable(gammahatDisplay, align="c", caption=header))
            }
        }


        ## Return significance test results
        sigtable <- output[[6]]
        header <- "Significance Test Results"
        html.output <- paste0(html.output, htmlTable(sigtable, align="c", caption=header))


        ## Return significance test results
        S.result <- output[[8]]
        if (is.matrix(S.result)) {
            header <- "Significance Test Results"
            html.output <- paste0(html.output, htmlTable(S.result, align="c", caption=header))
        }


        ## Print MVN test
        if (datatype == "rawdata") {
            MardiaSK <- output[[9]]
            if (deletion == "pairwise") {
                html.output <- paste0(html.output, htmlTable(MardiaSK[[1]], align="c", caption="Assessment of the Distribution of the Observed Marginals"))
                html.output <- paste0(html.output, htmlTable(MardiaSK[[2]], align="c", caption="Assessment of Multivariate Normality"))
            } else {
                html.output <- paste0(html.output, htmlTable(MardiaSK[[1]], align="c", caption="Assessment of Multivariate Skewness"))
                html.output <- paste0(html.output, htmlTable(MardiaSK[[2]], align="c", caption="Assessment of Multivariate Kurtosis"))
            }
        }

        HTML(html.output)

    })

    observeEvent(input$runButton, {
        updateTabsetPanel(session, "inTabset", 'out')
    })

    output$finaloutput <- renderUI({ 
        wbsctOutput()
    })

    

})
