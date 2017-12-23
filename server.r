shinyServer(function(input, output, session) {
   options(shiny.maxRequestSize = 30 * 1024^2 )
 
    output$estimationmethodInput <- renderUI({
        html_ui <- " "
        if (input$datatype == "rawdata") {
            html_ui <- paste0(radioButtons("estimationmethod", "Estimation method:", c("GLS" = "GLS", "TSGLS" = "TSGLS", "ADF" = "ADF", "TSADF" = "TSADF")))
        } else {
            html_ui <- paste0(radioButtons("estimationmethod", "Estimation method:", c("GLS" ="GLS", "TSGLS" = "TSGLS")))
        }
        HTML(html_ui)
    })

    output$wbsctOutput <- eventReactive(input$runButton, {

		# Import data files
        validate(need(input$datafile, ""))
        data <- input$datafile

		# Check that R can read the data file as a .csv
		result = tryCatch({
		  read.csv(file=input$datafile[[4]], head=FALSE, sep=",")
		}, warning = function(w) {
		  'problem'
		}, error = function(e) {
		  'problem'
		}, finally = {
		})
		if ('problem' %in% result) {
		  return(capture.output(cat('<br>Error: There was an problem reading data file #', i, '; it may not be a .csv file.', sep="")))
		} else { # If so import it as a matrix
			data <- as.matrix(read.csv(file=input$datafile[[4]], head=FALSE, sep=","))
		}

		if (ncol(data) > 16) {
		  return(capture.output(cat('<br>Error: WBCORR does not support more than 16 variables.', sep="")))
		}

		# Import N (calculate if raw data) for each group
        NList <- c()
        if (input$datatype == 'correlation') {
			    N <- as.numeric(input$N)
        } else {
			    N <- nrow(data)
        }

		# Check that hypothesis file is readable as a .csv
        validate(need(input$hypothesisfile, ""))
        result = tryCatch({
          read.csv(file=input$hypothesisfile[[4]], head=FALSE, sep=",")
        }, warning = function(w) {
          'problem'
        }, error = function(e) {
          'problem'
        }, finally = {
        })
        if ('problem' %in% result) {
          return(capture.output(cat('<br>Error: There was an problem reading the hypothesis file; it may not be a .csv file.', sep="")))
        } else {
          hypothesis <- as.matrix(read.csv(file=input$hypothesisfile[[4]], head=FALSE, sep=","))
        }

        if (input$datatype == 'rawdata') {
          if (input$estimationmethod %in% c('ADF','TSADF')) {
              deletion <- input$deletion2
          } else if (input$estimationmethod %in% c('GLS','TSGLS')) {
              deletion <- input$deletion1
          }
        } else {
          deletion <- 'no'
        }

        ComputeWBCorrChiSquare <- dget("ztest.R")
        capture.output(ComputeWBCorrChiSquare(data, N, hypothesis, input$datatype, input$estimationmethod, deletion))
    
    })

    observeEvent(input$runButton, {
      updateTabsetPanel(session, "inTabset", 'out')
    })
    

})
