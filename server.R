#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)

options(shiny.maxRequestSize=100*1024^2)

# Define server logic required to draw plots
shinyServer(function(input, output) {
  # update dropdown menus based on user input expression data
  output$myMenu <- renderUI({
    if (is.null(input$inputFile))
      return()
    myGenes <- getGeneList()
    return(wellPanel(selectInput(inputId="gene", 
                          label="Gene", 
                          choices=myGenes),
              selectInput(inputId="condition",
                          label="Condition",
                          choices=c("Ctrl","Diabetes")))
    )
  })
  # check if file is uploaded
  output$fileUploaded <- reactive({
    return(!is.null(getExpData()))
  })
  outputOptions(output, "fileUploaded", suspendWhenHidden=FALSE)
  # load expression data if available
  getExpData <- reactive({
    inputFile <- input$inputFile
    if (is.null(inputFile))
      return(NULL)
    return(read.table(inputFile$datapath, header=T, check.names=F, row.names=1, sep="\t"))
  })
  # load available genes
  getGeneList <- reactive({
    # get expression data
    expData <- getExpData()
    if (is.null(expData))
      return(NULL)
    return(rownames(expData))
  })
  # output$distPlot <- renderPlot({
  #   
  #   # generate bins based on input$bins from ui.R
  #   x    <- faithful[, 2] 
  #   bins <- seq(min(x), max(x), length.out = input$bins + 1)
  #   
  #   # draw the histogram with the specified number of bins
  #   hist(x, breaks = bins, col = 'darkgray', border = 'white')
  #   
  # })
  
})
