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
    if (is.null(input$expFile) || is.null(input$clustFile))
      return()
    myGenes <- getGeneList()
    myClusters <- getClusterList()
    return(wellPanel(selectInput(inputId="gene", 
                          label="Gene", 
                          choices=myGenes),
              selectInput(inputId="condition",
                          label="Condition",
                          choices=myClusters))
    )
  })
  # check if file is uploaded
  output$fileUploaded <- reactive({
    return(!is.null(mergeData()))
  })
  outputOptions(output, "fileUploaded", suspendWhenHidden=FALSE)
  # load expression data if available
  getExpData <- reactive({
    expFile <- input$expFile
    if (is.null(expFile))
      return(NULL)
    return(read.table(expFile$datapath, header=T, check.names=F, row.names=1, sep="\t"))
  })
  # load cluster info data if available
  getClustData <- reactive({
    clustFile <- input$clustFile
    if (is.null(clustFile))
      return(NULL)
    return(read.table(clustFile$datapath, header=T, check.names=F, row.names=1, sep="\t"))
  })
  # merge expression data with clust info data
  mergeData <- reactive({
    # get expression data
    expData <- getExpData()
    if (is.null(expData))
      return(NULL)
    # get clust info data
    clustData <- getClustData()
    if (is.null(clustData))
      return(NULL)
    # merge two tables
    combinedData <- merge(t(expData), clustData, by=0, all=T)
    rownames(combinedData) <- combinedData$Row.names
    return(combinedData <- combinedData[,-1])
  })
  # load available genes
  getGeneList <- reactive({
    # get expression data
    expData <- getExpData()
    if (is.null(expData))
      return(NULL)
    return(rownames(expData))
  })
  # load clusters
  getClusterList <- reactive({
    # get clust info data
    clustData <- getClustData()
    if (is.null(clustData))
      return(NULL)
    return(unique(clustData$ident))
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
