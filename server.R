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
                          choices=myGenes,
                          selected=intersect(c("Gapdh","GAPDH"),myGenes)),
              selectInput(inputId="cluster",
                          label="Cluster",
                          choices=myClusters),
                          selected="All"
              )
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
    ##combinedData <- combinedData[,-1]
    return(combinedData)
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
    return(c("All",unique(as.vector(clustData$ident))))
  })
  # get selected gene
  getGene <- reactive({
    # get gene list
    myGenes <- getGeneList()
    if (is.null(myGenes))
      return(NULL)
    if (is.null(input$gene))
      return(intersect(c("Gapdh","GAPDH"),myGenes))
    return(input$gene)
  })
  # extract data columns for plotting
  prepareData <- reactive({
    # get combined data table
    combinedData <- mergeData()
    if (is.null(combinedData))
      return(NULL)
    # get current selected gene
    curGene <- getGene()
    # extract needed columns
    dataToPlot <- combinedData[,c("Row.names","ident","tSNE_1","tSNE_2",curGene)]
    colnames(dataToPlot) <- c("cellID","ident","tSNE_1","tSNE_2","gene")
    return(dataToPlot)
  })
  # draw reference tSNE plot
  output$tsneRef <- renderPlot({
    # get data for plotting
    dataToPlot <- prepareData()
    if (is.null(dataToPlot))
      return(NULL)
    # draw reference tSNE plot
    g <- ggplot(dataToPlot, aes(x=tSNE_1,y=tSNE_2,color=ident))
    g <- g + geom_point(shape=19, size=3, alpha=.8)
    g <- g + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
    g <- g + theme(legend.justification=c(0,0), legend.title=element_blank())
    return(g)
  })
  # get expression upperbound
  getMaxExp <- reactive({
    # get data for plotting
    dataToPlot <- prepareData()
    if (is.null(dataToPlot))
      return(NULL)
    # find upperbound and round up
    maxexp <- max(dataToPlot$gene)
    roundexp <- round(maxexp)
    up <- ceiling(maxexp)
    bottom <- floor(maxexp)
    hcut <- bottom + 0.5
    if(roundexp == up){
      # > .5
      hcut <- up
    }
    return(hcut)
  })
  # draw expression tSNE plot
  output$tsneExp <- renderPlot({
    # get data for plotting
    dataToPlot <- prepareData()
    if (is.null(dataToPlot))
      return(NULL)
    # get expression upperbound
    hcut <- getMaxExp()
    # draw expression tSNE plot
    g <- ggplot(dataToPlot, aes(x=tSNE_1,y=tSNE_2,color=gene))
    g <- g + geom_point(shape=19, size=3, alpha=.8)
    #g <- g + coord_cartesian(xlim=c(-40, 35), ylim=c(-45,35))
    g <- g + scale_color_gradient(low="grey", high="red", limits=c(0,hcut))
    g <- g + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
    g <- g + theme(legend.justification=c(0,0), legend.title=element_blank())
    g <- g + theme(legend.key=element_blank()) + theme(legend.text=element_text(size=18))
    g <- g + theme(axis.text=element_text(size=24), axis.title=element_text(size=26,face="bold"))
    return(g)
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
