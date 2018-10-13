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
    myConditions <- getNonNumericCols()
    return(wellPanel(selectInput(inputId="gene", 
                                 label="Gene",
                                 choices=myGenes,
                                 selected=intersect(c("Gapdh","GAPDH"),myGenes)
                    ),
                    selectInput(inputId="cluster",
                                label="Cluster",
                                choices=c("All",myClusters),
                                selected="All"
                    ),
                    selectInput(inputId="condition",
                                label="Condition",
                                choices=c("None",myConditions),
                                selected="None"
                    )
    ))
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
    return(unique(as.vector(clustData$ident)))
  })
  # fetch current selected gene
  getGene <- reactive({
    # get gene list
    myGenes <- getGeneList()
    if (is.null(myGenes))
      return(NULL)
    if (is.null(input$gene))
      return(intersect(c("Gapdh","GAPDH"),myGenes))
    return(input$gene)
  })
  # fetch current selected cluster
  getCluster <- reactive({
    # get cluster list
    myClusters <- getClusterList()
    if (is.null(myClusters))
      return(NULL)
    if (is.null(input$cluster) || input$cluster == "All")
      return(myClusters)
    else
      return(input$cluster)
  })
  # fetch current selected condition
  getCondition <- reactive({
    # get condition name list
    myConditions <- getNonNumericCols()
    if (is.null(myConditions))
      return(NULL)
    if (is.null(input$condition))
      return("None")
    return(input$condition)
  })
  # extract non-numeric columns
  getNonNumericCols <- reactive({
    # get combined data table
    combinedData <- mergeData()
    if (is.null(combinedData))
      return(NULL)
    # get non-numeric columns
    selCols <- colnames(combinedData)[!sapply(combinedData, is.numeric)]
    # remove columns "ident", "orig.ident"
    return(selCols[!selCols %in% c("ident","orig.ident","Row.names")])
  })
  # extract data columns for plotting
  prepareData <- reactive({
    # get combined data table
    combinedData <- mergeData()
    if (is.null(combinedData))
      return(NULL)
    # get current selected gene
    curGene <- getGene()
    # get current selected condition column name
    curCondition <- getCondition()
    # extract needed columns and
    # add a new column by combining the "ident" column and the selected condition column
    if (curCondition != "None"){
      dataToPlot <- combinedData[,c("Row.names","ident","tSNE_1","tSNE_2",curGene,curCondition)]
      colnames(dataToPlot) <- c("cellID","ident","tSNE_1","tSNE_2","gene","condition")
      dataToPlot$cohort <- paste(dataToPlot$ident, dataToPlot$condition, sep=" - ")
    }
    else{
      dataToPlot <- combinedData[,c("Row.names","ident","tSNE_1","tSNE_2",curGene)]
      colnames(dataToPlot) <- c("cellID","ident","tSNE_1","tSNE_2","gene")
      dataToPlot$cohort <- dataToPlot$ident
    }
    return(dataToPlot)
  })
  # order dataToPlot by cohort (for barplot)
  orderDataByCohort <- reactive({
    # get data for plotting
    dataToPlot <- prepareData()
    if (is.null(dataToPlot))
      return(NULL)
    # order by cohort
    dataToPlotOrdered <- dataToPlot[with(dataToPlot, order(cohort)),]
    dataToPlotOrdered$cohort <- factor(dataToPlotOrdered$cohort)
    ##dataToPlotOrdered$cohort <- factor(dataToPlotOrdered$cohort, levels=c("Endothelial cells - None","Endothelial cells - Diabetes","Mesangial cells - None","Mesangial cells - Diabetes","Podocyte - None","Podocyte - Diabetes","Immune cells - None","Immune cells - Diabetes","Tubular cells - None","Tubular cells - Diabetes"))
    return(dataToPlotOrdered)
  })
  # order dataToPlot by gene expression (for tSNE)
  orderDataByExpression <- reactive({
    # get data for plotting
    dataToPlot <- prepareData()
    if (is.null(dataToPlot))
      return(NULL)
    # order by cohort
    dataToPlotOrdered <- dataToPlot[with(dataToPlot, order(gene)),]
    return(dataToPlotOrdered)
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
  # get lower and upper bound of tSNE coordinates
  getAxisRange <- reactive({
    # get data for plotting
    dataToPlot <- prepareData()
    if (is.null(dataToPlot))
      return(NULL)
    # find minimum and maximum for both axes
    min.x <- min(dataToPlot$tSNE_1)
    max.x <- max(dataToPlot$tSNE_1)
    min.y <- min(dataToPlot$tSNE_2)
    max.y <- max(dataToPlot$tSNE_2)
    # round up with base 5
    myBase <- 5
    lower.x <- round(min.x / myBase) * myBase
    if (lower.x > min.x)
      lower.x <- min.x - myBase
    upper.x <- round(max.x / myBase) * myBase
    if (upper.x < max.x)
      upper.x <- max.x + myBase
    lower.y <- round(min.y / myBase) * myBase
    if (lower.y > min.y)
      lower.y <- min.y - myBase
    upper.y <- round(max.y / myBase) * myBase
    if (upper.y < max.y)
      upper.y <- max.y + myBase
    return(c(lower.x,upper.x,lower.y,upper.y))
  })
  # draw reference tSNE plot
  output$tsneRef <- renderPlot({
    # get ordered data for plotting
    dataToPlotOrdered <- orderDataByExpression()
    if (is.null(dataToPlotOrdered))
      return(NULL)
    # get axis range
    axis.range <- getAxisRange()
    lower.x <- axis.range[1]
    upper.x <- axis.range[2]
    lower.y <- axis.range[3]
    upper.y <- axis.range[4]
    # get current selected cluster
    curCluster <- getCluster()
    # draw reference tSNE plot
    g <- ggplot(dataToPlotOrdered[dataToPlotOrdered$ident %in% curCluster, ], aes(x=tSNE_1,y=tSNE_2,color=ident))
    g <- g + geom_point(shape=19, size=3, alpha=.8)
    g <- g + coord_cartesian(xlim=c(lower.x, upper.x), ylim=c(lower.y, upper.y))
    g <- g + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
    g <- g + theme(legend.justification=c(0,0), legend.title=element_blank())
    return(g)
  })
  # draw expression tSNE plot
  output$tsneExp <- renderPlot({
    # get ordered data for plotting
    dataToPlotOrdered <- orderDataByExpression()
    if (is.null(dataToPlotOrdered))
      return(NULL)
    # get expression upperbound
    hcut <- getMaxExp()
    # get axis range
    axis.range <- getAxisRange()
    lower.x <- axis.range[1]
    upper.x <- axis.range[2]
    lower.y <- axis.range[3]
    upper.y <- axis.range[4]
    # get current selected cluster
    curCluster <- getCluster()
    # draw expression tSNE plot
    g <- ggplot(dataToPlotOrdered[dataToPlotOrdered$ident %in% curCluster, ], aes(x=tSNE_1,y=tSNE_2,color=gene))
    g <- g + geom_point(shape=19, size=3, alpha=.8)
    g <- g + coord_cartesian(xlim=c(lower.x, upper.x), ylim=c(lower.y, upper.y))
    g <- g + scale_color_gradient(low="grey", high="red", limits=c(0,hcut))
    g <- g + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
    g <- g + theme(legend.justification=c(0,0), legend.title=element_blank())
    g <- g + theme(legend.key=element_blank()) + theme(legend.text=element_text(size=18))
    g <- g + theme(axis.text=element_text(size=24), axis.title=element_text(size=26,face="bold"))
    return(g)
  })
  # draw expression bar plot
  output$barExp <- renderPlot({
    # get ordered data for plotting
    dataToPlotOrdered <- orderDataByCohort()
    if (is.null(dataToPlotOrdered))
      return(NULL)
    # generate plot
    g <- ggplot(dataToPlotOrdered, aes(x=factor(cellID), y=gene, fill=factor(ident)))
    g <- g + geom_bar(stat="identity")
    g <- g + facet_grid(. ~ cohort, scales="free_x", space="free_x", switch="x")
    g <- g + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.justification=c(0,0), legend.title=element_blank())
    g <- g + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
    g <- g + theme(strip.background=element_rect(color="white", fill=NA), panel.border=element_rect(color="lightgray", fill=NA))
    g <- g + theme(axis.line = element_line(color = "black")) + theme(strip.text.x = element_text(angle = 90))
    g <- g + theme(legend.key=element_blank()) + theme(legend.text=element_text(size=12)) + theme(axis.text=element_text(size=18), axis.title=element_text(size=18,face="bold"))
    ##g <- g + ggtitle(paste(gene,"expression per cell",sep=" - ")) + theme(plot.title = element_text(size=14, face="bold", hjust = 0.5))
    g <- g + ylab("Expression (log-scale)") + theme(axis.title=element_text(size=12,face="bold"),axis.title.x=element_blank(),axis.line.x=element_blank())
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
