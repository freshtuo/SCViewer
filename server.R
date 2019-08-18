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
library(ggrepel)
library(dplyr)
library(Matrix)

options(shiny.maxRequestSize=300*1024^2)

# load expression matrix
myLoadExpression <- function(texpFile, tgeneFile, tcellFile){
  #texp <- as(readMM(gzfile(paste(tdir, "expression.mtx.gz", sep="/"))), "dgCMatrix")
  texp <- as(readMM(gzfile(texpFile)), "dgCMatrix")
  rownames(texp) <- as.vector(read.table(file=tgeneFile, header=F, check.names=F)$V1)
  colnames(texp) <- as.vector(read.table(file=tcellFile, header=F, check.names=F)$V1)
  return(texp)
}

# Define server logic required to draw plots
shinyServer(function(input, output) {
  # update dropdown menus based on user input expression data
  output$myMenu <- renderUI({
    if (is.null(input$expFile) || is.null(input$clustFile))
      return()
    myGenes <- getGeneList()
    myClusters <- getClusterList()
    myConditions <- getNonNumericCols()
    if (input$tab == "tSNE")
      return(wellPanel(selectInput(inputId="gene", 
                                  label="Gene",
                                  choices=myGenes,
                                  selected=intersect(c("Gapdh","GAPDH"),myGenes)
                      ),
                      selectInput(inputId="cluster",
                                  label="Cluster",
                                  choices=c("All",myClusters),
                                  selected="All"
                      )
      ))
    else if (input$tab == "Bar" || input$tab == "Violin")
      return(wellPanel(selectInput(inputId="gene", 
                                   label="Gene",
                                   choices=myGenes,
                                   selected=intersect(c("Gapdh","GAPDH"),myGenes)
                      ),
                      selectInput(inputId="condition",
                                  label="Condition",
                                  choices=c("None",myConditions),
                                  selected="None"
                      )
      ))
    else if (input$tab == "Pie")
      return(wellPanel(selectInput(inputId="condition",
                                   label="Condition",
                                   choices=c("None",myConditions),
                                   selected="None")
                       )
             )
    else if (input$tab == "Density")
      return(wellPanel(selectInput(inputId="gene", 
                                   label="Gene",
                                   choices=myGenes,
                                   selected=intersect(c("Gapdh","GAPDH"),myGenes)
                      ),
                      selectInput(inputId="condition",
                                   label="Condition",
                                   choices=c("None",myConditions),
                                   selected="None")
                      )
             )
    else if (input$tab == "BarLine")
      return(wellPanel(selectInput(inputId="gene",
                                   label="Gene",
                                   choices=myGenes,
                                   selected=intersect(c("Gapdh","GAPDH"),myGenes)
                      ),
                      selectInput(inputId="condition",
                                  label="Condition",
                                  choices=c("None",myConditions),
                                  selected="None")
                      )
             )
    else
      return()
  })
  # check if file is uploaded
  output$fileUploaded <- reactive({
    return(!is.null(prepareData()))
  })
  outputOptions(output, "fileUploaded", suspendWhenHidden=FALSE)
  # track the current tab selection
  #observe({
  #  print(input$tab)
  #})
  # update descriptions based on current selected tab
  output$myDescription <- renderUI({
    if (input$tab == "tSNE")
      return(wellPanel(helpText("t-SNE plots showing different types of cells (right)",
                                "and highlighting expressions of a given gene (left).",br(),
                                "User can choose a gene and show its expression in cells of a given cell type.")))
    else if (input$tab == "Bar")
      return(wellPanel(helpText("Bar plot showing per-cell expressions.",
                                "Each bar represents a cell and its height indicates expression.",br(),
                                "By default, cells are grouped according to cell type clusters.",br(),
                                "User can choose an additional condition factor based on which",
                                "cells will be arranged in refined groups.")))
    else if (input$tab == "Violin")
      return(wellPanel(helpText("Violin plot showing expression distribution within a cell group.",br(),
                                "By default, cells are grouped according to cell type clusters.",br(),
                                "User can choose an additional condition factor based on which",
                                "cells will be arranged in refined groups.")))
    else if (input$tab == "Pie")
      return(wellPanel(helpText("Pie chart showing percentage of cells in each cluster.",br(),
                                "User can choose an additional condition factor based on which",
                                "pie chart will be generated per condition.")))
    else if (input$tab == "Density")
      return(wellPanel(helpText("Density plot showing distribution of expression for a given gene.",br(),
                                "User can choose an additional condition factor based on which",
                                "Density plot will be generated highlighting differences among conditions.")))
    else if (input$tab == "BarLine")
      return(wellPanel(helpText("Bar&Line plot showing percentage of cells expressing a given gene and their median expression.",br(),
                                "User can choose an additional condition factor based on which",
                                "Bar&Line plot will be generated highlighting differences among conditions.")))
    else
      return()
  })
  # load expression data if available
  getExpData <- reactive({
    expFile <- input$expFile
    geneFile <- input$geneFile
    cellFile <- input$cellFile
    if (is.null(expFile) | is.null(geneFile) | is.null(cellFile))
      return(NULL)
    #return(read.table(expFile$datapath, header=T, check.names=F, row.names=1, sep="\t"))
    return(myLoadExpression(expFile$datapath, geneFile$datapath, cellFile$datapath))
  })
  # load cluster info data if available
  getClustData <- reactive({
    clustFile <- input$clustFile
    if (is.null(clustFile))
      return(NULL)
    return(read.table(clustFile$datapath, header=T, check.names=F, row.names=1, sep="\t"))
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
    tclusters <- unique(as.vector(clustData$ident))
    # reorder clusters and return
    return(tclusters[order(tclusters)])
  })
  # fetch current selected gene
  getGene <- reactive({
    # get gene list
    myGenes <- getGeneList()
    if (is.null(myGenes))
      return(NULL)
    if (is.null(input$gene)){
      #return(intersect(c("Gapdh","GAPDH"),myGenes))
      return(NULL)
    }
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
    # get clust info data
    clustData <- getClustData()
    if (is.null(clustData))
      return(NULL)
    # get non-numeric columns
    selCols <- colnames(clustData)[!sapply(clustData, is.numeric)]
    # remove columns "ident", "orig.ident"
    return(selCols[!selCols %in% c("ident","orig.ident","Row.names")])
  })
  # extract data columns for plotting
  prepareData <- reactive({
    # get expression data
    expData <- getExpData()
    if (is.null(expData))
      return(NULL)
    # get clust info data
    clustData <- getClustData()
    if (is.null(clustData))
      return(NULL)
    # get current selected gene
    curGene <- getGene()
    # get current selected condition column name
    curCondition <- getCondition()
    # initial data
    combinedData <- clustData
    # add expression of the selected gene
    if (is.null(curGene)){
      combinedData$gene <- rep(0, nrow(clustData))
    } else {
      combinedData$gene <- as.vector(expData[curGene, rownames(clustData)])
    }
    combinedData$Row.names <- rownames(combinedData)
    # extract needed columns and
    # add a new column by combining the "ident" column and the selected condition column
    if (curCondition != "None"){
      dataToPlot <- combinedData[,c("Row.names","ident","tSNE_1","tSNE_2","gene",curCondition)]
      colnames(dataToPlot) <- c("cellID","ident","tSNE_1","tSNE_2","gene","condition")
      dataToPlot$cohort <- paste(dataToPlot$ident, dataToPlot$condition, sep=" - ")
    }
    else{
      dataToPlot <- combinedData[,c("Row.names","ident","tSNE_1","tSNE_2","gene")]
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
    # factorize ident
    dataToPlot$ident <- factor(dataToPlot$ident)
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
    g <- g + theme(legend.justification=c(0,0), legend.title=element_blank(), legend.text=element_text(size=12))
    g <- g + theme(axis.text=element_text(size=15), axis.title=element_text(size=16,face="bold"))
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
    # in case a non-expressed gene, use gray color
    g <- ggplot(dataToPlotOrdered[dataToPlotOrdered$ident %in% curCluster, ], aes(x=tSNE_1,y=tSNE_2))
    g <- g + geom_point(shape=19, size=3, alpha=.8, color="grey")
    # otherwise
    if (hcut > 0){
      g <- ggplot(dataToPlotOrdered[dataToPlotOrdered$ident %in% curCluster, ], aes(x=tSNE_1,y=tSNE_2,color=gene))
      g <- g + geom_point(shape=19, size=3, alpha=.8)
      g <- g + scale_color_gradient(low="grey", high="red", limits=c(0,hcut))
    }
    g <- g + coord_cartesian(xlim=c(lower.x, upper.x), ylim=c(lower.y, upper.y))
    g <- g + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
    g <- g + theme(legend.justification=c(0,0), legend.title=element_blank())
    g <- g + theme(legend.key=element_blank()) + theme(legend.text=element_text(size=15))
    g <- g + theme(axis.text=element_text(size=15), axis.title=element_text(size=16,face="bold"))
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
    g <- g + theme(axis.line = element_line(color = "black")) + theme(strip.text.x=element_text(angle=90, size=12, face="bold"))
    g <- g + theme(legend.key=element_blank()) + theme(legend.text=element_text(size=10)) + theme(axis.text=element_text(size=15), axis.title=element_text(size=16,face="bold"))
    g <- g + ylab("Expression (log-scale)") + theme(axis.title.x=element_blank(),axis.line.x=element_blank())
    return(g)
  })
  # draw expression violin plot
  output$violinExp <- renderPlot({
    # get ordered data for plotting
    dataToPlotOrdered <- orderDataByCohort()
    if (is.null(dataToPlotOrdered))
      return(NULL)
    # generate plot
    g <- ggplot(dataToPlotOrdered, aes(x=cohort, y=gene, fill=cohort))
    g <- g + geom_violin() + stat_summary(fun.data="median_hilow", geom="pointrange")
    g <- g + theme(legend.position="none") + theme(axis.text.x=element_text(angle=90, size=12, face="bold"))
    g <- g + ylab("Expression (log-scale)") 
    g <- g + theme(axis.title=element_text(size=16,face="bold"), axis.text=element_text(size=15), axis.title.x=element_blank())
    return(g)
  })
  # draw cell number pie chart
  output$pieCell <- renderPlot({
    # make sure expression file is also loaded though it is not necessary for pie chart
    if (is.null(input$expFile))
      return(NULL)
    # get cluster info
    clustData <- getClustData()
    if (is.null(clustData))
      return(NULL)
    # get current selected condition column name
    curCondition <- getCondition()
    # make pie chart
    if (curCondition != "None"){
      dataToPlot <- as.data.frame(table(clustData[,c("ident",curCondition)]))
      colnames(dataToPlot) <- c("Cluster","Condition","Freq")
      # calculate percentage
      dataToPlot.per <- as.data.frame(group_by(dataToPlot, Condition) %>% 
                                        mutate(Percent=round(Freq/sum(Freq)*100,digits=1),
                                               Pos=100-cumsum(Percent)+Percent/2,
                                               Label=paste(Percent,"%",sep="")))
      g <- ggplot(dataToPlot.per, aes(x="", y=Percent, fill=Cluster))
      g <- g + geom_bar(width=1, stat="identity") + facet_wrap(~Condition)
      g <- g + coord_polar("y", start=0, direction=-1) + theme_minimal()
      g <- g + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
      g <- g + theme(panel.border=element_blank(), panel.grid=element_blank())
      g <- g + theme(axis.ticks=element_blank(), axis.text.x=element_blank())
      g <- g + theme(strip.text=element_text(size=16, face="bold"))
      ##g <- g + scale_y_continuous(breaks=cumsum(dataToPlot.per$Freq) - dataToPlot.per$Freq/2, labels=dataToPlot.per$Percent)
      g <- g + geom_text_repel(aes(y=Pos, label=Label), size=5)
      return(g)
    }
    else{
      dataToPlot <- as.data.frame(table(clustData[,c("ident")]))
      colnames(dataToPlot) <- c("Cluster","Freq")
      # calculate percentage
      dataToPlot$Percent <- round(dataToPlot$Freq/sum(dataToPlot$Freq)*100,digits=1)
      dataToPlot$Pos <- 100-cumsum(dataToPlot$Percent)+dataToPlot$Percent/2
      dataToPlot$Label <- paste(dataToPlot$Percent,"%",sep="")
      # plot
      g <- ggplot(dataToPlot, aes(x="", y=Percent, fill=Cluster))
      g <- g + geom_bar(width=1, stat="identity")
      g <- g + coord_polar("y", start=0, direction=-1) + theme_minimal()
      g <- g + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
      g <- g + theme(panel.border=element_blank(), panel.grid=element_blank())
      g <- g + theme(axis.ticks=element_blank(), axis.text.x=element_blank())
      g <- g + theme(strip.text=element_text(size=16, face="bold"))
      g <- g + geom_text_repel(aes(y=Pos, label=Label), size=5)
      return(g)
    }
  })
  # draw expression density plot
  output$densityExp <- renderPlot({
    # get data for plotting
    dataToPlot <- prepareData()
    if (is.null(dataToPlot))
      return(NULL)
    # get current selected condition column name
    curCondition <- getCondition()
    # make density plot
    if (curCondition != "None"){
      g <- ggplot(dataToPlot, aes(x=gene, fill=condition)) + geom_density(alpha=0.5)
      g <- g + facet_wrap(~ident)
      g <- g + theme_classic()
      g <- g + xlab("Expression (log-scale)")
      g <- g + theme(strip.background=element_rect(color="white", fill=NA), panel.border=element_rect(color="lightgray", fill=NA))
      g <- g + theme(axis.line = element_line(color = "black")) + theme(strip.text.x=element_text(size=12, face="bold"))
      return(g)
    }
    else{
      g <- ggplot(dataToPlot, aes(x=gene, fill=ident)) + geom_density(alpha=0.5)
      g <- g + theme_classic()
      g <- g + xlab("Expression (log-scale)") 
      return(g)
    }
  })
  # draw bar&line plot
  output$barLineExp <- renderPlot({
    # get data for plotting
    dataToPlot <- prepareData()
    if (is.null(dataToPlot))
      return(NULL)
    # get current selected condition column name
    curCondition <- getCondition()
    # get current selected gene
    curGene <- getGene()
    # make density plot
    if (curCondition != "None"){
      # calculate percentage of cells expressing a gene for each cluster+condition
      tpercent <- aggregate(gene ~ ident+condition, data=dataToPlot, FUN=function(x) sum(x>0)/length(x))
      colnames(tpercent) <- c("ident","condition","percent")
      # calculate median expression value for the above sets of cells
      texp <- aggregate(gene ~ ident+condition, data=dataToPlot, FUN=function(x) median(x[x>0]))
      colnames(texp) <- c("ident","condition","medExp")
      # fix NA items (i.e. no cells in the given cluster&condition express the given gene)
      texp[is.na(texp)] <- 0
      # merge results
      tdata <- merge(tpercent, texp, by=c("ident","condition"), all=T)
      # draw plot
      tfactor <- ceiling(max(tdata$medExp)/max(tdata$percent))
      g <- ggplot(tdata)
      g <- g + geom_bar(aes(x=condition, y=percent), stat="identity", fill="lightsalmon", width=0.6)
      g <- g + geom_line(aes(x=as.numeric(condition), y=medExp/tfactor), size=1, color="steelblue2")
      g <- g + geom_point(aes(x=condition, y=medExp/tfactor), shape=19, size=2, color="steelblue3")
      g <- g + facet_wrap(~ident)
      g <- g + scale_y_continuous(name=paste("Fraction of cells expressing",curGene,sep=" "), sec.axis=sec_axis(~.*tfactor, name="Median expression in log-scale"))
      g <- g + theme_classic()
      g <- g + theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12, face="bold", angle=90))
      g <- g + theme(axis.title.y=element_text(size=15, face="bold"), axis.text.y=element_text(size=12, face="bold"))
      g <- g + theme(strip.background=element_rect(color="white", fill=NA), panel.border=element_rect(color="lightgray", fill=NA))
      return(g)
    }
    else{
      # calculate percentage of cells expressing a gene for each cluster
      tpercent <- aggregate(gene ~ ident, data=dataToPlot, FUN=function(x) sum(x>0)/length(x))
      colnames(tpercent) <- c("ident","percent")
      # calculate median expression value for the above sets of cells
      texp <- aggregate(gene ~ ident, data=dataToPlot, FUN=function(x) median(x[x>0]))
      colnames(texp) <- c("ident","medExp")
      # fix NA items (i.e. no cells in the given cluster express the given gene)
      texp[is.na(texp)] <- 0
      # merge results
      tdata <- merge(tpercent, texp, by="ident", all=T)
      # draw plot
      tfactor <- ceiling(max(tdata$medExp)/max(tdata$percent))
      g <- ggplot(tdata)
      g <- g + geom_bar(aes(x=ident, y=percent), stat="identity", fill="lightsalmon", width=0.6)
      g <- g + geom_line(aes(x=as.numeric(ident), y=medExp/tfactor), size=1, color="steelblue2")
      g <- g + geom_point(aes(x=ident, y=medExp/tfactor), shape=19, size=2, color="steelblue3")
      g <- g + scale_y_continuous(name=paste("Fraction of cells expressing",curGene,sep=" "), sec.axis=sec_axis(~.*tfactor, name="Median expression in log-scale"))
      g <- g + theme_classic()
      g <- g + theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12, face="bold", angle=90))
      g <- g + theme(axis.title.y=element_text(size=15, face="bold"), axis.text.y=element_text(size=12, face="bold"))
      return(g)
    }
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
