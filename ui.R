#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that make plots for single cell data
shinyUI(fluidPage(
  verticalLayout(
    # Application title
    titlePanel("Single Cell Data Viewer"),
    
    tabsetPanel(id="tab", type="tabs",
                tabPanel(title="Pie", plotOutput("pieCell")),
                tabPanel(title="tSNE", fluidRow(column(6,plotOutput("tsneExp", width="95%", height="400px")),
                                                column(6,plotOutput("tsneRef", width="100%", height="400px")))),
                tabPanel(title="UMAP", fluidRow(column(6,plotOutput("umapExp", width="95%", height="400px")),
                                                column(6,plotOutput("umapRef", width="100%", height="400px")))),
                tabPanel(title="Bar", plotOutput("barExp")),
                tabPanel(title="Violin", plotOutput("violinExp")),
                tabPanel(title="Density", plotOutput("densityExp")),
                tabPanel(title="BarLine", plotOutput("barLineExp")),
                tabPanel(title="About", 
                         h4("Thank you for using Single Cell Data Viewer!"),
                         HTML('Please address questions and comments to <a href="mailto:taz2008@med.cornell.edu" >Tuo Zhang</a>.'))
                )
    ),

    hr(),
  
    # Bottom bar with file uploader and parameter selectors 
    fluidRow(
      column(4,
             wellPanel(fileInput(inputId="expFile",
                                 label="upload expression matrix:",
                                 multiple=F),
                       fileInput(inputId="geneFile",
                                 label="upload gene list:",
                                 multiple=F),
                       fileInput(inputId="cellFile",
                                 label="upload cell list:",
                                 multiple=F),
                       fileInput(inputId="clustFile",
                                 label="upload annotations:", multiple=F)
             )
      ),
      column(4,
             conditionalPanel(condition="output.fileUploaded", uiOutput("myMenu"))
      ),
      column(4,
             ##h3("Test")
             uiOutput("myDescription")
      )
    )
  )
)
