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
                tabPanel(title="tSNE", fluidRow(column(6,plotOutput("tsneExp")),
                                                column(6,plotOutput("tsneRef")))),
                tabPanel(title="Bar", plotOutput("barExp")),
                tabPanel(title="Violin", plotOutput("violinExp")),
                tabPanel(title="About", h4("This panel is intentionally left blank"), h3("test"))
                )
    ),

    hr(),
  
    # Bottom bar with file uploader and parameter selectors 
    fluidRow(
      column(4,
             wellPanel(fileInput(inputId="expFile",
                                 label="upload expression data:",
                                 multiple=F),
                       fileInput(inputId="clustFile",
                                 label="upload cluster data:", multiple=F)
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
