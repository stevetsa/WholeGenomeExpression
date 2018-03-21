library(shiny)
library(DT)
library(shinyjs)
library(plotly)

options(shiny.maxRequestSize = 500*1024^2)

shinyUI(
  fluidPage(
    useShinyjs(),
      fluidRow(align='Top',
      column(2,
             # fileInput("SRRcode1", label= h6("SRR"), width="300px", multiple = F)
             textInput('SRRcode1', label=h6('SRR'))
      ),
      column(2,
             # fileInput("SRRcode2", label = h6("SRR"), width = "300px", multiple = F)
             textInput('SRRcode2', label = h6('SRR'))
             ),
      column(2,
             fileInput("listOfGenes", label=h6("Select List of Genes of Interest"),multiple =F)
      ),
      br(),
      br(),
      actionButton(inputId="SRAbutton", label="Go")
    ),
    br(),
    br(),
    br(),
    br(),
    mainPanel(
        navbarPage(title = '',
             navbarMenu(title = 'Visualization',
                tabPanel('Normalized Gene Expression', DT::dataTableOutput("sra")),
                tabPanel('Dot Plot', textInput("geneName", h6("Gene name"), value='', width = "200px"),
                  actionButton(inputId = "dotPlotButton", label = "Go"),
                  plotlyOutput("dotPlot",height = "400px",width = "400px")),
                tabPanel('FastQC', htmlOutput('inc'))))
    )
  )
)
 

 




    