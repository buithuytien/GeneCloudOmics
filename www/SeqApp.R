library(shiny)
library(shinycssloaders)
library(shinyFiles)
library(ggplot2)
library(UniprotR)

uiseq <- fluidPage(sidebarLayout(
  sidebarPanel(
    fileInput(
      'file1',
      "Upload protein's accessions",
      accept = c(
        'text/csv',
        'text/comma-separated-values',
        'text/tab-separated-values',
        'text/plain',
        '.csv'
      )
    ),
    br(),
    selectInput(
      "select",
      label = h3("Select choice"),
      choices = list(
        "Charge Plot" = 1,
        "Gravy index Plot" = 2,
        "Physicochemical properties" = 3
      ),
      selected = 1
    ),
    actionButton("SubmitSeq", "Submit"),
    
  ),
  mainPanel(h3("Protein Sequences"),
            uiOutput("help_text_prot_seq"), 
            withSpinner(plotOutput("plot")))
))

ServerSeq <- function(input, output)
{
  output$help_text_prot_seq <- renderUI({
    HTML(
      "<h3><b>This page retrieves the full protein sequences from <a href ='https://www.uniprot.org/'>UniProt.org</a> of a given set of UniProt accessions.</b></h3>"
    )
  })
  SeqObj <- NULL
  GetAccs <- function()
  {
    Accessions <-  read.csv(input$file1$datapath)
    Accessions <- unique(as.character(Accessions[, 1]))
    return(Accessions)
  }
 output$plot <- renderPlot(
   {
     SeqObj <- GetSequences(GetAccs())
     
     if(input$select == 1)
     {
       PlotCharge(SeqObj)
     }
     else if (input$select == 2)
     {
      PlotGravy(SeqObj) 
     }
     else {
       PlotPhysicochemical(SeqObj)
     }
   }
 )
}
shinyApp(ui = uiseq, server = ServerSeq)
