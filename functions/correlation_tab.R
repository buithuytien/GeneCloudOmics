correlation_tab <- tabPanel('  Correlation  ',
         sidebarPanel(
           radioButtons('cor_method',"Method:",
                        c('Pearson correlation','Spearman correlation')),
           conditionalPanel(
             condition = "input.cor_tabs == 'Correlation heatmap'",
             downloadButton("downloadcorrplot", "Download as PDF")
           ),
           conditionalPanel(
             condition = "input.cor_tabs == 'Correlation plot'",
             downloadButton("downloadcorrplot2", "Download as PDF")
           ),
           conditionalPanel(
             condition = "input.cor_tabs == 'Correlation matrix'",
             downloadButton("downloadcorrmat","Download as CSV")
           )
         ),
         mainPanel(
           conditionalPanel(
             condition = "input.cor_method=='Pearson correlation'",
             h3('Pearson correlation')
           ),
           conditionalPanel(
             condition = "input.cor_method=='Spearman correlation'",
             h3('Spearman correlation')
           ),
           tabsetPanel(type = "tabs",id="cor_tabs",
                       tabPanel("Correlation heatmap", plotOutput('corr.plot')),
                       tabPanel("Correlation plot", plotOutput('corr.plot2')),
                       tabPanel("Correlation matrix", div(tableOutput('corr.matrix'), style = "font-size:80%"))
           )
         )
)