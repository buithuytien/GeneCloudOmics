scatter_tab <- tabPanel('    Scatter    ',
         sidebarPanel(
           selectInput(inputId = 'scatter.x',label = 'X-axis',choices = ""),
           selectInput(inputId = 'scatter.y',label = 'Y-axis',choices = ""),
           radioButtons('trans',"Transformation:",
                        c('None','Natural log','log2','log10')),
           downloadButton("downloadscatter", "Download as PDF"),
           h6('Download all pairs of samples in one PDF (this may take some time to run) :'),
           downloadButton("downloadscatter_collage","Download collage")),
         mainPanel(
           h3('Heatscatter'),
           plotOutput('scatter.plot')
         ))