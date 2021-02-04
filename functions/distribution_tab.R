distribution_tab <- tabPanel('Distribution Fit',
         sidebarPanel(
           conditionalPanel(
             condition= "input.dist_tabs=='Distribution Fit'",
             selectInput(inputId = 'dist.var',label = 'Choose a column',choices = colnames('dataset')),
             checkboxGroupInput("distributions", "Distributions:",
                                choices = c("Log-normal","Log-logistic","Pareto","Burr","Weibull","Gamma"),selected = c("Log-normal","Pareto")),
             radioButtons('dist_zoom',"Zoom to see fit",c('slider','text input')),
             conditionalPanel(
               condition = "input.dist_zoom=='slider'",
               sliderInput("dist_range", "Range:",
                           min = 0.1, max = 1000,step=1,
                           value = c(0.1,1000))
             ),
             conditionalPanel(
               condition = "input.dist_zoom=='text input'",
               textOutput('dist_range_allowed'),
               numericInput('dist_range_min',"min",value=0.1,min=0.1,max=1000),
               numericInput('dist_range_max',"max",value=1000,min=0.1,max=1000)
             ),
             downloadButton("downloaddist", "Download as PDF")
           ),
           conditionalPanel(
             condition = "input.dist_tabs=='AIC table'",
             downloadButton("downloaddistaic", "Download as CSV")
           )
         ),
         mainPanel(
           tabsetPanel(type = "tabs",id="dist_tabs",
                       tabPanel("Distribution Fit", plotOutput("dist.plot")),
                       tabPanel("AIC table",
                                conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                 div(img(src="load.gif",width=240,height=180),
                                                     h4("Processing ... Please wait"),style="text-align: center;")
                                ), 
                                conditionalPanel(condition="!$('html').hasClass('shiny-busy')",
                                                 div(tableOutput('dist.aic'), style = "font-size:80%")
                                ))
           )
         ))