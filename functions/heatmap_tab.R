heat_map_tab <- tabPanel('Heatmap',
         sidebarPanel(
           conditionalPanel(
             condition = "input.heatmap_tabs=='Heatmap'",
             
             radioButtons("heatmap_de_ind",label="Choose data",choices=c("Indenpendent"="ind","DE result"="de")),
             numericInput('numOfCluster',"Number of clusters on rows",value=2,min=2,max=30,step=1),
             conditionalPanel(
               condition = "input.heatmap_de_ind == 'ind' ",
               # selectInput('numOfGeno',"Number of genotypes (mutants)",choices=c(1)),
               splitLayout(
                 numericInput('fold',"Fold change",value=2,min=1,step=1),
                 numericInput("fold_ncol", "min. column",value=2,min=1,step=1)
               )
               
               # uiOutput("refGeno"),
               # radioButtons('heatmap_value',"Values",
               #              c('Fold change','Log fold change'))
             ),
             
             downloadButton("downloadheatmap","Download as PDF"),
             actionButton('heatmap_plot',"Plot",width='65px',style="color: #fff; background-color: #337ab7; border-color: #337ab7;float:right")
             
             # conditionalPanel(
             #   condition = "input.heatmap_de_ind == 'ind' ",
             #   h5('Specify names of the genotypes'),
             #   uiOutput("expand_genonames")
             # )
           ),
           
           conditionalPanel(
             condition = "input.heatmap_tabs=='Gene clusters'",
             uiOutput("heatmap_display"),
             conditionalPanel(
               condition = "input.display_cluster=='ALL'",
               downloadButton("downloadclusters","Download as CSV")
             )
           )
           
         ),
         mainPanel(
           tabsetPanel(type = "tabs",id="heatmap_tabs",
                       tabPanel("Heatmap", 
                                conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                 div(img(src="load.gif",width=240,height=180),
                                                     h4("Processing ... Please wait"),style="text-align: center;")
                                ), 
                                conditionalPanel(condition="!$('html').hasClass('shiny-busy')",
                                                 plotOutput("heatmap.plot")
                                )),
                       tabPanel("Gene clusters", dataTableOutput('cluster.info'))
           )
         ))