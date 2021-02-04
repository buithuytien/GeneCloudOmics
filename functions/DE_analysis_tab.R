DE_analysis_tab <- tabPanel("DE Analysis",
                            # useShinyjs(),
                            sidebarPanel(
                              radioButtons("n_rep","Replicates?",choices=c("Multiple"=1,"Single"=0)),
                              conditionalPanel(
                                condition="input.n_rep=='1'",
                                radioButtons("de_method1","DE Method",choices=c("EdgeR","DESeq2","NOISeq"))
                              ),
                              conditionalPanel(
                                condition="input.n_rep=='0'",
                                radioButtons("de_method0","DE Method",choices=c("NOISeq"))
                              ),
                              h5("Choose 2 experiment conditions for DE analysis"),
                              selectInput("f1","Condition 1",choices = ""),
                              selectInput("f2","Condition 2",choices = ""),
                              
                              h5("DE criteria"),
                              splitLayout(
                                numericInput("p_val","FDR",min=0.01,max=1,value=0.05,step=0.01),
                                numericInput("fc","Fold Change",min=1,value=2,step=0.1)
                              ),
                              fluidRow(
                                column(4,
                                       actionButton("submit_DE","Submit")
                                ),
                                column(6,
                                       conditionalPanel(
                                         condition = "input.DE_tabs=='DE genes' ",
                                         downloadButton("download_de_table","Download table (csv)")
                                       ),
                                       conditionalPanel(
                                         condition = "input.DE_tabs=='Volcano plot' ",
                                         downloadButton("download_volcano","Download plot (PDF)")
                                       ),
                                       conditionalPanel(
                                         condition = "input.DE_tabs=='Dispersion plot' ",
                                         downloadButton("download_dispersion","Download plot (PDF)")
                                       )
                                       # conditionalPanel(
                                       #   condition = "input.DE_tabs=='Heatmap plot' ",
                                       #   downloadButton("download_heatmap","Download plot")
                                       # )
                                )
                              )
                              
                            ),
                            mainPanel(
                              tabsetPanel(type = "tabs", id= "DE_tabs",
                                          tabPanel("DE genes",
                                                   # h3("Differential Expression Analysis"),
                                                   conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                                    div(img(src="load.gif",width=240,height=180),
                                                                        h4("Processing ... Please wait"),style="text-align: center;")
                                                   ), 
                                                   conditionalPanel(condition="!$('html').hasClass('shiny-busy')",
                                                                    DT::dataTableOutput("DE_table")
                                                   )
                                          ),
                                          tabPanel("Volcano plot",    # for DESeq and edgeR
                                                   h6("Volcano plot is only available for edgeR and DESeq2 methods"),
                                                   conditionalPanel(
                                                     condition = "input.n_rep=='1' && input.method1!='NOISeq'",
                                                     conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                                      div(img(src="load.gif",width=240,height=180),
                                                                          h4("Processing ... Please wait"),style="text-align: center;")
                                                     ), 
                                                     conditionalPanel(condition="!$('html').hasClass('shiny-busy')",
                                                                      plotOutput("volcano_plot") 
                                                     )
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.method0=='NOISeq' || input.method1=='NOISeq'",
                                                     h6("Volcano Plot is only applicable to DESeq2 and edgeR")
                                                   )
                                          ),
                                          tabPanel("Dispersion plot", # for edgeR
                                                   h6("Dispersion plot is only available for edgeR and DESeq2 methods"),
                                                   conditionalPanel(
                                                     condition = "input.n_rep=='1' && input.method1!='NOISeq'",
                                                     # h3("Dispersion plot"),
                                                     conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                                      div(img(src="load.gif",width=240,height=180),
                                                                          h4("Processing ... Please wait"),style="text-align: center;")
                                                     ), 
                                                     conditionalPanel(condition="!$('html').hasClass('shiny-busy')",
                                                                      plotOutput("dispersion_plot")
                                                     )
                                                   ),
                                                   conditionalPanel(
                                                     condition = "input.method0=='NOISeq' || input.method1=='NOISeq'",
                                                     h6("Dispersion Plot is only applicable to DESeq2 and edgeR")
                                                   )
                                          )
                              )
                            )
)