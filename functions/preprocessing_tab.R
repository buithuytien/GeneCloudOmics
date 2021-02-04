preprocessing_tab <- tabPanel('Preprocessing',
                              sidebarPanel(
                                h4("Filtering"),
                                splitLayout(
                                  numericInput("min_val","Min. value", min=0.1,step=0.1,value=1.0),
                                  numericInput("min_col","Min. columns", min=1, value=2)
                                ),
                                conditionalPanel(
                                  condition = "input.file_type=='raw'",
                                  radioButtons('norm_method',"Normalisation method",
                                               c("None (Black)"="None",
                                                 'RPKM (Blue)'='RPKM','FPKM (Dark cyan)'='FPKM',
                                                 'TPM (Dark green)'='TPM',
                                                 "RUV (Brown)"='RUV'))
                                ),
                                actionButton("submit_preprocessing","Submit"),
                                conditionalPanel(
                                  condition = "input.preprocessing_tabs == 'Data table' ",
                                  br(),
                                  br(),
                                  downloadButton("download_norm_data", "Download table (csv)")
                                )
                              ),
                              mainPanel(
                                tabsetPanel(type = "tabs",id="preprocessing_tabs",
                                            tabPanel("RLE plot",
                                                     conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                                      div(img(src="load.gif",width=240,height=180),
                                                                          h4("Processing ... Please wait"),style="text-align: center;")
                                                     ), 
                                                     conditionalPanel(condition="!$('html').hasClass('shiny-busy')",
                                                                      plotOutput("RLE.plot2")
                                                     ),
                                                     
                                                     conditionalPanel(
                                                       condition = "input.file_type=='raw'",
                                                       conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                                        div(img(src="load.gif",width=240,height=180),
                                                                            h4("Processing ... Please wait"),style="text-align: center;")
                                                       ), 
                                                       conditionalPanel(condition="!$('html').hasClass('shiny-busy')",
                                                                        plotOutput("RLE.plot")
                                                       )
                                                     )
                                            ),
                                            tabPanel("Data table",
                                                     h3("Normalized data"),
                                                     DT::dataTableOutput("norm_table")
                                            ),
                                            tabPanel("Description table",
                                                     h3("Data description"),
                                                     DT::dataTableOutput("meta_table")
                                            )
                                )
                              )
)