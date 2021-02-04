PCA_tab <- tabPanel('PCA',
                    sidebarPanel(
                      conditionalPanel(
                        condition = "input.pca_tabs == 'PCA-2D plot'",
                        selectInput(inputId = 'pca.x',label = 'X-axis',choices = ""),
                        selectInput(inputId = 'pca.y',label = 'Y-axis',choices = "")
                      ),
                      selectInput(inputId = 'gene_size',label = 'Gene sample size',choices = ""),
                      radioButtons('gene_order',"Gene sample order (wrt column 1)",
                                   c('Descending (highest to lowest)'='Descending','Ascending (lowest to highest)'='Ascending','Random')),
                      conditionalPanel(
                        condition = "input.pca_tabs == 'PCA-2D plot' || input.pca_tabs == 'PCA-3D plot'",
                        checkboxInput('pca_cluster',strong('Kmeans clustering on columns'),FALSE),
                        conditionalPanel(
                          condition = "input.pca_cluster == true",
                          sliderInput("pca_cluster_num","Number of clusters:",value=1,min=1,max=1,step=1),
                          checkboxInput('pca_text',strong('Display sample name'),FALSE)
                        )
                      ),
                      conditionalPanel(
                        condition = "input.gene_order=='Random'",
                        helpText('* Click multiple times to resample'),
                        actionButton('pca_refresh',"Resample",style="background-color: #337ab7;border-color:#337ab7"),
                        br(),br()
                      ),
                      conditionalPanel(
                        condition = "input.pca_tabs == 'PCA variance'",
                        downloadButton("downloadpcavar", "Download as PNG")
                      ),
                      conditionalPanel(
                        condition = "input.pca_tabs == 'PCA-2D plot'",
                        downloadButton("downloadpca2d","Download as PNG")
                      ),
                      conditionalPanel(
                        condition = "input.pca_tabs == 'PCA-3D plot'",
                        downloadButton("downloadpca3d","Download as PNG")
                      )
                    ),
                    mainPanel(
                      tabsetPanel(type = "tabs",id="pca_tabs",
                                  tabPanel("PCA variance", plotlyOutput("pcavar.plot")),
                                  tabPanel("PCA-2D plot", plotlyOutput("pca2d.plot")),
                                  tabPanel("PCA-3D plot",plotlyOutput("pca3d.plot"))
                      )
                    ))