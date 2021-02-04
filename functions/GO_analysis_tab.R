GO_analysis_tab <- ################## GO analysis ###################
##################################################
tabPanel("GO Analysis",
         # useShinyjs(),
         sidebarPanel(
           conditionalPanel(
             condition = "input.go_tab == 'go_table' || input.go_tab == 'go_pie' ",
             helpText("* One-column csv file"),
             fileInput("filego","Upload list of DE genes"),
             fileInput("filebg","List of background genes"),
             selectInput("go_method","Select GO package", choices = c("clusterProfiler","GOstats","enrichR")), #new
             conditionalPanel(
               condition = "input.go_method=='clusterProfiler' || input.go_method=='GOstats'",
               # update if enrichR
               selectInput('go_species',"Select species",selected="org.EcK12.eg.db",choices=c("Homo sapiens"='org.Hs.eg.db',"Mus musculus"='org.Mm.eg.db',"Rattus norvegicus"='org.Rn.eg.db',"Gallus gallus"='org.Gg.eg.db',"Danio rerio"='org.Dr.eg.db',"Drosophila melanogaster"='org.Dm.eg.db',"Caenorhabditis elegans"='org.Ce.eg.db',"Saccharomyces cereviasiae"='org.Sc.sgd.db',"Arabidopsis thaliana"='org.At.tair.db',"Escherichia coli (strain K12)"='org.EcK12.eg.db',"Escherichia coli (strain Sakai)"='org.EcSakai.eg.db',"Anopheles gambiae"='org.Ag.eg.db',"Bos taurus"='org.Bt.eg.db',"Canis familiaris"='org.Cf.eg.db',"Macaca mulatta"='org.Mmu.eg.db',"Plasmodium falciparum"='org.Pf.plasmo.db',"Pan troglodytes"='org.Pt.eg.db',"Sus scrofa"='org.Ss.eg.db',"Xenopus tropicalis"='org.Xl.eg.db')), # species.choices
               selectInput('go_geneidtype',"Select identifier",choices=NULL),
               selectInput('subontology',"Select subontology",selected="BP",choices = c("biological process"="BP","molecular function"="MF","cellular component"="CC")),
               numericInput("go_max_p","Adjusted p cutoff",value=0.05,min=0,max=1, step=0.05)
               # selectInput('go_level',"Select level",choices = c(2:10))
             ), conditionalPanel(
               condition = "input.go_method == 'enrichR'",
               selectInput('enrichR_dbs',"Select database",choices='enrichRdbs')
             ),
             numericInput("go_min_no","Min. number of genes",value=3,min=1,step=1),
             br(),
             fluidRow(
               column(4,
                      actionButton("submit_go","Submit")
               ),
               conditionalPanel(
                 condition = "input.go_tab == 'go_table' ",
                 column(7,
                        downloadButton("download_go_table","Save CSV")
                 )
               ),
               conditionalPanel(
                 condition = "input.go_tab == 'go_pie' ",
                 column(7,
                        downloadButton("download_go_pie","Save PNG")
                 )
               )
             )
           ),
           conditionalPanel(
             condition = "input.go_tab == 'go_graph'",
             selectInput("go_term_slect","Select GO terms",choices="",multiple=T),
             checkboxInput("show_gene_names","Show gene names", value = F),
             actionButton("submit_go_graph","Submit"),
             br(),
             fluidRow(
               column(5,
                      br(),
                      radioButtons("download_go_graph_type","File Type",choices=c("png","pdf"))
               ),
               column(6,
                      br(),
                      br(),
                      downloadButton("download_go_graph","Save Graph")
               )
             )
           )
         ),
         mainPanel(
           tabsetPanel(type = "tabs", id = "go_tab",
                       tabPanel("go_table",
                                conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                 div(img(src="load.gif",width=240,height=180),
                                                     h4("Processing ... Please wait"),style="text-align: center;")
                                ), 
                                conditionalPanel(condition="!$('html').hasClass('shiny-busy')",
                                                 DT::dataTableOutput("go_table")
                                )
                       ),
                       tabPanel("go_pie",
                                conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                 div(img(src="load.gif",width=240,height=180),
                                                     h4("Processing ... Please wait"),style="text-align: center;")
                                ), 
                                conditionalPanel(condition="!$('html').hasClass('shiny-busy')",
                                                 h6("Pie chart is only available for clusterProfiler or GOstats method"),
                                                 h4("Relative size of GO terms level 2"),
                                                 plotlyOutput("go_pie")
                                )
                       ),
                       tabPanel("go_graph",
                                h6("Graph visualization is only available for clusterProfiler method"),
                                conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                 div(img(src="load.gif",width=240,height=180),
                                                     h4("Processing ... Please wait"),style="text-align: center;")
                                ), 
                                conditionalPanel(condition="!$('html').hasClass('shiny-busy')",
                                                 plotOutput("go_graph") 
                                )
                       )
           )
         )
)