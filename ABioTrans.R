#Sys.setenv("plotly_username"=" your_plotly_username")
#Sys.setenv("plotly_api_key"="your_api_key")
## test repo

print("start loading")
start.load <- Sys.time()   ### time

if(length(find.package(package = 'shiny',quiet = T))>0){
  library(shiny)
}else{
  print("Package shiny not installed")
  install.packages("shiny")
  print("Package shiny installed")
  library(shiny)
}


if(length(find.package(package = 'shinythemes',quiet = T))>0){
  library(shinythemes)
}else{
  print("Package shinythemes not installed")
  install.packages("shinythemes")
  print("Package shinythemes installed")
  library(shinythemes)
}

if(length(find.package(package = 'rstudioapi',quiet = T))>0){
  library(rstudioapi)
}else{
  install.packages("rstudioapi")
  library(rstudioapi)
}

wd <- dirname(rstudioapi::getActiveDocumentContext()$path)  #set wd as the current folder
print(wd == getwd())
print(wd)
print(getwd())
if(! wd == getwd()){
  setwd(wd)
}

# 
# ## sourcing util files
source(paste0("./www/utils.R"))
# source("ui.R")
# 
loadPkg()

species.choices <<- c("Homo sapiens"='org.Hs.eg.db',"Mus musculus"='org.Mm.eg.db',"Rattus norvegicus"='org.Rn.eg.db',"Gallus gallus"='org.Gg.eg.db',"Danio rerio"='org.Dr.eg.db',"Drosophila melanogaster"='org.Dm.eg.db',"Caenorhabditis elegans"='org.Ce.eg.db',"Saccharomyces cereviasiae"='org.Sc.sgd.db',"Arabidopsis thaliana"='org.At.tair.db',"Escherichia coli (strain K12)"='org.EcK12.eg.db',"Escherichia coli (strain Sakai)"='org.EcSakai.eg.db',"Anopheles gambiae"='org.Ag.eg.db',"Bos taurus"='org.Bt.eg.db',"Canis familiaris"='org.Cf.eg.db',"Macaca mulatta"='org.Mmu.eg.db',"Plasmodium falciparum"='org.Pf.plasmo.db',"Pan troglodytes"='org.Pt.eg.db',"Sus scrofa"='org.Ss.eg.db',"Xenopus tropicalis"='org.Xl.eg.db')
DBS <<- list('org.Hs.eg.db'=org.Hs.eg.db,'org.Mm.eg.db'=org.Mm.eg.db,'org.Rn.eg.db'=org.Rn.eg.db,"org.Gg.eg.db"=org.Gg.eg.db,"org.Dr.eg.db"=org.Dr.eg.db,"org.Dm.eg.db"=org.Dm.eg.db,"org.Ce.eg.db"=org.Ce.eg.db,"org.Sc.sgd.db"=org.Sc.sgd.db,"org.At.tair.db"=org.At.tair.db,"org.EcK12.eg.db"=org.EcK12.eg.db,"org.EcSakai.eg.db"=org.EcSakai.eg.db,"org.Ag.eg.db"=org.Ag.eg.db,"org.Bt.eg.db"=org.Bt.eg.db,"org.Cf.eg.db"=org.Cf.eg.db,"org.Mmu.eg.db"=org.Mmu.eg.db,"org.Pf.plasmo.db"=org.Pf.plasmo.db,"org.Pt.eg.db"=org.Pt.eg.db,"org.Ss.eg.db"=org.Ss.eg.db,"org.Xl.eg.db"=org.Xl.eg.db)
enrichRdbs <- as.character(read.csv(paste0(wd,"/www/enrichRdbs.csv"))[,1])

end.load <- Sys.time()
print("loading time")
print(end.load-start.load)
##### UI from here ###########
ui <- navbarPage(id = "navbar",
  theme = shinytheme("flatly"),
  title = 'ABioTrans',
  tabPanel('Home',
           # useShinyjs(),
           sidebarPanel(
             radioButtons('file_type',"Choose File Type",
                          c('Raw file (read count)'='raw','Normalised file'='norm')),
             conditionalPanel(
               condition = "input.file_type=='raw'",  # raw
               p("Example ",a("here", href="https://github.com/buithuytien/ABioTrans/blob/master/Test%20data/Eg_raw.png")),  # ADD EXAMPLE
               fileInput('file1','Choose Raw Counts'),
               # radioButtons('norm_method',"Normalisation method",
               #              c('RPKM','FPKM','TPM')),
               p("Example ",a("here", href = "https://github.com/buithuytien/ABioTrans/blob/master/Test%20data/Eg_gene_length.png")),  # ADD EXAMPLE
               fileInput('length1','Choose Gene Length'), #gene id + length
               p("Example ",a("here", href = "https://github.com/buithuytien/ABioTrans/blob/master/Test%20data/Eg_negative_control_genes.png")),  # ADD EXAMPLE
               fileInput('spikes1','Choose Negative Control Genes')
               # helpText("* Format requirement: CSV file. The first column contains gene names; the read counts of each genotype (conditions: wildtype, mutants, replicates, etc.) are in the following columns.Each genotype column should have a column name. ")
             ),
             conditionalPanel(
               condition = "input.file_type=='norm'", # normalized
               p("Example ",a("here", href = "https://github.com/buithuytien/ABioTrans/blob/master/Test%20data/Eg_normalised.png")),  # ADD EXAMPLE
               fileInput('file2','Choose Normalized Expression')
               # helpText("* Format requirement: CSV file. Gene names in rows and genotypes in columns, following the usual format of files deposited in the GEO database.")
             ),
             p("Example ",a("here", href="https://github.com/buithuytien/ABioTrans/blob/master/Test%20data/Eg_metadata.png")),  # ADD EXAMPLE
             fileInput('metafile1','Choose Meta Data File'),
             actionButton("submit_input","Submit")
           ),
           mainPanel(
             h3('Welcome to ABioTrans --'),
             h3('A Biostatistical tool for Transcriptomics Analysis'),
             img(src="Abiotrans-logo.png",
                 width = 570,height = 370)
           )
  ),
  tabPanel('Preprocessing',
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
  ),
  tabPanel('    Scatter    ',
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
           )),
  tabPanel('Distribution Fit',
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
           )),
  tabPanel('  Correlation  ',
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
  ),
  tabPanel('PCA',
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
           )),
  tabPanel("DE Analysis",
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
  ),
  tabPanel('Heatmap',
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
           )),
  
  ######## NOISE ######
  #############################################
  tabPanel('Noise',
           sidebarPanel(
             radioButtons('noise_situation',"Select desired noise plot between",choices = c('replicates'='a','genotypes (average of replicates)'='b','genotypes (no replicate)'='c')),
             conditionalPanel(
               condition = "input.noise_situation=='a' | input.noise_situation=='b' ",
               textInput('noise_numOfRep',"Number of replicates",value=1),
               helpText("* Please order the sample columns in input file properly. Replicates of the same genotype should be put in adjacent columns.")
             ),
             conditionalPanel(
               condition = "input.noise_situation=='b'",
               uiOutput("noise_anchor_choices")
             ),
             conditionalPanel(
               condition = "input.noise_situation=='c'",
               selectInput('noise_anchor_c',"Anchor genotype",choices = "")
             ),
             radioButtons('noise_graph_type',"Graph type:",
                          c('Bar chart','Line chart')),
             downloadButton("downloadnoise","Download as PNG"),
             actionButton('noise_plot',"Plot",width='65px',style="color: #fff; background-color: #337ab7; border-color:#337ab7;float:right"),
             conditionalPanel(
               condition = "input.noise_situation=='a' | input.noise_situation=='b' ",
               h5("Specify names of the genotypes"),
               uiOutput("expand_genonames_noise")
             )
           ),
           mainPanel(
             conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                              div(img(src="load.gif",width=240,height=180),h4("Processing ... Please wait"), style="text-align: center;")
             ), 
             conditionalPanel(condition="!$('html').hasClass('shiny-busy')",
                              plotlyOutput('noise.plot')
             )
           )),
  
  
  ###### ENTROPY #############
  #########################################
  tabPanel('Entropy',
           sidebarPanel(
             checkboxInput('tsflag',strong("Time series data"),FALSE),
             conditionalPanel(
               condition = "input.tsflag==true",
               textInput('entropy_timepoints',"Number of time points"),
               helpText("* Please order the sample columns in input file properly. Time series data of the same genotype should be put in adjacent columns.")
             ),
             radioButtons('entropy_graph_type',"Graph type:",
                          c('Bar chart','Line chart')),
             downloadButton("downloadentropy","Download as PNG"),
             conditionalPanel(
               condition = "input.tsflag==true",
               h5("Specify names of the genotypes"),
               uiOutput("expand_genonames_entropy")
             )
           ),
           mainPanel(
             h3('Shannon entropy'),
             plotlyOutput('entropy.plot')
           )
  ),
  
  
  ################## GO analysis ###################
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
                 selectInput('enrichR_dbs',"Select database",choices=enrichRdbs)
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
)
  ####################################################

server <- function(input,output,session){
  
  ########################################
  ##### get variable names for input #####
  ########################################
  
  observe({
    type <- input$file_type
    if(type=='norm'){
      DS <- df_norm()
    }else if(type=='raw'){
      DS <- df_raw()
    }
    nms <- colnames(DS)
    updateSelectInput(session, "scatter.x", choices = nms,selected = nms[1])
    updateSelectInput(session, "scatter.y", choices = nms,selected = nms[2])
    updateSelectInput(session, "dist.var", choices = nms)
    col_num <- ncol(DS)
    updateSliderInput(session,"pca_cluster_num",max=col_num-1)
    genotype_num <- NULL
    if(is.null(DS)==FALSE){
      for(i in 2:col_num){ 
        if(col_num %% i == 0)
          genotype_num <- c(genotype_num,i)
      }
    }
    updateSelectInput(session,"numOfGeno",choices=genotype_num)
    updateSelectInput(session,"noise_anchor_c",choices = nms)
    
    ### preprocessing tab
    f <- group_names()
    f <- unique(as.character(f))
    if(is.null(f)){
      hideTab(inputId="preprocessing_tabs", target="Description table")
      # hideTab(inputId="preprocessing_tabs", target="Description table")
    } else {
      showTab(inputId="preprocessing_tabs", target="Description table")
      updateSelectInput(session,"f1",choices=f,selected =f[1])
      updateSelectInput(session,"f2",choices=f,selected =f[2])
    }
    
    ### gene expression range for distribution fit ###
    if(is.null(DS)==FALSE){
      DS_dist <- distfit_df()
      range_min <- min(DS_dist)
      range_max <- max(DS_dist)
      updateSliderInput(session,"dist_range",max=round(range_max),value = c(0.1,range_max))
      updateNumericInput(session,"dist_range_min",min=0.000001,max=round(range_max),value = 0.1)
      updateNumericInput(session,"dist_range_max",min=0.000001,max=round(range_max),value = round(range_max))
    }
    
    ### gene sample size choices for PCA ###
    # print("line 647 check input$submit_preprocessing")
    # v=input$submit_preprocessing
    if(input$submit_preprocessing > 0){
      if(type=='norm'){
        DS_filt <- df_shiny()
      }else if(type=='raw'){
        DS_filt <- df_raw_shiny()
      }
    } else{
      DS_filt <- DS
    }
    
    i <- 1
    min_size <- 25
    samplesize <- NULL
    while(i*min_size<length(DS_filt[,1])){
      samplesize <- c(samplesize,i*min_size)
      i <- i*2
    }
    if(is.null(samplesize)){
      samplesize <- c(samplesize,length(DS_filt[,1]))
    }else if(samplesize[length(samplesize)]!=length(DS_filt[,1])){
      samplesize <- c(samplesize,length(DS_filt[,1]))
    }
    updateSelectInput(session,"gene_size", choices = samplesize,selected = samplesize[length(samplesize)])

    ### pca choices for PCA-2D ###
    pcchoices <- NULL
    if(is.null(DS)==FALSE)
      for (i in 1:ncol(DS)){
        pcchoices <- c(pcchoices,paste("PC",i,sep=""))
      }
    updateSelectInput(session,"pca.x",choices = pcchoices,selected = pcchoices[1])
    updateSelectInput(session,"pca.y",choices = pcchoices,selected = pcchoices[2])
    
    ### Gene ontology ####
    # database update
    go_method <- input$go_method
    dbs_name <- input$go_species
    dbs <- DBS[[dbs_name]]
    # slect_id <- input$go_geneidtype
    if( ! is.null (dbs)){
      id_choices <- AnnotationDbi::keytypes(dbs)
      slt <- input$go_geneidtype
      if (! slt %in% id_choices ){
        slt = id_choices[1]
      }
      updateSelectInput(session,"go_geneidtype", choices=id_choices, selected = slt)
    } else{
      id_choices <- NULL
      updateSelectInput(session,"go_geneidtype", choices=id_choices)
    }
    
    # go terms on GO graph
    go_method <- input$go_method
    res <- go_res_filt()
    if (go_method == "clusterProfiler" & !is.null(res) ){
      go_terms <- as.character(res$Description)
    } else {
      go_terms <- NULL
    }
    updateSelectInput(session,"go_term_slect",choices=go_terms)
  })
  
  
  observeEvent(input$submit_input, {
    type <- input$file_type
    if(type=='norm'){
      DS <- df_norm()
      lengths <- 0
    }else if(type=='raw'){
      DS <- df_raw()
      lengths <- gene_length()
      # if( length(intersect(rownames(lengths), rownames(DS))) < 1000 )
      #   length <- NULL
    }
    f <- group_names()
    spikes <- neg_control()
    
    # if any NULL value, throw error. TO CHANGE TO BE MORE SPECIFIC
    input.list <- list(DS, f)
    input.null <- sapply(input.list, is.null)
    names(input.null) <- c("Expression/Counts","Meta Data")
    
    if( any(input.null) ){
      index.null <- which(input.null)
      errors <- paste(names(input.null)[index.null],collapse = ', ')
      # print(errors)
      showModal(modalDialog(
        type = "Error",
        paste("Please check these input:",errors,"and try again!")
      ))
    } else{
      updateNavbarPage(session, inputId="navbar",selected="Preprocessing")
    }
    
    # update input
    updateNumericInput(session,"min_col",max=ncol(DS))   # update max column nunmber in filtering
    if(is.null(spikes)){
      updateRadioButtons(session,"norm_method",choices = c("None (Black)"="None",
                                                           'RPKM (Blue)'='RPKM','FPKM (Dark cyan)'='FPKM',
                                                           'TPM (Dark green)'='TPM',
                                                           "Upper Quartile (Brown)"='RUV') ) 
      #c("None",'RPKM','FPKM','TPM',"Upper Quartile"="RUV")
    } else {
      updateRadioButtons(session,"norm_method",choices = c("None (Black)"="None",
                                                           'RPKM (Blue)'='RPKM','FPKM (Dark cyan)'='FPKM',
                                                           'TPM (Dark green)'='TPM',
                                                           "RUV (Brown)"='RUV'))
    }
    if(is.null(lengths) & !(is.null(spikes)) ){
      updateRadioButtons(session,"norm_method",choices = c("None (Black)"="None","RUV (Brown)"="RUV"))
    } else if(is.null(lengths) & (is.null(spikes)) ){
      updateRadioButtons(session,"norm_method",choices = c("None (Black)"="None","Upper Quartile (Brown)"="RUV"))
    }
    
    if(is.null(f)){
      hideTab(inputId="navbar",target="DE Analysis")
    } else{
      showTab(inputId="navbar",target="DE Analysis")
    }
    # if(is.null(f)){
    #   hideTab(inputId="preprocessing_tabs", target="Description table")
    # } else {
    #   showTab(inputId="preprocessing_tabs", target="Description table")
    # }
  })
  
  # observeEvent(input$submit_preprocessing, {
  #   type <- input$file_type
  #   if(type=='norm'){
  #     DS <- df_shiny()
  #   }else if(type=='raw'){
  #     DS <- df_raw_shiny()
  #   }
  #   ### gene sample size choices for PCA ###
  #   i <- 1
  #   min_size <- 25
  #   samplesize <- NULL
  #   while(i*min_size<length(DS[,1])){
  #     samplesize <- c(samplesize,i*min_size)
  #     i <- i*2
  #   }
  #   if(is.null(samplesize)){
  #     samplesize <- c(samplesize,length(DS[,1]))
  #   }else if(samplesize[length(samplesize)]!=length(DS[,1])){
  #     samplesize <- c(samplesize,length(DS[,1]))
  #   }
  #   updateSelectInput(session,"gene_size", choices = samplesize,selected = samplesize[length(samplesize)])
  # })
  
  ######################################
  ######### read in / get data #########
  ######################################
  
  #####################
  ## get data #########
  #####################
  
  # get normalized counts
  df_norm <- reactive({        # get normalized counts
    if (is.null(input$file2))
      return (NULL)
    parts <- strsplit(input$file2$datapath,".",fixed=TRUE)
    type <- parts[[1]][length(parts[[1]])]
    if(type!="csv"){
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return (NULL)
    }
    ds <- read.csv(input$file2$datapath)
    ds <- na.omit(ds)
    ds <- ds[!duplicated(ds[,1]),]   # remove duplicated gene names
    
    row_names <- ds[,1]
    DS <- data.frame(ds)
    if(ncol(DS)<=1){
      showModal(modalDialog(
        title = "Error",
        "Please check normalised data file format (Eg_normalised.png) and try again!"
      ))
      return(NULL)
    }
    DS <- DS[,-1]
    row.names(DS) <- row_names
    for (i in 1:ncol(DS)){
      if(class(DS[,i])!="numeric" & class(DS[,i])!="integer"){
        showModal(modalDialog(
          title = "Error",
          "Please check normalised data file format (Eg_normalised.png) and try again!"
        ))
        return(NULL)
      }
    }
    return(DS)
  })  
  
  # get raw counts
  df_raw <- reactive ({          
    if(is.null(input$file1))
      return(NULL)
    parts <- strsplit(input$file1$datapath,".",fixed=TRUE)
    type <- parts[[1]][length(parts[[1]])]
    if(type!="csv"){
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return (NULL)
    }
    raw_ds <- read.csv(input$file1$datapath)
    raw_ds <- na.omit(raw_ds)
    raw_ds <- raw_ds[!duplicated(raw_ds[,1]),]   # remove duplicated gene names
    
    # raw_ds <- as.data.frame(raw_ds)
    if(ncol(raw_ds)<=1){
      showModal(modalDialog(
        title = "Error",
        "Data file must contain at least 2 columns. Please check raw data format and try again!"
      ))
      return(NULL)
    }
    
    row_names <- raw_ds[,1]
    rownames(raw_ds) <- row_names
    raw_DS <- raw_ds[,-1]  # remove the first column, which is gene Id
    
    for (i in 1:ncol(raw_DS)){
      if(class(raw_DS[,i])!="numeric" & class(raw_DS[,i])!="integer"){
        showModal(modalDialog(
          title = "Error",
          "Raw counts must be integer. Please check raw data formate and try again!"
        ))
        return(NULL)
      }
    }
    return(raw_DS)
  })
  
  # get gene length
  gene_length <- reactive({
    if (is.null(input$length1))
      return (NULL)
    lengths_df <- read.csv(input$length1$datapath)
    lengths_df2 <- data.frame("len" = lengths_df[,2]); 
    rownames(lengths_df2) <- as.character(lengths_df[,1])
    return(lengths_df2)
  })
  
  # get spikes / negative control genes
  neg_control <- reactive({
    if(is.null(input$spikes1))
      return(NULL)
    spikes <- read.csv(input$spikes1$datapath,header=F)
    spikes <- as.character(spikes[,1])
    # print(spikes[1:10])
    return(spikes)
  })
  
  # get meta data table
  group_names <- reactive({
    # if no data
    if(is.null(input$metafile1))
      return(NULL)
    
    # read in group names (metadata)
    groups <- read.csv(input$metafile1$datapath)
    group_colnames <- as.character(groups[,1])
    
    type <- input$file_type
    if(type=='norm'){
      DS <- df_norm()
    }else if(type=='raw'){
      DS <- df_raw()
    }
    col_names <- colnames(DS)   # columm names of DS in order
    
    # check if groups and column names are similar
    if ( !all(col_names %in% group_colnames) || ncol(groups) < 2 ){
      showNotification(type = "error", "group names and DS column names not similar")
      return(NULL)
    }
    
    if(ncol(groups)==2){
      f <- groups[match(col_names,groups[,1]),] [,2]   # arrange f in the same order as col_names
    } else {
      f <- groups[match(col_names,groups[,1]),] [,2]
      for(i in 3:ncol(groups)){
        f <- paste0(f,"_",groups[,i])
      }
    }
    f <- as.factor(make.names(f))
    # return(as.factor(f))
    return(f)
  })
  
  ### Gene ontology
  
  gene_list <- reactive({
    if(is.null(input$filego))
      return (NULL)
    parts <- strsplit(input$filego$datapath,".",fixed=TRUE)
    type <- parts[[1]][length(parts[[1]])]
    if(type!="csv"){
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return (NULL)
    }
    ds <- read.csv(input$filego$datapath,header=FALSE)
    if(ncol(ds) >= 2){
      col1 <- ds[-1,1]
    } else if( ncol(ds) == 1){
      col1 <- ds[,1]
    } else {
      showModal(modalDialog(
        title = "Error",
        "No data found! Please check required data format and try again!"
      ))
      return (NULL)
    }
    gene_list <- as.character(col1)
    print("gene list from gene_list")
    print(head(gene_list))
    return(gene_list)
  })
  
  bg_list <- reactive({
    if(is.null(input$filebg))
      return (NULL)
    parts <- strsplit(input$filebg$datapath,".",fixed=TRUE)
    type <- parts[[1]][length(parts[[1]])]
    if(type!="csv"){
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return (NULL)
    }
    ds <- read.csv(input$filebg$datapath,header=FALSE)
    if(ncol(ds) > 1){
      col1 <- ds[-1,1]
    } else if( ncol(ds) == 1){
      col1 <- ds[,1]
    } else {
      showModal(modalDialog(
        title = "Error",
        "No data found! Please check required data format and try again!"
      ))
      return (NULL)
    }
    bg_list <- as.character(col1)
    return(bg_list)
  })
  
  ####################################
  ########## PREPROCESSING ###########
  ####################################
  
  # filter normalized counts
  df_shiny <- eventReactive(input$submit_preprocessing, {
    DS_norm <- df_norm()
    min_val <- input$min_val
    min_col <- input$min_col
    keep <- rowSums(DS_norm >= min_val) >= min_col
    DS <- DS_norm[keep,]
    # DS <- apply(DS_norm, 1, function(x) length(x[x>min_val])>=min_col) 
    return(DS)
  })
  
  # filter raw counts
  df_raw_filt <- eventReactive(input$submit_preprocessing, {
    DS_raw <- df_raw()
    min_val <- input$min_val
    min_col <- input$min_col
    keep <- rowSums(DS_raw >= min_val) >= min_col
    DS_filt <- DS_raw[keep,]
    # DS_filt <- apply(DS_raw, 1, function(x) length(x[x>min_val])>=min_col) 
    return(DS_filt)    
  })
  
  # normalizing raw counts
  df_raw_shiny <- reactive({
    raw_DS <- df_raw_filt()     # get filtered raw counts
    method <- input$norm_method
    
    if(method%in% c("TPM","RPKM","FPKM")){
      lengths_df <- gene_length()
      merge_DS <- merge(raw_DS,lengths_df,by="row.names")
      rownames(merge_DS) <- merge_DS[,1]; merge_DS <- merge_DS[,-1]; 
      raw_DS <- merge_DS[,-ncol(merge_DS)]
      lengths <- merge_DS[,ncol(merge_DS)]
      # print("length")
      # print(head(merge_DS))
    }
    # print("from line 981 df_raw_shiny")
    # print(method)
    # print("raw_DS")
    # print(head(raw_DS[,1:4]))
    # print("dimension of raw_DS")
    # print(dim(raw_DS))
    
    if(method=='TPM'){
      tpm.matrix<- apply(raw_DS, 2, function(x) tpm(x, lengths))
      tpm.df <- data.frame(tpm.matrix)
      return (tpm.df)
    }else if(method=='RPKM'){
      rpkm.matrix <- edgeR::rpkm(raw_DS,lengths)
      rpkm.df <- data.frame(rpkm.matrix)
      return (rpkm.df)
    }else if(method=='FPKM'){
      fpkm.matrix<- apply(raw_DS, 2, function(x) fpkm(x, lengths))
      fpkm.df <- data.frame(fpkm.matrix)
      return (fpkm.df)
    }else if(method=='None'){
      return (raw_DS)
    }else if(method=='RUV'){
      spikes <- neg_control()
      if (!is.null(spikes))
        spikes <- intersect(spikes,rownames(raw_DS))
      # f <- group_names()
      # if( is.null(spikes) )
      #   spikes <- getEmpirical(rawDS,f)
      set1 <- RUVg.apply(raw_DS,spikes)
      RUV.df <- as.data.frame(normCounts(set1))
      return (RUV.df)
    }
    
  })
  
  ### for distribution fitting
  distfit_df <- reactive({
    type <- input$file_type
    if(type=='norm'){
      DS <- df_shiny()
    }else if(type=='raw'){
      DS <- df_raw_shiny()
    }
    for (i in 1:ncol(DS)) {
      DS <- DS[which(DS[,i] > 0),]
      DS <- na.omit(DS) 
    }
    return(DS)
  })
  
  
  ######### ANALYSIS FROM HERE ############
  ######## RLEplot and Preprocessing ###########
  #############################################
  RLE.plot <- reactive({
    type <- input$file_type
    if(type=='norm'){
      DS <- df_shiny()
    }else if(type=='raw'){
      DS <- df_raw_shiny()
    }
    set1 <- newSeqExpressionSet(as.matrix(DS))
    norm_method_name <- input$norm_method
    colors <- c('RPKM'='blue','FPKM'='darkcyan','TPM'='darkgreen',"RUV"='Brown',"Upper Quartile"='Brown')
    if(norm_method_name!="None" &input$submit_preprocessing != 0){
      spikes <- neg_control()
      if(norm_method_name == "RUV" & is.null(spikes))
        norm_method_name <- "Upper Quartile"
      plotRLE(set1, ylim=c(-1.5,1.5),outline=FALSE, col=colors[norm_method_name],
              main= paste(norm_method_name,"Normalized"))
    }
  })
  
  output$RLE.plot <- renderPlot({
    RLE.plot()
  })
  
  output$RLE.plot2 <- renderPlot({   # for raw data
    start.rle <- Sys.time()
    type <- input$file_type
    if(type=='norm'){
      raw_DS <- df_shiny()
      main_title <- "Input data"
    }else if(type=='raw'){
      raw_DS <- df_raw()
      main_title <- "Raw data"
    }
    set1 <- newSeqExpressionSet(as.matrix(raw_DS))
    if(input$submit_preprocessing != 0)
      plotRLE(set1, ylim=c(-1.5,1.5),outline=FALSE, main=main_title)
    end.rle <- Sys.time()
    print("time for RLE plot and preprocessing")
    print(end.rle - start.rle) 
  })
  
  
  output$norm_table <- DT::renderDataTable({
    type <- input$file_type
    if(type=='norm'){
      DS <- df_shiny()
    }else if(type=='raw'){
      DS <- df_raw_shiny()
    }
    # if(input$submit_preprocessing != 0)
    DS   # with filtering and normalization
  })
  
  output$meta_table <- DT::renderDataTable({
    f <- group_names()
    type <- input$file_type
    if(type=='norm'){
      DS <- df_shiny()
    }else if(type=='raw'){
      DS <- df_raw_shiny()
    }
    if(! is.null(f)){
      meta_df <- data.frame("Column names"=colnames(DS),"Description"=f)
      meta_df
    }
  })
  
  output$download_norm_data <- downloadHandler(
    filename = function(){
      method <- input$norm_method
      paste(method,"normalized.csv")
    },
    content = function(file){
      type <- input$file_type
      if(type=='norm'){
        DS <- df_shiny()
      }else if(type=='raw'){
        DS <- df_raw_shiny()
      }
      write.csv(DS, file, row.names = F)
    }
  )
  
  ############################
  ######## scatter ###########
  ############################
  
  plotScatter <- reactive({
    scatter.start <- Sys.time()
    trans <- input$trans
    x <- input$scatter.x
    y <- input$scatter.y
    type <- input$file_type
    if(type=='norm'){
      DS <- df_shiny()
    }else if(type=='raw'){
      DS <- df_raw_shiny()
    }
    if(trans=='None'){
      scatter.data <- DS
    }else if(trans=='Natural log'){
      scatter.data <- log1p(DS)
    }else if(trans=='log2'){
      scatter.data <- log2(DS+1)
    }else if(trans=='log10'){
      scatter.data <- log10(DS+1)
    }
    scatter.end <- Sys.time()
    print("Scatter plot time")
    print(scatter.end - scatter.start)
    return (list(x,y,scatter.data))
  })
  
  scatterplot <- function(){
    li <- plotScatter()
    x <- li[[1]]
    y <- li[[2]]
    scatter.data <- li[[3]]
    d <- kde2d(scatter.data[,x],scatter.data[,y])
    ColorLevels <- round(seq(min(d$z), max(d$z), length=5),4)
    heatscatter(x=scatter.data[,x],y=scatter.data[,y],xlab = x, ylab=y, main="")
    legend("topleft", paste("R=",round(cor(scatter.data[,x],scatter.data[,y]),3)), bty="n")
    legend("bottomright",title="KDE",legend=ColorLevels, pch=19,col=LSD::colorpalette("heat"))
    if(x!=y){
      lines(lowess(scatter.data[,x],scatter.data[,y]),col="black")
    }
  }
  scatterplot_collage <- function(){
    li <- plotScatter()
    scatter.data <- li[[3]]
    par(mfrow=c(3,3))
    for(i in 1:ncol(scatter.data)){
      for(j in i:ncol(scatter.data)){
        d <- kde2d(scatter.data[,i],scatter.data[,j])
        ColorLevels <- round(seq(min(d$z), max(d$z), length=5),4)
        heatscatter(x=scatter.data[,i],y=scatter.data[,j],xlab = colnames(scatter.data)[i], ylab=colnames(scatter.data)[j], main="")
        legend("topleft", paste("R=",round(cor(scatter.data[,i],scatter.data[,j]),3)), bty="n")
        legend("bottomright",title="KDE",legend=ColorLevels, pch=19,col=LSD::colorpalette("heat"))
        if(i!=j){
          lines(lowess(scatter.data[,i],scatter.data[,j]),col="black")
        }
      }
    }
  }
  
  output$scatter.plot <- renderPlot({
    scatterplot()
  })
  
  output$downloadscatter_collage <- downloadHandler(
    filename = function(){
      paste("heatscatter_collage",".pdf",sep="")
    },
    content = function(file){
      pdf(file)
      scatterplot_collage()
      dev.off()
    }
  )
  
  output$downloadscatter <- downloadHandler(
    filename = function(){
      paste("heatscatter",".pdf",sep="")
    },
    content = function(file){
      pdf(file) 
      scatterplot()
      dev.off()
    }
  )
  
  ############################
  ######## distfit ###########
  ############################
  
  output$downloaddist <- downloadHandler(
    filename = function(){
      paste("distribution_fit",".pdf",sep="")
    },
    content = function(file){
      pdf(file) 
      distplot()
      dev.off()
    }
  )
  
  output$dist_range_allowed <- renderText({
    DS <- distfit_df()
    paste("Suggested range: ( 0"," ~ ",round(max(DS))," ]",sep="")
  })
  
  plotDist <- reactive({
    dist.start <- Sys.time()
    dis <- input$distributions
    var <- input$dist.var
    DS <- distfit_df()
    fits <- list()
    distrs <- NULL
    numcol <- c(0,0,0,0,0,0)
    dist_zoom <- input$dist_zoom
    if(dist_zoom=='slider'){
      fit_range <- input$dist_range
    }else if(dist_zoom=='text input'){
      fit_range <- c(input$dist_range_min,input$dist_range_max)
    }
    if("Log-normal" %in% dis){
      fit_ln <- fitdist(DS[,var], "lnorm")
      fits <- c(fits,list(fit_ln))
      distrs <- c(distrs,"Log-normal")
      numcol[1]=1
    }
    if("Log-logistic" %in% dis){ 
      fit_ll <- fitdist(DS[,var], "llogis", start = list(shape = 10, scale = 10),lower=c(0,0))
      fits <- c(fits,list(fit_ll))
      distrs <- c(distrs,"Log-logistic")
      numcol[2]=1
    }
    if("Pareto" %in% dis){
      fit_P <- fitdist(DS[,var], "pareto", start = list(shape = 10, scale = 10),lower=c(0,0))
      fits <- c(fits,list(fit_P))
      distrs <- c(distrs,"Pareto")
      numcol[3]=1
    }
    if("Burr" %in% dis){
      fit_B <- fitdist(DS[,var], "burr", start = list(shape1 = 0.3, shape2 = 1, rate = 1),lower=c(0,0,0))
      fits <- c(fits,list(fit_B))
      distrs <- c(distrs,"Burr")
      numcol[4]=1
    }
    if("Weibull" %in% dis){
      fit_W <- fitdist(DS[,var], "weibull",lower=c(0,0))
      fits <- c(fits,list(fit_W))
      distrs <- c(distrs,"Weibull")
      numcol[5]=1
    }
    if("Gamma" %in% dis){
      fit_G <- fitdist(DS[,var], "gamma",lower=c(0, 0),start=list(scale=1,shape=1))
      fits <- c(fits,list(fit_G))
      distrs <- c(distrs,"Gamma")
      numcol[6]=1
    }
    dist.end <- Sys.time()
    print("Distribution fitting time")
    print(dist.end - dist.start)
    return (list(fits,distrs,numcol,var,fit_range))
    
  })
  
  output$dist.plot <- renderPlot({
    distplot()
  })
  
  distaic <- reactive({
    dist.start <- Sys.time()
    DS <- distfit_df()
    AIC.df <- as.data.frame(matrix(nrow=ncol(DS),ncol=6))
    rownames(AIC.df) <- colnames(DS)
    colnames(AIC.df) <- c("Log-normal","Log-logistic","Pareto", "Burr", "Weibull", "Gamma")
    for(i in 1:nrow(AIC.df)){
      fit_ln <- fitdist(DS[,i], "lnorm")
      fit_ll <- fitdist(DS[,i], "llogis", start = list(shape = 10, scale = 10),lower=c(0,0))
      fit_P <- fitdist(DS[,i], "pareto", start = list(shape = 10, scale = 10),lower=c(0,0))
      fit_B <- fitdist(DS[,i], "burr", start = list(shape1 = 0.3, shape2 = 1, rate = 1),lower=c(0,0,0))
      fit_W <- fitdist(DS[,i], "weibull",lower=c(0,0))
      fit_G <- fitdist(DS[,i], "gamma",lower=c(0, 0),start=list(scale=1,shape=1))
      fits <- list(fit_ln,fit_ll,fit_P,fit_B,fit_W,fit_G)
      AIC.df[i,] <- gofstat(fits)$aic
    }
    for(i in 1:nrow(AIC.df)){
      AIC.df$min.AIC[i]<-colnames(AIC.df)[which.min(AIC.df[i,1:6])]
    }
    dist.end <- Sys.time()
    print('distribution fitting time')
    print(dist.end - dist.start)
    return (AIC.df)
  })
  
  output$dist.aic <- renderTable({
    distaic()
  },rownames=TRUE)
  
  distplot <- function(){
    li <- plotDist()
    fits <- li[[1]]
    distrs <- li[[2]]
    numcol <- li[[3]]
    var <- li[[4]]
    fit_range <- li[[5]]
    line_types <- c(1,2,3,4,5,6) #par lty
    if(length(fits)!=0)
      cdfcomp(fits, xlogscale = TRUE, ylogscale = TRUE,
              ylab = "CDF", xlab = "Expression levels (log)", xlim = c(fit_range[1], fit_range[2]),
              legendtext = distrs, cex = 0.5 ,lwd=2, main = var,fitcol=rainbow(6)[which(numcol==1)],fitlty = line_types[which(numcol==1)])
  } 
  
  output$downloaddistaic <- downloadHandler(
    filename = function(){
      paste("aic",".csv",sep="")
    },
    content = function(file){
      write.csv(distaic(),file,row.names = TRUE)
    }
  )
  
  ############################
  ####### correlation ########
  ############################
  
  COR <- function(d, i,myMethod){
    Result2 <- cor(x = d[,i], y = d[,i], method = myMethod)
    return(format(round(Result2, 5), nsmall = 5))
  }
  
  cor_df <- reactive({
    cor.start <- Sys.time()
    type <- input$file_type
    if(type=='norm'){
      DS <- df_shiny()
    }else if(type=='raw'){
      DS <- df_raw_shiny()
    }
    method <- input$cor_method
    if(method=="Pearson correlation"){
      Cor2 <- data.frame(COR((DS),1:length(DS),"pearson"))
    }else if(method=="Spearman correlation"){
      Cor2 <- data.frame(COR((DS),1:length(DS),"spearman"))
    }
    Cor2 <- na.omit(Cor2)
    cor.end <- Sys.time()
    print("correlation time")
    print(cor.end - cor.start)
    return (Cor2)
  })
  
  output$corr.plot <- renderPlot({
    corrplot1()
  })
  
  output$corr.plot2 <- renderPlot({
    corrplot2()
  })
  
  output$corr.matrix <- renderTable({
    cor_df()
  },rownames=TRUE)
  
  corrplot1 <- function(){
    corr <- as.matrix(cor_df())
    corr <- apply(corr,2,as.numeric)
    rownames(corr) <- rownames(cor_df())
    if(ncol(corr)<=20){
      fontsize <- 1
    }else{
      fontsize <- 20/ncol(corr)
    }
    corrplot(corr,method="shade",shade.col=NA,tl.col="black",cl.lim=c(min(corr),1),is.corr = FALSE,tl.cex = fontsize)
  }
  
  corrplot2 <- function(){
    corr <- as.matrix(cor_df())
    corr <- apply(corr,2,as.numeric)
    rownames(corr) <- rownames(cor_df())
    if(ncol(corr)<=20){
      fontsize <- 1
    }else{
      fontsize <- 20/ncol(corr)
    }
    corrplot(corr,type="upper",tl.col="black",cl.lim=c(min(corr),1),is.corr = FALSE,tl.cex = fontsize)
  }
  
  output$downloadcorrplot <- downloadHandler(
    filename = function(){
      paste("corrheatmap",".pdf",sep="")
    },
    content = function(file){
      pdf(file) 
      corrplot1()
      dev.off()
    }
  )
  
  output$downloadcorrplot2 <- downloadHandler(
    filename = function(){
      paste("corrplot",".pdf",sep="")
    },
    content = function(file){
      pdf(file)
      corrplot2()
      dev.off()
    }
  )
  
  output$downloadcorrmat <- downloadHandler(
    filename = function(){
      paste("correlation",".csv",sep="")
    },
    content = function(file){
      write.csv(cor_df(),file,row.names = TRUE)
    }
  )
  
  ############################
  #######     PCA     ########
  ############################
  
  refreshDS1 <- eventReactive(input$pca_refresh,{
    type <- input$file_type
    if(type=='norm'){
      DS <- df_shiny()
    }else if(type=='raw'){
      DS <- df_raw_shiny()
    }
    DS1 <- DS[sample(nrow(DS), nrow(DS), replace = FALSE), ]
    return (DS1)
  })
  
  plotPCA <- reactive({ #process and return data
    pca.start <- Sys.time()
    type <- input$file_type
    if(type=='norm'){
      DS <- df_shiny()
    }else if(type=='raw'){
      DS <- df_raw_shiny()
    }
    order <- input$gene_order
    size <- input$gene_size
    x <- input$pca.x
    y <- input$pca.y
    cluster_flag <- input$pca_cluster
    rindex <- as.numeric(substring(x,3))
    cindex <- as.numeric(substring(y,3))
    if(order=='Ascending'){
      DS1 <- DS[order(DS[,1]),]
    }else if(order=='Descending'){
      DS1 <- DS[rev(order(DS[,1])),]
    }else if(order=='Random'){
      DS1 <- refreshDS1()
    }
    
    DSample <- head(DS1, n = size)
    PR <- prcomp(t(DSample),center=TRUE)
    PCA.var <- PR$sdev^2
    PCA.var.per <- round(PCA.var/sum(PCA.var)*100,1)
    xlabel <- paste(colnames(PR$x)[rindex]," - ", PCA.var.per[rindex], "%", sep="")
    ylabel <- paste(colnames(PR$x)[cindex]," - ", PCA.var.per[cindex], "%", sep="")
    if(cluster_flag==TRUE){
      num <- as.numeric(input$pca_cluster_num)
      kmeans.data <- data.frame(x=PR$x[,x],y=PR$x[,y])
      kmeans.result <- kmeans(kmeans.data,num)
      return (list(PR,PCA.var,PCA.var.per,rindex,cindex,xlabel,ylabel,cluster_flag,kmeans.result))
    }
    pca.end <- Sys.time()
    print("pca time")
    print(pca.end - pca.start)
    return (list(PR,PCA.var,PCA.var.per,rindex,cindex,xlabel,ylabel,cluster_flag))
  })
  
  pcavarplot <- function(){
    li <- plotPCA()
    PCA.var.per <- li[[3]]/100
    type <- input$file_type
    if(type=='norm'){
      DS <- df_shiny()
    }else if(type=='raw'){
      DS <- df_raw_shiny()
    }
    pcchoices <- NULL
    for (i in 1:length(PCA.var.per)){
      pcchoices <- c(pcchoices,paste("PC",i,sep=""))
    }
    xform <- list(categoryorder = "array",
                  categoryarray = pcchoices)
    p <- plot_ly(
      x = pcchoices,
      y = PCA.var.per,
      name = "PCA variance",
      type = "bar"
    ) %>% layout(xaxis = xform)
    
    return (p)
  }
  
  pca2dplot <- function(){
    li <- plotPCA()
    PR <- li[[1]]
    rindex <- li[[4]]
    cindex <- li[[5]]
    xlabel <- li[[6]]
    ylabel <- li[[7]]
    cluster_flag <- li[[8]]
    if(cluster_flag==FALSE){
      p <- plot_ly(
        x=PR$x[,rindex],
        y=PR$x[,cindex],
        type = "scatter",
        mode="markers"
      ) %>% layout(xaxis = list(title = xlabel), yaxis = list(title = ylabel))
    }else if(cluster_flag==TRUE){
      kmeans.result <- li[[9]]
      text_flag <- input$pca_text
      if(text_flag==TRUE){
        p <- plot_ly(
          x=PR$x[,rindex],
          y=PR$x[,cindex],
          type = "scatter",
          color=as.character(kmeans.result$cluster),
          mode="markers",
          colors = "Set1"
        ) %>% hide_colorbar() %>% 
          add_trace(
            x=PR$x[,rindex],
            y=PR$x[,cindex],
            type = 'scatter',
            mode = 'text', 
            text = names(kmeans.result$cluster), 
            textposition = 'top right'
          ) %>% layout(xaxis = list(title = xlabel), yaxis = list(title = ylabel),showlegend=FALSE)
      }else if(text_flag==FALSE){
        p <- plot_ly(
          x=PR$x[,rindex],
          y=PR$x[,cindex],
          type = "scatter",
          color=as.character(kmeans.result$cluster),
          mode="markers",
          text = names(kmeans.result$cluster),
          colors = "Set1"
        ) %>% hide_colorbar() %>% layout(xaxis = list(title = xlabel), yaxis = list(title = ylabel),showlegend=FALSE)
      }
    }
    
    
  }
  
  pca3dplot <- function(){
    li <- plotPCA()
    PR <- li[[1]]
    xlabel <- "PC1"
    ylabel <- "PC2"
    zlabel <- "PC3"
    cluster_flag <- li[[8]]
    if(cluster_flag==FALSE){
      p <- plot_ly(
        x=PR$x[,1],
        y=PR$x[,2],
        z=PR$x[,3],
        type="scatter3d",
        mode="markers",
        marker=list(size=5)
      ) %>% layout(scene=list(xaxis = list(title = xlabel), yaxis = list(title = ylabel),zaxis=list(title=zlabel)))
    }else if(cluster_flag==TRUE){
      kmeans.result <- li[[9]]
      text_flag <- input$pca_text
      if(text_flag==TRUE){
        p <- plot_ly(
          x=PR$x[,1],
          y=PR$x[,2],
          z=PR$x[,3],
          type = 'scatter3d',
          mode = 'text', 
          text = names(kmeans.result$cluster), 
          color = as.character(kmeans.result$cluster),
          textfont = list(size=10),
          textposition = 'top right'
        ) %>% layout(scene=list(xaxis = list(title = xlabel), yaxis = list(title = ylabel),zaxis=list(title=zlabel)),showlegend=FALSE)
      }else if(text_flag==FALSE){
        p <- plot_ly(
          x=PR$x[,1],
          y=PR$x[,2],
          z=PR$x[,3],
          type = "scatter3d",
          color=as.character(kmeans.result$cluster),
          mode="markers",
          marker=list(size=5),
          text = names(kmeans.result$cluster),
          colors = "Set1"
        ) %>% hide_colorbar() %>% layout(scene=list(xaxis = list(title = xlabel), yaxis = list(title = ylabel),zaxis=list(title=zlabel)),showlegend=FALSE)
      }
    }
  }
  
  output$pcavar.plot <- renderPlotly({
    pcavarplot()
  })
  
  output$pca2d.plot <- renderPlotly({
    pca2dplot()
  })
  
  output$pca3d.plot <- renderPlotly({
    pca3dplot()
  })
  
  output$downloadpcavar <- downloadHandler(
    filename = function(){
      paste("pca_variance",".png",sep="")
    },
    content = function(file){
      p <- pcavarplot()
      orca(p, file = "pca_variance.png")
    }
  )
  
  output$downloadpca2d <- downloadHandler(
    filename = function(){
      paste("pca2d",".png",sep="")
    },
    content = function(file){
      p <- pca2dplot()
      orca(p, file = "pca2d.png")
    }
  )
  
  output$downloadpca3d <- downloadHandler(
    filename = function(){
      paste("pca3d",".png",sep="")
    },
    content = function(file){
      p <- pca3dplot()
      plotly_IMAGE(p,format = "png",out_file = "pca3d.png") 
    }
  )
  
  ############################
  ######## DE analysis #######
  ############################
  group_names_de <- reactive({
    f <- group_names()
    f1 <- input$f1
    f2 <- input$f2
    if(is.null(f) || is.null(f1) || is.null(f2) )
      return(NULL)
    
    type <- input$file_type
    if(type=='norm'){
      raw_DS <- df_shiny()
    }else if(type=='raw'){
      raw_DS <- df_raw_filt()
    }
    f.df <- data.frame("f"=f); rownames(f.df) <- colnames(raw_DS)
    # f.df.slect <- subset(f.df,f %in% c(f1,f2) )
    # f.df.slect <- rbind( subset(f.df,f %in% f1),  subset(f.df,f %in% f2) )
    # f.df.slect2 <- f.df.slect; 
    # f_new <- f.df.slect[,1]
    # f.df.slect2[,1] <- droplevels(f_new, except = levels(f_new)%in%f_new)
    return(f.df)
  })
  
  df_raw_de <- reactive({
    f_de <- group_names_de()
    if(is.null(f_de))
      return(NULL)
    type <- input$file_type
    rep_number <- input$n_rep
    if(rep_number == 1) {
      de_type <- input$de_method1
    } else {
      de_type <- input$de_method0
    }
    
    if(type=='norm'){
      raw_DS <- df_shiny()         # filtered and normalized
    }else if(type=='raw'){
      if( de_type == "NOISeq") {
        raw_DS <- df_raw_shiny()   # filtered and normalized
      } else {
        raw_DS <- df_raw_filt()    # filtered and UN-NORMALIZED
      }
    }
    # raw_DS_de <- raw_DS[,rownames(f_de)]
    return(raw_DS)
  })
  
  # all helper functions are in utils.R file
  de_no_filt<- eventReactive(input$submit_DE, {       # return as table object
    start.de.table <- Sys.time()
    f_de <- group_names_de()  # for edgeR >> f=f_de[,1]; for the rest factors=f_de
    if(is.null(f_de) )
      return (NULL)
    DS_de <- df_raw_de()      # with only 2 conditions for DE analysis
    p_val <- 1 # input$p_val
    fc <- 1 # input$fc
    f1 <- input$f1
    f2 <- input$f2
    rep_number <- input$n_rep # either 0 = no replicates or 1 = have replicates
    
    
    spikes <- neg_control()
    norm_method <- input$norm_method
    if(! is.null(spikes) & norm_method=="RUV" ){
      set1 <- RUVg.apply(DS_de,spikes)
      W_1 <- pData(set1)$W_1
    } else{
      W_1 <- NULL
    }
    # print("from de_no_filt")
    # print(pData(set1))
    # print(W_1)
    
    if(rep_number == 1){  # have replicates
      de_type <- input$de_method1
      if(de_type == "EdgeR") {
        res <- edgerApply(DS=DS_de,f = f_de[,1],W_1=W_1,f1=f1,f2=f2)        # edgeR, return edgeR object
        res.df <- edgerFilter(res, FC=fc, p_val=p_val)  # fitler edgeR object result
      } else if(de_type == "DESeq2") {
        res <- deseqApply(DS=DS_de,f.df = f_de,W_1=W_1,f1=f1,f2=f2)            # DESeq, return DESeq object
        res.df <- deseqFilter(res, FC=fc, p_val=p_val)  # fitler DESeq object result
      } else if(de_type == "NOISeq") {
        res <- noiseqbioApply(DS=DS_de,f.df = f_de,f1=f1,f2=f2)           # NOISeqbio, return NOIseq object
        res.df <- noiseqbioFilter(res, FC=fc, p_val=p_val) # filter return NOIseq object
      }
    } else { # no replicates
      de_type <- input$de_method0 # NOISeq
      res <- noiseqsimApply(DS=DS_de,f.df = f_de,f1=f1,f2=f2)           # NOISeqbio, return NOIseq object
      res.df <- noiseqsimFilter(res, FC=fc)
    }
    end.de.table <- Sys.time()
    print("de table time")
    print(end.de.table - start.de.table)
    return(res.df)
  })
  
  
  de_filt <- function(res.df,p_val,fc,rep_number ) {
    # res.df <- de_no_filt()
    if(is.null(res.df))
      return(NULL)
    # p_val <- input$p_val
    # fc <- input$fc
    # rep_number <- input$n_rep
    
    if(rep_number == 1){
      res.df.filt <- filter(res.df, FDR<=p_val, log2FCabs>=log2(fc) )
    } else{
      res.df.filt <- filter(res.df, FDR<=p_val, log2FCabs>=log2(fc) )
    }
    return(res.df.filt)
  }
  
  output$DE_table <- DT::renderDataTable({
    res.df <- de_no_filt()
    p_val <- input$p_val
    fc <- input$fc
    rep_number <- input$n_rep
    if(input$submit_DE > 0){
      res.df.filt <- de_filt(res.df,p_val,fc,rep_number)
      res.df.filt
    }
  })
  
  ##### volcano plot ######
  volcano_plot <- eventReactive(input$submit_DE, {
    volcano.start.time <- Sys.time()
    rep_number <- input$n_rep
    if( rep_number == 0)
      return (NULL)
    de_type <- input$de_method1
    if(de_type == "NOISeq")
      return (NULL)
    
    res <- de_no_filt()   # de result, no filter
    if(is.null (res))
      return (NULL)
    
    p_val <- input$p_val
    fc <- input$fc
    res$Gene <- rownames(res)
    res <- na.omit(res)
    # plot
    ymax <- max(-log10(res$PValue)); if(ymax > 5) ymax = 5
    # print("from volcano plot - range(res$PValue)")
    # print(range(res$PValue))
    with(res, plot(log2FC, -log10(PValue), pch=20, main="Volcano plot",xlim=c(-2.5,2.5),ylim=c(0,ymax))) #xlim=c(-5,5),ylim=c(0,ymax)
    # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
    with(subset(res, FDR< p_val), points(log2FC, -log10(PValue), pch=20, col="red"))
    with(subset(res, abs(log2FC)>log2(fc)), points(log2FC, -log10(PValue), pch=20, col="orange"))
    with(subset(res, FDR< p_val & abs(log2FC)>log2(fc)), points(log2FC, -log10(PValue), pch=20, col="green"))
    legend("topleft",bty="n",col=c("red","orange","green","black"),pch=19,
           legend=c("FDR < FDR limit","FC > FC limit","Both","Other"))
    # library(calibrate)
    # with(subset(res, FDR<.05 & abs(log2FC)>1), textxy(log2FC, -log10(PValue), labs=Gene, cex=.8))
    volcano.end.time <- Sys.time()
    print("volcano time")
    print(volcano.end.time - volcano.start.time)
  })
  
  output$volcano_plot <- renderPlot({
    volcano_plot()
  })
  
  ##### dispersion plot ######
  dispersion_plot <- eventReactive(input$submit_DE, {
    dispersion.start.time <- Sys.time()
    f_de <- group_names_de()  # for edgeR >> f=f_de[,1]; for the rest factors=f_de
    if(is.null(f_de) )
      return (NULL)
    rep_number <- input$n_rep # either 0 or 1
    if(rep_number == 0)
      return(NULL)
    DS_de <- df_raw_de()      # with only 2 conditions for DE analysis
    p_val <- 1 # input$p_val
    fc <- 1 # input$fc
    de_type <- input$de_method1
    if(de_type == "EdgeR"){
      edgerDisp(DS_de, f_de[,1])
    } else if(de_type == "DESeq2"){
      deseqDisp(DS_de,f_de)
    }
    dispersion.end.time <- Sys.time()
    print("dispersion time")
    print(dispersion.end.time - dispersion.start.time)
  })
  
  output$dispersion_plot <- renderPlot({
    dispersion_plot()
  })
  
  ########## download buttons DE analysis ###########
  output$download_de_table <- downloadHandler(
    filename = function(){
      paste0("DE analysis",".csv")
    },
    content = function(file){
      res.df <- de_no_filt()
      p_val <- input$p_val
      fc <- input$fc
      rep_number <- input$n_rep
      res.df.filt <- de_filt(res.df,p_val,fc,rep_number)
      write.csv(res.df.filt, file, row.names = F)
    }
  )
  
  output$download_volcano <- downloadHandler(
    filename =function() {
      paste0("Volcano",".pdf")
    },
    content = function(file) {
      pdf(file)
      volcano_plot()
      dev.off()
    }
  )
  
  output$download_dispersion <- downloadHandler(
    filename=function(){
      paste0("Dispersion plot",".pdf")
    },
    content = function(file) {
      pdf(file)
      dispersion_plot()
      dev.off()
    }
  )
  
  # output$download_heatmap <- downloadHandler(
  #   filename = function(){
  #     paste0("Heatmap",".pdf")
  #   },
  #   content = function(file) {
  #     pdf(file)
  #     print(heatmap_plot())
  #     dev.off()
  #   }
  # )
  
  ############################
  ######### heatmap ##########
  ############################
  
  ####### heatmap renderUI commented ######### 
  # output$expand_genonames <- renderUI({
  #   type <- input$file_type
  #   if(type=='norm'){
  #     DS <- df_shiny()
  #   }else if(type=='raw'){
  #     DS <- df_raw_shiny()
  #   }
  #   if(ncol(DS)==input$numOfGeno){
  #     lapply(1:input$numOfGeno, function(i) {
  #       textInput(paste('type',i,sep=""), paste('Type',i,sep=" "),value = colnames(DS)[i])
  #     })
  #   }else{
  #     lapply(1:input$numOfGeno, function(i) {
  #       textInput(paste('type',i,sep=""), paste('Type',i,sep=" "))
  #     })
  #   }
  # })
  # 
  # output$refGeno <- renderUI({
  #   selectInput('heatmap_anchor',"Reference genotype",choices=c(1:input$numOfGeno))
  # })
  # 
  output$heatmap_display <- renderUI({
    display <- "ALL"
    for (i in 1:input$numOfCluster){
      display <- c(display,i)
    }
    selectInput('display_cluster',"Display cluster",choices=display)
  })
  ################
  
  setOneWithinFold <- function(arr){ #logFC
    fold <- as.numeric(input$fold)
    for (i in 1:length(arr)){
      if((arr[i]<=(fold)) & (arr[i]>=(1/fold))){
        arr[i]=1
      }
    }
    return(arr)
  }
  
  plotHeatmap <- eventReactive(input$heatmap_plot, {#process and return data 
    heatmap.start.time <- Sys.time()
    
    type <- input$file_type
    value <- input$heatmap_value
    de_type <- input$heatmap_de_ind
      
    if(type=='norm'){
      DS <- df_shiny()
    }else if(type=='raw'){
      DS <- df_raw_shiny()
    }
    names <- NULL
    # for(i in 1:input$numOfGeno){
    #   id <- paste('type',i,sep="")
    #   names <- c(names,input[[id]])
    # }
    
    # numOfGeno <- input$numOfGeno
    # ref <- as.numeric(input$heatmap_anchor)
    clusterNum <- input$numOfCluster
    
    if(de_type == "ind"){
      fold <- as.numeric(input$fold)
      fold_ncol <- input$fold_ncol
      DS2 <- deWithoutStats(DS, FC=fold, n_col=fold_ncol)
      de_genes <- rownames(DS2)
      # print("from heatmap YT version")
      # print("de genes")
      # print(head(de_genes))
      # print(paste(length(de_genes),"genes"))
      
    } else if(de_type == "de"){
      res.df <- de_no_filt()
      if(is.null (res.df) ) return(NULL)
      p_val <- input$p_val
      fc <- input$fc
      rep_number <- input$n_rep
      res.df.filt <- de_filt(res.df, p_val,fc,rep_number) 
      de_genes <- res.df.filt$Gene
      # print("from line heatmap de result from DE analysis")
      # print("res.df.filt")
      # print(head(res.df.filt))
    }
      
    de_genes_exp <- DS[rownames(DS)%in%de_genes,]
    DS3 <- t(scale(t(de_genes_exp)))
    DS3 <- na.omit(DS3)
    # print("from line 1894 - heatmap de type")
    # print("DS3")
    # print(head(DS3))
    
    set.seed(110)
    a <- ComplexHeatmap::Heatmap(DS3, name="Normalized expression",
                                 col = colorRamp2(c(min(DS3),0,max(DS3)), c("red","black", "green")),
                                 row_names_gp = gpar(fontsize = 1), 
                                 row_dend_gp = gpar(fontsize = 1),
                                 row_title_gp = gpar(fontsize = 10),
                                 cluster_columns = FALSE,
                                 row_dend_width = unit(3, "cm"),
                                 split = clusterNum, clustering_distance_rows = "pearson",                       
                                 show_heatmap_legend = TRUE,
                                 show_row_names = FALSE, show_column_names = T,
                                 heatmap_legend_param = list(title = "Normalized expression") )
    set.seed(110)
    rcl.list <- row_order(a); 
    DS3.1 <- as.matrix(rownames(DS3))

    # Cluster <- NULL
    # for(i in 1:length(rcl.list)){
    #   for(j in 1:length(rcl.list[[i]])){
    #     pair <- c(i,DS3.1[rcl.list[[i]][j]])
    #     Cluster <- rbind(Cluster,pair)
    #   }
    # }
    # Cluster <- data.frame(Cluster,row.names = NULL)
    # colnames(Cluster) <- c("cluster","GeneID")
    
    rcl.list2 <- rcl.list
    for(i in 1:length(rcl.list)){
      rcl.list2[[i]] <- rownames(DS3)[rcl.list[[i]] ]
    }

    for(i in 1:length(rcl.list2)){
      genes <- rcl.list2[[i]]
      group_name <- rep(i, length(genes) )
      Cluster_i <- data.frame("GeneID"=genes, "Cluster"=i)
      if( i ==1 ){
        Cluster <- Cluster_i
      } else{
        Cluster <- rbind(Cluster, Cluster_i)
      }
    }
    
    print("line 2089, Cluster")
    print(head(Cluster))
    print(paste0("de_type = ",de_type))
    print("length of input DS")
    print(dim(DS3))
    
    # end heat map analysis
    heatmap.end.time <- Sys.time()
    print("heat map time")
    print(heatmap.end.time - heatmap.start.time)
    return (list(a,DS3,Cluster))
  })
  
  getCluster <- eventReactive(input$heatmap_plot, {
    set.seed(110) 
    ll <- plotHeatmap(); a <- ll[[1]]; DS3 <- ll[[2]]
    rcl.list <- row_order(a)
    DS3.1 <- as.matrix(rownames(DS3))
    
    Cluster <- NULL
    for(i in 1:length(rcl.list)){
      for(j in 1:length(rcl.list[[i]])){
        pair <- c(i,DS3.1[rcl.list[[i]][j]])
        Cluster <- rbind(Cluster,pair)
      }
    }
    Cluster <- data.frame(Cluster,row.names = NULL)
    colnames(Cluster) <- c("cluster","gene.id")
    
    return (Cluster[,c("gene.id", "cluster")])
  })
  
  mapPlot <- function(){
    myHeatmap <- plotHeatmap()[[1]]
    myHeatmap <- draw(myHeatmap)
  }
  
  output$heatmap.plot <- renderPlot({
    mapPlot()
  })
  
  output$downloadheatmap <- downloadHandler(
    filename = function(){
      paste("heatmap",".pdf",sep="")
    },
    content = function(file){
      pdf(file)
      p <- mapPlot()
      dev.off()
    }
  )
  
  
  output$cluster.info <- DT::renderDataTable({
    clusternum <- input$display_cluster
    gl <- plotHeatmap()[[3]]  #getCluster()
    if(! is.null(gl) ) {
      if(clusternum=="ALL"){
        gl
      }else{
        clusternum <- as.numeric(clusternum)
        dplyr::filter(gl,cluster==clusternum)
      }
    }
  })
  
  output$downloadclusters <- downloadHandler(
    filename = function(){
      paste("genelist",".csv",sep="")
    },
    content = function(file){
      gl <- plotHeatmap()[[3]]
      write.csv(gl,file, row.names = FALSE)
    }
  )
  
  ############################
  ########## noise ###########
  ############################
  
  SQCO <- function(MT){
    if(ncol(MT)==1){
      res <- matrix(0)
      return(res)
    }
    temp <- NULL
    for(i in 1:nrow(MT)){
      m <- sum(MT[i,])/length(MT[i,])
      v <- stats::var(MT[i,])
      if(m!=0)
        temp <- c(temp,v/(m*m))
    }
    res <- matrix(sum(temp)/length(temp))
    return (res)
  }
  
  output$expand_genonames_noise <- renderUI({
    type <- input$file_type
    if(type=='norm'){
      DS <- df_shiny()
    }else if(type=='raw'){
      DS <- df_raw_shiny()
    }
    numOfRep <- as.numeric(input$noise_numOfRep)
    numOfGeno <- ncol(DS)/numOfRep
    
    lapply(1:numOfGeno, function(i) {
      textInput(paste('noisetype',i,sep=""), paste('Type',i,sep=" "),value=colnames(DS)[(i-1)*numOfRep+1])
    })
  })
  
  output$noise_anchor_choices <- renderUI({
    type <- input$file_type
    if(type=='norm'){
      DS <- df_shiny()
    }else if(type=='raw'){
      DS <- df_raw_shiny()
    }
    numOfRep <- as.numeric(input$noise_numOfRep)
    numOfGeno <- ncol(DS)/numOfRep
    names <- NULL
    for(i in 1:numOfGeno){
      id <- paste('noisetype',i,sep="")
      names <- c(names,input[[id]])
    }
    selectInput('noise_anchor_b',"Anchor genotype",choices = names)
  })
  
  noisePlot <- eventReactive(input$noise_plot, {
    noise.start.time <- Sys.time()
    type <- input$file_type
    if(type=='norm'){
      DS <- df_shiny()
    }else if(type=='raw'){
      DS <- df_raw_shiny()
    }
    numOfRep <- as.numeric(input$noise_numOfRep)
    numOfGeno <- ncol(DS)/numOfRep
    graph <- input$noise_graph_type
    names <- NULL
    for(i in 1:numOfGeno){
      id <- paste('noisetype',i,sep="")
      names <- c(names,input[[id]])
    }
    
    situation <- input$noise_situation
    if(situation=='a'){
      DS1 <- list()
      for(j in 1:numOfGeno){
        DS1[[j]] <- as.matrix(DS[,((j-1)*numOfRep+1):(j*numOfRep)]) 
      } 
      Noise <- NULL
      for(y in 1:numOfGeno){
        Noise <- c(Noise,SQCO(DS1[[y]]))
      }
      xform <- list(categoryorder = "array",
                    categoryarray = names)
      if(graph=="Bar chart"){
        p <- plot_ly(
          x = names,
          y = Noise,
          type = "bar"
        ) %>% layout(xaxis = xform)
      }else if(graph=="Line chart"){
        p <- plot_ly(
          x = names,
          y = Noise,
          type="scatter",
          mode="lines+markers"
        ) %>% layout(xaxis = xform,yaxis=list(range = c(0, max(Noise)+0.001)))
      }
    }else if(situation=='b'){
      DS_ave <- NULL
      for(j in 1:numOfGeno){
        part_DS <- as.matrix(DS[,((j-1)*numOfRep+1):(j*numOfRep)]) 
        DS_ave <- cbind(DS,data.frame(matrixStats::rowMeans2(part_DS)))
      }
      anchor <- input$noise_anchor_b
      names <- NULL
      for(i in 1:numOfGeno){
        id <- paste('noisetype',i,sep="")
        names <- c(names,input[[id]])
      }
      anchor_index <- match(anchor,names)
      Noise <- NULL
      for(i in 1:numOfGeno){
        if(i!=anchor_index){
          Noise <- c(Noise,SQCO(cbind(DS_ave[,anchor_index],DS_ave[,i])))
        }
      }
      names_wo_anchor <- NULL
      for(i in 1:numOfGeno){
        if(i!=anchor_index){
          id <- paste('noisetype',i,sep="")
          names_wo_anchor <- c(names_wo_anchor,input[[id]])
        }
      }
      xform <- list(categoryorder = "array",
                    categoryarray = names_wo_anchor)
      if(graph=="Bar chart"){
        p <- plot_ly(
          x = names_wo_anchor,
          y = Noise,
          type = "bar"
        ) %>% layout(xaxis = xform)
      }else if(graph=="Line chart"){
        p <- plot_ly(
          x = names_wo_anchor,
          y = Noise,
          type="scatter",
          mode="lines+markers"
        ) %>% layout(xaxis = xform,yaxis=list(range = c(0, max(Noise)+0.001)))
      }
    }else if(situation=='c'){
      anchor <- input$noise_anchor_c
      names <- colnames(DS)
      anchor_index <- match(anchor,names)
      Noise <- NULL
      for(i in 1:ncol(DS)){
        if(i!=anchor_index){
          Noise <- c(Noise,SQCO(cbind(DS[,anchor_index],DS[,i])))
        }
      }
      names_wo_anchor <- names[-anchor_index]
      xform <- list(categoryorder = "array",
                    categoryarray = names_wo_anchor)
      if(graph=="Bar chart"){
        p <- plot_ly(
          x = names_wo_anchor,
          y = Noise,
          type = "bar"
        ) %>% layout(xaxis = xform)
      }else if(graph=="Line chart"){
        p <- plot_ly(
          x = names_wo_anchor,
          y = Noise,
          type="scatter",
          mode="lines+markers"
        ) %>% layout(xaxis = xform,yaxis=list(range = c(0, max(Noise)+0.001)))
      }
    }
    noise.end.time <- Sys.time()
    print("noise time")
    print(noise.end.time - noise.start.time)
    return (p)
  })
  
  output$noise.plot <- renderPlotly({
    noisePlot()
  })
  
  output$downloadnoise <- downloadHandler(
    filename = function(){
      paste("noise",".png",sep="")
    },
    content = function(file){
      p <- noisePlot()
      export(p, file = "noise.png")
    }
  )
  
  ############################
  ######### entropy ##########
  ############################
  
  computeBin <- function(arr){  # Doane's rule
    n <- length(arr)
    gx <- moments::skewness(arr)
    sigmag <- sqrt(6*(n-2)/((n+1)*n+3))
    bin <- 1+log2(n)+log2(1+abs(gx)/sigmag)
    return (bin)
  }
  
  getBinCounts <- function(arr){
    vec <- entropy::discretize(arr,computeBin(arr),r=range(arr))
    return (vec)
  }
  
  output$expand_genonames_entropy <- renderUI({
    type <- input$file_type
    if(type=='norm'){
      DS <- df_shiny()
    }else if(type=='raw'){
      DS <- df_raw_shiny()
    }
    tp <- as.numeric(input$entropy_timepoints)
    numOfGeno <- ncol(DS)/tp
    
    lapply(1:numOfGeno, function(i) {
      textInput(paste('entropytype',i,sep=""), paste('Type',i,sep=" "),value=colnames(DS)[(i-1)*tp+1])
    })
  })
  
  entropyPlot <- reactive({
    entropy.start.time <- Sys.time()
    type <- input$file_type
    if(type=='norm'){
      DS <- df_shiny()
    }else if(type=='raw'){
      DS <- df_raw_shiny()
    }
    if(is.null(DS)==FALSE){
      tsflag <- input$tsflag
      graph <- input$entropy_graph_type
      names <- colnames(DS)
      xform <- list(categoryorder = "array",
                    categoryarray = names)
      entropy.vector <- NULL #entropy of each column
      for (i in 1:length(DS)){
        binCount <- getBinCounts(DS[,i])
        entropy <- entropy.empirical(binCount,unit="log2")
        entropy.vector <- c(entropy.vector,entropy)
      }
      if(tsflag==FALSE){
        if(graph=="Bar chart"){
          p <- plot_ly(
            x = names,
            y = entropy.vector,
            type = "bar"
          ) %>% layout(xaxis = xform)
        }else if(graph=="Line chart"){
          p <- plot_ly(
            x = names,
            y = entropy.vector,
            type="scatter",
            mode="lines+markers"
          ) %>% layout(xaxis = xform)
          #yaxis=list(range = c(0, max(ent)+0.002))
        }
      }else if(tsflag==TRUE){
        tp <- as.numeric(input$entropy_timepoints)
        numOfGeno <- ncol(DS)/tp
        names <- NULL
        for(i in 1:numOfGeno){
          id <- paste('entropytype',i,sep="")
          names <- c(names,input[[id]])
        }
        time_index <- c(1:tp)
        ent <- data.frame(time_index)
        for(j in 1:numOfGeno){
          part_ent <- entropy.vector[(tp*j-(tp-1)):(tp*j)]
          ent <- cbind(ent,part_ent)
        }
        if(graph=="Bar chart"){
          p <- plot_ly(x=ent[,1],y=ent[,2],name=names[1],type='bar')
          for(i in 1:(numOfGeno-1)){
            p <- add_trace(p,y=ent[,i+2],name=names[i+1],type='bar')
          }
          p <- layout(p,xaxis = list(title = 'Time'),yaxis=list(title='Entropy'))
        }else if(graph=="Line chart"){
          p <- plot_ly(x=ent[,1],y=ent[,2],name=names[1],type='scatter',mode='lines+markers')
          for(i in 1:(numOfGeno-1)){
            p <- add_trace(p,y=ent[,i+2],name=names[i+1],type='scatter',mode='lines+markers')
          }
          p <- layout(p,xaxis = list(title = 'Time'),yaxis=list(title='Entropy'))
        }
      }
      entropy.end.time <- Sys.time()
      print("entropy time")
      print(entropy.end.time - entropy.start.time)
      return (p)
    }
  })
  
  output$entropy.plot <- renderPlotly({
    entropyPlot()
  })
  
  output$downloadentropy <- downloadHandler(
    filename = function(){
      paste("entropy",".png",sep="")
    },
    content = function(file){
      p <- entropyPlot()
      export(p, file = "entropy.png")
    }
  )
  
  ############################
  ###### gene ontology #######
  ############################
  
  go_res <- eventReactive(input$submit_go, {
    go.start <- Sys.time()
    genes <- gene_list()
    if(is.null(genes)){
      showModal(modalDialog(
        title = "Error","Please check your input DE genes!"
      ))
      return (NULL)
    }
    go_method <- input$go_method
    if(go_method %in% c("clusterProfiler","GOstats")){
      dbs_name <- input$go_species
      orgDb <- DBS[[dbs_name]]
      keyType <- input$go_geneidtype
      ont <- input$subontology
      bg <- bg_list()
      if(go_method == "clusterProfiler"){
        temp <- goPrep(fg=genes,bg=bg,keyType=keyType,orgDb=orgDb) # list(fg_mapped, bg_mapped, fg_unmapped)
        if(is.null (temp) )
          return(NULL)
        res <- enrichgoApply(gene_list=genes, keyType, orgDb,ont=ont, pvalueCutoff=1, qvalueCutoff=1)
      } else if(go_method == "GOstats"){
        temp <- goPrep(fg=genes,bg=bg,keyType=keyType,orgDb=orgDb) # list(fg_mapped, bg_mapped, fg_unmapped)
        if(is.null (temp) )
          return(NULL)
        fg_mapped <- temp[[1]] ; bg_mapped <- temp[[2]]
        if(dbs_name == "org.Sc.sgd.db"){
          primary_id <-  "ORF"
        } else if(dbs_name == "org.At.tair.db"){
          primary_id <- "TAIR"
        } else{
          primary_id <- "ENTREZID"
        }
        res <- gostatsApply(fg_mapped, bg_mapped, keyType, orgDb, ont=ont,primary_id=primary_id, pvalueCutoff=1)
      }
    } else if(go_method == "enrichR"){
      if (!havingIP()){
        showModal(modalDialog(
          title = "Error","enrichR requires internet connection. Please check your internet connection and try again!"
        ))
        return(NULL)
      }
      dbs <- input$enrichR_dbs
      res <- enrichrApply(genes,dbs)
      if (nrow(res) == 0){
        showModal(modalDialog(
          title = "Warning","No term matched! Please ensure all input DE genes are in SYMBOL format!"
        ))
      }
    }
    
    if(is.null(res)){
      showModal(modalDialog(
        title = "Warning","NULL result given! Organism or identifier might not be correct! Or try another method"
      ))
    }
    go.end <- Sys.time()
    print("Go time")
    print(go.end - go.start)
    return(res)   # no filter by min counts in each go term
  })
  
  go_res_filt <- function(){
    res <- go_res()
    if(! is.null (res)){
      if(! is.data.frame(res) ) res <- as.data.frame(res)
      go_method <- input$go_method
      min_no <- input$go_min_no
      max_p <- input$go_max_p
      res_filt <- filter(res, Count>=min_no)   # filter by gene counts in go terms
      if(go_method == "clusterProfiler"){
        res_filt <- filter(res, p.adjust <= max_p)
      } else if(go_method == "GOstats"){
        res_filt <- filter(res, Pvalue <= max_p)
      }
      # res_filt
    } else{
      res_filt = NULL
    }
    return(res_filt)
  }
  
  output$go_table <- DT::renderDataTable({
    go_res_filt()
  })
  
  go_pie_res <- eventReactive(input$submit_go, {
    genes <- gene_list()
    if(is.null(genes))
      return (NULL)
    go_method <- input$go_method
    if(go_method == "enrichR")
      return (NULL)
    dbs_name <- input$go_species
    orgDb <- DBS[[dbs_name]]
    keyType <- input$go_geneidtype
    ont <- input$subontology
    bg <- bg_list()

    temp <- goPrep(fg=genes,bg=bg,keyType=keyType,orgDb=orgDb) # list(fg_mapped, bg_mapped, fg_unmapped)
    if(is.null (temp) )
      return(NULL)
    GOresult <- groupGO(genes, keyType= keyType, OrgDb = orgDb, ont = ont, level = 2,
                        readable = FALSE)
    res <- dplyr::filter(GOresult@result,Count!=0)
    res$Term <- paste(res$Description," (",res$ID,")",sep="")
    res <- dplyr::arrange(res,desc(Count))
    res2 <- res[,c("Term","Count")]
    res2 <- rbind(res2,c("Unclassified",length(temp[[3]])))
    return(res2)
  })
  
  goPie <- function(){
    res <- go_pie_res()
    print("from goPie line 2543 print res")
    print(res)
    if(! is.null(res) ){
      p <- plot_ly(res, labels = ~Term, values = ~Count, type = 'pie',
                   textposition = 'inside',
                   textinfo = 'label+percent',
                   insidetextfont = list(color = '#FFFFFF'),
                   #hoverinfo = 'text',
                   #text = ~paste('$', X1960, ' billions'),
                   marker = list(line = list(color = '#FFFFFF', width = 1)),
                   #The 'pull' attribute can also be used to create space between the sectors
                   showlegend = FALSE) %>%
        layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
      # title = 'Gene ontology',
      return (p)
    }
  }
  
  output$go_pie <- renderPlotly({
    goPie()
  })
  
  goList <- eventReactive(input$submit_go_graph, {
    go_method <- input$go_method
    go_terms <- input$go_term_slect
    res <- go_res_filt()
    # print("from goList")
    # print(go_method)
    # print(head(res))
    
    if(go_method=="clusterProfiler" & !is.null(res) ){
      res_filt <- res[res$Description%in%go_terms,]
      go_list <- goToList(res_filt)
      # print(head(go_list)) ##
      return(go_list)
    } else {
      return(NULL)
    }
  })
  
  go_graph_plot <- function() {
    go_list <-goList()
    show_gene_names <- input$show_gene_names
    if(! is.null (go_list)){
      temp <- graphGene(go_list); p <- temp[[1]]; q <- temp[[2]]
      if(show_gene_names){
        p <- p + geom_node_text(aes_(label=~name), repel=TRUE)
      }
    } else {
      return(NULL)
    }
    return(p)
  }
  
  output$go_graph <- renderPlot({
    # go_list <-goList()
    # if(! is.null (go_list)){
    #   p <- graphGene(go_list)
    #   p
    # }
    go_graph_plot()
  })
  
  output$download_go_table <- downloadHandler(
    filename = function() {
      paste0("GO analysis",".csv")
    },
    content = function(file) {
      res <- go_res()
      min_no <- input$go_min_no
      res_filt <- filter(res, Count>=min_no) 
      write.csv(res_filt,file, row.names = FALSE)
    }
  )
  
  output$download_go_pie <- downloadHandler(
    filename = function() {
      paste0("GO term lvl 2",".png")
    },
    content = function(file) {
      p <- goPie()
      export(p, file = paste("GO term lvl 2",".png",sep=""))
    }
  )
  
  output$download_go_graph <- downloadHandler(
    filename = function () {
      paste0("GO graph",".", input$download_go_graph_type )
    },
    content = function(file) {
      if(input$download_go_graph_type == "pdf"){
        pdf(file, width=15, height=7)
        print(go_graph_plot())
        dev.off()
      } else if(input$download_go_graph_type == "png" ) {
        png(file,  width=10, height=7)
        print(go_graph_plot())
        dev.off()
      }
    }
    
  )
  
  #session$onSessionEnded(stopApp)
}

shinyApp(ui,server)
