# Sys.setenv("plotly_username"=" your_plotly_username")
# Sys.setenv("plotly_api_key"="your_api_key")
## test repo

# wd <- dirname(rstudioapi::getActiveDocumentContext()$path) # set wd as the current folder
# print(wd == getwd())
# print(wd)
# print(getwd())
# if (!wd == getwd()) {
#   setwd(wd)
# }

print("start loading")
start.load <- Sys.time() ### time

if (length(find.package(package = "shiny", quiet = T)) > 0) {
  library(shiny)
} else {
  print("Package shiny not installed")
  install.packages("shiny")
  print("Package shiny installed")
  library(shiny)
}

if (length(find.package(package = "cyjShiny", quiet = T)) > 0) {
  library(cyjShiny)
} else {
  remotes::install_github("cytoscape/cyjShiny")
  library(cyjShiny)
}

if (length(find.package(package = "shinythemes", quiet = T)) > 0) {
  library(shinythemes)
} else {
  print("Package shinythemes not installed")
  install.packages("shinythemes")
  print("Package shinythemes installed")
  library(shinythemes)
}

if (length(find.package(package = "rstudioapi", quiet = T)) > 0) {
  library(rstudioapi)
} else {
  install.packages("rstudioapi")
  library(rstudioapi)
}

if (!length(find.package(package = "rlang", quiet = T)) > 0) {
  install.packages("rlang")
}

#################################
if (length(find.package(package = "RColorBrewer", quiet = T)) > 0) {
  library(RColorBrewer)
} else {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

if (length(find.package(package = "kohonen", quiet = T)) > 0) {
  library(kohonen)
} else {
  install.packages("kohonen")
  library(kohonen)
}

if (length(find.package(package = "sparsepca", quiet = T)) > 0) {
  library(sparsepca)
} else {
  install.packages("sparsepca")
  library(sparsepca)
}

if (length(find.package(package = "randomForest", quiet = T)) > 0) {
  library(randomForest)
} else {
  install.packages("randomForest")
  library(randomForest)
}

if (length(find.package(package = "cluster", quiet = T)) > 0) {
  library(cluster)
} else {
  install.packages("cluster")
  library(cluster)
}


## for t-sne
if (length(find.package(package = "reticulate", quiet = T)) > 0) {
  library(reticulate)
} else {
  install.packages("reticulate")
  library(reticulate)
}

if (length(find.package(package = "Rtsne", quiet = T)) > 0) {
  library(Rtsne)
} else {
  install.packages("Rtsne")
  library(Rtsne)
}


####################### Dependencies For RAFSIL ###################################
if (length(find.package(package = "RAFSIL", quiet = T)) > 0) {
  library(RAFSIL)
}

if (length(find.package(package = "gridGraphics", quiet = T)) > 0) {
  library(gridGraphics)
} else {
  install.packages("gridGraphics")
  library(gridGraphics)
}

if (length(find.package(package = "gridExtra", quiet = T)) > 0) {
  library(gridExtra)
} else {
  install.packages("gridExtra")
  library(gridExtra)
}

if (length(find.package(package = "tidyverse", quiet = T)) > 0) {
  library(tidyverse)
} else {
  install.packages("tidyverse", dependencies = TRUE)
  library(tidyverse)
}

if (length(find.package(package = "ggpubr", quiet = T)) > 0) {
  library(ggpubr)
} else {
  install.packages("ggpubr")
  library(ggpubr)
}

if (length(find.package(package = "networkD3", quiet = T)) > 0) {
  library(networkD3)
} else {
  install.packages("https://cran.r-project.org/src/contrib/Archive/networkD3/networkD3_0.2.10.tar.gz", repo=NULL, type="source")
  library(networkD3)
}

if (length(find.package(package = "data.tree", quiet = T)) > 0) {
  library(data.tree)
} else {
  install.packages("data.tree")
  library(data.tree)
}


if (length(find.package(package = "bubbles", quiet = T)) > 0) {
  library(bubbles)
} else {
  devtools::install_github("jcheng5/bubbles", upgrade = FALSE)
  library(bubbles)
}

###################################################################################
####################### Dependencies For Uniprot ###################################

if (length(find.package(package = "UniprotR", quiet = T)) > 0) {
  library(UniprotR)
} else {
  install.packages("UniprotR")
  library(UniprotR)
}

if (length(find.package(package = "scales", quiet = T)) > 0) {
  library(scales)
} else {
  install.packages("scales")
  library(scales)
}

###################################################################################

####################### Dependencies For Pathway Enrichment ###################################

if (length(find.package(package = "gprofiler2", quiet = T)) > 0) {
  library(gprofiler2)
} else {
  install.packages("gprofiler2")
  library(gprofiler2)
}

###################################################################################

####################### Dependencies For Protein Interactions ###################################

if (length(find.package(package = "httr", quiet = T)) > 0) {
  library(httr)
} else {
  install.packages("httr")
  library(httr)
}

if (length(find.package(package = "curl", quiet = T)) > 0) {
  library(curl)
} else {
  install.packages("curl")
  library(curl)
}

if (length(find.package(package = "later", quiet = T)) > 0) {
  library(later)
} else {
  install.packages("later")
  library(later)
}

if (length(find.package(package = "qdapTools", quiet = T)) > 0) {
  library(qdapTools)
} else {
  install.packages("qdapTools")
  library(qdapTools)
}

if (length(find.package(package = "alakazam", quiet = T)) > 0) {
  library(alakazam)
} else {
  install.packages("https://cran.r-project.org/src/contrib/Archive/alakazam/alakazam_1.0.0.tar.gz", repo=NULL, type="source")
  library(alakazam)
}

if (length(find.package(package = "msa", quiet = T)) > 0) {
  library(msa)
} else {
  BiocManager::install("msa", update = FALSE)
  library(msa)
}

if (length(find.package(package = "ape", quiet = T)) > 0) {
  library(ape)
} else {
  install.packages("ape")
  library(ape)
}

if (length(find.package(package = "seqinr", quiet = T)) > 0) {
  library(seqinr)
} else {
  install.packages("seqinr")
  library(seqinr)
}

if (length(find.package(package = "qdapRegex", quiet = T)) > 0) {
  library(qdapRegex)
} else {
  install.packages("qdapRegex")
  library(qdapRegex)
}

###################################################################################

####################### Dependencies For Co-expression ###################################

if (length(find.package(package = "shinyjs", quiet = T)) > 0) {
  library(shinyjs)
} else {
  install.packages("shinyjs")
  library(shinyjs)
}

###################################################################################

####################### Dependencies For Microarray ###################################

if (length(find.package(package = "devtools", quiet = T)) > 0) {
  library(devtools)
} else {
  install.packages("devtools")
  library(devtools)
}

if (length(find.package(package = "remotes", quiet = T)) > 0) {
  library(remotes)
} else {
  devtools::install_github("r-lib/remotes")
  library(remotes)
}

if (length(find.package(package = "maEndToEnd", quiet = T)) > 0) {
  suppressPackageStartupMessages({library("maEndToEnd")})
} 
# else {
#   remotes::install_github("b-klaus/maEndToEnd", ref="master")
#   suppressPackageStartupMessages({library("maEndToEnd")})
# }

if (length(find.package(package = "oligoClasses", quiet = T)) > 0) {
  library(moments)
} else {
  print("Package oligoClasses not installed")
  install.packages("oligoClasses")
  print("Package oligoClasses installed")
  library(oligoClasses)
}

if (length(find.package(package = "ArrayExpress", quiet = T)) > 0) {
  library(moments)
} else {
  print("Package ArrayExpress not installed")
  install.packages("ArrayExpress")
  print("Package ArrayExpress installed")
  library(ArrayExpress)
}

if (length(find.package(package = "pd.hugene.1.0.st.v1", quiet = T)) > 0) {
  library(moments)
} else {
  print("Package pd.hugene.1.0.st.v1 not installed")
  install.packages("pd.hugene.1.0.st.v1")
  print("Package pd.hugene.1.0.st.v1 installed")
  library(pd.hugene.1.0.st.v1)
}

if (length(find.package(package = "hugene10sttranscriptcluster.db", quiet = T)) > 0) {
  library(moments)
} else {
  print("Package hugene10sttranscriptcluster.db not installed")
  install.packages("hugene10sttranscriptcluster.db")
  print("Package hugene10sttranscriptcluster.db installed")
  library(hugene10sttranscriptcluster.db)
}

if (length(find.package(package = "arrayQualityMetrics", quiet = T)) > 0) {
  library(moments)
} else {
  print("Package arrayQualityMetrics not installed")
  install.packages("arrayQualityMetrics")
  print("Package arrayQualityMetrics installed")
  library(arrayQualityMetrics)
}

if (length(find.package(package = "limma", quiet = T)) > 0) {
  library(moments)
} else {
  print("Package limma not installed")
  install.packages("limma")
  print("Package limma installed")
  library(limma)
}

if (length(find.package(package = "topGO", quiet = T)) > 0) {
  library(moments)
} else {
  print("Package topGO not installed")
  install.packages("topGO")
  print("Package topGO installed")
  library(topGO)
}

if (length(find.package(package = "ReactomePA", quiet = T)) > 0) {
  library(moments)
} else {
  print("Package ReactomePA not installed")
  install.packages("ReactomePA")
  print("Package ReactomePA installed")
  library(ReactomePA)
}

###################################################################################
########################### Style files for Cytoscape.js ################

styles <- c(
  "generic style"="./www/style/basicStyle.js",
  "red-yellow"="./www/style/red-yellow.js",
  "red-pink" = "./www/style/red-pink.js",
  "green-blue"="./www/style/green-blue.js",
  "green-blue(ppi)"="./www/style/green-blue(ppi).js")

##########################################################################

#
# ## sourcing util files
source("./www/utils.R")
source("./www/PhyscochemicalSep.R")

#
loadPkg()

id_to_name <- read.csv(paste0("./www/TransTable_Human.csv"))


#################### Complex Enrichment ##########################

complexes <- load(paste0("./www/allComplexes.RData")) #allComplexes is masked under complexes
up_corum_mapping <- read.csv(paste0("./www/UniProt_CORUM_Mapping.csv"))

##################################################################


end.load <- Sys.time()
print("loading time")
print(end.load - start.load)

##### UI from here ###########
ui <- tagList(
  shinyjs::useShinyjs(),
  navbarPage(
    id = "navbar",
    theme = shinytheme("flatly"),
    title = "",
    tabPanel(
      "GeneCloudOmics",
      br(),
      sidebarLayout(
        sidebarPanel(
          img(
            src = "GeneCloudOmics-logo.png",
            width = "100%", height = "100%"
          )
        ),
        mainPanel(
          h2("Welcome to GeneCloudOmics", align = "center",style = "color:#73C6B6;font-weight: bold;"),
          p("The Biostatistical Tool for Gene Expression Data Analysis", align = "center"),
          h4(span("GeneCloudOmics", style = "color:#73C6B6;font-weight: bold;")," is a web server for transcriptome data analysis and visualization. It supports the analysis of 
      microarray and RNASeq data and performs ten different bio-statistical analyses that cover the common analytics for gene expression data. Furthermore, 
      it gives the users access to several bioinformatics tools to perform 12 different bioinformatics analyses on gene/protein datasets."),
          h4(span("GeneCloudOmics", style = "color:#73C6B6;font-weight: bold;"),"  is designed as a one-stop server that helps the users perform all tasks through an intuitive graphical 
      user interface (GUI) that waves the hassle of coding, installing tools,  packages or libraries and dealing with operating systems compatibility and versioning issues, some of 
      the complications that make data analysis tasks more challenging for biologists. GeneCloudOmics is an open-source tool and the website is free and open to all users 
      and there is no login requirement."),
          h4(span("Supported Transcriptome data:", style = "color:#73C6B6;font-weight: bold;"), " RNA-Seq and Microarray "),
          h4(span("Data Preprocessing:", style = "color:#73C6B6;font-weight: bold;"), " GeneCloudOmics performs raw data normalization using four normalization methods RPKM, 
      FPKM, TPM and RUV. The raw vs. normalized data are visualized as boxplots and violin plots."),
          h4(span("Differential Gene Expression (DGE) Analysis:", style = "color:#73C6B6;font-weight: bold;"), " GDE using five methods EdgeR, DESeq2, NOISeq, and LIMMA."),
          h4(span("Bio-statistical Analysis:", style = "color:#73C6B6;font-weight: bold;"), " GeneCloudOmics provides the user with the following bio-statistical analyses: 
      Pearson and Spearman rank correlations, PCA, k-means and hierarchical clustering, 
      Shannon entropy and noise (square of the coefficient of variation), t-SNE, random forest and SOM analyses. All analyses include proper high-resolution visualization."),
          h4(span("Bioinformatics Analysis of Gene and Protein sets:", style = "color:#73C6B6;font-weight: bold;"), " For the differential expressed genes (DEG), GeneCloudOmics provides 
      the users with multiple bioinformatics tools to investigate their 
      gene/protein list including gene ontology (GO), pathway enrichment analysis, PPI, co-expression, gene/protein function, subcellular localization, complex enrichment, protein domains, tissue expression, sequence properties (acidity, hydrophobicity and charge),
       evolutionary analysis (gene tree, phylogenetic tree and species/chromosome location)  and pathological analysis (diseases that these genes/proteins are involved in). The analyses include proper high-resolution visualization, when applicable."),
        )
      )
    ),
    navbarMenu('Preprocessing',
               tabPanel(
                 "RnaSeq Data",
                 value = "active_tab_rnaseq",
                 tabsetPanel(
                   id = "Rnaseq_pre",
                   tabPanel(
                     "Upload data",
                     sidebarPanel(
                       radioButtons(
                         "file_type", "Choose File Type",
                         c("Raw file (read count)" = "raw", "Normalised file" = "norm")
                       ),
                       conditionalPanel(
                         condition = "input.file_type=='raw'", # raw
                         withTags({
                           div(class="header", checked=NA,
                               p("Example ", a(href="https://github.com/buithuytien/GeneCloudOmics/blob/master/Test%20data/Eg_raw.png", "here"))
                           )
                         }),
                         fileInput("file1", "Choose Raw Counts"),
                         
                         withTags({
                           div(class="header", checked=NA,
                               p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/master/Test%20data/Eg_gene_length.png")), # ADD EXAMPLE
                           )
                         }),
                         fileInput("length1", "Choose Gene Length"), # gene id + length
                         
                         withTags({
                           div(class="header", checked=NA,
                               p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/master/Test%20data/Eg_negative_control_genes.png")), # ADD EXAMPLE
                           )
                         }),
                         fileInput("spikes1", "Choose Negative Control Genes")
                       ),
                       conditionalPanel(
                         condition = "input.file_type=='norm'", # normalized
                         withTags({
                           div(class = "header",
                               p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/master/Test%20data/Eg_normalised.png")), # ADD EXAMPLE
                           )
                         }),
                         fileInput("file2", "Choose Normalized Expression")
                         # helpText("* Format requirement: CSV file. Gene names in rows and genotypes in columns, following the usual format of files deposited in the GEO database.")
                       ),
                       withTags({
                         div(class = "header",
                             p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/master/Test%20data/Eg_metadata.png")), # ADD EXAMPLE
                         )
                       }),
                       fileInput("metafile1", "Choose Meta Data File"),
                       actionButton("submit_input", "Submit")
                       
                       
                     ),
                     mainPanel(
                       # h3("Welcome to GeneCloudOmics --"),
                       # h3("A Biostatistical tool for Transcriptomics Analysis"),
                       # img(
                       #   src = "GeneCloudOmics-logo.png",
                       #   width = 570, height = 370
                       # )
                     )
                   ),
                   tabPanel(
                     "Preprocessing",
                     sidebarPanel(
                       h4("Filtering"),
                       splitLayout(
                         numericInput("min_val", "Min. value", min = 0.1, step = 0.1, value = 1.0),
                         numericInput("min_col", "Min. columns", min = 1, value = 2)
                       ),
                       conditionalPanel(
                         condition = "input.file_type=='raw'",
                         radioButtons(
                           "norm_method", "Normalisation method",
                           c(
                             "None (Black)" = "None",
                             "RPKM (Blue)" = "RPKM", "FPKM (Dark cyan)" = "FPKM",
                             "TPM (Dark green)" = "TPM",
                             "RUV (Brown)" = "RUV"
                           )
                         )
                       ),
                       actionButton("submit_preprocessing", "Submit"),
                       conditionalPanel(
                         condition = "input.preprocessing_tabs == 'Data table' ",
                         br(),
                         br(),
                         downloadButton("download_norm_data", "Download table (csv)")
                       )
                     ),
                     mainPanel(
                       h3("Preprocessing Rnaseq data"),
                       tabsetPanel(
                         type = "tabs", id = "preprocessing_tabs",
                         tabPanel(
                           "RLE plot",
                           conditionalPanel(
                             condition = "$('html').hasClass('shiny-busy')",
                             div(img(src = "load.gif", width = 240, height = 180),
                                 h4("Processing ... Please wait"),
                                 style = "text-align: center;"
                             )
                           ),
                           conditionalPanel(
                             condition = "!$('html').hasClass('shiny-busy')",
                             plotOutput("RLE.plot2")
                           ),
                           
                           conditionalPanel(
                             condition = "input.file_type=='raw'",
                             conditionalPanel(
                               condition = "$('html').hasClass('shiny-busy')",
                               div(img(src = "load.gif", width = 240, height = 180),
                                   h4("Processing ... Please wait"),
                                   style = "text-align: center;"
                               )
                             ),
                             conditionalPanel(
                               condition = "!$('html').hasClass('shiny-busy')",
                               plotOutput("RLE.plot")
                             )
                           )
                         ),
                         tabPanel(
                           "Violin Plot",
                           conditionalPanel(
                             condition = "$('html').hasClass('shiny-busy')",
                             div(img(src = "load.gif", width = 240, height = 180),
                                 h4("Processing ... Please wait"),
                                 style = "text-align: center;"
                             )
                           ),
                           conditionalPanel(
                             condition = "!$('html').hasClass('shiny-busy')",
                             plotlyOutput("violin_plot2")
                           ),
                           conditionalPanel(
                             condition = "!$('html').hasClass('shiny-busy')",
                             plotlyOutput("violin_plot")
                           )
                         ),
                         tabPanel(
                           "Data table",
                           h3("Normalized data"),
                           DT::dataTableOutput("norm_table")
                         ),
                         tabPanel(
                           "Description table",
                           h3("Data description"),
                           DT::dataTableOutput("meta_table")
                         )
                       )
                     )
                   )
                 )
               ),
               tabPanel(
                 "Microarray Data",
                 value = "active_tab_micro",
                 sidebarPanel(
                   withTags({
                     div(class = "header",
                         p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/master/Test%20data/Eg_raw.png")), # ADD EXAMPLE ( have to change )
                     )
                   }),
                   fileInput("file_micro", "Choose Microarray Data"),
                   downloadButton("downloadMicroRaw", "Download Raw Data as CSV"),
                   br(), br(),
                   downloadButton("downloadMicroMeta", "Download Meta Data as CSV")
                 ),
                 
                 mainPanel(
                   h3("Preprocessing Microarray Data"),
                   conditionalPanel(
                     condition = "$('html').hasClass('shiny-busy')",
                     div(img(src = "load.gif", width = 240, height = 180),
                         h4("Processing ... Please wait"),
                         style = "text-align: center;"
                     )
                   )
                 )
               )
               
    ),
    navbarMenu('Transcriptome Analysis',
               tabPanel(
                 "    Scatter    ",
                 sidebarPanel(
                   selectInput(inputId = "scatter.x", label = "X-axis", choices = ""),
                   selectInput(inputId = "scatter.y", label = "Y-axis", choices = ""),
                   radioButtons(
                     "trans", "Transformation:",
                     c("None", "Natural log", "log2", "log10")
                   ),
                   checkboxInput("regline", "Display regression line", value = FALSE),
                   downloadButton("downloadscatter", "Download as PNG"),
                   h6("Download all pairs of samples in one PDF (this may take some time to run) :"),
                   downloadButton("downloadscatter_collage", "Download collage")
                 ),
                 mainPanel(
                   h3("Heatscatter"),
                   uiOutput("help_text_scatter"),
                   plotlyOutput("scatter.plot")
                 )
               ),
               tabPanel(
                 "Distribution Fit",
                 sidebarPanel(
                   conditionalPanel(
                     condition = "input.dist_tabs=='Distribution Fit'",
                     selectInput(inputId = "dist.var", label = "Choose a column", choices = colnames("dataset")),
                     checkboxGroupInput("distributions", "Distributions:",
                                        choices = c("Log-normal", "Log-logistic", "Pareto", "Burr", "Weibull", "Gamma"), selected = c("Log-normal", "Pareto")
                     ),
                     radioButtons("dist_zoom", "Zoom to see fit", c("slider", "text input")),
                     conditionalPanel(
                       condition = "input.dist_zoom=='slider'",
                       sliderInput("dist_range", "Range:",
                                   min = 0.1, max = 1000, step = 1,
                                   value = c(0.1, 1000)
                       )
                     ),
                     conditionalPanel(
                       condition = "input.dist_zoom=='text input'",
                       textOutput("dist_range_allowed"),
                       numericInput("dist_range_min", "min", value = 0.1, min = 0.1, max = 1000),
                       numericInput("dist_range_max", "max", value = 1000, min = 0.1, max = 1000)
                     ),
                     downloadButton("downloaddist", "Download as PDF")
                   ),
                   conditionalPanel(
                     condition = "input.dist_tabs=='AIC table'",
                     downloadButton("downloaddistaic", "Download as CSV")
                   )
                 ),
                 mainPanel(
                   h3("Distribution Fit"),
                   tabsetPanel(
                     type = "tabs", id = "dist_tabs",
                     tabPanel(
                       "Distribution Fit",
                       uiOutput("help_text_dis_fit"),
                       plotOutput("dist.plot")),
                     tabPanel(
                       "AIC table",
                       conditionalPanel(
                         condition = "$('html').hasClass('shiny-busy')",
                         div(img(src = "load.gif", width = 240, height = 180),
                             h4("Processing ... Please wait"),
                             style = "text-align: center;"
                         )
                       ),
                       conditionalPanel(
                         condition = "!$('html').hasClass('shiny-busy')",
                         div(tableOutput("dist.aic"), style = "font-size:80%")
                       )
                     )
                   )
                 )
               ),
               tabPanel(
                 "  Correlation  ",
                 sidebarPanel(
                   radioButtons(
                     "cor_method", "Method:",
                     c("Pearson correlation", "Spearman correlation")
                   ),
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
                     downloadButton("downloadcorrmat", "Download as CSV")
                   )
                 ),
                 mainPanel(
                   conditionalPanel(
                     condition = "input.cor_method=='Pearson correlation'",
                     h3("Pearson correlation")
                   ),
                   conditionalPanel(
                     condition = "input.cor_method=='Spearman correlation'",
                     h3("Spearman correlation")
                   ),
                   tabsetPanel(
                     type = "tabs", id = "cor_tabs",
                     tabPanel(
                       "Correlation heatmap",
                       uiOutput("help_text_correlation"),
                       plotOutput("corr.plot")),
                     tabPanel("Correlation plot", plotOutput("corr.plot2")),
                     tabPanel("Correlation matrix", div(tableOutput("corr.matrix"), style = "font-size:80%"))
                   )
                 )
               ),
               tabPanel(
                 "PCA",
                 sidebarPanel(
                   conditionalPanel(
                     condition = "input.pca_tabs == 'PCA-2D plot'",
                     selectInput(inputId = "pca.x", label = "X-axis", choices = ""),
                     selectInput(inputId = "pca.y", label = "Y-axis", choices = "")
                   ),
                   selectInput(inputId = "gene_size", label = "Gene sample size", choices = ""),
                   radioButtons(
                     "gene_order", "Gene sample order (wrt column 1)",
                     c("Descending (highest to lowest)" = "Descending", "Ascending (lowest to highest)" = "Ascending", "Random")
                   ),
                   conditionalPanel(
                     condition = "input.pca_tabs == 'PCA-2D plot' || input.pca_tabs == 'PCA-3D plot'",
                     checkboxInput("pca_cluster", strong("Kmeans clustering on columns"), FALSE),
                     conditionalPanel(
                       condition = "input.pca_cluster == true",
                       sliderInput("pca_cluster_num", "Number of clusters:", value = 1, min = 1, max = 1, step = 1),
                       checkboxInput("pca_text", strong("Display sample name"), FALSE)
                     )
                   ),
                   ######################################
                   radioButtons(
                     "pca_type", "Type of PCA",
                     c("PCA" = "PCA", "Sparse PCA" = "SPCA")
                   ),
                   ######################################
                   conditionalPanel(
                     condition = "input.gene_order=='Random'",
                     helpText("* Click multiple times to resample"),
                     actionButton("pca_refresh", "Resample", style = "background-color: #337ab7;border-color:#337ab7"),
                     br(), br()
                   ),
                   conditionalPanel(
                     condition = "input.pca_tabs == 'PCA variance'",
                     downloadButton("downloadpcavar", "Download as PNG")
                   ),
                   conditionalPanel(
                     condition = "input.pca_tabs == 'PCA-2D plot'",
                     downloadButton("downloadpca2d", "Download as PNG")
                   ),
                   conditionalPanel(
                     condition = "input.pca_tabs == 'PCA-3D plot'",
                     downloadButton("downloadpca3d", "Download as PNG")
                   )
                 ),
                 mainPanel(
                   h3("PCA"),
                   tabsetPanel(
                     type = "tabs", id = "pca_tabs",
                     tabPanel("PCA variance", 
                              uiOutput("help_text_PCA"),
                              plotlyOutput("pcavar.plot")),
                     tabPanel("PCA-2D plot", plotlyOutput("pca2d.plot")),
                     tabPanel("PCA-3D plot", plotlyOutput("pca3d.plot"))
                   )
                 )
               ),
               tabPanel(
                 "DE Analysis",
                 sidebarPanel(
                   radioButtons("n_rep", "Replicates?", choices = c("Multiple" = 1, "Single" = 0)),
                   conditionalPanel(
                     condition = "input.n_rep=='1'",
                     radioButtons("de_method1", "DE Method", choices = c("EdgeR", "DESeq2", "NOISeq"))
                   ),
                   conditionalPanel(
                     condition = "input.n_rep=='0'",
                     radioButtons("de_method0", "DE Method", choices = c("NOISeq"))
                   ),
                   h5("Choose 2 experiment conditions for DE analysis"),
                   selectInput("f1", "Condition 1", choices = ""),
                   selectInput("f2", "Condition 2", choices = ""),
                   
                   h5("DE criteria"),
                   splitLayout(
                     numericInput("p_val", "FDR", min = 0.01, max = 1, value = 0.05, step = 0.01),
                     numericInput("fc", "Fold Change", min = 1, value = 2, step = 0.1)
                   ),
                   fluidRow(
                     column(
                       4,
                       actionButton("submit_DE", "Submit")
                     ),
                     column(
                       6,
                       conditionalPanel(
                         condition = "input.DE_tabs=='DE genes' ",
                         downloadButton("download_de_table", "Download table (csv)")
                       ),
                       conditionalPanel(
                         condition = "input.DE_tabs=='Volcano plot' ",
                         downloadButton("download_volcano", "Download plot (PDF)")
                       ),
                       conditionalPanel(
                         condition = "input.DE_tabs=='Dispersion plot' ",
                         downloadButton("download_dispersion", "Download plot (PDF)")
                       )
                       # conditionalPanel(
                       #   condition = "input.DE_tabs=='Heatmap plot' ",
                       #   downloadButton("download_heatmap","Download plot")
                       # )
                     )
                   )
                 ),
                 mainPanel(
                   h3("DE Analysis"),
                   tabsetPanel(
                     type = "tabs", id = "DE_tabs",
                     tabPanel(
                       "DE genes",
                       uiOutput("help_text_DE_anal"),
                       # h3("Differential Expression Analysis"),
                       conditionalPanel(
                         condition = "$('html').hasClass('shiny-busy')",
                         div(img(src = "load.gif", width = 240, height = 180),
                             h4("Processing ... Please wait"),
                             style = "text-align: center;"
                         )
                       ),
                       conditionalPanel(
                         condition = "!$('html').hasClass('shiny-busy')",
                         DT::dataTableOutput("DE_table")
                       )
                     ),
                     tabPanel(
                       "Volcano plot", # for DESeq and edgeR
                       h6("Volcano plot is only available for edgeR and DESeq2 methods"),
                       conditionalPanel(
                         condition = "input.n_rep=='1' && input.method1!='NOISeq'",
                         conditionalPanel(
                           condition = "$('html').hasClass('shiny-busy')",
                           div(img(src = "load.gif", width = 240, height = 180),
                               h4("Processing ... Please wait"),
                               style = "text-align: center;"
                           )
                         ),
                         conditionalPanel(
                           condition = "!$('html').hasClass('shiny-busy')",
                           plotOutput("volcano_plot")
                         )
                       ),
                       conditionalPanel(
                         condition = "input.method0=='NOISeq' || input.method1=='NOISeq'",
                         h6("Volcano Plot is only applicable to DESeq2 and edgeR")
                       )
                     ),
                     tabPanel(
                       "Dispersion plot", # for edgeR
                       h6("Dispersion plot is only available for edgeR and DESeq2 methods"),
                       conditionalPanel(
                         condition = "input.n_rep=='1' && input.method1!='NOISeq'",
                         # h3("Dispersion plot"),
                         conditionalPanel(
                           condition = "$('html').hasClass('shiny-busy')",
                           div(img(src = "load.gif", width = 240, height = 180),
                               h4("Processing ... Please wait"),
                               style = "text-align: center;"
                           )
                         ),
                         conditionalPanel(
                           condition = "!$('html').hasClass('shiny-busy')",
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
               tabPanel(
                 "Heatmap",
                 sidebarPanel(
                   conditionalPanel(
                     condition = "input.heatmap_tabs=='Heatmap'",
                     
                     radioButtons("heatmap_de_ind", label = "Choose data", choices = c("Indenpendent" = "ind", "DE result" = "de")),
                     numericInput("numOfCluster", "Number of clusters on rows", value = 2, min = 2, max = 30, step = 1),
                     conditionalPanel(
                       condition = "input.heatmap_de_ind == 'ind' ",
                       # selectInput('numOfGeno',"Number of genotypes (mutants)",choices=c(1)),
                       splitLayout(
                         numericInput("fold", "Fold change", value = 2, min = 1, step = 1),
                         numericInput("fold_ncol", "min. column", value = 2, min = 1, step = 1)
                       )
                       
                       # uiOutput("refGeno"),
                       # radioButtons('heatmap_value',"Values",
                       #              c('Fold change','Log fold change'))
                     ),
                     
                     downloadButton("downloadheatmap", "Download as PDF"),
                     actionButton("heatmap_plot", "Plot", width = "65px", style = "color: #fff; background-color: #337ab7; border-color: #337ab7;float:right")
                     
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
                       downloadButton("downloadclusters", "Download as CSV")
                     )
                   )
                 ),
                 mainPanel(
                   h3("Heatmap"),
                   tabsetPanel(
                     type = "tabs", id = "heatmap_tabs",
                     tabPanel(
                       "Heatmap",
                       uiOutput("help_text_heatmap"),
                       conditionalPanel(
                         condition = "$('html').hasClass('shiny-busy')",
                         div(img(src = "load.gif", width = 240, height = 180),
                             h4("Processing ... Please wait"),
                             style = "text-align: center;"
                         )
                       ),
                       conditionalPanel(
                         condition = "!$('html').hasClass('shiny-busy')",
                         plotOutput("heatmap.plot")
                       )
                     ),
                     tabPanel("Gene clusters", dataTableOutput("cluster.info"))
                   )
                 )
               ),
               
               ######## NOISE ######
               #############################################
               tabPanel(
                 "Noise",
                 sidebarPanel(
                   radioButtons("noise_situation", "Select desired noise plot between", choices = c("replicates" = "a", "genotypes (average of replicates)" = "b", "genotypes (no replicate)" = "c")),
                   conditionalPanel(
                     condition = "input.noise_situation=='a' | input.noise_situation=='b' ",
                     textInput("noise_numOfRep", "Number of replicates", value = 1),
                     helpText("* Please order the sample columns in input file properly. Replicates of the same genotype should be put in adjacent columns.")
                   ),
                   conditionalPanel(
                     condition = "input.noise_situation=='b'",
                     uiOutput("noise_anchor_choices")
                   ),
                   conditionalPanel(
                     condition = "input.noise_situation=='c'",
                     selectInput("noise_anchor_c", "Anchor genotype", choices = "")
                   ),
                   radioButtons(
                     "noise_graph_type", "Graph type:",
                     c("Bar chart", "Line chart")
                   ),
                   downloadButton("downloadnoise", "Download as PNG"),
                   actionButton("noise_plot", "Plot", width = "65px", style = "color: #fff; background-color: #337ab7; border-color:#337ab7;float:right"),
                   conditionalPanel(
                     condition = "input.noise_situation=='a' | input.noise_situation=='b' ",
                     h5("Specify names of the genotypes"),
                     uiOutput("expand_genonames_noise")
                   )
                 ),
                 mainPanel(
                   h3("Noise"),
                   uiOutput("help_text_Noise"),
                   conditionalPanel(
                     condition = "$('html').hasClass('shiny-busy')",
                     div(img(src = "load.gif", width = 240, height = 180), h4("Processing ... Please wait"), style = "text-align: center;")
                   ),
                   conditionalPanel(
                     condition = "!$('html').hasClass('shiny-busy')",
                     plotlyOutput("noise.plot")
                   )
                 )
               ),
               
               
               ###### ENTROPY #############
               #########################################
               tabPanel(
                 "Entropy",
                 sidebarPanel(
                   checkboxInput("tsflag", strong("Time series data"), FALSE),
                   conditionalPanel(
                     condition = "input.tsflag==true",
                     textInput("entropy_timepoints", "Number of time points"),
                     helpText("* Please order the sample columns in input file properly. Time series data of the same genotype should be put in adjacent columns.")
                   ),
                   radioButtons(
                     "entropy_graph_type", "Graph type:",
                     c("Bar chart", "Line chart")
                   ),
                   downloadButton("downloadentropy", "Download as PNG"),
                   conditionalPanel(
                     condition = "input.tsflag==true",
                     h5("Specify names of the genotypes"),
                     uiOutput("expand_genonames_entropy")
                   )
                 ),
                 mainPanel(
                   h3("Shannon entropy"),
                   uiOutput("help_text_Entropy"),
                   plotlyOutput("entropy.plot")
                 )
               ),
               
               ################### Support Vector Machine ##################### 
               
               # tabPanel('SVM',
               #          sidebarPanel(
               #            actionButton("submit_svm","Submit")),
               #          mainPanel(
               #            h3('SVM Plot'),
               #            tabsetPanel(
               #              type = "tabs", id = "SVM_tabs",
               #              tabPanel("SVM Classification",
               #                       uiOutput("help_text_SVM"),
               #                       plotOutput('svm_plot')),
               #              tabPanel("Raw Plot", plotOutput("svm_df_plot"))
               #            )
               #          )),
               
               #####################################################################
               
               tabPanel(
                 't-SNE',
                 sidebarPanel(
                   splitLayout(
                     numericInput("perplexity_value","Perplexity value", min=1, value=2),
                     numericInput("no_of_pca","No. of PCs", min=1, value=30)
                     # numericInput("no_of_clusters","No. of clusters", min=2, value=2)
                   ),
                   radioButtons('tsne2_trans',"Transformation:",
                                c('None', 'log10')),
                   checkboxInput("tsne_cluster", strong("Kmeans clustering on columns"), FALSE),
                   conditionalPanel(
                     condition = "input.tsne_cluster == true",
                     numericInput("tsne_cluster_num", "Number of clusters:", min = 1, value = 2)
                   ),
                   checkboxInput("tsne_text", strong("Display sample name"), FALSE),
                   actionButton("submit_tsne2","Submit"),
                   
                   conditionalPanel(
                     condition = "input.tsne_tabs=='t-SNE table'",
                     downloadButton("download_tsne", "Download as CSV")
                   )
                   
                 ),
                 mainPanel(
                   h3('t-SNE Plot'),
                   
                   tabsetPanel(
                     type = "tabs", id = "tsne_tabs",
                     tabPanel("t-SNE plot", 
                              uiOutput("help_text_tsne"),
                              plotlyOutput('tsne2.plot')),
                     tabPanel("t-SNE table", 
                              DT::dataTableOutput("tsne_table") )
                   )
                 )),
               
               tabPanel(
                 "Random Forest",
                 sidebarPanel(
                   radioButtons(
                     "analysis_type", "Choose Analysis Type",
                     c("RF clustering" = "rf", "RAFSIL" = "rafsil")
                   ),
                   conditionalPanel(
                     condition = "input.analysis_type=='rf'",  #rf
                     splitLayout(
                       numericInput("num_trees", "No. of trees", min = 1, value = 25),
                       numericInput("num_clusters", "No. of clusters", min = 1, value = 2)
                     ),
                     radioButtons(
                       "rf_trans", "Transformation:",
                       c("None", "log10")
                     ),
                     actionButton("submit_rf", "Submit")
                   ),
                   conditionalPanel(
                     condition = "input.analysis_type=='rafsil'",  #rafsil
                     actionButton("submit_rafsil", "Submit")
                   )
                   # conditionalPanel(
                   #          condition = "input.rf_tabs == 'RF plot'",
                   #          downloadButton("downloadrfplot", "Download as PDF")
                   #        ),
                   #        conditionalPanel(
                   #          condition = "input.rf_tabs == 'RF matrix'",
                   #          downloadButton("downloadrfmatrix", "Download as PDF")
                   #        )
                 ),
                 mainPanel(
                   h3("Clustering With Random Forest"),
                   tabsetPanel(type = "tabs",id="rf_tabs",
                               tabPanel("RF plot", 
                                        uiOutput("help_text_rf"),
                                        plotlyOutput("rf.plot")),
                               tabPanel("RAFSIL plot", plotOutput("RAFSIL.plot")),
                               tabPanel("RF matrix", div(tableOutput('rf.matrix'), style = "font-size:80%"))
                   )
                 )
               ),
               
               tabPanel(
                 "SOM",
                 sidebarPanel(
                   selectInput(inputId = "som_samples", label = "Samples used", choices = ""),
                   splitLayout(
                     numericInput("som_grid_h", "No. of horizontal grids", min = 1, value = 2),
                     numericInput("som_grid_v", "No. of vertical grids", min = 1, value = 2)
                   ),
                   numericInput("som_cluster_size", "No. of clusters (for cluster plot)", min = 2, value = 2),
                   radioButtons(
                     "som_trans", "Transformation:",
                     c("None", "log10")
                   ),
                   actionButton("submit_som", "Submit"),
                   br(),
                   br(),
                   conditionalPanel(
                     condition = "input.som_tabs == 'Property plot'",
                     downloadButton("downloadProperty", "Download as PDF")
                   ),
                   conditionalPanel(
                     condition = "input.som_tabs == 'Count plot'",
                     downloadButton("downloadCount", "Download as PDF")
                   ),
                   conditionalPanel(
                     condition = "input.som_tabs == 'Codes plot'",
                     downloadButton("downloadCodes","Download as PDF")
                   ),
                   conditionalPanel(
                     condition = "input.som_tabs == 'Distance plot'",
                     downloadButton("downloadDistance", "Download as PDF")
                   ),
                   conditionalPanel(
                     condition = "input.som_tabs == 'Cluster plot'",
                     downloadButton("downloadCluster","Download as PDF")
                   )
                 ),
                 mainPanel(
                   h3("SOM Analysis"),
                   tabsetPanel(
                     type = "tabs", id = "som_tabs",
                     tabPanel("Property plot", 
                              uiOutput("help_text_SOM"),
                              plotOutput("som_property.plot")),
                     tabPanel("Count plot", plotOutput("som_count.plot")),
                     tabPanel("Codes plot", plotOutput("som_codes.plot")),
                     tabPanel("Distance plot", plotOutput("som_dist.plot")),
                     tabPanel("Cluster plot", plotOutput("som_cluster.plot"))
                   )
                 ),
                 tags$style(type = 'text/css', 
                            '.navbar { font-size: 17px;}'
                 )
               )
    ),
    ###############################################
    ###############################################
    ###############################################
    navbarMenu(
      "Gene Set Analysis",
      
      
      ########## Pathway Enrichment ##############
      #########################################
      tabPanel(
        "Gene Pathways Enrichment",
        tags$head(tags$style("#path_enri_visu{height:95vh !important;}")),
        sidebarLayout(
          sidebarPanel(
              withTags({
                div(class = "header",
                    p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/online-version/Test%20data/gPro_gene_names.csv")),
                )
              }),
            fileInput("file_path_enri_gene", "Upload genes CSV file"),
            textInput("text_path_enri_gene", "Enter gene id"),
            actionButton("submit_path_enri_gene", "Submit"),br(),br(),
            selectInput("loadStyleFile_path_gene", "Select Style: ", choices=styles),
            # selectInput(inputId = "overlap_min", label = "Minimum Overlap", choices = ""),
            sliderInput("overlap_min_gene", "Minimum Overlap",
                        min = 0, max = 100,
                        value = 15),
            selectInput("doLayout_path_gene", "Select Layout:",
                        choices=c("",
                                  "cose",
                                  "cola",
                                  "circle",
                                  "concentric",
                                  "breadthfirst",
                                  "grid",
                                  "random",
                                  "dagre",
                                  "cose-bilkent")),
            actionButton("sfn_path_gene", "Select First Neighbor"),
            br(),br(),
            actionButton("fit_path_gene", "Fit Graph"),br(),br(),
            actionButton("fitSelected_path_gene", "Fit Selected"),br(),br(),
            actionButton("clearSelection_path_gene", "Clear Selection"), br(),br(),
            actionButton("removeGraphButton_path_gene", "Remove Graph"), br(),br(),
            actionButton("addRandomGraphFromDataFramesButton_path_gene", "Add Random Graph"),br(),br(),
            actionButton("getSelectedNodes_path_gene", "Get Selected Nodes"), br(),br(),
            htmlOutput("selectedNodesDisplay_path_gene"),
           
          ),
          mainPanel(
            h3("Pathways Enrichment"),
            tabsetPanel(
              type = "tabs", id = "path_enri_tab_gene",
              tabPanel("Plot",
                       uiOutput("help_text_path_enri_gene"),
                       conditionalPanel(
                         condition = "$('html').hasClass('shiny-busy')",
                         div(img(src = "load.gif", width = 240, height = 180),
                             h4("Processing ... Please wait"),
                             style = "text-align: center;"
                         )
                       ),
                       conditionalPanel(
                         condition = "!$('html').hasClass('shiny-busy')",
                         plotlyOutput("path_enri.plot_gene")
                       ), 
              ),
              tabPanel(
                "Visualization",
                conditionalPanel(
                  condition = "$('html').hasClass('shiny-busy')",
                  div(img(src = "load.gif", width = 240, height = 180),
                      h4("Processing ... Please wait"),
                      style = "text-align: center;"
                  )
                ),
                conditionalPanel(
                  condition = "!$('html').hasClass('shiny-busy')",
                  cyjShinyOutput('path_enri_visu_gene', height=350)
                ),
              )
            )
          )
        )),
      
      
      ###### Tissue Expression #############
      #########################################
      tabPanel(
        "Tissue Expression",
        sidebarPanel(
          withTags({
            div(class = "header",
                p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/online-version/Test%20data/gene_id.csv")),
            )
          }),
          fileInput("file_prot_expr", "Upload UniProt accession CSV file"),
          textInput("text_prot_expr","Enter Uniprot accession numbers"),
          actionButton("submit_prot_expr", "Submit"),br(),br(),
          downloadButton("prot_expr_download", "Download as CSV")
        ),
        mainPanel(
          h3("Tissue Expression"),
          uiOutput("help_text_prot_exp"),
          conditionalPanel(
            condition = "$('html').hasClass('shiny-busy')",
            div(img(src = "load.gif", width = 240, height = 180),
                h4("Processing ... Please wait"),
                style = "text-align: center;"
            )
          ),
          conditionalPanel(
            condition = "!$('html').hasClass('shiny-busy')",
            DT::dataTableOutput("prot_expr_table")
          ), 
        )),
      
      
      ########## Co-expression #############
      #########################################
      tabPanel(
        "Co-expression",
        sidebarPanel(
          selectInput("organismID", "Choose organism:", 
                      choices= list(
                        "Arabidopsis_thaliana",
                        "Caenorhabditis_elegans",
                        "Danio_rerio",
                        "Drosophila_melanogaster",
                        "Escherichia_coli",
                        "Homo_sapiens",
                        "Mus_musculus",
                        "Rattus_norvegicus",
                        "Saccharomyces_cerevisiae"
                      )),
          withTags({
            div(class = "header",
                p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/online-version/Test%20data/gene_names.csv")),
            )
          }),
          fileInput("file_gene", "Upload genes CSV file"),
          textInput("text_gene","Enter gene ids"),
          actionButton("genemania_submit", "Submit")
        ),
        mainPanel(
          h3("Co-expression"),
          uiOutput("help_text_gene_mania"),
          conditionalPanel(
            condition = "$('html').hasClass('shiny-busy')",
            div(img(src = "load.gif", width = 240, height = 180),
                h4("Processing ... Please wait"),
                style = "text-align: center;"
            )
          ),
          conditionalPanel(
            condition = "!$('html').hasClass('shiny-busy')",
            div(id = "hide_link", 
                p("Please click", htmlOutput("linkCo")))
            %>% shinyjs::hidden()
          ),
        ))
      #########################################
      
      #########################################
      
    ),
    navbarMenu(
      "Protein Set Analysis",
      
      ########## Gene Ontology #############
      #########################################
      tabPanel(
        "Gene ontology",
        sidebarPanel(
          withTags({
            div(class = "header",
                p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/online-version/Test%20data/gene_id.csv")),
            )
          }),
          fileInput("file_uniprot", "Upload UniProt accession CSV file"),
          textInput("text_uniprot", "Enter UniProt accession csv file"),
          actionButton("submit_uniprot", "Submit"),br(),br(),
          conditionalPanel(
            condition = "input.uniprot_tabs == 'Biological process'",
            downloadButton("download_bio_plot", "Download Plot"),br(),br(),
            downloadButton("download_bio_pro", "Download Table")
          ),
          conditionalPanel(
            condition = "input.uniprot_tabs == 'Molecular function'",
            downloadButton("download_mole_plot", "Download Plot"),br(),br(),
            downloadButton("download_mole_func", "Download Table")
          ),
          conditionalPanel(
            condition = "input.uniprot_tabs == 'Cellular component'",
            downloadButton("download_cell_plot", "Download Plot"),br(),br(),
            downloadButton("download_cell_comp", "Download Table")
          )
        ),
        mainPanel(
          h3("Gene ontology"),
          tabsetPanel(
            type = "tabs", id = "uniprot_tabs",
            tabPanel("Biological process",
                     uiOutput("help_text_bio_pr"),
                     conditionalPanel(
                       condition = "$('html').hasClass('shiny-busy')",
                       div(img(src = "load.gif", width = 240, height = 180),
                           h4("Processing ... Please wait"),
                           style = "text-align: center;"
                       )
                     ),
                     conditionalPanel(
                       condition = "!$('html').hasClass('shiny-busy')",
                       plotOutput("uniprotbioplot"),
                       shiny::dataTableOutput("uniprot_biotable")),
            ),
            tabPanel("Molecular function",
                     conditionalPanel(
                       condition = "$('html').hasClass('shiny-busy')",
                       div(img(src = "load.gif", width = 240, height = 180),
                           h4("Processing ... Please wait"),
                           style = "text-align: center;"
                       )
                     ),
                     conditionalPanel(
                       condition = "!$('html').hasClass('shiny-busy')",
                       plotOutput("uniprot_molcplot"),
                       shiny::dataTableOutput("uniprot_molctable")
                     ),
            ),
            tabPanel("Cellular component",
                     conditionalPanel(
                       condition = "$('html').hasClass('shiny-busy')",
                       div(img(src = "load.gif", width = 240, height = 180),
                           h4("Processing ... Please wait"),
                           style = "text-align: center;"
                       )
                     ),
                     conditionalPanel(
                       condition = "!$('html').hasClass('shiny-busy')",
                       plotOutput("uniprot_celplot"),
                       shiny::dataTableOutput("uniprot_celtable")
                     ),
            )
          )
        )), ##### Gene ontlogy closing 
      
      
      ###### Protein Interaction #############
      #########################################
      tabPanel(
        "P-P Interactions",
        tags$head(tags$style("#cyjShiny{height:95vh !important;}")),
        sidebarLayout(
          sidebarPanel(
            withTags({
              div(class = "header",
                  p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/online-version/Test%20data/gene_id.csv")),
              )
            }),
            fileInput("file_prot_Int", "Upload UniProt accession CSV file"),
            textInput("text_prot_Int","Enter UniProt accession numbers"),
            actionButton("submit_prot_Int", "Submit"),br(),br(),
            selectInput("loadStyleFile", "Select Style: ", choices=styles),
            selectInput("doLayout", "Select Layout:",
                        choices=c("",
                                  "cose",
                                  "cola",
                                  "circle",
                                  "concentric",
                                  "breadthfirst",
                                  "grid",
                                  "random",
                                  "dagre",
                                  "cose-bilkent")),
            # selectInput("showCondition", "Select Condition:", choices=rownames(output$tbl.lfc)),
            # selectInput("selectName", "Select Node by ID:", choices = c("", sort(tbl.nodes$id))),
            actionButton("sfn", "Select First Neighbor"),
            br(),br(),
            actionButton("fit", "Fit Graph"),br(),br(),
            actionButton("fitSelected", "Fit Selected"),br(),br(),
            actionButton("clearSelection", "Clear Selection"), br(),br(),
            actionButton("removeGraphButton", "Remove Graph"), br(),br(),
            actionButton("addRandomGraphFromDataFramesButton", "Add Random Graph"),br(),br(),
            actionButton("getSelectedNodes", "Get Selected Nodes"), br(),br(),
            htmlOutput("selectedNodesDisplay"),
          ),
          mainPanel(
            h3("Protein-Protein Interactions"),
            tabsetPanel(
              type = "tabs", id = "prot_inte_tab",
              tabPanel("Visualization",
                       uiOutput("help_text_p_inte"),
                       conditionalPanel(
                         condition = "$('html').hasClass('shiny-busy')",
                         div(img(src = "load.gif", width = 240, height = 180),
                             h4("Processing ... Please wait"),
                             style = "text-align: center;"
                         )
                       ),
                       conditionalPanel(
                         condition = "!$('html').hasClass('shiny-busy')",
                         cyjShinyOutput('cyjShiny', height=350)
                       )),
              tabPanel("Protein Interaction",
                       conditionalPanel(
                         condition = "$('html').hasClass('shiny-busy')",
                         div(img(src = "load.gif", width = 240, height = 180),
                             h4("Processing ... Please wait"),
                             style = "text-align: center;"
                         )
                       ),
                       conditionalPanel(
                         condition = "!$('html').hasClass('shiny-busy')",
                         DT::dataTableOutput("prot_int_table")
                       )),
              tabPanel("Protein Name",
                       conditionalPanel(
                         condition = "$('html').hasClass('shiny-busy')",
                         div(img(src = "load.gif", width = 240, height = 180),
                             h4("Processing ... Please wait"),
                             style = "text-align: center;"
                         )
                       ),
                       conditionalPanel(
                         condition = "!$('html').hasClass('shiny-busy')",
                         DT::dataTableOutput("prot_name_table")
                       ))
            ),
           
          )
          #  mainPanel(cyjShinyOutput('cyjShiny', height=400),
          #           width=10,
          #           tabPanel("Protein Name",
          #           DT::dataTableOutput("prot_name_table"))
          #  )
        )),   #### PPI Closing
      
      
      ###### Protein Function #############
      #########################################
      tabPanel(
        "Protein Function",
        sidebarPanel(
          withTags({
            div(class = "header",
                p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/online-version/Test%20data/gene_id.csv")),
            )
          }),
          fileInput("file_prot_func", "Upload UniProt accession CSV file"),
          textInput("text_prot_func","Enter UniProt accession numbers"),
          actionButton("submit_prot_func", "Submit"),br(),br(),
          downloadButton("prot_func_download", "Download as CSV")
        ),
        mainPanel(
          h3("Protein Function"),
          uiOutput("help_text_prot_fn"),
          conditionalPanel(
            condition = "$('html').hasClass('shiny-busy')",
            div(img(src = "load.gif", width = 240, height = 180),
                h4("Processing ... Please wait"),
                style = "text-align: center;"
            )
          ),
          conditionalPanel(
            condition = "!$('html').hasClass('shiny-busy')",
            DT::dataTableOutput("prot_func_table")
          ),
        )),  ### Protein function closing 
      
      
      ###### Subcellular Localization #############
      #########################################
      tabPanel(
        "Subcellular Localization",
        sidebarPanel(
          withTags({
            div(class = "header",
                p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/online-version/Test%20data/gene_id.csv")),
            )
          }),
          fileInput("file_prot_local", "Upload UniProt accession CSV file"),
          textInput("text_prot_local","Enter UniProt accession numbers"),
          actionButton("submit_prot_local", "Submit"),br(),br(),
          downloadButton("prot_local_download", "Download as CSV")
        ),
        mainPanel(
          h3("Subcellular Localization"),
          uiOutput("help_text_sub_loc"),
          conditionalPanel(
            condition = "$('html').hasClass('shiny-busy')",
            div(img(src = "load.gif", width = 240, height = 180),
                h4("Processing ... Please wait"),
                style = "text-align: center;"
            )
          ),
          conditionalPanel(
            condition = "!$('html').hasClass('shiny-busy')",
            DT::dataTableOutput("prot_local_table")
          ), 
        )),  ### Subcellular closing 
      
      ########## Protein Domains ##############
      #########################################
      tabPanel(
        "Protein Domains",
        sidebarPanel(
          withTags({
            div(class = "header",
                p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/online-version/Test%20data/gene_id.csv")),
            )
          }),
          fileInput("file_prot_domain", "Upload UniProt accession CSV file"),
          textInput("text_prot_domain","Enter UniProt accession numbers"),
          actionButton("submit_prot_domain", "Submit"),br(),br(),
          downloadButton("prot_domain_download", "Download as CSV")
        ),
        mainPanel(
          h3("Protein Domains"),
          uiOutput("help_text_pro_dom"),
          conditionalPanel(
            condition = "$('html').hasClass('shiny-busy')",
            div(img(src = "load.gif", width = 240, height = 180),
                h4("Processing ... Please wait"),
                style = "text-align: center;"
            )
          ),
          conditionalPanel(
            condition = "!$('html').hasClass('shiny-busy')",
            DT::dataTableOutput("prot_domain_table")
          ), 
        )),   ### Protein domain closing 
      
      
      ######### Protein Sequences #############
      #########################################
      tabPanel(
        "Protein properties",
        sidebarPanel(
          fileInput("file_prot_seq", "Upload UniProt accession CSV file"),
          textInput("text_prot_seq","Enter UniProt accession numbers"),
          actionButton("submit_prot_Seq", "Submit"),br(),br(),
          shinyjs::hidden(downloadButton('downloadData', 'Download Sequence FASTA')),
          
        ),
        mainPanel(
          h3("Protein Sequences"),
          tabsetPanel(type = "tabs",
                      tabPanel("Sequence information",uiOutput("help_text_prot_seq")),
                      tabPanel(
                        "Sequence charge", 
                        conditionalPanel(
                          condition = "$('html').hasClass('shiny-busy')",
                          div(img(src = "load.gif", width = 240, height = 180),
                              h4("Processing ... Please wait"),
                              style = "text-align: center;"
                          )
                        ),
                        conditionalPanel(
                          condition = "!$('html').hasClass('shiny-busy')",
                          plotOutput("ChargePlot") #  , width="750px",height="750px"
                        ),
                      ),
                      tabPanel(
                        "Sequence acidity", 
                        conditionalPanel(
                          condition = "$('html').hasClass('shiny-busy')",
                          div(img(src = "load.gif", width = 240, height = 180),
                              h4("Processing ... Please wait"),
                              style = "text-align: center;"
                          )
                        ),
                        conditionalPanel(
                          condition = "!$('html').hasClass('shiny-busy')",
                          plotOutput("AcidityPlot") #  , width="750px",height="750px"
                        ),
                      ),
                      tabPanel(
                        "Sequence gravy index",
                        conditionalPanel(
                          condition = "$('html').hasClass('shiny-busy')",
                          div(img(src = "load.gif", width = 240, height = 180),
                              h4("Processing ... Please wait"),
                              style = "text-align: center;"
                          )
                        ),
                        conditionalPanel(
                          condition = "!$('html').hasClass('shiny-busy')",
                          plotOutput("GravyPlot") #  , width="750px",height="750px"
                        ),
                      ),
                      tabPanel(
                        "All physicochemical properties", 
                        conditionalPanel(
                          condition = "$('html').hasClass('shiny-busy')",
                          div(img(src = "load.gif", width = 240, height = 180),
                              h4("Processing ... Please wait"),
                              style = "text-align: center;"
                          )
                        ),
                        conditionalPanel(
                          condition = "!$('html').hasClass('shiny-busy')",
                          plotOutput("SequencePlot" , width="900px",height="750px")
                        ),
                      )
                      #################################################################################
                      
          ),
          
        )), # For protein properties tab panel 
      
      #Here we start evolution tab panel 
      tabPanel(
        "Evolutionary Analysis",
        sidebarPanel(
          fileInput("file_prot_seq_evol", "Upload UniProt accession CSV file"),
          textInput("text_prot_seq_evol","Enter UniProt accession numbers"),
          actionButton("submit_prot_seq_evol","Submit")
          
        ),
        mainPanel(
          h3("Protein Evolutionary analysis"),
          tabsetPanel(type = "tabs",
                      tabPanel("Evolutionary information",uiOutput("help_text_prot_seq_evol")),
                      tabPanel(
                        "Protein's genes tree", 
                        conditionalPanel(
                          condition = "$('html').hasClass('shiny-busy')",
                          div(img(src = "load.gif", width = 240, height = 180),
                              h4("Processing ... Please wait"),
                              style = "text-align: center;"
                          )
                        ),
                        conditionalPanel(
                          condition = "!$('html').hasClass('shiny-busy')",
                          radialNetworkOutput("GenePlot", width="900px",height="750px") #  , width="750px",height="750px"
                        ),
                      ),
                      tabPanel(
                        "Protein's chromosomal location", 
                        conditionalPanel(
                          condition = "$('html').hasClass('shiny-busy')",
                          div(img(src = "load.gif", width = 240, height = 180),
                              h4("Processing ... Please wait"),
                              style = "text-align: center;"
                          )
                        ),
                        conditionalPanel(
                          condition = "!$('html').hasClass('shiny-busy')",
                          plotOutput("Chromo", width="1200px",height="1500px") #  , width="750px",height="750px"
                        ),
                      ),
                      tabPanel(
                        "Evolutionary analysis", 
                        conditionalPanel(
                          condition = "$('html').hasClass('shiny-busy')",
                          div(img(src = "load.gif", width = 240, height = 180),
                              h4("Processing ... Please wait"),
                              style = "text-align: center;"
                          )
                        ),
                        conditionalPanel(
                          condition = "!$('html').hasClass('shiny-busy')",
                          plotOutput("Phylogenetic", width="900px",height="750px") # , width="900px",height="750px"
                        ),
                      )
                      
          ),
          
        )), #For protein evolution tab panel 
      
      tabPanel(
        "Pathological Analysis",
        sidebarPanel(
          fileInput("file_prot_seq_Patho", "Upload UniProt accession CSV file"),
          textInput("text_prot_seq_Patho","Enter UniProt accession numbers"),
          actionButton("submit_prot_seq_Patho","Submit")
        ),
        mainPanel(
          h3("Protein pathological analysis"),
          tabsetPanel(type = "tabs",
                      tabPanel("Pathological information",uiOutput("help_text_prot_seq_Patho")),
                      tabPanel(
                        "Protein's disease role", 
                        conditionalPanel(
                          condition = "$('html').hasClass('shiny-busy')",
                          div(img(src = "load.gif", width = 240, height = 180),
                              h4("Processing ... Please wait"),
                              style = "text-align: center;"
                          )
                        ),
                        conditionalPanel(
                          condition = "!$('html').hasClass('shiny-busy')",
                          dataTableOutput("DisaeseTable") #  , width="750px",height="750px"
                        ),
                      ),
                      tabPanel(
                        "Protein's disease distribution", 
                        conditionalPanel(
                          condition = "$('html').hasClass('shiny-busy')",
                          div(img(src = "load.gif", width = 240, height = 180),
                              h4("Processing ... Please wait"),
                              style = "text-align: center;"
                          )
                        ),
                        conditionalPanel(
                          condition = "!$('html').hasClass('shiny-busy')",
                          bubblesOutput("DiseasePlot") #  , width="750px",height="750px"
                        ),
                      )
                      
          ),
          
        )), #For protein pathology tab panel 
      
      
      #########################################
      tabPanel(
        "Protein Complex Enrichment",
        sidebarPanel(
          withTags({
            div(class = "header",
                p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/online-version/Test%20data/gene_id.csv")),
            )
          }),
          fileInput("file_complex_prot", "Upload UniProt accession CSV file"),
          textInput("text_complex_prot","Enter UniProt accession numbers"),
          actionButton("submit_complex_prot", "Submit"),br(),br(),
          downloadButton("complex_download_prot", "Download as CSV")
        ),
        mainPanel(
          h3("Complex Enrichment"),
          uiOutput("help_text_complex_en_prot"),
          conditionalPanel(
            condition = "$('html').hasClass('shiny-busy')",
            div(img(src = "load.gif", width = 240, height = 180),
                h4("Processing ... Please wait"),
                style = "text-align: center;"
            )
          ),
          conditionalPanel(
            condition = "!$('html').hasClass('shiny-busy')",
            shiny::dataTableOutput("complex_table_prot")
          ),
        )), ###Complex enrichment
      
      
      ########## Pathway Enrichment ##############
      #########################################
      tabPanel(
        "Pathways Enrichment",
        tags$head(tags$style("#path_enri_visu{height:95vh !important;}")),
        sidebarLayout(
          sidebarPanel(
            withTags({
              div(class = "header",
                  p("Example ", a("here", href = "https://github.com/buithuytien/GeneCloudOmics/blob/online-version/Test%20data/gPro_gene_names.csv")),
              )
            }),
            fileInput("file_path_enri_prot", "Upload genes CSV file"),
            textInput("text_path_enri_prot", "Enter gene id"),
            actionButton("submit_path_enri_prot", "Submit"),br(),br(),
            selectInput("loadStyleFile_path_prot", "Select Style: ", choices=styles),
            # selectInput(inputId = "overlap_min", label = "Minimum Overlap", choices = ""),
            sliderInput("overlap_min_prot", "Minimum Overlap",
                        min = 0, max = 100,
                        value = 15),
            selectInput("doLayout_path_prot", "Select Layout:",
                        choices=c("",
                                  "cose",
                                  "cola",
                                  "circle",
                                  "concentric",
                                  "breadthfirst",
                                  "grid",
                                  "random",
                                  "dagre",
                                  "cose-bilkent")),
            actionButton("sfn_path_prot", "Select First Neighbor"),
            br(),br(),
            actionButton("fit_path_prot", "Fit Graph"),br(),br(),
            actionButton("fitSelected_path_prot", "Fit Selected"),br(),br(),
            actionButton("clearSelection_path_prot", "Clear Selection"), br(),br(),
            actionButton("removeGraphButton_path_prot", "Remove Graph"), br(),br(),
            actionButton("addRandomGraphFromDataFramesButton_path_prot", "Add Random Graph"),br(),br(),
            actionButton("getSelectedNodes_path_prot", "Get Selected Nodes"), br(),br(),
            htmlOutput("selectedNodesDisplay_path_prot")
          ),
          mainPanel(
            h3("Pathways Enrichment"),
            tabsetPanel(
              type = "tabs", id = "path_enri_tab_prot",
              tabPanel("Plot",
                       uiOutput("help_text_path_enri_prot"),
                       conditionalPanel(
                         condition = "$('html').hasClass('shiny-busy')",
                         div(img(src = "load.gif", width = 240, height = 180),
                             h4("Processing ... Please wait"),
                             style = "text-align: center;"
                         )
                       ),
                       conditionalPanel(
                         condition = "!$('html').hasClass('shiny-busy')",
                         plotlyOutput("path_enri.plot_prot")
                       ), 
              ),
              tabPanel(
                "Visualization",
                conditionalPanel(
                  condition = "$('html').hasClass('shiny-busy')",
                  div(img(src = "load.gif", width = 240, height = 180),
                      h4("Processing ... Please wait"),
                      style = "text-align: center;"
                  )
                ),
                conditionalPanel(
                  condition = "!$('html').hasClass('shiny-busy')",
                  cyjShinyOutput('path_enri_visu_prot', height=350)
                ),
              )
            )
          )
        )) #Pathway enrichemnt analysis 
      
    ) ####Nav bar closing 
  )
)
####################################################

server <- function(input, output, session) {
  
  
  gene_mania_link <- reactiveVal("https://genemania.org")
  count_fasta <- reactiveVal(0)
  count_id <- reactiveVal(0)
  
  ########################################
  ##### Increases the Upload Limit #######
  ########################################
  
  options(shiny.maxRequestSize=250*1024^2)
  
  ########################################
  ##### get variable names for input #####
  ########################################
  
  observe({
    type <- input$file_type
    if (type == "norm") {
      DS <- df_norm()
    } else if (type == "raw") {
      DS <- df_raw()
    }
    nms <- colnames(DS)
    updateSelectInput(session, "scatter.x", choices = nms, selected = nms[1])
    updateSelectInput(session, "scatter.y", choices = nms, selected = nms[2])
    updateSelectInput(session, "dist.var", choices = nms)
    col_num <- ncol(DS)
    updateSliderInput(session, "pca_cluster_num", max = col_num - 1)
    genotype_num <- NULL
    if (is.null(DS) == FALSE) {
      for (i in 2:col_num) {
        if (col_num %% i == 0) {
          genotype_num <- c(genotype_num, i)
        }
      }
    }
    updateSelectInput(session, "numOfGeno", choices = genotype_num)
    updateSelectInput(session, "noise_anchor_c", choices = nms)
    
    ##############################################
    ##############################################
    updateSelectInput(session, "som_samples", choices = c("All", nms), selected = "All")
    updateSelectInput(session, "overlay.x1", choices = nms, selected = nms[1])
    updateSelectInput(session, "overlay.y1", choices = nms, selected = nms[2])
    updateSelectInput(session, "overlay.x2", choices = nms, selected = nms[3])
    updateSelectInput(session, "overlay.y2", choices = nms, selected = nms[4])
    ##############################################
    ##############################################
    
    observeEvent(input$start_rnaseq, {
      updateNavbarPage(session, inputId = "navbar", selected = "active_tab_rnaseq")
    })
    
    observeEvent(input$start_micro, {
      updateNavbarPage(session, inputId = "navbar", selected = "active_tab_micro")
    })
    
    ### preprocessing tab
    f <- group_names()
    f <- unique(as.character(f))
    if (is.null(f)) {
      hideTab(inputId = "preprocessing_tabs", target = "Description table")
      # hideTab(inputId="preprocessing_tabs", target="Description table")
    } else {
      showTab(inputId = "preprocessing_tabs", target = "Description table")
      updateSelectInput(session, "f1", choices = f, selected = f[1])
      updateSelectInput(session, "f2", choices = f, selected = f[2])
    }
    
    ### gene expression range for distribution fit ###
    if (is.null(DS) == FALSE) {
      DS_dist <- distfit_df()
      range_min <- min(DS_dist)
      range_max <- max(DS_dist)
      updateSliderInput(session, "dist_range", max = round(range_max), value = c(0.1, range_max))
      updateNumericInput(session, "dist_range_min", min = 0.000001, max = round(range_max), value = 0.1)
      updateNumericInput(session, "dist_range_max", min = 0.000001, max = round(range_max), value = round(range_max))
    }
    
    ### gene sample size choices for PCA ###
    # print("line 647 check input$submit_preprocessing")
    # v=input$submit_preprocessing
    if (input$submit_preprocessing > 0) {
      if (type == "norm") {
        DS_filt <- df_shiny()
      } else if (type == "raw") {
        DS_filt <- df_raw_shiny()
      }
    } else {
      DS_filt <- DS
    }
    
    i <- 1
    min_size <- 25
    samplesize <- NULL
    while (i * min_size < length(DS_filt[, 1])) {
      samplesize <- c(samplesize, i * min_size)
      i <- i * 2
    }
    if (is.null(samplesize)) {
      samplesize <- c(samplesize, length(DS_filt[, 1]))
    } else if (samplesize[length(samplesize)] != length(DS_filt[, 1])) {
      samplesize <- c(samplesize, length(DS_filt[, 1]))
    }
    updateSelectInput(session, "gene_size", choices = samplesize, selected = samplesize[length(samplesize)])
    
    ### pca choices for PCA-2D ###
    pcchoices <- NULL
    if (is.null(DS) == FALSE) {
      for (i in 1:ncol(DS)) {
        pcchoices <- c(pcchoices, paste("PC", i, sep = ""))
      }
    }
    updateSelectInput(session, "pca.x", choices = pcchoices, selected = pcchoices[1])
    updateSelectInput(session, "pca.y", choices = pcchoices, selected = pcchoices[2])
    
  })
  
  
  observeEvent(input$submit_input, {
    type <- input$file_type
    if (type == "norm") {
      DS <- df_norm()
      lengths <- 0
    } else if (type == "raw") {
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
    names(input.null) <- c("Expression/Counts", "Meta Data")
    
    if (any(input.null)) {
      index.null <- which(input.null)
      errors <- paste(names(input.null)[index.null], collapse = ", ")
      # print(errors)
      showModal(modalDialog(
        type = "Error",
        paste("Please check these input:", errors, "and try again!")
      ))
    } else {
      updateTabsetPanel(session, inputId = "Rnaseq_pre", selected = "Preprocessing")
    }
    
    # update input
    updateNumericInput(session, "min_col", max = ncol(DS)) # update max column nunmber in filtering
    if (is.null(spikes)) {
      updateRadioButtons(session, "norm_method", choices = c(
        "None (Black)" = "None",
        "RPKM (Blue)" = "RPKM", "FPKM (Dark cyan)" = "FPKM",
        "TPM (Dark green)" = "TPM",
        "Upper Quartile (Brown)" = "RUV"
      ))
      # c("None",'RPKM','FPKM','TPM',"Upper Quartile"="RUV")
    } else {
      updateRadioButtons(session, "norm_method", choices = c(
        "None (Black)" = "None",
        "RPKM (Blue)" = "RPKM", "FPKM (Dark cyan)" = "FPKM",
        "TPM (Dark green)" = "TPM",
        "RUV (Brown)" = "RUV"
      ))
    }
    if (is.null(lengths) & !(is.null(spikes))) {
      updateRadioButtons(session, "norm_method", choices = c("None (Black)" = "None", "RUV (Brown)" = "RUV"))
    } else if (is.null(lengths) & (is.null(spikes))) {
      updateRadioButtons(session, "norm_method", choices = c("None (Black)" = "None", "Upper Quartile (Brown)" = "RUV"))
    }
    
    if (is.null(f)) {
      hideTab(inputId = "navbar", target = "DE Analysis")
    } else {
      showTab(inputId = "navbar", target = "DE Analysis")
    }
    # if(is.null(f)){
    #   hideTab(inputId="preprocessing_tabs", target="Description table")
    # } else {
    #   showTab(inputId="preprocessing_tabs", target="Description table")
    # }
    
    hide("help_text_scatter")
    hide("help_text_dis_fit")
    hide("help_text_correlation")
    hide("help_text_PCA")
    hide("help_text_DE_anal")
    hide("help_text_heatmap")
    hide("help_text_Noise")
    hide("help_text_Entropy")
    # hide("help_text_SVM")
    hide("help_text_tsne")
    hide("help_text_rf")
    hide("help_text_SOM")
    
  })
  
  
  ############################ Micro array ###############################
  
  
  output$downloadMicroRaw <- downloadHandler(
    filename = function() {
      paste("Microarray_Raw", ".csv", sep = "")
    },
    content = function(file) {
      raw_data <- df_micro()
      
      micro_raw_data <- as.data.frame(raw_data@assayData[["exprs"]][1:2000,])
      micro_raw_data <- tibble::rowid_to_column(micro_raw_data, "Geneid")
      
      n <- ncol(micro_raw_data)
      for(i in 2:n)
      {
        colnames(micro_raw_data)[i] <- paste0('new_',colnames(micro_raw_data)[i])
      }
      
      write.csv(micro_raw_data, file, row.names = FALSE)
    }
  )
  
  output$downloadMicroMeta <- downloadHandler(
    filename = function() {
      paste("Microarray_Meta", ".csv", sep = "")
    },
    content = function(file) {
      raw_data <- df_micro()
      
      micro_raw_data <- as.data.frame(raw_data@assayData[["exprs"]][1:2000,])
      micro_meta_data <- as.data.frame(raw_data@phenoData@data[["Factor.Value.disease."]])
      micro_meta_data <- add_column(micro_meta_data, as.data.frame(colnames(micro_raw_data)), .after = 0)
      
      names(micro_meta_data)[1] <- "Id"
      names(micro_meta_data)[2] <- "Types"
      
      n <- nrow(micro_meta_data)
      micro_meta_data[,1] <- as.character(micro_meta_data[,1])
      for(i in 1:n)
      {
        micro_meta_data[i,1] <- paste0("new_",micro_meta_data[i,1])
      }
      
      write.csv(micro_meta_data, file, row.names = FALSE)
    }
  )

  
  ####################################################################
  
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
  df_norm <- reactive({ # get normalized counts
    if (is.null(input$file2)) {
      return(NULL)
    }
    parts <- strsplit(input$file2$datapath, ".", fixed = TRUE)
    type <- parts[[1]][length(parts[[1]])]
    if (type != "csv") {
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return(NULL)
    }
    ds <- read.csv(input$file2$datapath)
    ds <- na.omit(ds)
    ds <- ds[!duplicated(ds[, 1]), ] # remove duplicated gene names
    
    row_names <- ds[, 1]
    DS <- data.frame(ds)
    if (ncol(DS) <= 1) {
      showModal(modalDialog(
        title = "Error",
        "Please check normalised data file format (Eg_normalised.png) and try again!"
      ))
      return(NULL)
    }
    DS <- DS[, -1]
    row.names(DS) <- row_names
    for (i in 1:ncol(DS)) {
      if (class(DS[, i]) != "numeric" & class(DS[, i]) != "integer") {
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
  df_raw <- reactive({
    if (is.null(input$file1)) {
      return(NULL)
    }
    parts <- strsplit(input$file1$datapath, ".", fixed = TRUE)
    type <- parts[[1]][length(parts[[1]])]
    if (type != "csv") {
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return(NULL)
    }
    raw_ds <- read.csv(input$file1$datapath)
    raw_ds <- na.omit(raw_ds)
    raw_ds <- raw_ds[!duplicated(raw_ds[, 1]), ] # remove duplicated gene names
    
    # raw_ds <- as.data.frame(raw_ds)
    if (ncol(raw_ds) <= 1) {
      showModal(modalDialog(
        title = "Error",
        "Data file must contain at least 2 columns. Please check raw data format and try again!"
      ))
      return(NULL)
    }
    
    row_names <- raw_ds[, 1]
    rownames(raw_ds) <- row_names
    raw_DS <- raw_ds[, -1] # remove the first column, which is gene Id
    
    for (i in 1:ncol(raw_DS)) {
      if (class(raw_DS[, i]) != "numeric" & class(raw_DS[, i]) != "integer") {
        showModal(modalDialog(
          title = "Error",
          "Raw counts must be integer. Please check raw data formate and try again!"
        ))
        return(NULL)
      }
    }
    return(raw_DS)
  })
  
  ######################## Micro array ##########################
  
  
  df_micro <- reactive({
    print("running")
    if (is.null(input$file_micro)) {
      return(NULL)
    }
    parts <- strsplit(input$file_micro$datapath, ".", fixed = TRUE)
    type <- parts[[1]][length(parts[[1]])]
    
    if (type != "zip") {
      showModal(modalDialog(
        title = "Error",
        "Please input a zip file!"
      ))
      return(NULL)
    }
    
    unzip(input$file_micro$datapath,exdir = parts[[1]][1])
    fol_name <- print(list.files(parts[[1]][1]))
    micro_data_dir <- paste0(parts[[1]][1],"/",fol_name)
    print(micro_data_dir)
    
    sdrf_location <- file.path(micro_data_dir, "E-MTAB-2967.sdrf.txt")
    SDRF <- read.delim(sdrf_location)
    
    rownames(SDRF) <- SDRF$Array.Data.File
    SDRF <- AnnotatedDataFrame(SDRF)
    
    raw_data <- oligo::read.celfiles(filenames = file.path(micro_data_dir, 
                                                           SDRF$Array.Data.File),
                                     verbose = FALSE, phenoData = SDRF)
    
    print(raw_data)
    
    return(raw_data)
    
  })
  
  
  ###############################################################
  
  # get gene length
  gene_length <- reactive({
    if (is.null(input$length1)) {
      return(NULL)
    }
    lengths_df <- read.csv(input$length1$datapath)
    lengths_df2 <- data.frame("len" = lengths_df[, 2])
    rownames(lengths_df2) <- as.character(lengths_df[, 1])
    return(lengths_df2)
  })
  
  # get spikes / negative control genes
  neg_control <- reactive({
    if (is.null(input$spikes1)) {
      return(NULL)
    }
    spikes <- read.csv(input$spikes1$datapath, header = F)
    spikes <- as.character(spikes[, 1])
    # print(spikes[1:10])
    return(spikes)
  })
  
  # get meta data table
  group_names <- reactive({
    # if no data
    if (is.null(input$metafile1)) {
      return(NULL)
    }
    
    # read in group names (metadata)
    groups <- read.csv(input$metafile1$datapath)
    group_colnames <- as.character(groups[, 1])
    
    type <- input$file_type
    if (type == "norm") {
      DS <- df_norm()
    } else if (type == "raw") {
      DS <- df_raw()
    }
    col_names <- colnames(DS) # columm names of DS in order
    
    # check if groups and column names are similar
    if (!all(col_names %in% group_colnames) || ncol(groups) < 2) {
      showNotification(type = "error", "group names and DS column names not similar")
      return(NULL)
    }
    
    if (ncol(groups) == 2) {
      f <- groups[match(col_names, groups[, 1]), ] [, 2] # arrange f in the same order as col_names
    } else {
      f <- groups[match(col_names, groups[, 1]), ] [, 2]
      for (i in 3:ncol(groups)) {
        f <- paste0(f, "_", groups[, i])
      }
    }
    f <- as.factor(make.names(f))
    # return(as.factor(f))
    return(f)
  })
  
  ### Gene ontology
  
  gene_list <- reactive({
    if (is.null(input$filego)) {
      return(NULL)
    }
    parts <- strsplit(input$filego$datapath, ".", fixed = TRUE)
    type <- parts[[1]][length(parts[[1]])]
    if (type != "csv") {
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return(NULL)
    }
    ds <- read.csv(input$filego$datapath, header = FALSE)
    if (ncol(ds) >= 2) {
      col1 <- ds[-1, 1]
    } else if (ncol(ds) == 1) {
      col1 <- ds[, 1]
    } else {
      showModal(modalDialog(
        title = "Error",
        "No data found! Please check required data format and try again!"
      ))
      return(NULL)
    }
    gene_list <- as.character(col1)
    print("gene list from gene_list")
    print(head(gene_list))
    return(gene_list)
  })
  
  bg_list <- reactive({
    if (is.null(input$filebg)) {
      return(NULL)
    }
    parts <- strsplit(input$filebg$datapath, ".", fixed = TRUE)
    type <- parts[[1]][length(parts[[1]])]
    if (type != "csv") {
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return(NULL)
    }
    ds <- read.csv(input$filebg$datapath, header = FALSE)
    if (ncol(ds) > 1) {
      col1 <- ds[-1, 1]
    } else if (ncol(ds) == 1) {
      col1 <- ds[, 1]
    } else {
      showModal(modalDialog(
        title = "Error",
        "No data found! Please check required data format and try again!"
      ))
      return(NULL)
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
    DS <- DS_norm[keep, ]
    # DS <- apply(DS_norm, 1, function(x) length(x[x>min_val])>=min_col)
    return(DS)
  })
  
  # filter raw counts
  df_raw_filt <- eventReactive(input$submit_preprocessing, {
    DS_raw <- df_raw()
    min_val <- input$min_val
    min_col <- input$min_col
    keep <- rowSums(DS_raw >= min_val) >= min_col
    DS_filt <- DS_raw[keep, ]
    # DS_filt <- apply(DS_raw, 1, function(x) length(x[x>min_val])>=min_col)
    return(DS_filt)
  })
  
  # normalizing raw counts
  df_raw_shiny <- reactive({
    raw_DS <- df_raw_filt() # get filtered raw counts
    method <- input$norm_method
    
    if (method %in% c("TPM", "RPKM", "FPKM")) {
      lengths_df <- gene_length()
      merge_DS <- merge(raw_DS, lengths_df, by = "row.names")
      rownames(merge_DS) <- merge_DS[, 1]
      merge_DS <- merge_DS[, -1]
      raw_DS <- merge_DS[, -ncol(merge_DS)]
      lengths <- merge_DS[, ncol(merge_DS)]
      # print("length")
      # print(head(merge_DS))
    }
    # print("from line 981 df_raw_shiny")
    # print(method)
    # print("raw_DS")
    # print(head(raw_DS[,1:4]))
    # print("dimension of raw_DS")
    # print(dim(raw_DS))
    
    if (method == "TPM") {
      tpm.matrix <- apply(raw_DS, 2, function(x) tpm(x, lengths))
      tpm.df <- data.frame(tpm.matrix)
      return(tpm.df)
    } else if (method == "RPKM") {
      rpkm.matrix <- edgeR::rpkm(raw_DS, lengths)
      rpkm.df <- data.frame(rpkm.matrix)
      return(rpkm.df)
    } else if (method == "FPKM") {
      fpkm.matrix <- apply(raw_DS, 2, function(x) fpkm(x, lengths))
      fpkm.df <- data.frame(fpkm.matrix)
      return(fpkm.df)
    } else if (method == "None") {
      return(raw_DS)
    } else if (method == "RUV") {
      spikes <- neg_control()
      if (!is.null(spikes)) {
        spikes <- intersect(spikes, rownames(raw_DS))
      }
      # f <- group_names()
      # if( is.null(spikes) )
      #   spikes <- getEmpirical(rawDS,f)
      set1 <- RUVg.apply(raw_DS, spikes)
      RUV.df <- as.data.frame(normCounts(set1))
      return(RUV.df)
    }
  })
  
  ### for distribution fitting
  distfit_df <- reactive({
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    for (i in 1:ncol(DS)) {
      DS <- DS[which(DS[, i] > 0), ]
      DS <- na.omit(DS)
    }
    return(DS)
  })
  
  
  ######### ANALYSIS FROM HERE ############
  ######## RLEplot and Preprocessing ###########
  #############################################
  
  
  RLE.plot <- reactive({
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    set1 <- newSeqExpressionSet(as.matrix(DS))
    norm_method_name <- input$norm_method
    colors <- c("RPKM" = "blue", "FPKM" = "darkcyan", "TPM" = "darkgreen", "RUV" = "Brown", "Upper Quartile" = "Brown")
    if (norm_method_name != "None" & input$submit_preprocessing != 0) {
      spikes <- neg_control()
      if (norm_method_name == "RUV" & is.null(spikes)) {
        norm_method_name <- "Upper Quartile"
      }
      par(mar = c(7, 4, 4, 4) + 1.2)
      plotRLE(set1,
              ylim = c(-1.5, 1.5), outline = FALSE, col = colors[norm_method_name],
              las = 2,
              hjust = 1,
              main = paste(norm_method_name, "Normalized")
      )
    }
  })
  
  
  violin_plot <- reactive({
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    
    norm_method_name <- input$norm_method
    
    if (norm_method_name != "None" & input$submit_preprocessing != 0) {
      
      df <- as.data.frame(DS)
      df <- setNames(stack(df),c("norm_type","Genotype"))
      df$norm_type <- log(df$norm_type+1)
      
      p <- ggplot(df, aes(x=Genotype, y=norm_type)) + 
        geom_violin(trim=FALSE) + 
        scale_color_brewer(palette="Dark2") +
        labs(title=paste(norm_method_name,"Normalized",sep = ' '), y = paste("log(",norm_method_name,"+1)",sep = '') )+
        stat_summary(fun.data=mean_sdl, mult=1, 
                     geom="pointrange", color="red") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
      p
      
      ggplotly(p, tooltip = c("text"))
      
    }
  })
  
  
  output$RLE.plot <- renderPlot({
    RLE.plot()
  })
  
  output$violin_plot <- renderPlotly({
    violin_plot()
  })
  
  output$RLE.plot2 <- renderPlot({ # for raw data
    start.rle <- Sys.time()
    type <- input$file_type
    if (type == "norm") {
      raw_DS <- df_shiny()
      main_title <- "Input data"
    } else if (type == "raw") {
      raw_DS <- df_raw()
      main_title <- "Raw data"
    }
    
    set1 <- newSeqExpressionSet(as.matrix(raw_DS))
    if (input$submit_preprocessing != 0) {
      par(mar = c(7, 4, 4, 4) + 1.2)
      plotRLE(set1, ylim = c(-1.5, 1.5), outline = FALSE, main = main_title, las = 2)
    }
    end.rle <- Sys.time()
    print("time for RLE plot and preprocessing")
    print(end.rle - start.rle)
  })
  
  violin_plot2 <- reactive({
    type <- input$file_type
    if (type == "norm") {
      raw_DS <- df_shiny()
      main_title <- "Input data"
    } else if (type == "raw") {
      raw_DS <- df_raw()
      main_title <- "Raw data"
    }
    
    df <- as.data.frame(raw_DS)
    df <- setNames(stack(df),c("norm_type","Genotype"))
    df$norm_type <- log(df$norm_type+1)
    
    p <- ggplot(df, aes(x=Genotype, y=norm_type)) + 
      geom_violin(trim=FALSE) + 
      scale_color_brewer(palette="Dark2") +
      labs(title="Raw Data", y = paste("log(Raw Data+1)",sep = '') ) +
      stat_summary(fun.data=mean_sdl, mult=1, 
                   geom="pointrange", color="red") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    p
    
    ggplotly(p, tooltip = c("text"))
  })
  
  output$violin_plot2 <- renderPlotly({
    violin_plot2()
  })
  
  output$norm_table <- DT::renderDataTable({
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    # if(input$submit_preprocessing != 0)
    DS # with filtering and normalization
  })
  
  output$meta_table <- DT::renderDataTable({
    f <- group_names()
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    if (!is.null(f)) {
      meta_df <- data.frame("Column names" = colnames(DS), "Description" = f)
      meta_df
    }
  })
  
  output$download_norm_data <- downloadHandler(
    filename = function() {
      method <- input$norm_method
      paste(method, "normalized.csv")
    },
    content = function(file) {
      type <- input$file_type
      if (type == "norm") {
        DS <- df_shiny()
      } else if (type == "raw") {
        DS <- df_raw_shiny()
      }
      write.csv(DS, file, row.names = F)
    }
  )
  
  ############################
  ######## scatter ###########
  ############################
  
  # input scatter data
  plotScatter <- reactive({
    scatter.start <- Sys.time()
    trans <- input$trans
    x <- input$scatter.x
    y <- input$scatter.y
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    if (trans == "None") {
      scatter.data <- DS
    } else if (trans == "Natural log") {
      scatter.data <- log1p(DS)
    } else if (trans == "log2") {
      scatter.data <- log2(DS + 1)
    } else if (trans == "log10") {
      scatter.data <- log10(DS + 1)
    }
    scatter.end <- Sys.time()
    print("Scatter plot time")
    print(scatter.end - scatter.start)
    return(list(x, y, scatter.data))
  })
  
  
  get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }
  
  scatterplot <- function() {
    # get values from list
    li <- plotScatter()
    xval <- li[[1]]
    yval <- li[[2]]
    scatter.data <- li[[3]]
    
    # data frame needed for ggplot
    df <- data.frame(t1 = scatter.data[, xval], t2 = scatter.data[, yval])
    
    # get 2d kernel density estimate
    df$density <- get_density(df$t1, df$t2, n = 100)
    
    # plot heat scatter w/ ggplot
    p <- ggplot(df, aes(x = t1, y = t2, color = density, text = paste(xval, ": ", round(t1, 4), "\n", yval, ": ", round(t2, 4), sep = ""), group = 1)) +
      geom_point(shape = 19, size = 0.25) +
      scale_color_viridis()
    
    # modify label and fill defaults
    p <- p + xlab(xval) + ylab(yval) + labs(color = "KDE", title = paste("R=", round(cor(scatter.data[, xval], scatter.data[, yval]), 3)))
    
    # if checkbox is ticked, display regression line
    if (input$regline == TRUE) {
      p <- p + geom_smooth(method = lm, se = FALSE, size = 0.5, color = "blue")
    }
    p
    
    # add interactivity w/ plotly
    ggplotly(p, tooltip = c("text"))
  }
  
  scatterplot_collage <- function() {
    li <- plotScatter()
    scatter.data <- li[[3]]
    par(mfrow = c(3, 3))
    for (i in 1:ncol(scatter.data)) {
      for (j in i:ncol(scatter.data)) {
        d <- kde2d(scatter.data[, i], scatter.data[, j])
        ColorLevels <- round(seq(min(d$z), max(d$z), length = 5), 4)
        heatscatter(x = scatter.data[, i], y = scatter.data[, j], xlab = colnames(scatter.data)[i], ylab = colnames(scatter.data)[j], main = "")
        legend("topleft", paste("R=", round(cor(scatter.data[, i], scatter.data[, j]), 3)), bty = "n")
        legend("bottomright", title = "KDE", legend = ColorLevels, pch = 19, col = LSD::colorpalette("heat"))
        if (i != j) {
          lines(lowess(scatter.data[, i], scatter.data[, j]), col = "black")
        }
      }
    }
  }
  
  output$scatter.plot <- renderPlotly({
    scatterplot()
  })
  
  output$downloadscatter_collage <- downloadHandler(
    filename = function() {
      paste("heatscatter_collage", ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      scatterplot_collage()
      dev.off()
    }
  )
  
  output$downloadscatter <- downloadHandler(
    filename = function() {
      paste("heatscatter", ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      scatterplot()
      dev.off()
    }
  )
  
  output$help_text_scatter <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          Scatter plot compares global gene expression
          between 2 conditions.<br>Color density is calculated based on 2D Gaussian kernel density.
          </b>
        </p>
      </center>
    ")
  })
  
  ############################
  ######## distfit ###########
  ############################
  
  output$downloaddist <- downloadHandler(
    filename = function() {
      paste("distribution_fit", ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      distplot()
      dev.off()
    }
  )
  
  output$dist_range_allowed <- renderText({
    DS <- distfit_df()
    paste("Suggested range: ( 0", " ~ ", round(max(DS)), " ]", sep = "")
  })
  
  plotDist <- reactive({
    dist.start <- Sys.time()
    dis <- input$distributions
    var <- input$dist.var
    DS <- distfit_df()
    fits <- list()
    distrs <- NULL
    numcol <- c(0, 0, 0, 0, 0, 0)
    dist_zoom <- input$dist_zoom
    if (dist_zoom == "slider") {
      fit_range <- input$dist_range
    } else if (dist_zoom == "text input") {
      fit_range <- c(input$dist_range_min, input$dist_range_max)
    }
    if ("Log-normal" %in% dis) {
      fit_ln <- fitdist(DS[, var], "lnorm")
      fits <- c(fits, list(fit_ln))
      distrs <- c(distrs, "Log-normal")
      numcol[1] <- 1
    }
    if ("Log-logistic" %in% dis) {
      fit_ll <- fitdist(DS[, var], "llogis", start = list(shape = 10, scale = 10), lower = c(0, 0))
      fits <- c(fits, list(fit_ll))
      distrs <- c(distrs, "Log-logistic")
      numcol[2] <- 1
    }
    if ("Pareto" %in% dis) {
      fit_P <- fitdist(DS[, var], "pareto", start = list(shape = 10, scale = 10), lower = c(0, 0))
      fits <- c(fits, list(fit_P))
      distrs <- c(distrs, "Pareto")
      numcol[3] <- 1
    }
    if ("Burr" %in% dis) {
      fit_B <- fitdist(DS[, var], "burr", start = list(shape1 = 0.3, shape2 = 1, rate = 1), lower = c(0, 0, 0))
      fits <- c(fits, list(fit_B))
      distrs <- c(distrs, "Burr")
      numcol[4] <- 1
    }
    if ("Weibull" %in% dis) {
      fit_W <- fitdist(DS[, var], "weibull", lower = c(0, 0))
      fits <- c(fits, list(fit_W))
      distrs <- c(distrs, "Weibull")
      numcol[5] <- 1
    }
    if ("Gamma" %in% dis) {
      fit_G <- fitdist(DS[, var], "gamma", lower = c(0, 0), start = list(scale = 1, shape = 1))
      fits <- c(fits, list(fit_G))
      distrs <- c(distrs, "Gamma")
      numcol[6] <- 1
    }
    dist.end <- Sys.time()
    print("Distribution fitting time")
    print(dist.end - dist.start)
    return(list(fits, distrs, numcol, var, fit_range))
  })
  
  output$dist.plot <- renderPlot({
    distplot()
  })
  
  distaic <- reactive({
    dist.start <- Sys.time()
    DS <- distfit_df()
    AIC.df <- as.data.frame(matrix(nrow = ncol(DS), ncol = 6))
    rownames(AIC.df) <- colnames(DS)
    colnames(AIC.df) <- c("Log-normal", "Log-logistic", "Pareto", "Burr", "Weibull", "Gamma")
    for (i in 1:nrow(AIC.df)) {
      fit_ln <- fitdist(DS[, i], "lnorm")
      fit_ll <- fitdist(DS[, i], "llogis", start = list(shape = 10, scale = 10), lower = c(0, 0))
      fit_P <- fitdist(DS[, i], "pareto", start = list(shape = 10, scale = 10), lower = c(0, 0))
      fit_B <- fitdist(DS[, i], "burr", start = list(shape1 = 0.3, shape2 = 1, rate = 1), lower = c(0, 0, 0))
      fit_W <- fitdist(DS[, i], "weibull", lower = c(0, 0))
      fit_G <- fitdist(DS[, i], "gamma", lower = c(0, 0), start = list(scale = 1, shape = 1))
      fits <- list(fit_ln, fit_ll, fit_P, fit_B, fit_W, fit_G)
      AIC.df[i, ] <- gofstat(fits)$aic
    }
    for (i in 1:nrow(AIC.df)) {
      AIC.df$min.AIC[i] <- colnames(AIC.df)[which.min(AIC.df[i, 1:6])]
    }
    dist.end <- Sys.time()
    print("distribution fitting time")
    print(dist.end - dist.start)
    return(AIC.df)
  })
  
  output$dist.aic <- renderTable(
    {
      distaic()
    },
    rownames = TRUE
  )
  
  distplot <- function() {
    li <- plotDist()
    fits <- li[[1]]
    distrs <- li[[2]]
    numcol <- li[[3]]
    var <- li[[4]]
    fit_range <- li[[5]]
    line_types <- c(1, 2, 3, 4, 5, 6) # par lty
    if (length(fits) != 0) {
      cdfcomp(fits,
              xlogscale = TRUE, ylogscale = TRUE,
              ylab = "CDF", xlab = "Expression levels (log)", xlim = c(fit_range[1], fit_range[2]), ylim = c(10^-3, 1),
              legendtext = distrs, cex = 0.5, main = var, fitcol = rainbow(6)[which(numcol == 1)], fitlty = line_types[which(numcol == 1)]
      )
    }
  }
  
  output$downloaddistaic <- downloadHandler(
    filename = function() {
      paste("aic", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(distaic(), file, row.names = TRUE)
    }
  )
  
  output$help_text_dis_fit <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          Fitting the selected probability distribution(s) to transcriptome-wide data of the selected sample.<br>
           Cumulative distribution function will be showned, with black lines being the empirical (transcriptomic) data, 
           and colored lines being the probability distribution(s) with best fitted parameters. <br>
           The AIC (Akaike information criterion) table provides a comparison in goodness-of-fit to transcriptomic data
           among selected probability distribution(s).
          </b>
        </p>
      </center>
    ")
  })
  
  ############################
  ####### correlation ########
  ############################
  
  COR <- function(d, i, myMethod) {
    Result2 <- cor(x = d[, i], y = d[, i], method = myMethod)
    return(format(round(Result2, 5), nsmall = 5))
  }
  
  cor_df <- reactive({
    cor.start <- Sys.time()
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    method <- input$cor_method
    if (method == "Pearson correlation") {
      Cor2 <- data.frame(COR((DS), 1:length(DS), "pearson"))
    } else if (method == "Spearman correlation") {
      Cor2 <- data.frame(COR((DS), 1:length(DS), "spearman"))
    }
    Cor2 <- na.omit(Cor2)
    cor.end <- Sys.time()
    print("correlation time")
    print(cor.end - cor.start)
    return(Cor2)
  })
  
  output$corr.plot <- renderPlot({
    corrplot1()
  })
  
  output$corr.plot2 <- renderPlot({
    corrplot2()
  })
  
  output$corr.matrix <- renderTable(
    {
      cor_df()
    },
    rownames = TRUE
  )
  
  corrplot1 <- function() {
    corr <- as.matrix(cor_df())
    corr <- apply(corr, 2, as.numeric)
    rownames(corr) <- rownames(cor_df())
    if (ncol(corr) <= 20) {
      fontsize <- 1
    } else {
      fontsize <- 20 / ncol(corr)
    }
    corrplot(corr, method = "shade", shade.col = NA, tl.col = "black", cl.lim = c(min(corr), 1), is.corr = FALSE, tl.cex = fontsize)
  }
  
  corrplot2 <- function() {
    corr <- as.matrix(cor_df())
    corr <- apply(corr, 2, as.numeric)
    rownames(corr) <- rownames(cor_df())
    if (ncol(corr) <= 20) {
      fontsize <- 1
    } else {
      fontsize <- 20 / ncol(corr)
    }
    corrplot(corr, type = "upper", tl.col = "black", cl.lim = c(min(corr), 1), is.corr = FALSE, tl.cex = fontsize)
  }
  
  output$downloadcorrplot <- downloadHandler(
    filename = function() {
      paste("corrheatmap", ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      corrplot1()
      dev.off()
    }
  )
  
  output$downloadcorrplot2 <- downloadHandler(
    filename = function() {
      paste("corrplot", ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      corrplot2()
      dev.off()
    }
  )
  
  output$downloadcorrmat <- downloadHandler(
    filename = function() {
      paste("correlation", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(cor_df(), file, row.names = TRUE)
    }
  )
  
  output$help_text_correlation <- renderUI({
    HTML("
    <br>
    <br>
      <center>
      <ol>
        <li>
        <p>
          <b>
          <h4>Pearson correlation</h4><br>
          Pearson correlation measures linear relationship between two vectors, where r = 1 if the two vectors are identical, 
          and r = 0 if there are no linear relationships between the vectors.<br>
          The correlation coefficient r between two vectors (e.g. transcriptome in two different samples), 
          containing n observations (e.g. gene expression values), is defined by (for large n):<br>
          <img src='https://www.statisticssolutions.com/wp-content/uploads/2019/09/ehtsht.png' alt='Pearson-correlation' border='0'><br>
          where x<sub>i</sub> and y<sub>i</sub> are the ith observation in the vectors X and Y, respectively.
          </b>
        </p>
        </li>
        <br>
        <li>
        <p>
          <b>
          <h4>Spearman correlation</h4><br>
          Spearman correlation is a non-parametric test that measrues the degree of association 
          between two vectors  (e.g. transcriptome in two different samples).  The Spearman rank correlation 
          test does not carry any assumptions about the distribution of the data and is the appropriate correlation 
          analysis when the variables are measured on a scale that is at least ordinal.
          The following formula is used to calculate the Spearman rank correlation:<br>
          <img src='https://i.ibb.co/rkjbg1d/Spearman-correlation.png' alt='Spearman-correlation' border='0'><br>
          ρ = Spearman rank correlation<br>
          r<sub>x,i</sub>, r<sub>y,i</sub> = ranks of corresponding variables (or gene)<br>
          n = number of observations
          </b>
        </p>
        </li>
        </ol>
      </center>
    ")
  })
  
  ############################
  #######     PCA     ########
  ############################
  
  refreshDS1 <- eventReactive(input$pca_refresh, {
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    DS1 <- DS[sample(nrow(DS), nrow(DS), replace = FALSE), ]
    return(DS1)
  })
  
  plotPCA <- reactive({ # process and return data
    pca.start <- Sys.time()
    type <- input$file_type
    ############################
    pca_type <- input$pca_type
    ############################
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    order <- input$gene_order
    size <- input$gene_size
    x <- input$pca.x
    y <- input$pca.y
    cluster_flag <- input$pca_cluster
    rindex <- as.numeric(substring(x, 3))
    cindex <- as.numeric(substring(y, 3))
    if (order == "Ascending") {
      DS1 <- DS[order(DS[, 1]), ]
    } else if (order == "Descending") {
      DS1 <- DS[rev(order(DS[, 1])), ]
    } else if (order == "Random") {
      DS1 <- refreshDS1()
    }
    
    DSample <- head(DS1, n = size)
    
    ##### PCA & Sparse PCA #####
    if (pca_type == "PCA") {
      PR <- prcomp(t(DSample), center = TRUE)
      print("Normal PCA selected")
    }
    else if (pca_type == "SPCA") {
      PR <- spca(t(DSample), scale = FALSE, center = TRUE, max_iter = 10)
      PR$x <- PR$scores
      print("Sparse PCA selected")
    }
    
    col_val_x <- as.numeric(gsub("[^[:digit:]]", "", x))
    col_val_y <- as.numeric(gsub("[^[:digit:]]", "", y))
    #####################
    
    
    PCA.var <- PR$sdev^2
    PCA.var.per <- round(PCA.var / sum(PCA.var) * 100, 1)
    xlabel <- paste(colnames(PR$x)[rindex], " - ", PCA.var.per[rindex], "%", sep = "")
    ylabel <- paste(colnames(PR$x)[cindex], " - ", PCA.var.per[cindex], "%", sep = "")
    if (cluster_flag == TRUE) {
      num <- as.numeric(input$pca_cluster_num)
      ####################################################################################
      kmeans.data <- data.frame(x = PR$x[, col_val_x], y = PR$x[, col_val_y])
      print(kmeans.data)
      ####################################################################################
      kmeans.result <- kmeans(kmeans.data, num)
      return(list(PR, PCA.var, PCA.var.per, rindex, cindex, xlabel, ylabel, cluster_flag, kmeans.result))
    }
    pca.end <- Sys.time()
    print("pca time")
    print(pca.end - pca.start)
    return(list(PR, PCA.var, PCA.var.per, rindex, cindex, xlabel, ylabel, cluster_flag))
  })
  
  pcavarplot <- function() {
    li <- plotPCA()
    PCA.var.per <- li[[3]] / 100
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    pcchoices <- NULL
    for (i in 1:length(PCA.var.per)) {
      pcchoices <- c(pcchoices, paste("PC", i, sep = ""))
    }
    xform <- list(
      categoryorder = "array",
      categoryarray = pcchoices
    )
    p <- plot_ly(
      x = pcchoices,
      y = PCA.var.per,
      name = "PCA variance",
      type = "bar"
    ) %>% layout(xaxis = xform)
    
    return(p)
  }
  
  pca2dplot <- function() {
    li <- plotPCA()
    PR <- li[[1]]
    rindex <- li[[4]]
    cindex <- li[[5]]
    xlabel <- li[[6]]
    ylabel <- li[[7]]
    cluster_flag <- li[[8]]
    if (cluster_flag == FALSE) {
      p <- plot_ly(
        x = PR$x[, rindex],
        y = PR$x[, cindex],
        type = "scatter",
        mode = "markers"
      ) %>% layout(xaxis = list(title = xlabel), yaxis = list(title = ylabel))
    } else if (cluster_flag == TRUE) {
      kmeans.result <- li[[9]]
      text_flag <- input$pca_text
      if (text_flag == TRUE) {
        p <- plot_ly(
          x = PR$x[, rindex],
          y = PR$x[, cindex],
          type = "scatter",
          color = as.character(kmeans.result$cluster),
          mode = "markers",
          colors = "Set1"
        ) %>%
          hide_colorbar() %>%
          add_trace(
            x = PR$x[, rindex],
            y = PR$x[, cindex],
            type = "scatter",
            mode = "text",
            text = names(kmeans.result$cluster),
            textposition = "top right"
          ) %>%
          layout(xaxis = list(title = xlabel), yaxis = list(title = ylabel), showlegend = FALSE)
      } else if (text_flag == FALSE) {
        p <- plot_ly(
          x = PR$x[, rindex],
          y = PR$x[, cindex],
          type = "scatter",
          color = as.character(kmeans.result$cluster),
          mode = "markers",
          text = names(kmeans.result$cluster),
          colors = "Set1"
        ) %>%
          hide_colorbar() %>%
          layout(xaxis = list(title = xlabel), yaxis = list(title = ylabel), showlegend = FALSE)
      }
    }
  }
  
  pca3dplot <- function() {
    li <- plotPCA()
    PR <- li[[1]]
    xlabel <- "PC1"
    ylabel <- "PC2"
    zlabel <- "PC3"
    cluster_flag <- li[[8]]
    if (cluster_flag == FALSE) {
      p <- plot_ly(
        x = PR$x[, 1],
        y = PR$x[, 2],
        z = PR$x[, 3],
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 5)
      ) %>% layout(scene = list(xaxis = list(title = xlabel), yaxis = list(title = ylabel), zaxis = list(title = zlabel)))
    } else if (cluster_flag == TRUE) {
      kmeans.result <- li[[9]]
      text_flag <- input$pca_text
      if (text_flag == TRUE) {
        p <- plot_ly(
          x = PR$x[, 1],
          y = PR$x[, 2],
          z = PR$x[, 3],
          type = "scatter3d",
          mode = "text",
          text = names(kmeans.result$cluster),
          color = as.character(kmeans.result$cluster),
          textfont = list(size = 10),
          textposition = "top right"
        ) %>% layout(scene = list(xaxis = list(title = xlabel), yaxis = list(title = ylabel), zaxis = list(title = zlabel)), showlegend = FALSE)
      } else if (text_flag == FALSE) {
        p <- plot_ly(
          x = PR$x[, 1],
          y = PR$x[, 2],
          z = PR$x[, 3],
          type = "scatter3d",
          color = as.character(kmeans.result$cluster),
          mode = "markers",
          marker = list(size = 5),
          text = names(kmeans.result$cluster),
          colors = "Set1"
        ) %>%
          hide_colorbar() %>%
          layout(scene = list(xaxis = list(title = xlabel), yaxis = list(title = ylabel), zaxis = list(title = zlabel)), showlegend = FALSE)
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
    filename = function() {
      paste("pca_variance", ".png", sep = "")
    },
    content = function(file) {
      p <- pcavarplot()
      orca(p, file = "pca_variance.png")
    }
  )
  
  output$downloadpca2d <- downloadHandler(
    filename = function() {
      paste("pca2d", ".png", sep = "")
    },
    content = function(file) {
      p <- pca2dplot()
      orca(p, file = "pca2d.png")
    }
  )
  
  output$downloadpca3d <- downloadHandler(
    filename = function() {
      paste("pca3d", ".png", sep = "")
    },
    content = function(file) {
      p <- pca3dplot()
      plotly_IMAGE(p, format = "png", out_file = "pca3d.png")
    }
  )
  
  output$help_text_PCA <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          Principal Components Analysis (PCA) is an multivariate statistical technique for simplifying high-dimentional 
          data sets (Basilevsky 1994). Given m observations on n variables, the goal of PCA is to reduce the dimensionality 
          of the data matrix by finding r new variables, where r is less than n. Termed principal components, these r new 
          variables together account for as much of the variance in the original n variables as possible while remaining 
          mutually uncorrelated and orthogonal. Each principal component is a linear combination of the original variables, 
          and so it is often possible to ascribe meaning to what the components represent.
          A PCA analysis of transcriptomic data consider the genes as variables, creating a set of “principal gene components” that indicate the features of genes that best explain the experimental responses they produce.
          </b>
        </p>
        <p>
          <b>
          A PCA analysis of transcriptomic data consider the genes as variables, creating a set of 
          “principal gene components” that indicate the features of genes that best explain the experimental responses they produce.
          </b>
        </p>
        <p>
          <b>
          To compute the principal components, the n eigenvalues and their corresponding eigenvectors are calculated 
          from the n×n covariance matrix of conditions. Each eigenvector defines a principal component. A component can be 
          viewed as a weighted sum of the conditions, where the coefficients of the eigenvectors are the weights. 
          The projection of gene i along the axis defined by the i<sup>th</sup> principal component is:<br>
          <img src='https://i.ibb.co/n6M8n79/Screenshot-from-2020-09-05-15-33-31.png' alt='Screenshot-from-2020-09-05-15-33-31' border='0'><br>
          Where <var>v<sub>tj</sub></var> is the t<sup>th</sup> coefficient for the i<sup>th</sup> principal component; <var>a<sub>it</sub></var> is the 
          expression measurement for gene i under the t<sup>th</sup> condition. 
          A′ is the data in terms of principal components. Since V is an orthonormal matrix, A′ is a rotation of the data from the original space of 
          observations to a new space with principal component axes.
          </b>
        </p>
        <p>
          <b>
          The variance accounted for by each of the components is its associated eigenvalue; 
          it is the variance of a component over all genes. Consequently, the eigenvectors with large eigenvalues 
          are the ones that contain most of the information; eigenvectors with small eigenvalues are uninformative.
          </b>
        </p>
      </center>
    ")
  })
  
  ############################
  ######## DE analysis #######
  ############################
  group_names_de <- reactive({
    f <- group_names()
    f1 <- input$f1
    f2 <- input$f2
    if (is.null(f) || is.null(f1) || is.null(f2)) {
      return(NULL)
    }
    
    type <- input$file_type
    if (type == "norm") {
      raw_DS <- df_shiny()
    } else if (type == "raw") {
      raw_DS <- df_raw_filt()
    }
    f.df <- data.frame("f" = f)
    rownames(f.df) <- colnames(raw_DS)
    # f.df.slect <- subset(f.df,f %in% c(f1,f2) )
    # f.df.slect <- rbind( subset(f.df,f %in% f1),  subset(f.df,f %in% f2) )
    # f.df.slect2 <- f.df.slect;
    # f_new <- f.df.slect[,1]
    # f.df.slect2[,1] <- droplevels(f_new, except = levels(f_new)%in%f_new)
    return(f.df)
  })
  
  df_raw_de <- reactive({
    f_de <- group_names_de()
    if (is.null(f_de)) {
      return(NULL)
    }
    type <- input$file_type
    rep_number <- input$n_rep
    if (rep_number == 1) {
      de_type <- input$de_method1
    } else {
      de_type <- input$de_method0
    }
    
    if (type == "norm") {
      raw_DS <- df_shiny() # filtered and normalized
    } else if (type == "raw") {
      if (de_type == "NOISeq") {
        raw_DS <- df_raw_shiny() # filtered and normalized
      } else {
        raw_DS <- df_raw_filt() # filtered and UN-NORMALIZED
      }
    }
    # raw_DS_de <- raw_DS[,rownames(f_de)]
    return(raw_DS)
  })
  
  # all helper functions are in utils.R file
  de_no_filt <- eventReactive(input$submit_DE, { # return as table object
    start.de.table <- Sys.time()
    f_de <- group_names_de() # for edgeR >> f=f_de[,1]; for the rest factors=f_de
    if (is.null(f_de)) {
      return(NULL)
    }
    DS_de <- df_raw_de() # with only 2 conditions for DE analysis
    p_val <- 1 # input$p_val
    fc <- 1 # input$fc
    f1 <- input$f1
    f2 <- input$f2
    rep_number <- input$n_rep # either 0 = no replicates or 1 = have replicates
    
    
    spikes <- neg_control()
    norm_method <- input$norm_method
    if (!is.null(spikes) & norm_method == "RUV") {
      set1 <- RUVg.apply(DS_de, spikes)
      W_1 <- pData(set1)$W_1
    } else {
      W_1 <- NULL
    }
    # print("from de_no_filt")
    # print(pData(set1))
    # print(W_1)
    
    if (rep_number == 1) { # have replicates
      de_type <- input$de_method1
      if (de_type == "EdgeR") {
        res <- edgerApply(DS = DS_de, f = f_de[, 1], W_1 = W_1, f1 = f1, f2 = f2) # edgeR, return edgeR object
        res.df <- edgerFilter(res, FC = fc, p_val = p_val) # fitler edgeR object result
      } else if (de_type == "DESeq2") {
        res <- deseqApply(DS = DS_de, f.df = f_de, W_1 = W_1, f1 = f1, f2 = f2) # DESeq, return DESeq object
        res.df <- deseqFilter(res, FC = fc, p_val = p_val) # fitler DESeq object result
      } else if (de_type == "NOISeq") {
        res <- noiseqbioApply(DS = DS_de, f.df = f_de, f1 = f1, f2 = f2) # NOISeqbio, return NOIseq object
        res.df <- noiseqbioFilter(res, FC = fc, p_val = p_val) # filter return NOIseq object
      }
    } else { # no replicates
      de_type <- input$de_method0 # NOISeq
      res <- noiseqsimApply(DS = DS_de, f.df = f_de, f1 = f1, f2 = f2) # NOISeqbio, return NOIseq object
      res.df <- noiseqsimFilter(res, FC = fc)
    }
    end.de.table <- Sys.time()
    print("de table time")
    print(end.de.table - start.de.table)
    return(res.df)
  })
  
  
  de_filt <- function(res.df, p_val, fc, rep_number) {
    # res.df <- de_no_filt()
    if (is.null(res.df)) {
      return(NULL)
    }
    # p_val <- input$p_val
    # fc <- input$fc
    # rep_number <- input$n_rep
    
    if (rep_number == 1) {
      res.df.filt <- filter(res.df, FDR <= p_val, log2FCabs >= log2(fc))
    } else {
      res.df.filt <- filter(res.df, FDR <= p_val, log2FCabs >= log2(fc))
    }
    return(res.df.filt)
  }
  
  output$DE_table <- DT::renderDataTable({
    res.df <- de_no_filt()
    p_val <- input$p_val
    fc <- input$fc
    rep_number <- input$n_rep
    if (input$submit_DE > 0) {
      res.df.filt <- de_filt(res.df, p_val, fc, rep_number)
      res.df.filt
    }
  })
  
  ##### volcano plot ######
  volcano_plot <- eventReactive(input$submit_DE, {
    volcano.start.time <- Sys.time()
    rep_number <- input$n_rep
    if (rep_number == 0) {
      return(NULL)
    }
    de_type <- input$de_method1
    if (de_type == "NOISeq") {
      return(NULL)
    }
    
    res <- de_no_filt() # de result, no filter
    if (is.null(res)) {
      return(NULL)
    }
    
    p_val <- input$p_val
    fc <- input$fc
    res$Gene <- rownames(res)
    res <- na.omit(res)
    # plot
    ymax <- quantile(-log10(res$PValue), c(0.98))
    xmax <- quantile(abs(res$log2FC), c(0.98))
    # if (ymax > 5) ymax <- 5
    # print("from volcano plot - range(res$PValue)")
    # print(range(res$PValue))
    with(res, plot(log2FC, -log10(PValue), pch = 20, main = "Volcano plot", 
                   xlim = c(-xmax, xmax), ylim = c(0, ymax))) # xlim=c(-5,5),ylim=c(0,ymax)
    # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
    with(subset(res, FDR < p_val), points(log2FC, -log10(PValue), pch = 20, col = "red"))
    with(subset(res, abs(log2FC) > log2(fc)), points(log2FC, -log10(PValue), pch = 20, col = "orange"))
    with(subset(res, FDR < p_val & abs(log2FC) > log2(fc)), points(log2FC, -log10(PValue), pch = 20, col = "green"))
    legend("topleft",
           bty = "n", col = c("red", "orange", "green", "black"), pch = 19,
           legend = c("FDR < FDR limit", "FC > FC limit", "Both", "Other")
    )
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
    f_de <- group_names_de() # for edgeR >> f=f_de[,1]; for the rest factors=f_de
    if (is.null(f_de)) {
      return(NULL)
    }
    rep_number <- input$n_rep # either 0 or 1
    if (rep_number == 0) {
      return(NULL)
    }
    DS_de <- df_raw_de() # with only 2 conditions for DE analysis
    p_val <- 1 # input$p_val
    fc <- 1 # input$fc
    de_type <- input$de_method1
    if (de_type == "EdgeR") {
      edgerDisp(DS_de, f_de[, 1])
    } else if (de_type == "DESeq2") {
      deseqDisp(DS_de, f_de)
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
    filename = function() {
      paste0("DE analysis", ".csv")
    },
    content = function(file) {
      res.df <- de_no_filt()
      p_val <- input$p_val
      fc <- input$fc
      rep_number <- input$n_rep
      res.df.filt <- de_filt(res.df, p_val, fc, rep_number)
      write.csv(res.df.filt, file, row.names = F)
    }
  )
  
  output$download_volcano <- downloadHandler(
    filename = function() {
      paste0("Volcano", ".pdf")
    },
    content = function(file) {
      pdf(file)
      volcano_plot()
      dev.off()
    }
  )
  
  output$download_dispersion <- downloadHandler(
    filename = function() {
      paste0("Dispersion plot", ".pdf")
    },
    content = function(file) {
      pdf(file)
      dispersion_plot()
      dev.off()
    }
  )
  
  output$help_text_DE_anal <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          DE analysis identifies the genes that are statistically different in expression levels between the 2 selected conditions. Two important threshold are:
          </b>
        </p>
        <p>
          <b>
            1. The lower bound of expression fold change  between the 2 selected conditions<br>
            2. The upper bound of hypothesis test p-value
          </b>
        </p>
        <p>
          <b>
            GeneCloudOmics implements 3 popular methods to identify DE genes:<br>
              1. <a href='https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8'>DESeq2</a><br>
              2. <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/'>EdgeR</a><br>
              3. <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4666377/'>NOISeq</a>
          </b>
        </p>
      </center>
    ")
  })
  
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
    for (i in 1:input$numOfCluster) {
      display <- c(display, i)
    }
    selectInput("display_cluster", "Display cluster", choices = display)
  })
  ################
  
  setOneWithinFold <- function(arr) { # logFC
    fold <- as.numeric(input$fold)
    for (i in 1:length(arr)) {
      if ((arr[i] <= (fold)) & (arr[i] >= (1 / fold))) {
        arr[i] <- 1
      }
    }
    return(arr)
  }
  
  plotHeatmap <- eventReactive(input$heatmap_plot, { # process and return data
    heatmap.start.time <- Sys.time()
    
    type <- input$file_type
    value <- input$heatmap_value
    de_type <- input$heatmap_de_ind
    
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
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
    
    if (de_type == "ind") {
      fold <- as.numeric(input$fold)
      fold_ncol <- input$fold_ncol
      DS2 <- deWithoutStats(DS, FC = fold, n_col = fold_ncol)
      de_genes <- rownames(DS2)
      # print("from heatmap YT version")
      # print("de genes")
      # print(head(de_genes))
      # print(paste(length(de_genes),"genes"))
    } else if (de_type == "de") {
      res.df <- de_no_filt()
      if (is.null(res.df)) {
        return(NULL)
      }
      p_val <- input$p_val
      fc <- input$fc
      rep_number <- input$n_rep
      res.df.filt <- de_filt(res.df, p_val, fc, rep_number)
      de_genes <- res.df.filt$Gene
      # print("from line heatmap de result from DE analysis")
      # print("res.df.filt")
      # print(head(res.df.filt))
    }
    
    de_genes_exp <- DS[rownames(DS) %in% de_genes, ]
    DS3 <- t(scale(t(de_genes_exp)))
    DS3 <- na.omit(DS3)
    # print("from line 1894 - heatmap de type")
    # print("DS3")
    # print(head(DS3))
    
    set.seed(110)
    a <- ComplexHeatmap::Heatmap(DS3,
                                 name = "Normalized expression",
                                 col = colorRamp2(c(min(DS3), 0, max(DS3)), c("red", "black", "green")),
                                 row_names_gp = gpar(fontsize = 1),
                                 row_dend_gp = gpar(fontsize = 1),
                                 row_title_gp = gpar(fontsize = 10),
                                 cluster_columns = FALSE,
                                 row_dend_width = unit(3, "cm"),
                                 split = clusterNum, clustering_distance_rows = "pearson",
                                 show_heatmap_legend = TRUE,
                                 show_row_names = FALSE, show_column_names = T,
                                 heatmap_legend_param = list(title = "Normalized expression")
    )
    set.seed(110)
    rcl.list <- row_order(a)
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
    for (i in 1:length(rcl.list)) {
      rcl.list2[[i]] <- rownames(DS3)[rcl.list[[i]]]
    }
    
    for (i in 1:length(rcl.list2)) {
      genes <- rcl.list2[[i]]
      group_name <- rep(i, length(genes))
      Cluster_i <- data.frame("GeneID" = genes, "Cluster" = i)
      if (i == 1) {
        Cluster <- Cluster_i
      } else {
        Cluster <- rbind(Cluster, Cluster_i)
      }
    }
    
    print("line 2089, Cluster")
    print(head(Cluster))
    print(paste0("de_type = ", de_type))
    print("length of input DS")
    print(dim(DS3))
    
    # end heat map analysis
    heatmap.end.time <- Sys.time()
    print("heat map time")
    print(heatmap.end.time - heatmap.start.time)
    return(list(a, DS3, Cluster))
  })
  
  getCluster <- eventReactive(input$heatmap_plot, {
    set.seed(110)
    ll <- plotHeatmap()
    a <- ll[[1]]
    DS3 <- ll[[2]]
    rcl.list <- row_order(a)
    DS3.1 <- as.matrix(rownames(DS3))
    
    Cluster <- NULL
    for (i in 1:length(rcl.list)) {
      for (j in 1:length(rcl.list[[i]])) {
        pair <- c(i, DS3.1[rcl.list[[i]][j]])
        Cluster <- rbind(Cluster, pair)
      }
    }
    Cluster <- data.frame(Cluster, row.names = NULL)
    colnames(Cluster) <- c("cluster", "gene.id")
    
    return(Cluster[, c("gene.id", "cluster")])
  })
  
  mapPlot <- function() {
    myHeatmap <- plotHeatmap()[[1]]
    myHeatmap <- draw(myHeatmap)
  }
  
  output$heatmap.plot <- renderPlot({
    mapPlot()
  })
  
  output$downloadheatmap <- downloadHandler(
    filename = function() {
      paste("heatmap", ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      p <- mapPlot()
      dev.off()
    }
  )
  
  
  output$cluster.info <- DT::renderDataTable({
    clusternum <- input$display_cluster
    gl <- plotHeatmap()[[3]] # getCluster()
    if (!is.null(gl)) {
      if (clusternum == "ALL") {
        gl
      } else {
        clusternum <- as.numeric(clusternum)
        dplyr::filter(gl, cluster == clusternum)
      }
    }
  })
  
  output$downloadclusters <- downloadHandler(
    filename = function() {
      paste("genelist", ".csv", sep = "")
    },
    content = function(file) {
      gl <- plotHeatmap()[[3]]
      write.csv(gl, file, row.names = FALSE)
    }
  )
  
  output$help_text_heatmap <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          Hierarchical clustering is used to find the groups of co-expressed genes. 
          The clustering is performed on normalized expressions of differentially expressed genes using Ward clustering method. 
          Normalized expression of the jth gene at time ti is defined as<br>
          <img src='https://i.ibb.co/tJgCMVD/Screenshot-from-2020-09-05-16-37-39.png' alt='Screenshot-from-2020-09-05-16-37-39' border='0'><br>
          where <var>x<sub>j</sub>(t<sub>i</sub>)</var> is the expression of the j<sup>th</sup> gene at time t<sub>i</sub>, 
          <var>x<sub>j</sub>&#772</var> is the mean expression across all time points, and <var>&#963<sub>j</sub></var> is the standard deviation.
          </b>
        </p>
      </center>
    ")
  })
  
  ############################
  ########## noise ###########
  ############################
  
  SQCO <- function(MT) {
    if (ncol(MT) == 1) {
      res <- matrix(0)
      return(res)
    }
    temp <- NULL
    for (i in 1:nrow(MT)) {
      m <- sum(MT[i, ]) / length(MT[i, ])
      v <- stats::var(MT[i, ])
      if (m != 0) {
        temp <- c(temp, v / (m * m))
      }
    }
    res <- matrix(sum(temp) / length(temp))
    return(res)
  }
  
  output$expand_genonames_noise <- renderUI({
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    numOfRep <- as.numeric(input$noise_numOfRep)
    numOfGeno <- ncol(DS) / numOfRep
    
    lapply(1:numOfGeno, function(i) {
      textInput(paste("noisetype", i, sep = ""), paste("Type", i, sep = " "), value = colnames(DS)[(i - 1) * numOfRep + 1])
    })
  })
  
  output$noise_anchor_choices <- renderUI({
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    numOfRep <- as.numeric(input$noise_numOfRep)
    numOfGeno <- ncol(DS) / numOfRep
    names <- NULL
    for (i in 1:numOfGeno) {
      id <- paste("noisetype", i, sep = "")
      names <- c(names, input[[id]])
    }
    selectInput("noise_anchor_b", "Anchor genotype", choices = names)
  })
  
  noisePlot <- eventReactive(input$noise_plot, {
    noise.start.time <- Sys.time()
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    numOfRep <- as.numeric(input$noise_numOfRep)
    numOfGeno <- ncol(DS) / numOfRep
    graph <- input$noise_graph_type
    names <- NULL
    for (i in 1:numOfGeno) {
      id <- paste("noisetype", i, sep = "")
      names <- c(names, input[[id]])
    }
    
    situation <- input$noise_situation
    if (situation == "a") {
      DS1 <- list()
      for (j in 1:numOfGeno) {
        DS1[[j]] <- as.matrix(DS[, ((j - 1) * numOfRep + 1):(j * numOfRep)])
      }
      Noise <- NULL
      for (y in 1:numOfGeno) {
        Noise <- c(Noise, SQCO(DS1[[y]]))
      }
      xform <- list(
        categoryorder = "array",
        categoryarray = names
      )
      if (graph == "Bar chart") {
        p <- plot_ly(
          x = names,
          y = Noise,
          type = "bar"
        ) %>% layout(xaxis = xform)
      } else if (graph == "Line chart") {
        p <- plot_ly(
          x = names,
          y = Noise,
          type = "scatter",
          mode = "lines+markers"
        ) %>% layout(xaxis = xform, yaxis = list(range = c(0, max(Noise) + 0.001)))
      }
    } else if (situation == "b") {
      DS_ave <- NULL
      for (j in 1:numOfGeno) {
        part_DS <- as.matrix(DS[, ((j - 1) * numOfRep + 1):(j * numOfRep)])
        DS_ave <- cbind(DS, data.frame(matrixStats::rowMeans2(part_DS)))
      }
      anchor <- input$noise_anchor_b
      names <- NULL
      for (i in 1:numOfGeno) {
        id <- paste("noisetype", i, sep = "")
        names <- c(names, input[[id]])
      }
      anchor_index <- match(anchor, names)
      Noise <- NULL
      for (i in 1:numOfGeno) {
        if (i != anchor_index) {
          Noise <- c(Noise, SQCO(cbind(DS_ave[, anchor_index], DS_ave[, i])))
        }
      }
      names_wo_anchor <- NULL
      for (i in 1:numOfGeno) {
        if (i != anchor_index) {
          id <- paste("noisetype", i, sep = "")
          names_wo_anchor <- c(names_wo_anchor, input[[id]])
        }
      }
      xform <- list(
        categoryorder = "array",
        categoryarray = names_wo_anchor
      )
      if (graph == "Bar chart") {
        p <- plot_ly(
          x = names_wo_anchor,
          y = Noise,
          type = "bar"
        ) %>% layout(xaxis = xform)
      } else if (graph == "Line chart") {
        p <- plot_ly(
          x = names_wo_anchor,
          y = Noise,
          type = "scatter",
          mode = "lines+markers"
        ) %>% layout(xaxis = xform, yaxis = list(range = c(0, max(Noise) + 0.001)))
      }
    } else if (situation == "c") {
      anchor <- input$noise_anchor_c
      names <- colnames(DS)
      anchor_index <- match(anchor, names)
      Noise <- NULL
      for (i in 1:ncol(DS)) {
        if (i != anchor_index) {
          Noise <- c(Noise, SQCO(cbind(DS[, anchor_index], DS[, i])))
        }
      }
      names_wo_anchor <- names[-anchor_index]
      xform <- list(
        categoryorder = "array",
        categoryarray = names_wo_anchor
      )
      if (graph == "Bar chart") {
        p <- plot_ly(
          x = names_wo_anchor,
          y = Noise,
          type = "bar"
        ) %>% layout(xaxis = xform)
      } else if (graph == "Line chart") {
        p <- plot_ly(
          x = names_wo_anchor,
          y = Noise,
          type = "scatter",
          mode = "lines+markers"
        ) %>% layout(xaxis = xform, yaxis = list(range = c(0, max(Noise) + 0.001)))
      }
    }
    noise.end.time <- Sys.time()
    print("noise time")
    print(noise.end.time - noise.start.time)
    return(p)
  })
  
  output$noise.plot <- renderPlotly({
    noisePlot()
  })
  
  output$downloadnoise <- downloadHandler(
    filename = function() {
      paste("noise", ".png", sep = "")
    },
    content = function(file) {
      p <- noisePlot()
      export(p, file = "noise.png")
    }
  )
  
  output$help_text_Noise <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          To quantify between gene expressions scatter of all replicates in one experimental condition, we 
          computed transcriptome-wide average noise for each cell type, defined as<br>
          <img src='https://i.ibb.co/kcxrCzv/Screenshot-from-2020-09-05-16-14-27.png' alt='Screenshot-from-2020-09-05-16-14-27' border='0'><br>
          where <var>n</var> is the number of genes and <var>n<sub>i</sub><sup>2</sup></var> is the pairwise noise of the i<sup>th</sup>
          gene (variability between any two replicates), defined as<br>
          <img src='https://i.ibb.co/dp21hsK/Screenshot-from-2020-09-05-16-14-57.png' alt='Screenshot-from-2020-09-05-16-14-57' border='0'><br>
          where <var>m</var> is the number of replicates in each condition and <var>n<sub>ijk</sub><sup>2</sup></var> is the expression noise of the i<sup>th</sup> gene, 
          defined by the variance divided by the squared mean expression in the pair of replicates (j,k).<br>
          Citation: <a href='https://www.nature.com/articles/srep07137'>https://www.nature.com/articles/srep07137</a> (Kumar’s embryonic development paper)
          </b>
        </p>
      </center>
    ")
  })
  
  ############################
  ######### entropy ##########
  ############################
  
  computeBin <- function(arr) { # Doane's rule
    n <- length(arr)
    gx <- moments::skewness(arr)
    sigmag <- sqrt(6 * (n - 2) / ((n + 1) * n + 3))
    bin <- 1 + log2(n) + log2(1 + abs(gx) / sigmag)
    return(bin)
  }
  
  getBinCounts <- function(arr) {
    vec <- entropy::discretize(arr, computeBin(arr), r = range(arr))
    return(vec)
  }
  
  output$expand_genonames_entropy <- renderUI({
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    tp <- as.numeric(input$entropy_timepoints)
    numOfGeno <- ncol(DS) / tp
    
    lapply(1:numOfGeno, function(i) {
      textInput(paste("entropytype", i, sep = ""), paste("Type", i, sep = " "), value = colnames(DS)[(i - 1) * tp + 1])
    })
  })
  
  entropyPlot <- reactive({
    entropy.start.time <- Sys.time()
    type <- input$file_type
    if (type == "norm") {
      DS <- df_shiny()
    } else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    if (is.null(DS) == FALSE) {
      tsflag <- input$tsflag
      graph <- input$entropy_graph_type
      names <- colnames(DS)
      xform <- list(
        categoryorder = "array",
        categoryarray = names
      )
      entropy.vector <- NULL # entropy of each column
      for (i in 1:length(DS)) {
        binCount <- getBinCounts(DS[, i])
        entropy <- entropy.empirical(binCount, unit = "log2")
        entropy.vector <- c(entropy.vector, entropy)
      }
      if (tsflag == FALSE) {
        if (graph == "Bar chart") {
          p <- plot_ly(
            x = names,
            y = entropy.vector,
            type = "bar"
          ) %>% layout(xaxis = xform)
        } else if (graph == "Line chart") {
          p <- plot_ly(
            x = names,
            y = entropy.vector,
            type = "scatter",
            mode = "lines+markers"
          ) %>% layout(xaxis = xform)
          # yaxis=list(range = c(0, max(ent)+0.002))
        }
      } else if (tsflag == TRUE) {
        tp <- as.numeric(input$entropy_timepoints)
        numOfGeno <- ncol(DS) / tp
        names <- NULL
        for (i in 1:numOfGeno) {
          id <- paste("entropytype", i, sep = "")
          names <- c(names, input[[id]])
        }
        time_index <- c(1:tp)
        ent <- data.frame(time_index)
        for (j in 1:numOfGeno) {
          part_ent <- entropy.vector[(tp * j - (tp - 1)):(tp * j)]
          ent <- cbind(ent, part_ent)
        }
        if (graph == "Bar chart") {
          p <- plot_ly(x = ent[, 1], y = ent[, 2], name = names[1], type = "bar")
          for (i in 1:(numOfGeno - 1)) {
            p <- add_trace(p, y = ent[, i + 2], name = names[i + 1], type = "bar")
          }
          p <- layout(p, xaxis = list(title = "Time"), yaxis = list(title = "Entropy"))
        } else if (graph == "Line chart") {
          p <- plot_ly(x = ent[, 1], y = ent[, 2], name = names[1], type = "scatter", mode = "lines+markers")
          for (i in 1:(numOfGeno - 1)) {
            p <- add_trace(p, y = ent[, i + 2], name = names[i + 1], type = "scatter", mode = "lines+markers")
          }
          p <- layout(p, xaxis = list(title = "Time"), yaxis = list(title = "Entropy"))
        }
      }
      entropy.end.time <- Sys.time()
      print("entropy time")
      print(entropy.end.time - entropy.start.time)
      return(p)
    }
  })
  
  output$entropy.plot <- renderPlotly({
    entropyPlot()
  })
  
  output$downloadentropy <- downloadHandler(
    filename = function() {
      paste("entropy", ".png", sep = "")
    },
    content = function(file) {
      p <- entropyPlot()
      export(p, file = "entropy.png")
    }
  )
  
  output$help_text_Entropy <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          Shannon entropy (Shannon, 1948) measures the disorder of a high-dimensional
          system, where higher values indicate increasing disorder.  Entropy of each transcriptome, X, is defined as<br>
          <img src='https://i.ibb.co/5W0KwMP/Screenshot-from-2020-09-05-16-31-03.png' alt='Screenshot-from-2020-09-05-16-31-03' border='0'><br>
          where p(xi) is the probability of gene expression value x=xi.
          </b>
        </p>
      </center>
    ")
  })
  
  ###################################
  ########## New-Features ###########
  ###################################
  ###################################
  
  
  ###################################
  ###################################
  ############   SVM     ############
  ###################################
  ###################################
  
  # data_svm <- eventReactive(input$submit_svm, {
  #   print("Running...")
  #   svm.start <- Sys.time()
  #   type <- input$file_type
  #   
  #   if (type == "norm") {
  #     DS <- df_shiny()
  #   }
  #   else if (type == "raw") {
  #     DS <- df_raw_shiny()
  #   }
  #   
  #   print("1")
  #   
  #   x1.1 <- as.data.frame(DS[,1:3])
  #   x1.11 <- as.data.frame(DS[1:100,1:3])
  #   x2.2 <- as.data.frame(DS[,4:6])
  #   x2.22 <- as.data.frame(DS[1:100,4:6])
  #   
  #   print("2")
  #   
  #   x1.1 <-  setNames(stack(x1.1),c("x1.1","colName"))
  #   x1.11 <-  setNames(stack(x1.11),c("x1.1","colName"))
  #   x2.2 <-  setNames(stack(x2.2),c("x2.2","colName"))
  #   x2.22 <-  setNames(stack(x2.22),c("x2.2","colName"))
  #   
  #   mut_type <- as.data.frame(x1.1[,2])
  #   mut_type1 <- as.data.frame(x1.11[,2])
  #   
  #   dat <- data.frame(matrix(ncol = 3, nrow = nrow(x1.1)))
  #   dat1 <- data.frame(matrix(ncol = 3, nrow = nrow(x1.11)))
  #   x <- c("x1.1", "x2.2", "y")
  #   colnames(dat) <- x
  #   colnames(dat1) <- x
  #   
  #   print("3")
  #   
  #   dat[,1] <- x1.1[,1]
  #   dat[,2] <- x2.2[,1]
  #   dat[,3] <- mut_type[,1]
  #   
  #   print("4")
  #   
  #   dat1[,1] <- x1.11[,1]
  #   dat1[,2] <- x2.22[,1]
  #   dat1[,3] <- mut_type1[,1]
  #   
  #   dat[,3] <- as.factor(as.numeric(dat[,3]))
  #   dat1[,3] <- as.factor(as.numeric(dat1[,3]))
  #   
  #   return(dat1)
  #   
  # })
  # 
  # plotSVM <- function() {
  #   
  #   dat <- data_svm()
  #   
  #   print("training")
  #   svmfit <- svm(y~., data = dat, kernel = "radial", cost = 10, gamma = 1)
  #   print("done")
  #   plot(svmfit , dat )
  #   
  # }
  # 
  # plotSVM_df <- function() {
  #   
  #   dat <- data_svm()
  #   
  #   ggplot(data = dat, aes(x = x2.2, y = x1.1, color = y, shape = y)) + 
  #     geom_point(size = 2) +
  #     scale_color_manual(values=c("#000000","#FF0000","#00BA00")) +
  #     theme(legend.position = "none")
  #   
  # }
  # 
  # output$svm_plot <- renderPlot({
  #   plotSVM()
  # })
  # 
  # output$svm_df_plot <- renderPlot({
  #   plotSVM_df()
  # })
  # 
  # output$help_text_SVM <- renderUI({
  #   HTML("<h3><b>To be implemented</b></h3>")
  # })
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  
  
  ###################################
  ###################################
  ############   t-SNE    ###########
  ############# Python ##############
  ###################################
  # data for t-sne
  plotTSNE2 <- eventReactive(input$submit_tsne2, {
    
    tsne2_trans <- input$tsne2_trans
    type <- input$file_type
    perplexity_value <- input$perplexity_value
    no_of_pca <- input$no_of_pca
    # no_of_clusters <- input$no_of_clusters
    
    
    if(type=='norm'){
      DS <- df_shiny()
    }else if(type=='raw'){ 
      DS <- df_raw_shiny()
    }
    if(tsne2_trans=='None'){
      tsne2.data <- t(apply(DS, MARGIN = 1, scale)); 
      colnames(tsne2.data) <- colnames(DS)
    }
    else if(tsne2_trans=='log10'){
      tsne2.data <- log10(DS+1)
    }
    
    tsne_cluster_flag <- input$tsne_cluster # ture or false
   
    return (list(tsne2.data, perplexity_value, no_of_pca, tsne_cluster_flag)) #, no_of_clusters
  })
  
  
  tsne2plot <- function(){
    set.seed(13)
    tsne2.start <- Sys.time()
    # get data 
    li <- plotTSNE2()
    tsne2.data <- li[[1]]
    perplexity_value <- li[[2]]
    no_of_pca <- li[[3]]
    tsne_cluster_flag <- li[[4]]
    tsne_text_flag <- input$tsne_text # display sample name or not
    # no_of_clusters <- li[[4]] 
    
    # get tsne value
    tsne_val <- Rtsne(t(tsne2.data), 
                      dims = 2,
                      initial_dims = no_of_pca,
                      perplexity = perplexity_value,
                      theta = 0.0)
    tsne_df <- data.frame(
      TSNE1 = tsne_val$Y[, 1],
      TSNE2 = tsne_val$Y[, 2],
      Sample = colnames(tsne2.data)
    )
    
    if(!tsne_cluster_flag){
      # plotting
      p <- plot_ly(data = tsne_df, x = ~TSNE1, y = ~TSNE2, text = ~Sample) %>% 
        add_trace(type = "scatter", mode = 'markers', opacity = 0.5)

    } else { # tsne_cluster_flag == TRUE
      set.seed(13)
      tsne_cluster_num <- as.numeric(input$tsne_cluster_num)
      tsne_kmeans_result <- kmeans(tsne_df[,1:2], tsne_cluster_num)
      tsne_df$cluster <- factor(tsne_kmeans_result$cluster, levels = 1:max(tsne_kmeans_result$cluster) )
      
      # plotting
      p <- plot_ly(data = tsne_df, x = ~TSNE1, y = ~TSNE2, text = ~Sample, color = ~cluster ) %>%
        add_trace(type = "scatter", mode = 'markers', opacity = 0.5)
    }
    
    if(tsne_text_flag){
      p <- p %>% hide_colorbar() %>%
        add_trace(type = "scatter", mode = 'text', textposition = "top right", showlegend = FALSE)
    }
    
    tsne2.end <- Sys.time()
    print("t-SNE plot time")
    print(tsne2.end - tsne2.start)
    
    return(list(p, tsne_df))
  }
  
  output$tsne2.plot <- renderPlotly({
    li <- tsne2plot()
    p <- li[[1]]
    tsne_table <- li[[2]]
    p
  })
  
  output$tsne_table <- DT::renderDataTable({
    tsne_table <- tsne2plot()[[2]] # get table
    tsne_table
  })
  
  output$download_tsne <- downloadHandler(
    filename = function() {
      paste("tsne_list", ".csv", sep = "")
    },
    content = function(file) {
      gl <- tsne2plot()[[2]]
      write.csv(gl, file, row.names = FALSE)
    }
  )
  
  output$help_text_tsne <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          t-SNE (t-distributed stochastic neighbor embedding) (Hinton, 2008) uses the local relationships between points 
          to create a low-dimensional mapping. This allows it to capture non-linear structure. </b>
          t-SNE creates a probability distribution using the Gaussian distribution that defines 
          the relationships between the points in high-dimensional space. t-SNE uses the Student 
          t-distribution to recreate the probability distribution in low-dimensional space. This 
          prevents the crowding problem, where points tend to get crowded in low-dimensional 
          space due to the curse of dimensionality. </b>
          t-SNE optimizes the embeddings directly using gradient descent. The cost function is 
          non-convex, thus there is the risk of getting stuck in local minima. t-SNE uses multiple 
          tricks to try to avoid this problem.
          </b>
        </p>
      </center>
    ")
  })
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  
  
  ###################################
  ###################################
  #######   Random Forest    ########
  ###################################
  ###################################
  # data for random forest
  plotRF <- eventReactive(input$submit_rf, {
    rf.start <- Sys.time()
    rf_trans <- input$rf_trans
    type <- input$file_type
    num_trees <- input$num_trees
    num_clusters <- input$num_clusters
    
    if (type == "norm") {
      DS <- df_shiny()
    }
    else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    if (rf_trans == "None") {
      rf.data <- DS
    }
    else if (rf_trans == "log10") {
      rf.data <- log10(DS + 1)
    }
    rf.end <- Sys.time()
    print("Random forest plot time")
    print(rf.end - rf.start)
    return(list(rf.data, num_trees, num_clusters))
  })
  
  
  plotRAFSIL <- eventReactive(input$submit_rafsil, {
    rf.start <- Sys.time()
    rf_trans <- input$rf_trans
    type <- input$file_type
    
    if (type == "norm") {
      DS <- df_shiny()
    }
    else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    if (rf_trans == "None") {
      rf.data <- DS
    }
    else if (rf_trans == "log10") {
      rf.data <- log10(DS + 1)
    }
    f <- group_names()
    if (!is.null(f)) {
      meta_df <- data.frame("Column names" = colnames(DS), "Description" = f)
      meta_df
    }
    
    meta_df <- meta_df %>% remove_rownames %>% column_to_rownames(var="Column.names")
    meta_df$Description <- as.numeric(as.factor(meta_df$Description))
    
    meta_df <- as.matrix(meta_df)
    DS <- as.matrix(DS)
    
    rf.end <- Sys.time()
    print("RFSIL plot time")
    print(rf.end - rf.start)
    return(list(meta_df,DS))
  })
  
  rafsilplot <- function() {
    
    tryCatch({
      # get data
      t_list <- plotRAFSIL()
      ord = order(t_list[[1]]) ; t_list[[2]]=t_list[[2]][,ord] ; t_list[[1]] = t_list[[1]][ord] ; rm(ord)
      
      #- run RAFSIL1 with 50 forests
      res.r1 = RAFSIL(t(t_list[[2]]),nrep = 50, method="RAFSIL1")
      res.r2 = RAFSIL(t(t_list[[2]]),           method="RAFSIL2")
      
      #- retriev the dissimilarities
      dis.r1  = res.r1$D
      dis.r2  = res.r2$D
      dis.cor = sqrt((1 - cor(t_list[[2]],method="spearman"))/2)
      
      par(mfrow=c(1,2))
      par(mai=c(.1,.1,.5,.1))
      plotTSNE(dis.r1,labels=t_list[[1]],is_distance=FALSE,verbose=TRUE,perplexity=5)
      mtext("rafsil-1 / embedding", line=1)
      plotTSNE(dis.r2,labels=t_list[[1]],is_distance=FALSE,verbose=TRUE,perplexity=5)
      mtext("rafsil-2 / embedding", line=1)
    }, error = function(error_condition) {
      plot_exception("RAFSIL cannot be applied on this dataset.\nPlease use random forest clustering instead")
    }) 
    
  }
  
  ####################################################################################
  
  plot_exception <-function(
    ...,
    sep=" ",
    type=c("message","warning","cat","print"),
    color="auto",
    console=TRUE,
    size = 6){      
    type=match.arg(type)
    txt = paste(...,collapse=sep)
    if(console){
      if(type == "message") message(txt)
      if(type == "warning") warning(txt)
      if(type == "cat") cat(txt)
      if(type == "print") print(txt)
    }
    if(color =="auto") color <- if(type == "cat") "black" else "red"
    if(txt == "warning") txt <- paste("warning:",txt)
    print(ggplot2::ggplot() +
            ggplot2::geom_text(ggplot2::aes(x=0,y=0,label=txt),color=color,size=size) + 
            ggplot2::theme_void())
    invisible(NULL)
  }
  
  ####################################################################################
  
  
  rfplot <- function() {
    # get data
    li <- plotRF()
    rf.data <- li[[1]]
    num_trees <- li[[2]]
    num_clusters <- li[[3]]
    
    # unsupervised random forest on data
    print("Running random forest...")
    rf.data <- t(rf.data)
    rf_out <- randomForest(rf.data, type = unsupervised, ntree = num_trees, proximity = TRUE)
    print("Done!")
    
    mds_out <- cmdscale(1 - rf_out$proximity, eig = TRUE, k = 2)
    clusters_pam <- pam(1 - rf_out$proximity, k = num_clusters, diss = TRUE)
    shape_lvl <- c(1:num_clusters)
    shape_legend <- factor(clusters_pam$clustering, levels = shape_lvl)
    
    # print the proximity matrix
    print(table(clusters_pam$clustering, shape_legend))
    print(str(mds_out))
    print(mds_out$points)
    print(rownames(mds_out$points))
    
    # plot the graph
    df <- data.frame(x = mds_out$points[, 1], y = mds_out$points[, 2], color = shape_legend, shape = shape_legend)
    p <- ggplot(data = df, aes(x = x, y = y, color = color, shape = shape, text = paste("x: ", round(x, 4), "\n", "y: ", round(y, 4), "\n", "Name: ", rownames(mds_out$points), "\n", "Cluster: ", shape, sep = ""), group = 1)) +
      geom_point(size = 1.40) + theme_bw()
    p <- p + theme(legend.position = "none")
    p
    
    # add interactivity w/ plotly
    ggplotly(p, tooltip = c("text"))
  }
  
  
  rf_matrix <- reactive({
    li <- plotRF()
    rf.data <- li[[1]]
    num_trees <- li[[2]]
    num_clusters <- li[[3]]
    
    # unsupervised random forest on data
    print("Running random forest...")
    rf.data <- t(rf.data)
    rf_out <- randomForest(rf.data, type = unsupervised, ntree = num_trees, proximity = TRUE)
    print("Done!")
    
    mds_out <- cmdscale(1 - rf_out$proximity, eig = TRUE, k = 2)
    
    return(mds_out$points)
  })
  
  
  output$rf.plot <- renderPlotly({
    rfplot()
  })
  
  output$RAFSIL.plot <- renderPlot({
    rafsilplot()
  })
  
  output$rf.matrix <- renderTable({
    rf_matrix()
  },rownames=TRUE)
  
  output$help_text_rf <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          Clustering belongs to unsupervised learning, in which each sample is clustered 
          into different classes, based on their similarity (usually based on Euclidean distance).
          Random forest algorithm is used to generate a proximity matrix - a rough estimate of the 
          distance between samples based on the proportion of times the samples end up in the same 
          leaf node of the decision tree. The proximity matrix is converted to a dist matrix which 
          is then input to the hierarchical clustering algorithm.
          </b>
        </p>
        <p>
          <b>
          Implementation adapted from - <a href ='https://nishanthu.github.io/articles/ClusteringUsingRandomForest.html'>
          https://nishanthu.github.io/articles/ClusteringUsingRandomForest.html</a>
          </b>
        </p>
      </center>
    ")
  })
  
  # output$downloadrfplot <- downloadHandler(
  #   filename = function(){
  #     paste("randomforestplot",".pdf",sep="")
  #   },
  #   content = function(file){
  #     pdf(file) 
  #     rfplot()
  #     dev.off()
  #   }
  # )
  
  # output$downloadrfmatrix <- downloadHandler(
  #   filename = function(){
  #     paste("randomforestmatrix",".pdf",sep="")
  #   },
  #   content = function(file){
  #     pdf(file) 
  #     rf_matrix()
  #     dev.off()
  #   }
  # )
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  
  ###################################
  ###################################
  ############   SOM    #############
  ###################################
  ###################################
  # data for SOM
  plotSOM <- eventReactive(input$submit_som, {
    som.start <- Sys.time()
    som_trans <- input$som_trans
    sample_choice <- input$som_samples
    grid_h <- input$som_grid_h
    grid_v <- input$som_grid_v
    plot_type <- input$som_plot_type
    cluster_size <- input$som_cluster_size
    type <- input$file_type
    
    if (type == "norm") {
      DS <- df_shiny()
    }
    else if (type == "raw") {
      DS <- df_raw_shiny()
    }
    if (som_trans == "None") {
      som.data <- DS
    }
    else if (som_trans == "log10") {
      som.data <- log10(DS + 1)
    }
    
    # Use all samples or individual
    if (sample_choice == "All") {
      som.data <- som.data
    }
    else {
      som.data <- som.data[, sample_choice]
    }
    
    # some parameters
    som.data <- as.matrix(som.data)
    som_grid <- somgrid(xdim = grid_h, ydim = grid_v, topo = "hexagonal")
    som_model <- som(som.data, grid = som_grid)
    
    som.end <- Sys.time()
    print("SOM plot time")
    print(som.end - som.start)
    return(list(som_model, cluster_size))
  })
  
  
  sompropertyplot <- function() {
    # get data
    li <- plotSOM()
    som_model <- li[[1]]
    
    # plot type: property
    colors <- function(n, alpha = 'Set1') {
      rev(brewer.pal(n, alpha))
    }
    # use codes vectors (weight) for property plot
    plot(som_model, type = "property", property = getCodes(som_model), main = "Property", palette.name = colors)
  }
  
  somcountplot <- function() {
    li <- plotSOM()
    som_model <- li[[1]]
    
    # plot type: count
    colors <- function(n, alpha = 'Set2') {
      rev(brewer.pal(n, alpha))
    }
    
    # show how many genes are mapped to each node
    plot(som_model, type = "count", main = "Count", palette.name = colors)
  }
  
  somcodesplot <- function() {
    li <- plotSOM()
    som_model <- li[[1]]
    
    # plot type: codes
    # shows codebook vectors of genes
    plot(som_model, type = "codes", main = "Codes")
  }
  
  somdistplot <- function() {
    li <- plotSOM()
    som_model <- li[[1]]
    
    # plot type: distance
    colors <- function(n, alpha = 'Set3') {
      rev(brewer.pal(n, alpha))
    }
    
    # show how close genes are from each other when they are mapped
    plot(som_model, type = "dist.neighbours", main = "Distance", palette.name = colors)
  }
  
  somclusterplot <- function() {
    li <- plotSOM()
    som_model <- li[[1]]
    cluster_size <- li[[2]]
    
    # plot type: cluster
    colors <- function(n, alpha = 1) {
      rev(heat.colors(n, alpha))
    }
    
    # define colors from RColorBrewer
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
    col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    # use hierarchical clustering to cluster the SOM
    som.hc <- cutree(hclust(object.distances(som_model, "codes")), cluster_size)
    plot(som_model, type = "mapping", bgcol = col_vector[som.hc], main = "Clusters")
    add.cluster.boundaries(som_model, som.hc)
  }
  
  
  output$som_property.plot <- renderPlot({
    sompropertyplot()
  })
  output$som_count.plot <- renderPlot({
    somcountplot()
  })
  output$som_codes.plot <- renderPlot({
    somcodesplot()
  })
  output$som_dist.plot <- renderPlot({
    somdistplot()
  })
  output$som_cluster.plot <- renderPlot({
    somclusterplot()
  })
  
  output$downloadProperty <- downloadHandler(
    filename = function(){
      paste("SOMProperty",".pdf",sep="")
    },
    content = function(file){
      pdf(file) 
      sompropertyplot()
      dev.off()
    }
  )
  
  output$downloadCount <- downloadHandler(
    filename = function(){
      paste("SOMCount",".pdf",sep="")
    },
    content = function(file){
      pdf(file) 
      somcountplot()
      dev.off()
    }
  )
  
  output$downloadCodes <- downloadHandler(
    filename = function(){
      paste("SOMCodes",".pdf",sep="")
    },
    content = function(file){
      pdf(file) 
      somcodesplot()
      dev.off()
    }
  )
  
  output$downloadDistance <- downloadHandler(
    filename = function(){
      paste("SOMDistance",".pdf",sep="")
    },
    content = function(file){
      pdf(file) 
      somdistplot()
      dev.off()
    }
  )
  
  output$downloadCluster <- downloadHandler(
    filename = function(){
      paste("SOMCluster",".pdf",sep="")
    },
    content = function(file){
      pdf(file) 
      somclusterplot()
      dev.off()
    }
  )
  
  output$help_text_SOM <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          A self-organizing map (SOM) produces a two-dimensional, discretized representation 
          of the high-dimensional gene expression matrix, and is therefore a dimensionality 
          reduction technique. Self-organizing maps apply uses a neighborhood function to 
          preserve the topological properties of the input gene expression matrix.
          </b>
        </p>
        <p>
          <b>
          Each data point (1 sample) in the input gene expression matrix recognizes 
          themselves by competeting for representation. SOM mapping steps starts 
          from initializing the weight vectors.From there a sample vector is 
          selected randomly and the map of weight vectors is searched to find 
          which weight best represents that sample. Each weight vector has neighboring 
          weights that are close to it. The weight that is chosen is rewarded by being able 
          to become more like that randomly selected sample vector. The neighbors of that 
          weight are also rewarded by being able to become more like the chosen sample vector. 
          This allows the map to grow and form different shapes. Most generally, they form square/rectangular/hexagonal/L shapes in 2D feature space.
          </b>
        </p>
        <p>
          <b>
          Citation: <a href ='https://doi.org/10.1016/S0925-2312(98)00037-X'>https://doi.org/10.1016/S0925-2312(98)00037-X</a>
          </b>
        </p>
      </center>
    ")
  })
  
  ###################################
  ######## Gene-Set Analysis ########
  ###################################
  ###################################
  
  
  
  ###################################
  ###################################
  ###### Complex Enrichment ########
  ###################################
  ###################################
  
  download_com_table <- reactiveVal(0)
  
  
  df_complex <- reactive({
    print("running...")
    if (is.null(input$file_complex_prot)&& is.null(input$text_complex_prot)) {
      return(NULL)
    }
    else if(!is.null(input$file_complex_prot)){
    parts <- strsplit(input$file_complex_prot$datapath, ".", fixed = TRUE)
    type <- parts[[1]][length(parts[[1]])]
    if (type != "csv") {
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return(NULL)
    }
    
    Accessions <- read.csv(input$file_complex_prot$datapath)
    Accessions <- na.omit(Accessions)
    Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
    }
    else{
      Acessions<-strsplit(input$text_complex_prot," ")
      Accessions <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        Accessions<-rbind(Accessions,Acessions[[1]][x])
      }
      
      print(Accessions)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
      
    }
    
    return(Accessions)
    
  })
  
  
  df_com_table <- function(){
    
    gene_id <- df_com_id()
    
    corum_table <- data.frame()
    n <- 1
    
    for (id in gene_id) {
      check <- lookup(id, as.data.frame(up_corum_mapping))
      print(class(check))
      if(!is.na(check))
      {
        for (c_id in as.matrix(check)) {
          row_name <- paste0(id," (",as.character(lookup(as.character(id), as.data.frame(id_to_name), missing="No Match"))," )")
          c_row <- data.frame(Uniprot_id = sprintf('<a href="https://www.uniprot.org/uniprot/%s" class="btn btn-primary">%s</a>',id,row_name),
                              Corum_id = c_id,
                              Complex_Name = as.character(allComplexes[paste0(c_id),"Complex_Name"]),
                              Complex_comment = allComplexes[paste0(c_id),"Complex_comment"],
                              row.names = n)
          corum_table <- rbind(corum_table, c_row)
          n = n + 1
        }
        
      } else {
        c_row <- data.frame(Uniprot_id = paste0(id," (",as.character(lookup(as.character(id), as.data.frame(id_to_name), missing="No Match"))," )"),
                            Corum_id = "No Match",
                            Complex_Name = "No Match",
                            Complex_comment = "No Match",
                            row.names = n)
        corum_table <- rbind(corum_table, c_row)
        n = n + 1
      }
    }
    download_com_table(corum_table)
    return(corum_table)
    
  }
  
  output$complex_table_prot <- shiny::renderDataTable({
    df_com_table()
  }, escape = FALSE)
  
  
  df_com_id <- eventReactive(input$submit_complex_prot, {
    hide("help_text_complex_en")
    df <- df_complex()
    return(df)
  })
  
  
  output$complex_download_prot <- downloadHandler(
    filename = function() {
      paste("complex", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(download_com_table(), file, row.names = FALSE)
    }
  )
  
  output$help_text_complex_en_prot <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page performs a Complex enrichment through the <a href ='http://mips.helmholtz-muenchen.de/corum/'>CORUM database</a> 
          of a given set of UniProt accessions and links the results to <a href ='http://UniProt.org'>UniProt.org</a>.
          </b>
        </p>
      </center>
    ")
  })
  
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  
  ###################################
  ###################################
  ######## Protein Function #########
  ###################################
  ###################################
  
  
  download_prot_func <- reactiveVal(0)
  
  df_prot_func <- reactive({
    print("running...")
    if (is.null(input$file_prot_func) && is.null(input$text_prot_func)) {
      return(NULL)
    }
    else if(!is.null(input$file_prot_func)){
    parts <- strsplit(input$file_prot_func$datapath, ".", fixed = TRUE)
    type <- parts[[1]][length(parts[[1]])]
    if (type != "csv") {
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return(NULL)
    }
    
    Accessions <- read.csv(input$file_prot_func$datapath)
    Accessions <- na.omit(Accessions)
    Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
    }
    else{
      Acessions<-strsplit(input$text_prot_func," ")
      Accessions <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        Accessions<-rbind(Accessions,Acessions[[1]][x])
      }
      
      print(Accessions)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
      
    }
    return(Accessions)
    
  })
  
  df_func_table <- function() {
    
    Accessions <- df_func_id()
    
    print("fetching...")
    df <- GetProteinFunction(Accessions)
    
    count <- 1
    for(id in row.names(df))
    {
      row_name <- paste0(id," (",as.character(lookup(as.character(id), as.data.frame(id_to_name), missing="No Match"))," )")
      row.names(df)[count] <- sprintf('<a href="https://www.uniprot.org/uniprot/%s" class="btn btn-primary">%s</a>',id,row_name)
      count <- count + 1
    }
    
    print("fetched...")
    output_table <- data.frame()
    output_table <- data.frame(
      "ID" = row.names(df),
      "Function" = df[,"Function..CC."]
    )
    download_prot_func(output_table)
    return(output_table)
    
  }
  
  output$prot_func_table <- DT::renderDataTable({
    df_func_table()
  }, escape = FALSE)
  
  df_func_id <- eventReactive(input$submit_prot_func, {
    hide("help_text_prot_fn")
    df <- df_prot_func()
    return(df)
  })
  
  output$prot_func_download <- downloadHandler(
    filename = function() {
      paste("Protein-Function", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(download_prot_func(), file, row.names = FALSE)
    }
  )
  
  output$help_text_prot_fn <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page retrieves the Protein function information 
          from <a href ='http://UniProt.org'>UniProt.org</a> of a given set of UniProt accessions.
          </b>
        </p>
      </center>
    ")
  })
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  ###################################
  ###################################
  ####### Tissue Expression ########
  ###################################
  ###################################
  
  download_prot_expr <- reactiveVal(0)
  
  df_prot_expr <- reactive({
    print("running...")
    if (is.null(input$file_prot_expr)&& is.null(input$text_prot_expr)) {
      return(NULL)
    }else if(!is.null(input$file_prot_expr)){
    parts <- strsplit(input$file_prot_expr$datapath, ".", fixed = TRUE)
    type <- parts[[1]][length(parts[[1]])]
    if (type != "csv") {
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return(NULL)
    }
    
    Accessions <- read.csv(input$file_prot_expr$datapath)
    Accessions <- na.omit(Accessions)
    Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
    }else{
      Acessions<-strsplit(input$text_prot_expr," ")
      print(Acessions)
      Accessions <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        print(x)
        Accessions<-rbind(Accessions,Acessions[[1]][x])
      }
      
      print(Accessions)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
      
    }
    
    return(Accessions)
    
  })
  
  df_expr_table <- function() {
    
    Accessions <- df_expr_id()
    print("fetching...")
    df <- GetExpression(Accessions)
    
    count <- 1
    for(id in row.names(df))
    {
      row_name <- paste0(id," (",as.character(lookup(as.character(id), as.data.frame(id_to_name), missing="No Match"))," )")
      row.names(df)[count] <- sprintf('<a href="https://www.uniprot.org/uniprot/%s" class="btn btn-primary">%s</a>',id,row_name)
      count <- count + 1
    }
    
    print("fetched...")
    output_table <- data.frame()
    output_table <- data.frame(
      "ID" = row.names(df),
      "Tissue Specificity" = df[,"Tissue.specificity"]
    )
    
    download_prot_expr(output_table)                            
    return(output_table)
    
  }
  
  output$prot_expr_table <- DT::renderDataTable({
    df_expr_table()
  }, escape = FALSE)
  
  df_expr_id <- eventReactive(input$submit_prot_expr, {
    hide("help_text_prot_exp")
    df <- df_prot_expr()
    return(df)
  })
  
  output$prot_expr_download <- downloadHandler(
    filename = function() {
      paste("Protein-Expression", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(download_prot_expr(), file, row.names = FALSE)
    }
  )
  
  output$help_text_prot_exp <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page retrieves the Tissue Expression information 
          from <a href ='http://UniProt.org'>UniProt.org</a> of a given set of UniProt accessions.
          </b>
        </p>
      </center>
    ")
  })
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  ###################################
  ###################################
  #### Subcellular Localization #####
  ###################################
  ###################################
  
  
  download_prot_local <- reactiveVal(0)
  
  df_prot_local <- reactive({
    print("running...")
    if (is.null(input$file_prot_local) && is.null(input$text_prot_local)) {
      return(NULL)
    }
    else if(!is.null(input$file_prot_local)){
    parts <- strsplit(input$file_prot_local$datapath, ".", fixed = TRUE)
    type <- parts[[1]][length(parts[[1]])]
    if (type != "csv") {
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return(NULL)
    }
    
    Accessions <- read.csv(input$file_prot_local$datapath)
    Accessions <- na.omit(Accessions)
    Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
    }
    else{
      Acessions<-strsplit(input$text_prot_local," ")
      Accessions <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        Accessions<-rbind(Accessions,Acessions[[1]][x])
      }
      
      print(Accessions)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
      
    }
    return(Accessions)
    
  })
  
  df_local_table <- function() {
    
    Accessions <- df_local_id()
    print("fetching...")
    df <- GetSubcellular_location(Accessions)
    
    count <- 1
    for(id in row.names(df))
    {
      row_name <- paste0(id," (",as.character(lookup(as.character(id), as.data.frame(id_to_name), missing="No Match"))," )")
      row.names(df)[count] <- sprintf('<a href="https://www.uniprot.org/uniprot/%s" class="btn btn-primary">%s</a>',id,row_name)
      count <- count + 1
    }
    
    print("fetched...")
    output_table <- data.frame()
    output_table <- data.frame(
      "ID" = row.names(df),
      "Subcellular Location" = df[,"Subcellular.location..CC."]
    )
    download_prot_local(output_table)                            
    return(output_table)
    
  }
  
  output$prot_local_table <- DT::renderDataTable({
    df_local_table()
  }, escape = FALSE)
  
  df_local_id <- eventReactive(input$submit_prot_local, {
    hide("help_text_sub_loc")
    df <- df_prot_local()
    return(df)
  })
  
  output$prot_local_download <- downloadHandler(
    filename = function() {
      paste("Subcellular-Localization", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(download_prot_local(), file, row.names = FALSE)
    }
  )
  
  output$help_text_sub_loc <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page retrieves the Subcellular LocalizationSubcellular Localization information 
          from <a href ='http://UniProt.org'>UniProt.org</a> of a given set of UniProt accessions.
          </b>
        </p>
      </center>
    ")
  })
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  ###################################
  ###################################
  ######## Protein Domains ##########
  ###################################
  ###################################
  
  
  download_prot_domain <- reactiveVal(0)
  
  df_prot_domain <- reactive({
    print("running...")
    if (is.null(input$file_prot_domain)&& is.null(input$text_prot_domain)) {
      return(NULL)
    }
    else if(!is.null(input$file_prot_domain)){
    parts <- strsplit(input$file_prot_domain$datapath, ".", fixed = TRUE)
    type <- parts[[1]][length(parts[[1]])]
    if (type != "csv") {
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return(NULL)
    }
    
    Accessions <- read.csv(input$file_prot_domain$datapath)
    Accessions <- na.omit(Accessions)
    Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
    }
    else{
      Acessions<-strsplit(input$text_prot_domain," ")
      Accessions <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        Accessions<-rbind(Accessions,Acessions[[1]][x])
      }
      
      print(Accessions)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
      
    }
    return(Accessions)
    
  })
  
  df_domain_table <- function() {
    
    
    Accessions <- df_domain_id()
    print("fetching...")
    df <- GetFamily_Domains(Accessions)
    
    count <- 1
    for(id in row.names(df))
    {
      row_name <- paste0(id," (",as.character(lookup(as.character(id), as.data.frame(id_to_name), missing="No Match"))," )")
      row.names(df)[count] <- sprintf('<a href="https://www.uniprot.org/uniprot/%s" class="btn btn-primary">%s</a>',id,row_name)
      count <- count + 1
    }
    
    print("fetched...")
    output_table <- data.frame()
    output_table <- data.frame(
      "ID" = row.names(df),
      "Protein Families" = df[,"Protein.families"],
      "Protein Domain" = df[,"Domain..FT."]
    )
    download_prot_domain(output_table)                           
    return(output_table)
    
  }
  
  output$prot_domain_table <- DT::renderDataTable({
    df_domain_table()
  }, escape = FALSE)
  
  df_domain_id <- eventReactive(input$submit_prot_domain, {
    hide("help_text_pro_dom")
    df <- df_prot_domain()
    return(df)
  })
  
  output$prot_domain_download <- downloadHandler(
    filename = function() {
      paste("Protein-Domains", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(download_prot_domain(), file, row.names = FALSE)
    }
  )
  
  output$help_text_pro_dom <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page retrieves the Protein Domains information 
          from <a href ='http://UniProt.org'>UniProt.org</a> of a given set of UniProt accessions.
          </b>
        </p>
      </center>
    ")
  })
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  
  ###################################
  ###################################
  ###### Pathways Enrichment ########
  ###################################
  ###################################
  
  
  pathway_enri_df <- reactiveVal(0)
  pathway_enri_nodes <- reactiveVal(0)
  
  df_path_enri_gene <- reactive({
    print("running pathway gene...")
    if (is.null(input$file_path_enri_gene) && is.null(input$text_path_enri_gene) ) {
      return(NULL)
    }else if(!is.null(input$file_path_enri_gene)){
      parts <- strsplit(input$file_path_enri_gene$datapath, ".", fixed = TRUE)
      type <- parts[[1]][length(parts[[1]])]
      if (type != "csv") {
        showModal(modalDialog(
          title = "Error",
          "Please input a csv file!"
        ))
        return(NULL)
      }
      
      Accessions <- read.csv(input$file_path_enri_gene$datapath)
      print(Accessions)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
    }
    else{
     
      Acessions<-strsplit(input$text_path_enri_gene," ")
      print(Acessions)
      Accessions <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        print(x)
        Accessions<-rbind(Accessions,Acessions[[1]][x])
      }
      
      print(Accessions)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
      
    }
    
      return(Accessions)
      })
  
  
  df_path_enri_prot <- reactive({
    print("running...")
    if (is.null(input$file_path_enri_prot)&& is.null(input$text_path_enri_prot)) {
      return(NULL)
    }
    else if(!is.null(input$file_path_enri_prot)){
    parts <- strsplit(input$file_path_enri_prot$datapath, ".", fixed = TRUE)
    type <- parts[[1]][length(parts[[1]])]
    if (type != "csv") {
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return(NULL)
    }
    
    Accessions <- read.csv(input$file_path_enri_prot$datapath)
    print(Accessions)
    Accessions <- na.omit(Accessions)
    Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
    }
    else{
      print("In else")
      Acessions<-strsplit(input$text_path_enri_prot," ")
      print(Acessions)
      Accessions <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        print(x)
        Accessions<-rbind(Accessions,Acessions[[1]][x])
      }
      
      print(Accessions)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
      
    }
    return(Accessions)
    
  })
  
  
  df_path_enri_id_gene <- eventReactive(input$submit_path_enri_gene,{
    print("running")
    
    hide("help_text_path_enri")
    df <- df_path_enri_gene()

    return(df)
  })
  
  df_path_enri_id_prot <- eventReactive(input$submit_path_enri_prot,{
    print("running")
    
    hide("help_text_path_enri")
    df <- df_path_enri_prot()
    return(df)
  })
  
  
  output$help_text_path_enri_gene <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page performs Pathways Enrichment from for a given set of genes using
          <a href ='https://biit.cs.ut.ee/gprofiler/gost'>g:Profiler</a>.
          </b>
        </p>
      </center>
    ")
  })
  
  output$help_text_path_enri_prot <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page performs Pathways Enrichment from for a given set of genes using
          <a href ='https://biit.cs.ut.ee/gprofiler/gost'>g:Profiler</a>.
          </b>
        </p>
      </center>
    ")
  })
  
  
  plot_path_enri_gene <- function() {
    df_path_enri_id_gene()
    gene_name <- as.data.frame(df_path_enri_id_gene())
    gene_name[,1] <- as.character(gene_name[,1])
    
    path_list <- gost(gene_name[,1],exclude_iea = TRUE,evcodes = TRUE ,sources = "GO:BP")
    path_df <- path_list[[1]]
    pathway_enri_df(path_df)
    prot_num <- data.frame()
    for(i in 1:nrow(path_df))
    {
      prot_num <- rbind(prot_num,nrow(as.data.frame(strsplit(path_df[i,"intersection"],","))))
    }
    
    path_enrich_df <- data.frame(
      "term_name" = path_df[,"term_name"],
      "intersection" = prot_num[,1]
    )
    
    pathway_enri_nodes(path_enrich_df)
    
    path_enrich_df <- path_enrich_df[order(path_enrich_df$intersection),]
    
    
    bar_plot <- ggplot(data=path_enrich_df, aes(x=reorder(path_enrich_df$term_name , path_enrich_df$intersection), y=path_enrich_df$intersection)) +
      geom_bar(stat="identity", fill="steelblue" , alpha = 0.7) + xlab("Molecular function") + ylab("Number of Genes") +
      geom_text(aes(label = path_enrich_df$intersection), vjust = -0.03) + theme(axis.text.x = element_text(angle = 90 , hjust = 1 , vjust = 0.2))+
      theme_minimal() +coord_flip() + theme_bw()+theme(text = element_text(size=12, face="bold", colour="black"),axis.text.x = element_text(vjust=2))
    
    
    return(bar_plot)
  }
  
  plot_path_enri_prot <- function() {
    df_path_enri_id_prot()
    gene_name <- as.data.frame(df_path_enri_id_prot())
    gene_name[,1] <- as.character(gene_name[,1])
    
    path_list <- gost(gene_name[,1],exclude_iea = TRUE,evcodes = TRUE ,sources = "GO:BP")
    path_df <- path_list[[1]]
    pathway_enri_df(path_df)
    prot_num <- data.frame()
    for(i in 1:nrow(path_df))
    {
      prot_num <- rbind(prot_num,nrow(as.data.frame(strsplit(path_df[i,"intersection"],","))))
    }
    
    path_enrich_df <- data.frame(
      "term_name" = path_df[,"term_name"],
      "intersection" = prot_num[,1]
    )
    
    pathway_enri_nodes(path_enrich_df)
    
    path_enrich_df <- path_enrich_df[order(path_enrich_df$intersection),]
    
    
    bar_plot <- ggplot(data=path_enrich_df, aes(x=reorder(path_enrich_df$term_name , path_enrich_df$intersection), y=path_enrich_df$intersection)) +
      geom_bar(stat="identity", fill="steelblue" , alpha = 0.7) + xlab("Molecular function") + ylab("Number of Genes") +
      geom_text(aes(label = path_enrich_df$intersection), vjust = -0.03) + theme(axis.text.x = element_text(angle = 90 , hjust = 1 , vjust = 0.2))+
      theme_minimal() +coord_flip() + theme_bw()+theme(text = element_text(size=12, face="bold", colour="black"),axis.text.x = element_text(vjust=2))
    
    
    return(bar_plot)
  }
  
  
  output$path_enri.plot_gene <- renderPlotly({
    df_path_enri_id_gene()
  gene_name <- as.data.frame(df_path_enri_id_gene())
    gene_name[,1] <- as.character(gene_name[,1])
    
    ggplotly(Pathway.Enr(gene_name[,1]), tooltip = c("text"))
  })
  
  output$path_enri.plot_prot <- renderPlotly({
    df_path_enri_id_prot()
    gene_name <- as.data.frame(df_path_enri_id_prot())
    gene_name[,1] <- as.character(gene_name[,1])
    
    ggplotly(Pathway.Enr(gene_name[,1]), tooltip = c("text"))
  })
  
  #visualization
  
  observeEvent(input$fit_path_gene, ignoreInit=TRUE, {
    fit(session, 80)
  })
  
  
  observeEvent(input$fit_path_prot, ignoreInit=TRUE, {
    fit(session, 80)
  })
  
  
  observeEvent(input$showCondition_gene, ignoreInit=TRUE, {
    condition.name <- isolate(input$showCondition_gene)
    values <- as.numeric(pathway_enri_nodes()[,2])
    node.names <- pathway_enri_nodes()[,1]
    print(values)
    setNodeAttributes(session, attributeName="lfc", nodes=node.names, values)
  })
  
  
  observeEvent(input$showCondition_prot, ignoreInit=TRUE, {
    condition.name <- isolate(input$showCondition_prot)
    values <- as.numeric(pathway_enri_nodes()[,2])
    node.names <- pathway_enri_nodes()[,1]
    print(values)
    setNodeAttributes(session, attributeName="lfc", nodes=node.names, values)
  })
  
  
  observeEvent(input$loadStyleFile_path_gene,  ignoreInit=TRUE, {
    if(input$loadStyleFile_path != ""){
      tryCatch({
        loadStyleFile(input$loadStyleFile_path_gene)
      }, error=function(e) {
        msg <- sprintf("ERROR in stylesheet file '%s': %s", input$loadStyleFile_path_gene, e$message)
        showNotification(msg, duration=NULL, type="error")
      })
      later(function() {updateSelectInput(session, "loadStyleFile", selected=character(0))}, 0.5)
    }
  })
  
  
  
  observeEvent(input$loadStyleFile_path_prot,  ignoreInit=TRUE, {
    if(input$loadStyleFile_path_prot != ""){
      tryCatch({
        loadStyleFile(input$loadStyleFile_path_prot)
      }, error=function(e) {
        msg <- sprintf("ERROR in stylesheet file '%s': %s", input$loadStyleFile_path_prot, e$message)
        showNotification(msg, duration=NULL, type="error")
      })
      later(function() {updateSelectInput(session, "loadStyleFile", selected=character(0))}, 0.5)
    }
  })
  
  
  
  observeEvent(input$doLayout_path_gene,  ignoreInit=TRUE,{
    if(input$doLayout_path_gene != ""){
      strategy <- input$doLayout_path_gene
      doLayout(session, strategy)
      later(function() {updateSelectInput(session, "doLayout", selected=character(0))}, 1)
    }
  })
  
  observeEvent(input$doLayout_path_prot,  ignoreInit=TRUE,{
    if(input$doLayout_path_prot != ""){
      strategy <- input$doLayout_path_prot
      doLayout(session, strategy)
      later(function() {updateSelectInput(session, "doLayout", selected=character(0))}, 1)
    }
  })
  
  
  
  observeEvent(input$sfn_path_gene,  ignoreInit=TRUE,{
    selectFirstNeighbors(session)
  })
  
  observeEvent(input$sfn_path_prot,  ignoreInit=TRUE,{
    selectFirstNeighbors(session)
  })
  
  
  
  observeEvent(input$fitSelected_path_gene,  ignoreInit=TRUE,{
    fitSelected(session, 100)
  })
  
  observeEvent(input$fitSelected_path_prot,  ignoreInit=TRUE,{
    fitSelected(session, 100)
  })
  
  
  
  observeEvent(input$getSelectedNodes_path_gene, ignoreInit=TRUE, {
    output$selectedNodesDisplay_path_gene <- renderText({" "})
    getSelectedNodes(session)
  })
  
  observeEvent(input$getSelectedNodes_path_prot, ignoreInit=TRUE, {
    output$selectedNodesDisplay_path_prot <- renderText({" "})
    getSelectedNodes(session)
  })
  
  
  
  observeEvent(input$clearSelection_path_gene,  ignoreInit=TRUE, {
    clearSelection(session)
  })  
  
  observeEvent(input$clearSelection_path_prot,  ignoreInit=TRUE, {
    clearSelection(session)
  })  
  
  
  
  observeEvent(input$removeGraphButton_path_gene, ignoreInit=TRUE, {
    removeGraph(session)
  })
  
  observeEvent(input$removeGraphButton_path_prot, ignoreInit=TRUE, {
    removeGraph(session)
  })
  
  
  observeEvent(input$addRandomGraphFromDataFramesButton_path_gene, ignoreInit=TRUE, {
    source.nodes <-  LETTERS[sample(1:5, 5)]
    target.nodes <-  LETTERS[sample(1:5, 5)]
    tbl.edges <- data.frame(source=source.nodes,
                            target=target.nodes,
                            interaction=rep("generic", length(source.nodes)),
                            stringsAsFactors=FALSE)
    all.nodes <- sort(unique(c(source.nodes, target.nodes, "orphan")))
    tbl.nodes <- data.frame(id=all.nodes,
                            type=rep("unspecified", length(all.nodes)),
                            stringsAsFactors=FALSE)
    addGraphFromDataFrame(session, tbl.edges, tbl.nodes)
  })
  
  observeEvent(input$addRandomGraphFromDataFramesButton_path_prot, ignoreInit=TRUE, {
    source.nodes <-  LETTERS[sample(1:5, 5)]
    target.nodes <-  LETTERS[sample(1:5, 5)]
    tbl.edges <- data.frame(source=source.nodes,
                            target=target.nodes,
                            interaction=rep("generic", length(source.nodes)),
                            stringsAsFactors=FALSE)
    all.nodes <- sort(unique(c(source.nodes, target.nodes, "orphan")))
    tbl.nodes <- data.frame(id=all.nodes,
                            type=rep("unspecified", length(all.nodes)),
                            stringsAsFactors=FALSE)
    addGraphFromDataFrame(session, tbl.edges, tbl.nodes)
  })
  
  
  
  # observeEvent(input$selectedNodes, {
  #       newNodes <- input$selectedNodes;
  #       output$selectedNodesDisplay <- renderText({
  #          paste(newNodes)
  #          })
  #       })
  
  pathway_overlap <- reactiveVal(0)
  
  new_source_var <- reactiveVal(0)
  new_target_var <- reactiveVal(0)
  overlap_wt <- reactiveVal(0)
  
  output$path_enri_visu_gene <- renderCyjShiny({
    
    print("visualization")
    df_path_enri_id_gene()
    Enrich <- gost(df_path_enri_id_gene(),evcodes = T, sources = c('KEGG', 'REAC'))
    Pathway <- Construct.COPathway(Enrich, input$overlap_min_gene)
    nodes_tot <- c(unique(Pathway[,1],unique(Pathway[,2])))
    
    
    path_enri.nodes <- data.frame(id=nodes_tot,
                                  type=nodes_tot,
                                  stringsAsFactors=FALSE)
    
    path_enri.edges <- data.frame(source=Pathway[,1],
                                  target=Pathway[,2],
                                  interaction=Pathway[,1],
                                  stringsAsFactors=FALSE)
    
    graph.json <- dataFramesToJSON(path_enri.edges, path_enri.nodes)
    cyjShiny(graph=graph.json, layoutName="cola", styleFile = "./www/style/basicStyle.js")
    
  })

  Construct.COPathway <- function(EnrichmentObject, threshold = 1)
{
  
  PathwayNetwork <- data.frame()
  PathwayDF <- EnrichmentObject[["result"]] 
  for (i in 1:nrow(PathwayDF))
  {
    if (dim(PathwayDF)[1] == 1)
      break
    Pathway_accessions <- strsplit(PathwayDF$intersection, ",") 
    for (accession in Pathway_accessions)
    {
      inx <- which(grepl(accession, PathwayDF$intersection) == T)
      inx <- inx[-i] 
      Source <- rep(PathwayDF$term_name[i], length(inx))
      Target <- PathwayDF$term_name[c(inx)]
      PathwayNetwork <- rbind(PathwayNetwork , cbind(Source, Target, accession))
    }
    PathwayDF <- PathwayDF[-i,]
    CoEnrichment <- setNames(aggregate(PathwayNetwork$accession, by = list(PathwayNetwork$Source, PathwayNetwork$Target),
                                       paste, collapse=","), c("Source", "Target", "Accesion"))
    
    CoEnrichment$ProteinCount <- str_count(CoEnrichment$Accesion, ",")
    CoEnrichment <- CoEnrichment[CoEnrichment$ProteinCount >= threshold,]
  }
  return(CoEnrichment)
}
  
  output$path_enri_visu_prot <- renderCyjShiny({
    
    print("visualization")
    df_path_enri_id_prot()
    Enrich <- gost(df_path_enri_id_prot(),evcodes = T, sources = c('KEGG', 'REAC'))
    Pathway <- Construct.COPathway(Enrich, input$overlap_min_prot)
    nodes_tot <- c(unique(Pathway[,1],unique(Pathway[,2])))
    
    
    path_enri.nodes <- data.frame(id=nodes_tot,
                                  type=nodes_tot,
                                  stringsAsFactors=FALSE)
    
    path_enri.edges <- data.frame(source=Pathway[,1],
                                  target=Pathway[,2],
                                  interaction=Pathway[,1],
                                  stringsAsFactors=FALSE)
    
    graph.json <- dataFramesToJSON(path_enri.edges, path_enri.nodes)
    cyjShiny(graph=graph.json, layoutName="cola", styleFile = "./www/style/basicStyle.js")
    
  })
  
  observeEvent(input$edge_wt_gene, ignoreInit=TRUE, {
    condition.name <- isolate(input$showCondition_gene)
    
    print(overlap_wt())
    setEdgeAttributes(session, attributeName="wt", sourceNodes=new_source_var(),
                      targetNodes=new_target_var(),
                      interactions=new_target_var(),
                      values=overlap_wt())
  })
  
  observeEvent(input$edge_wt_prot, ignoreInit=TRUE, {
    condition.name <- isolate(input$showCondition_prot)
    
    print(overlap_wt())
    setEdgeAttributes(session, attributeName="wt", sourceNodes=new_source_var(),
                      targetNodes=new_target_var(),
                      interactions=new_target_var(),
                      values=overlap_wt())
  })
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  
  
  ###################################
  ###################################
  ##########   Uniprot    ###########
  ###################################
  ###################################
  
  df_uniprot <- reactive({
    print("running")
    if (is.null(input$file_uniprot) && is.null(input$text_uniprot)) {
      return(NULL)
    }
   else if(!is.null(input$file_uniprot)){
     parts <- strsplit(input$file_uniprot$datapath, ".", fixed = TRUE)
     type <- parts[[1]][length(parts[[1]])]
     if (type != "csv") {
       showModal(modalDialog(
         title = "Error",
         "Please input a csv file!"
       ))
       return(NULL)
     }
     
     Accessions <- read.csv(input$file_uniprot$datapath)
     Accessions <- na.omit(Accessions)[,1]
     Accessions <- unique(Accessions)
     Accessions <- trimws(Accessions)
     print(Accessions)
   }else{
     Acessions<-strsplit(input$text_uniprot," ")
     Accessions <- data.frame(Acessions[[1]][1])
     for (x in 2:length(Acessions[[1]])) {
       Accessions<-rbind(Accessions,Acessions[[1]][x])
     }
     
     print(Accessions)
     Accessions <- na.omit(Accessions)
     Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
    
   }
    return(Accessions)
    
  })
  
  plotUniprot <-  eventReactive(input$submit_uniprot, {
    
    Accessions <- df_uniprot()
    print(Accessions)
    hide("help_text_bio_pr")
    #print("Please Wait... Fetching Taxa Object. It may take a while")
    #TaxaObj <- GetNamesTaxa(Accessions)
    print("Please Wait... Fetching Gene Ontology Object. It may take a while")
    GeneOntologyObj <- GetProteinGOInfo(Accessions) 
    print("Done") 
    return(GeneOntologyObj)
  })
  
  output$help_text_bio_pr <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page retrieves the Gene Ontology (GO) terms 
          from <a href ='http://UniProt.org'>UniProt.org</a> of a given set of UniProt accessions.
          </b>
        </p>
      </center>
    ")
  })
  
  plotCE <- function() {
    # get data
    GO_df <- plotUniprot()
    return(Plot.GOSubCellular(GO_df,20))
    # ggplotly(bar_plot, tooltip = c("text"))
  }
  
  output$download_cell_plot <- downloadHandler(
    filename = function(){paste("Cellular-Component",'.png',sep='')},
    content = function(file){
      ggsave(file,plot=plotCE())
    }
  )
  
  plotBIO <- function() {
    # get data
    GO_df <- plotUniprot()
    return(PlotGOBiological(GO_df,20))
    # ggplotly(bar_plot, tooltip = c("text"))
  }
  
  output$download_bio_plot <- downloadHandler(
    filename = function(){paste("Biological-Process",'.png',sep='')},
    content = function(file){
      ggsave(file,plot=plotBIO())
    }
  )
  
  plotMol <- function() {
    # get data
    GO_df <- plotUniprot()
    return(Plot.GOMolecular(GO_df, 20))
    # ggplotly(bar_plot, tooltip = c("text"))
  }
  
  output$download_mole_plot <- downloadHandler(
    filename = function(){paste("Molecular-Function",'.png',sep='')},
    content = function(file){
      ggsave(file,plot=plotMol())
    }
  )
  
  output$uniprot_celplot <- renderPlot({
    GO_df <- plotUniprot()
    Plot.GOSubCellular(GO_df,20)
    ##ggplotly(plotCE(), tooltip = c("text"))
  })
  
  output$uniprotbioplot <- renderPlot({
    #ggplotly(plotBIO(), tooltip = c("text"))
    GO_df <- plotUniprot()
    PlotGOBiological(GO_df,20)
  })
  
  output$uniprot_molcplot <- renderPlot({
    GO_df <- plotUniprot()
    Plot.GOMolecular(GO_df,20)
    #ggplotly(plotMol(), tooltip = c("text"))
  })
  
  download_cel_table <- NULL
  
  output$uniprot_celtable <- shiny::renderDataTable({
    
    GO_df <- plotUniprot()
    CellularDF <- Goparse(GO_df, 5)
    CellularDF <- na.omit(CellularDF)
    download_cel_table <- CellularDF
    CellularDF
    
    
  })
  
  output$download_cell_comp <- downloadHandler(
    filename = function() {
      paste("Cellular-Component", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(download_cel_table, file, row.names = FALSE)
    }
  )
  
  download_bio_table <- NULL
  
  output$uniprot_biotable <- shiny::renderDataTable({
    
    GO_df <- plotUniprot()
    BiologicalDF <- Goparse(GO_df, 3)
    BiologicalDF <- na.omit(BiologicalDF)
    download_bio_table <- BiologicalDF
    BiologicalDF
  })
  
  output$download_bio_pro <- downloadHandler(
    filename = function() {
      paste("Biological-Process", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(download_bio_table, file, row.names = FALSE)
    }
  )
  
  download_mol_table <- NULL
  
  output$uniprot_molctable <- shiny::renderDataTable({
    
    GO_df <- plotUniprot()
    MolecularDF <- Goparse(GO_df, 4)
    MolecularDF <- na.omit(MolecularDF)
    download_mol_table <- MolecularDF
    MolecularDF
    
  })
  
  output$download_mole_func <- downloadHandler(
    filename = function() {
      paste("Molecular-Function", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(download_mol_table, file, row.names = FALSE)
    }
  )
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  
  
  ###################################
  ###################################
  #######  P-P Interactions   #######
  ###################################
  ###################################
  
  
  ############ Initializing Variables ###############
  
  df_interaction <- reactiveVal(0)
  df_names <- reactiveVal(0)
  
  # tbl.nodes <- data.frame(id=c("A", "B", "C"),
  #                       type=c("kinase", "TF", "glycoprotein"),
  #                       lfc=c(1, 1, 1),
  #                       count=c(0, 0, 0),
  #                       stringsAsFactors=FALSE)
  
  # tbl.edges <- data.frame(source=c("A", "B", "C"),
  #                       target=c("B", "C", "A"),
  #                       interaction=c("phosphorylates", "synthetic lethal", "unknown"),
  #                       stringsAsFactors=FALSE)
  
  ####################################################
  
  df_prot_Int <- reactive({
    print("running")
    if (is.null(input$file_prot_Int)&& is.null(input$text_prot_Int)) {
      return(NULL)
    }
    else if(!is.null(input$file_prot_Int)){
    parts <- strsplit(input$file_prot_Int$datapath, ".", fixed = TRUE)
    type <- parts[[1]][length(parts[[1]])]
    if (type != "csv") {
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return(NULL)
    }
    
    Accessions <- read.csv(input$file_prot_Int$datapath)
    Accessions <- na.omit(Accessions)
    Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
    }else{
      Acessions<-strsplit(input$text_prot_Int," ")
      Accessions <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        Accessions<-rbind(Accessions,Acessions[[1]][x])
      }
      
      print(Accessions)
      Accessions <- na.omit(Accessions)
      Accessions <- Accessions[!duplicated(Accessions[, 1]), ]
    }
    
    return(Accessions)
    
  })
  
  
  observeEvent(input$fit, ignoreInit=TRUE, {
    fit(session, 80)
  })
  
  
  observeEvent(input$loadStyleFile,  ignoreInit=TRUE, {
    if(input$loadStyleFile != ""){
      tryCatch({
        loadStyleFile(input$loadStyleFile)
      }, error=function(e) {
        msg <- sprintf("ERROR in stylesheet file '%s': %s", input$loadStyleFile, e$message)
        showNotification(msg, duration=NULL, type="error")
      })
      later(function() {updateSelectInput(session, "loadStyleFile", selected=character(0))}, 0.5)
    }
  })
  
  
  observeEvent(input$doLayout,  ignoreInit=TRUE,{
    if(input$doLayout != ""){
      strategy <- input$doLayout
      doLayout(session, strategy)
      later(function() {updateSelectInput(session, "doLayout", selected=character(0))}, 1)
    }
  })
  
  
  observeEvent(input$selectName,  ignoreInit=TRUE,{
    selectNodes(session, input$selectName)
  })
  
  
  observeEvent(input$sfn,  ignoreInit=TRUE,{
    selectFirstNeighbors(session)
  })
  
  
  observeEvent(input$fitSelected,  ignoreInit=TRUE,{
    fitSelected(session, 100)
  })
  
  
  observeEvent(input$getSelectedNodes, ignoreInit=TRUE, {
    output$selectedNodesDisplay <- renderText({" "})
    getSelectedNodes(session)
  })
  
  
  observeEvent(input$clearSelection,  ignoreInit=TRUE, {
    clearSelection(session)
  })  
  
  
  observeEvent(input$removeGraphButton, ignoreInit=TRUE, {
    removeGraph(session)
  })
  
  
  observeEvent(input$addRandomGraphFromDataFramesButton, ignoreInit=TRUE, {
    source.nodes <-  LETTERS[sample(1:5, 5)]
    target.nodes <-  LETTERS[sample(1:5, 5)]
    tbl.edges <- data.frame(source=source.nodes,
                            target=target.nodes,
                            interaction=rep("generic", length(source.nodes)),
                            stringsAsFactors=FALSE)
    all.nodes <- sort(unique(c(source.nodes, target.nodes, "orphan")))
    tbl.nodes <- data.frame(id=all.nodes,
                            type=rep("unspecified", length(all.nodes)),
                            stringsAsFactors=FALSE)
    addGraphFromDataFrame(session, tbl.edges, tbl.nodes)
  })
  
  
  observeEvent(input$selectedNodes, {
    newNodes <- input$selectedNodes;
    output$selectedNodesDisplay <- renderText({
      paste(newNodes)
    })
  })
  
  
  output$cyjShiny <- renderCyjShiny({
    print(" renderCyjShiny invoked")
    print("graph.json:")
    
    
    print("running...")
    
    
    # tryCatch({
    
    Accessions <- df_prot_int_id()
    print("Please Wait... Fetching interaction data. It may take a while")
    protein_interaction_df <- getInteraction(Accessions)
    df_interaction(protein_interaction_df)
    print("Fetched...")
    
    #migrating rowId to first colunm 
    # protein_interaction_df <- cbind(ID = rownames(protein_interaction_df),protein_interaction_df)
    # rownames(protein_interaction_df) <- 1:nrow(protein_interaction_df)
    
    #making nodes
    nodes <- as.character(protein_interaction_df[,1])
    for (i in 1:nrow(protein_interaction_df))
    {
      if(!(is.na(protein_interaction_df[i,2])))
      {
        data_df <- strsplit(as.character(protein_interaction_df[i,2]),"; ")
        for(j in data_df)
        {
          nodes <- c(nodes,j)
        }
      }
    }
    
    print(nodes)
    
    print("Please Wait... Fetching Gene Names. It may take a while")
    protein_gene_name <- getGeneNames(nodes)
    df_names(protein_gene_name)
    print("........................")
    print(as.character(protein_gene_name[,1]))
    print("Fetched...")
    edge_source <- character()
    edge_target <- character()
    
    for (i in 1:nrow(protein_interaction_df))
    {
      if(!(is.na(protein_interaction_df[i,2])))
      {
        data_df <- strsplit(as.character(protein_interaction_df[i,2]),"; ")
        for(j in data_df)
        {
          edge_source <- c(edge_source,rep(as.character(protein_gene_name[as.character(protein_interaction_df[i,1]),1]),length(j)))
          print(as.character(protein_gene_name[j,1]))
          edge_target <- c(edge_target,as.character(protein_gene_name[j,1]))
        }
      }
    }
    
    tbl.nodes <- data.frame(id=as.character(protein_gene_name[,1]),
                            type=as.character(protein_gene_name[,1]),
                            stringsAsFactors=FALSE)
    
    
    tbl.edges <- data.frame(source=edge_source,
                            target=edge_target,
                            interaction=edge_target,
                            stringsAsFactors=FALSE)
    
    # }, error = function(error_condition) {
    #   print("using defauslt value")
    # })
    
    graph.json <- dataFramesToJSON(tbl.edges, tbl.nodes)
    
    print(fromJSON(graph.json))
    cyjShiny(graph=graph.json, layoutName="cola", styleFile = "./www/style/basicStyle.js")
  })
  
  
  # observeEvent(input$submit_prot_Int, {
  
  #   print("running...")
  #   Accessions <- df_prot_Int()
  #   print("Please Wait... Fetching interaction data. It may take a while")
  #   protein_interaction_df <- getInteraction(Accessions)
  #   df_interaction(protein_interaction_df)
  #   print("Fetched...")
  
  #migrating rowId to first colunm 
  # protein_interaction_df <- cbind(ID = rownames(protein_interaction_df),protein_interaction_df)
  # rownames(protein_interaction_df) <- 1:nrow(protein_interaction_df)
  
  #making nodes
  # nodes <- as.character(protein_interaction_df[,1])
  # for (i in 1:nrow(protein_interaction_df))
  # {
  #   if(!(is.na(protein_interaction_df[i,2])))
  #   {
  #     data_df <- strsplit(protein_interaction_df[i,2],"; ")
  #     for(j in data_df)
  #     {
  #       nodes <- c(nodes,j)
  #     }
  #   }
  # }
  
  # print("Please Wait... Fetching Gene Names. It may take a while")
  # protein_gene_name <- getGeneNames(nodes)
  # df_names(protein_gene_name)
  # print("Fetched...")
  
  # # print("Rendering Visualization using Cytoscape")
  
  # # g <- graphNEL(as.character(protein_gene_name[,1]), edgemode="undirected")
  # for (i in 1:nrow(protein_interaction_df))
  # {
  #   if(!(is.na(protein_interaction_df[i,2])))
  #   {
  #     data_df <- strsplit(protein_interaction_df[i,2],"; ")
  #     for(j in data_df)
  #     {
  #       g <- graph::addEdge(as.character(protein_gene_name[as.character(protein_interaction_df[i,1]),1]), as.character(protein_gene_name[j,1]), g)
  #     }
  #   }
  # }
  
  # nodeDataDefaults(g, attr="label") <- "undefined"
  # nodeDataDefaults(g, attr="type") <- "undefined"
  # nodeDataDefaults(g, attr="flux") <- 0
  # edgeDataDefaults(g, attr="edgeType") <- "undefined"
  
  # rcy <- RCyjs(title="RCyjs vignette")
  # setGraph(rcy, g)
  # print("set")
  # strategies <- getLayoutStrategies(rcy)
  # print(strategies)
  
  # RCyjs::layout(rcy, "grid")
  # # print("lay")
  # fit(rcy, padding=200)
  # print("fit")
  # setDefaultStyle(rcy)
  
  # })
  
  df_prot_int_id <- eventReactive(input$submit_prot_Int, {
    hide("help_text_p_inte")
    Accessions <- df_prot_Int()
    return(Accessions)
  })
  
  output$help_text_p_inte <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page retrieves the Protein-Protein Interactions 
          from <a href ='http://UniProt.org'>UniProt.org</a> of a given set of UniProt accessions.
          </b>
        </p>
      </center>
    ")
  })
  
  getInteraction <- function(ProteinAccList) {
    
    if(!has_internet())
    {
      message("Please connect to the internet as the package requires internect connection.")
      return()
    }
    protein_interaction_df = data.frame()
    baseUrl <- "http://www.uniprot.org/uniprot/"
    Colnames = "interactor"
    for (ProteinAcc in ProteinAccList)
    {
      #to see if Request == 200 or not
      Request <- tryCatch(
        {
          GET(paste0(baseUrl , ProteinAcc,".xml") , timeout(60))
        },error = function(cond)
        {
          message("Internet connection problem occurs and the function will return the original error")
          message(cond)
        }
      )
      #this link return information in tab formate (format = tab)
      ProteinName_url <- paste0("?query=accession:",ProteinAcc,"&format=tab&columns=",Colnames)
      RequestUrl <- paste0(baseUrl , ProteinName_url)
      RequestUrl <- URLencode(RequestUrl)
      
      if (Request$status_code == 200){
        # parse the information in DataFrame
        ProteinDataTable <- tryCatch(read.csv(RequestUrl, header = TRUE, sep = '\t'), error=function(e) NULL)
        if (!is.null(ProteinDataTable))
        {
          ProteinDataTable <- ProteinDataTable[1,]
          ProteinInfoParsed <- as.data.frame(ProteinDataTable,row.names = ProteinAcc)
          # add Dataframes together if more than one accession
          protein_interaction_df <- rbind(protein_interaction_df, ProteinInfoParsed)
          print(paste0(ProteinAcc," interactions Fetched.."))
        }
      }else {
        HandleBadRequests(Request$status_code)
      }
    }
    
    protein_interaction_df <- cbind(ID = rownames(protein_interaction_df),protein_interaction_df)
    rownames(protein_interaction_df) <- 1:nrow(protein_interaction_df)
    
    return(protein_interaction_df)
    
  }
  
  getGeneNames <- function(ProteinAccList) {
    
    # baseUrl <- "http://www.uniprot.org/uniprot/"
    # Colnames = "genes(PREFERRED)"
    
    # protein_gene_name = data.frame()
    # for (ProteinAcc in ProteinAccList)
    # {
    #   #to see if Request == 200 or not
    #   Request <- tryCatch(
    #     {
    #       GET(paste0(baseUrl , ProteinAcc,".xml") , timeout(10))
    #     },error = function(cond)
    #     {
    #       message("Internet connection problem occurs and the function will return the original error")
    #       message(cond)
    #     }
    #   ) 
    #   #this link return information in tab formate (format = tab)
    #   ProteinName_url <- paste0("?query=accession:",ProteinAcc,"&format=tab&columns=",Colnames)
    #   RequestUrl <- paste0(baseUrl , ProteinName_url)
    #   RequestUrl <- URLencode(RequestUrl)
    #   if (Request$status_code == 200){
    #     # parse the information in DataFrame
    #     ProteinDataTable <- tryCatch(read.csv(RequestUrl, header = TRUE, sep = '\t'), error=function(e) NULL)
    #     if (!is.null(ProteinDataTable))
    #     {
    #       ProteinDataTable <- ProteinDataTable[1,]
    #       ProteinInfoParsed <- as.data.frame(ProteinDataTable,row.names = ProteinAcc)
    #       # add Dataframes together if more than one accession
    #       protein_gene_name <- rbind(protein_gene_name, ProteinInfoParsed)
    #       print(paste0(ProteinAcc," name fetched"))
    #     }  else
    #   {
    #     ProteinDataTable <- as.character(ProteinAcc)
    #     ProteinInfoParsed <- as.data.frame(ProteinDataTable,row.names = ProteinAcc)
    #     print(ProteinInfoParsed)
    #     # add Dataframes together if more than one accession
    #     protein_gene_name <- rbind(protein_gene_name, ProteinInfoParsed)
    #   }
    
    #   }else {
    #     HandleBadRequests(Request$status_code)
    
    #       ProteinDataTable <- as.character(ProteinAcc)
    #       ProteinInfoParsed <- as.data.frame(ProteinDataTable,row.names = ProteinAcc)
    #       # add Dataframes together if more than one accession
    #       protein_gene_name <- rbind(protein_gene_name, ProteinInfoParsed)
    #   }
    # }
    
    # return(protein_gene_name)
    
    protein_gene_name = data.frame()
    # print(gene_names)
    # gene_names_df <- data.frame(
    #   key = gene_names[,1],
    #   pair = gene_names[,2]
    # )
    for (ProteinAcc in ProteinAccList)
    {
      ProteinDataTable <- as.character(lookup(ProteinAcc, as.data.frame(id_to_name), missing=ProteinAcc))
      ProteinInfoParsed <- as.data.frame(ProteinDataTable,row.names = ProteinAcc)
      # add Dataframes together if more than one accession
      protein_gene_name <- rbind(protein_gene_name, ProteinInfoParsed)
    }
    
    return(protein_gene_name)
    
    
  }
  
  output$prot_int_table <- DT::renderDataTable({
    
    
    protein_interaction_df <- df_interaction()
    protein_gene_name <- df_names()
    print(protein_interaction_df)
    print("here")
    print(class(protein_interaction_df))
    if(df_names() == 0)
    {
      
      p_int_formatted <- data.frame()
      
    } else {
      
      protein_interaction_df[,1] <- as.character(protein_interaction_df[,1])
      
      p_int_formatted <- data.frame()
      count = 0
      n = 1
      for ( id in protein_interaction_df[,1])
      {
        count = count + 1
        if(!is.null(protein_interaction_df[,2]))
        {
          a = strsplit(as.character(protein_interaction_df[,2]),"; ")
          
          for(int_with in a[[count]])
          {
            p_int_row <- data.frame(id = as.character(paste0(as.character(lookup(id, as.data.frame(id_to_name), missing="Not found"))," ( ", id," )")),
                                    Interacts_With = as.character(paste0(as.character(lookup(int_with, as.data.frame(id_to_name), missing="Not found"))," ( ", int_with," )")),
                                    row.names = n)
            p_int_formatted <- rbind(p_int_formatted,p_int_row)
            n = n + 1
          }
        }
      }
      
      # for(i in 1:nrow(protein_interaction_df))
      # {
      #     protein_interaction_df[i,1] <- paste0(protein_interaction_df[i,1],
      #                               ' (',
      #                               protein_gene_name[protein_interaction_df[i,1],1],
      #                               ')')
      # }
      # print(protein_interaction_df)
      # colnames(protein_interaction_df)[2] <- "Interacts With"
      
    }
    
    p_int_formatted
    
  })
  
  output$prot_name_table <- DT::renderDataTable({
    protein_gene_name <- df_names()
    if(protein_gene_name == 0)
    {
      protein_gene_name <- data.frame()
    } else {
      
      protein_gene_name <- cbind(ID = rownames(protein_gene_name),protein_gene_name)
      rownames(protein_gene_name) <- 1:nrow(protein_gene_name)
      colnames(protein_gene_name)[2] <- "Names"
      
    } 
    protein_gene_name
    
  })
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  
  
  ###################################
  ###################################
  #########  Gemne Mania  ###########
  ###################################
  ###################################
  
  df_genemania <- reactive({
    print("running")
    if (is.null(input$file_gene)&&is.null(input$text_gene)) {
      return(NULL)
    }else if(!is.null(input$file_gene)){
    parts <- strsplit(input$file_gene$datapath, ".", fixed = TRUE)
    type <- parts[[1]][length(parts[[1]])]
    if (type != "csv") {
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return(NULL)
    }
    
    gene_names <- read.csv(input$file_gene$datapath)
    gene_names <- na.omit(gene_names)
    gene_names <- gene_names[!duplicated(gene_names[, 1]), ]
    }else{
      gene_name<-strsplit(input$text_prot_Int," ")
      gene_names <- data.frame(gene_name[[1]][1])
      for (x in 2:length(gene_name[[1]])) {
        gene_names<-rbind(gene_names,gene_name[[1]][x])
      }
      
      print(gene_names)
      gene_names <- na.omit(gene_names)
      gene_names <- gene_names[!duplicated(gene_names[, 1]), ]
    }
    
    return(gene_names)
    
  })
  
  observeEvent(input$genemania_submit, {
    
    hide("help_text_gene_mania")
    print("running...")
    organism_id <- input$organismID
    gene_names <- df_genemania()
    base_url <- "http://genemania.org/search/"
    
    url <- paste0(base_url,organism_id)
    for ( names in as.character(gene_names) )
    {
      url <- paste0(url,"/",names)
    }
    # print(gene_mania_link())
    print(url)
    gene_mania_link(url)
    shinyjs::toggle("hide_link")
    
  })
  
  output$linkCo <- renderUI({
    tags$a(href = gene_mania_link(), "here", inline =TRUE)
  })
  
  output$help_text_gene_mania <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page submits a given gene list to 
          <a href ='http://genemania.org/'>GeneMania.org</a> to retrieve the Co-expression 
          </b>
        </p>
      </center>
    ")
  })
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  
  ###################################
  ###################################
  ######  Protein Sequences  ########
  ###################################
  ###################################
  Proteins <- NULL
  df_prot_seq <- eventReactive(input$submit_prot_Seq, {
    print("running")
    if (is.null(input$file_prot_seq)&& is.null(input$text_prot_seq)) {
      return(NULL)
    }
    else if(!is.null(input$file_prot_seq)){
    parts <- strsplit(input$file_prot_seq$datapath, ".", fixed = TRUE)
    type <- parts[[1]][length(parts[[1]])]
    if (type != "csv") {
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return(NULL)
    }
    
    protein_Id <- unique(as.character(na.omit(read.csv(input$file_prot_seq$datapath)[,1])))
    Proteins <<- protein_Id
    }
    else{
      Acessions<-strsplit(input$text_prot_seq," ")
      Proteins <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        Proteins<-rbind(Proteins,Acessions[[1]][x])
      }
      
      print(Proteins)
      Proteins <- na.omit(Proteins)
      Proteins <- Proteins[!duplicated(Proteins[, 1]), ]
      
    }
    
    shinyjs::show("downloadData")
    return(Proteins)
    
  })
  
  Seqdata <- NULL
  
  output$help_text_prot_seq <- renderUI({
    HTML("<br>
    <br>
      <center>
        <p>
          <b>This page retrieves the full protein sequences from <a href ='https://www.uniprot.org/'>UniProt.org</a> of a given set of UniProt accessions, Please upload accessions to start analysis.
          </b>
        </p>
      </center>
    ")
  })
  
  output$help_text_prot_seq_evol <- renderUI({
    HTML("<br>
    <br>
      <center>
        <p>
          <b>
          This page performs Evolutionary analysis of protein sequences retrieved from <a href ='https://www.uniprot.org/'>UniProt.org</a>, Please upload accessions to start analysis.
          </b>
        </p>
      </center>
    ")
  })
  
  output$help_text_prot_seq_Patho <- renderUI({
    HTML("
    <br>
    <br>
      <center>
        <p>
          <b>
          This page retrieves protein's pathological information from <a href ='https://www.uniprot.org/'>UniProt.org</a> of a given set of UniProt accessions, Please upload accessions to start analysis.
          </b>
        </p>
      </center>
    ")
    })
  output$SequencePlot <- renderPlot(
    {
      if (!is.null(df_prot_seq()))
      {
        hide("help_text_prot_seq")
        if(is.null(Seqdata))
        {
          Proteins <- df_prot_seq()
          Seqdata <<- GetSequences(Proteins)
        }
        PlotPhysicochemical(Seqdata)
      }
    }
    
  )
  output$GravyPlot <- renderPlot(
    {
      if (!is.null(df_prot_seq()))
      {
        hide("help_text_prot_seq")
        if(is.null(Seqdata))
        {
          Proteins <- df_prot_seq()
          Seqdata <<- GetSequences(Proteins)
        }
        PlotGravy(Seqdata)
      }
      
    }
  )
  output$ChargePlot <- renderPlot(
    {
      if (!is.null(df_prot_seq()))
      {
        hide("help_text_prot_seq")
        if(is.null(Seqdata))
        {
          Proteins <- df_prot_seq()
          Seqdata <<- GetSequences(Proteins)
        }
        PlotCharge(Seqdata)
      }
      
    }
  )
  output$AcidityPlot <- renderPlot(
    {
      if (!is.null(df_prot_seq()))
      {
        hide("help_text_prot_seq")
        if(is.null(Seqdata))
        {
          Proteins <- df_prot_seq()
          Seqdata <<- GetSequences(Proteins)
        }
        PlotAcidity(Seqdata)
      }
    }
  )
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("Sequences", ".FASTA")
    },    
    content = function(file) {
      
      Accessions <- df_prot_seq()
      for (Acc in Accessions)
      {
        Request <- tryCatch(
          {
            GET(paste0("https://www.uniprot.org/uniprot/" , Acc , ".Fasta") , timeout(10))
          },error = function(cond)
          {
            message("Internet connection problem occurs and the function will return the original error")
            message(cond)
          }
        )
        if (Request$status_code == 200)
        {
          OutNumber <<- OutNumber + 1
          Fastadata <- read.csv(paste0("https://www.uniprot.org/uniprot/" , Acc , ".Fasta") , header = F , sep = "\t")
          Sequences <- paste0(as.character(unlist(Fastadata)) , collapse = "\n")
          write.table(x = Sequences , file = paste0(FileName ,".fasta") , quote = F , row.names = F , col.names = F, append = T)
        }
        
      }
    }
  )
  
  ##################################
  ######## protein revolution ######
  
  df_prot_seq_evol <- eventReactive(input$submit_prot_seq_evol,{
    print("running")
    if (is.null(input$file_prot_seq_evol)&& is.null(input$text_prot_seq_evol)) {
      return(NULL)
    }
    else if(!is.null(input$file_prot_seq_evol)){
    parts <- strsplit(input$file_prot_seq_evol$datapath, ".", fixed = TRUE)
    type <- parts[[1]][length(parts[[1]])]
    if (type != "csv") {
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return(NULL)
    }
    
    protein_Id <- unique(as.character(na.omit(read.csv(input$file_prot_seq_evol$datapath)[,1])))
    Proteins <- protein_Id
    }
    else{
      Acessions<-strsplit(input$text_prot_seq_evol," ")
      Proteins <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        Proteins<-rbind(Proteins,Acessions[[1]][x])
      }
      
      print(Proteins)
      Proteins <- na.omit(Proteins)
      Proteins <- Proteins[!duplicated(Proteins[, 1]), ]
      
    }
    return(Proteins)
    
  })
  
  GenesObj <- NULL
  
  output$GenePlot <- renderRadialNetwork(
    {
      if (!is.null(df_prot_seq_evol()))
      {
        if (is.null(GenesObj))
        {
          Proteins <- df_prot_seq_evol()
          GenesObj <- GetNamesTaxa(Proteins)
        }
        ConstructGenes(GenesObj)
      }
    }
  )
  
  output$Chromo <- renderPlot(
    if (!is.null(df_prot_seq_evol()))
    {
      if(is.null(GenesObj))
      {
        Proteins <- df_prot_seq_evol()
        GenesObj <- GetNamesTaxa(Proteins)
      }
      PlotChromosomeInfo(GenesObj)
    }
  )
  
  output$Phylogenetic <- renderPlot(
    {
      if(!is.null(df_prot_seq_evol()))
        if(is.null(Seqdata))
        {
          Proteins <- df_prot_seq_evol()
          Seqdata <<- GetSequences(Proteins)
        }
      ConstructPhylogeny(Seqdata)
    }
  )
  
  ###################################
  ###################################
  ###################################
  ###################################
  
  #Pathogens
  df_prot_seq_Patho <- eventReactive(input$submit_prot_seq_Patho,{
    print("running")
    if (is.null(input$file_prot_seq_Patho)&& is.null(input$text_prot_seq_Patho)) {
      return(NULL)
    }
    else if(!is.null(input$file_prot_seq_Patho)){
    parts <- strsplit(input$file_prot_seq_Patho$datapath, ".", fixed = TRUE)
    type <- parts[[1]][length(parts[[1]])]
    if (type != "csv") {
      showModal(modalDialog(
        title = "Error",
        "Please input a csv file!"
      ))
      return(NULL)
    }
    
    protein_Id <- unique(as.character(na.omit(read.csv(input$file_prot_seq_Patho$datapath)[,1])))
    Proteins <- protein_Id
    }
    else{
      Acessions<-strsplit(input$text_prot_seq_Patho," ")
      Proteins <- data.frame(Acessions[[1]][1])
      for (x in 2:length(Acessions[[1]])) {
        Proteins<-rbind(Proteins,Acessions[[1]][x])
      }
      
      print(Proteins)
      Proteins <- na.omit(Proteins)
      Proteins <- Proteins[!duplicated(Proteins[, 1]), ]
      
    }
    
    return(Proteins)
    
  })
  
  Pathodata <- NULL
  DiseaseTable <- NULL
  
  output$DisaeseTable <- renderDataTable({
    if(!is.null(df_prot_seq_Patho()))
    {
      Proteins <- df_prot_seq_Patho()
      Pathodata <- GetPathology_Biotech(Proteins)
      DiseaseTable <- Get.diseases(Pathodata) 
    }
  }, escape = F)
  
  output$DiseasePlot <- renderBubbles({
    if(!is.null(df_prot_seq_Patho()))
    {
      if(!is.null(DiseaseTable))
      {
        Plot.NDiseases(DiseaseTable)
      }
      else {
        Proteins <- df_prot_seq_Patho()
        Pathodata <- GetPathology_Biotech(Proteins)
        DiseaseTable <- Get.diseases(Pathodata)
        Plot.NDiseases(DiseaseTable)
      }
    }
  })
  # session$onSessionEnded(stopApp)
}

app <- shinyApp(ui = ui, server = server)
app
