
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
  source(file.path("functions", "df_norm.R"), local = TRUE)$value
  df_norm
  # get raw counts
  source(file.path("functions","df_raw.R"), local = TRUE)$value
  df_raw
  
  # get gene length
  source(file.path("functions","gene_length.R"), local = TRUE)$value
  gene_length
  
  # get spikes / negative control genes
  source(file.path("functions","neg_control.R"), local = TRUE)$value
  neg_control
  
  # get meta data table
  source(file.path("functions","group_names.R"), local = TRUE)$value
  group_names
  
  ### Gene ontology
  source(file.path("functions","gene_ontology_functions.R"), local = TRUE)$value
  gene_list
  bg_list
  
  ####################################
  ########## PREPROCESSING ###########
  ####################################
  
  source(file.path("functions","preprocessing.R"), local = TRUE)$value
  df_shiny
  df_raw_filt
  df_raw_shiny
  distfit_df
  
  ######### ANALYSIS FROM HERE ############
  ######## RLEplot and Preprocessing ###########
  #############################################
  source(file.path("functions","RLE_plot.R"), local = TRUE)$value
  
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
  source(file.path("functions","scatterplot.R"), local = TRUE)$value
  plotScatter
  scatterplot
  scatterplot_collage
  
  
  ############################
  ######## distfit ###########
  ############################
  source(file.path("functions","distfit.R"), local = TRUE)$value
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
  plotDist
  
  output$dist.plot <- renderPlot({
    distplot()
  })
  
  distaic
  
  output$dist.aic <- renderTable({
    distaic()
  },rownames=TRUE)
  
  distplot
  
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
  source(file.path("functions","correlation.R"), local = TRUE)$value
  COR
  cor_df
  
  output$corr.plot <- renderPlot({
    corrplot1()
  })
  
  output$corr.plot2 <- renderPlot({
    corrplot2()
  })
  
  output$corr.matrix <- renderTable({
    cor_df()
  },rownames=TRUE)
  
  corrplot1
  corrplot2
  
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
  
  source(file.path("functions","PCA.R"), local = TRUE)$value
  refreshDS1
  plotPCA
  pcavarplot
  pca2dplot
  pca3dplot
  
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
  source(file.path("functions","DE analysis.R"), local = TRUE)$value
  group_names_de
  df_raw_de
  de_no_filt
  de_filt
  
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
  source(file.path("functions","volcano_plot.R"), local = TRUE)$value
  volcano_plot
  
  output$volcano_plot <- renderPlot({
    volcano_plot()
  })
  ##### dispersion plot ######
  source(file.path("functions","dispersion_plot.R"), local = TRUE)$value
  dispersion_plot
  
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
  ## Heat map
  source(file.path("functions","heatmap.R"), local = TRUE)$value
  plotHeatmap
  
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
  source(file.path("functions","noise.R"), local = TRUE)$value
  
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
  
  source(file.path("functions","entropy.R"), local = TRUE)$value
  
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
  
  source(file.path("functions","gene_ontology.R"), local = TRUE)$value
  output$go_pie <- renderPlotly({
    goPie()
  })
  goList
  go_graph_plot
  
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