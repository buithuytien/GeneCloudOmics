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