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