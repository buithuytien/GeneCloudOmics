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
