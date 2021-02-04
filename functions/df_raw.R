
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
