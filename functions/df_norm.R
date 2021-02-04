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