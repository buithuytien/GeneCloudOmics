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
