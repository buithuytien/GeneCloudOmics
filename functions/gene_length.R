gene_length <- reactive({
  if (is.null(input$length1))
    return (NULL)
  lengths_df <- read.csv(input$length1$datapath)
  lengths_df2 <- data.frame("len" = lengths_df[,2]); 
  rownames(lengths_df2) <- as.character(lengths_df[,1])
  return(lengths_df2)
})