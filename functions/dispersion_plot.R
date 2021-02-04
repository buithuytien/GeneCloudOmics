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
