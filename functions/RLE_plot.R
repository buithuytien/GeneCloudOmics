RLE.plot <- reactive({
  type <- input$file_type
  if(type=='norm'){
    DS <- df_shiny()
  }else if(type=='raw'){
    DS <- df_raw_shiny()
  }
  set1 <- newSeqExpressionSet(as.matrix(DS))
  norm_method_name <- input$norm_method
  colors <- c('RPKM'='blue','FPKM'='darkcyan','TPM'='darkgreen',"RUV"='Brown',"Upper Quartile"='Brown')
  if(norm_method_name!="None" &input$submit_preprocessing != 0){
    spikes <- neg_control()
    if(norm_method_name == "RUV" & is.null(spikes))
      norm_method_name <- "Upper Quartile"
    plotRLE(set1, ylim=c(-1.5,1.5),outline=FALSE, col=colors[norm_method_name],
            main= paste(norm_method_name,"Normalized"))
  }
})
