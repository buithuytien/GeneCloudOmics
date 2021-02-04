volcano_plot <- eventReactive(input$submit_DE, {
  volcano.start.time <- Sys.time()
  rep_number <- input$n_rep
  if( rep_number == 0)
    return (NULL)
  de_type <- input$de_method1
  if(de_type == "NOISeq")
    return (NULL)
  
  res <- de_no_filt()   # de result, no filter
  if(is.null (res))
    return (NULL)
  
  p_val <- input$p_val
  fc <- input$fc
  res$Gene <- rownames(res)
  res <- na.omit(res)
  # plot
  ymax <- max(-log10(res$PValue)); if(ymax > 5) ymax = 5
  # print("from volcano plot - range(res$PValue)")
  # print(range(res$PValue))
  with(res, plot(log2FC, -log10(PValue), pch=20, main="Volcano plot",xlim=c(-2.5,2.5),ylim=c(0,ymax))) #xlim=c(-5,5),ylim=c(0,ymax)
  # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
  with(subset(res, FDR< p_val), points(log2FC, -log10(PValue), pch=20, col="red"))
  with(subset(res, abs(log2FC)>log2(fc)), points(log2FC, -log10(PValue), pch=20, col="orange"))
  with(subset(res, FDR< p_val & abs(log2FC)>log2(fc)), points(log2FC, -log10(PValue), pch=20, col="green"))
  legend("topleft",bty="n",col=c("red","orange","green","black"),pch=19,
         legend=c("FDR < FDR limit","FC > FC limit","Both","Other"))
  # library(calibrate)
  # with(subset(res, FDR<.05 & abs(log2FC)>1), textxy(log2FC, -log10(PValue), labs=Gene, cex=.8))
  volcano.end.time <- Sys.time()
  print("volcano time")
  print(volcano.end.time - volcano.start.time)
})
