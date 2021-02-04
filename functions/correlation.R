COR <- function(d, i,myMethod){
  Result2 <- cor(x = d[,i], y = d[,i], method = myMethod)
  return(format(round(Result2, 5), nsmall = 5))
}

cor_df <- reactive({
  cor.start <- Sys.time()
  type <- input$file_type
  if(type=='norm'){
    DS <- df_shiny()
  }else if(type=='raw'){
    DS <- df_raw_shiny()
  }
  method <- input$cor_method
  if(method=="Pearson correlation"){
    Cor2 <- data.frame(COR((DS),1:length(DS),"pearson"))
  }else if(method=="Spearman correlation"){
    Cor2 <- data.frame(COR((DS),1:length(DS),"spearman"))
  }
  Cor2 <- na.omit(Cor2)
  cor.end <- Sys.time()
  print("correlation time")
  print(cor.end - cor.start)
  return (Cor2)
})

corrplot1 <- function(){
  corr <- as.matrix(cor_df())
  corr <- apply(corr,2,as.numeric)
  rownames(corr) <- rownames(cor_df())
  if(ncol(corr)<=20){
    fontsize <- 1
  }else{
    fontsize <- 20/ncol(corr)
  }
  corrplot(corr,method="shade",shade.col=NA,tl.col="black",cl.lim=c(min(corr),1),is.corr = FALSE,tl.cex = fontsize)
}

corrplot2 <- function(){
  corr <- as.matrix(cor_df())
  corr <- apply(corr,2,as.numeric)
  rownames(corr) <- rownames(cor_df())
  if(ncol(corr)<=20){
    fontsize <- 1
  }else{
    fontsize <- 20/ncol(corr)
  }
  corrplot(corr,type="upper",tl.col="black",cl.lim=c(min(corr),1),is.corr = FALSE,tl.cex = fontsize)
}
