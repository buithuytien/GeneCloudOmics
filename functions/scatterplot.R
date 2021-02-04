

plotScatter <- reactive({
  scatter.start <- Sys.time()
  trans <- input$trans
  x <- input$scatter.x
  y <- input$scatter.y
  type <- input$file_type
  if(type=='norm'){
    DS <- df_shiny()
  }else if(type=='raw'){
    DS <- df_raw_shiny()
  }
  if(trans=='None'){
    scatter.data <- DS
  }else if(trans=='Natural log'){
    scatter.data <- log1p(DS)
  }else if(trans=='log2'){
    scatter.data <- log2(DS+1)
  }else if(trans=='log10'){
    scatter.data <- log10(DS+1)
  }
  scatter.end <- Sys.time()
  print("Scatter plot time")
  print(scatter.end - scatter.start)
  return (list(x,y,scatter.data))
})

scatterplot <- function(){
  li <- plotScatter()
  x <- li[[1]]
  y <- li[[2]]
  scatter.data <- li[[3]]
  d <- kde2d(scatter.data[,x],scatter.data[,y])
  ColorLevels <- round(seq(min(d$z), max(d$z), length=5),4)
  heatscatter(x=scatter.data[,x],y=scatter.data[,y],xlab = x, ylab=y, main="")
  legend("topleft", paste("R=",round(cor(scatter.data[,x],scatter.data[,y]),3)), bty="n")
  legend("bottomright",title="KDE",legend=ColorLevels, pch=19,col=LSD::colorpalette("heat"))
  if(x!=y){
    lines(lowess(scatter.data[,x],scatter.data[,y]),col="black")
  }
}

scatterplot_collage <- function(){
  li <- plotScatter()
  scatter.data <- li[[3]]
  par(mfrow=c(3,3))
  for(i in 1:ncol(scatter.data)){
    for(j in i:ncol(scatter.data)){
      d <- kde2d(scatter.data[,i],scatter.data[,j])
      ColorLevels <- round(seq(min(d$z), max(d$z), length=5),4)
      heatscatter(x=scatter.data[,i],y=scatter.data[,j],xlab = colnames(scatter.data)[i], ylab=colnames(scatter.data)[j], main="")
      legend("topleft", paste("R=",round(cor(scatter.data[,i],scatter.data[,j]),3)), bty="n")
      legend("bottomright",title="KDE",legend=ColorLevels, pch=19,col=LSD::colorpalette("heat"))
      if(i!=j){
        lines(lowess(scatter.data[,i],scatter.data[,j]),col="black")
      }
    }
  }
}

