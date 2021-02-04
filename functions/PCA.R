refreshDS1 <- eventReactive(input$pca_refresh,{
  type <- input$file_type
  if(type=='norm'){
    DS <- df_shiny()
  }else if(type=='raw'){
    DS <- df_raw_shiny()
  }
  DS1 <- DS[sample(nrow(DS), nrow(DS), replace = FALSE), ]
  return (DS1)
})

plotPCA <- reactive({ #process and return data
  pca.start <- Sys.time()
  type <- input$file_type
  if(type=='norm'){
    DS <- df_shiny()
  }else if(type=='raw'){
    DS <- df_raw_shiny()
  }
  order <- input$gene_order
  size <- input$gene_size
  x <- input$pca.x
  y <- input$pca.y
  cluster_flag <- input$pca_cluster
  rindex <- as.numeric(substring(x,3))
  cindex <- as.numeric(substring(y,3))
  if(order=='Ascending'){
    DS1 <- DS[order(DS[,1]),]
  }else if(order=='Descending'){
    DS1 <- DS[rev(order(DS[,1])),]
  }else if(order=='Random'){
    DS1 <- refreshDS1()
  }
  
  DSample <- head(DS1, n = size)
  PR <- prcomp(t(DSample),center=TRUE)
  PCA.var <- PR$sdev^2
  PCA.var.per <- round(PCA.var/sum(PCA.var)*100,1)
  xlabel <- paste(colnames(PR$x)[rindex]," - ", PCA.var.per[rindex], "%", sep="")
  ylabel <- paste(colnames(PR$x)[cindex]," - ", PCA.var.per[cindex], "%", sep="")
  if(cluster_flag==TRUE){
    num <- as.numeric(input$pca_cluster_num)
    kmeans.data <- data.frame(x=PR$x[,x],y=PR$x[,y])
    kmeans.result <- kmeans(kmeans.data,num)
    return (list(PR,PCA.var,PCA.var.per,rindex,cindex,xlabel,ylabel,cluster_flag,kmeans.result))
  }
  pca.end <- Sys.time()
  print("pca time")
  print(pca.end - pca.start)
  return (list(PR,PCA.var,PCA.var.per,rindex,cindex,xlabel,ylabel,cluster_flag))
})

pcavarplot <- function(){
  li <- plotPCA()
  PCA.var.per <- li[[3]]/100
  type <- input$file_type
  if(type=='norm'){
    DS <- df_shiny()
  }else if(type=='raw'){
    DS <- df_raw_shiny()
  }
  pcchoices <- NULL
  for (i in 1:length(PCA.var.per)){
    pcchoices <- c(pcchoices,paste("PC",i,sep=""))
  }
  xform <- list(categoryorder = "array",
                categoryarray = pcchoices)
  p <- plot_ly(
    x = pcchoices,
    y = PCA.var.per,
    name = "PCA variance",
    type = "bar"
  ) %>% layout(xaxis = xform)
  
  return (p)
}

pca2dplot <- function(){
  li <- plotPCA()
  PR <- li[[1]]
  rindex <- li[[4]]
  cindex <- li[[5]]
  xlabel <- li[[6]]
  ylabel <- li[[7]]
  cluster_flag <- li[[8]]
  if(cluster_flag==FALSE){
    p <- plot_ly(
      x=PR$x[,rindex],
      y=PR$x[,cindex],
      type = "scatter",
      mode="markers"
    ) %>% layout(xaxis = list(title = xlabel), yaxis = list(title = ylabel))
  }else if(cluster_flag==TRUE){
    kmeans.result <- li[[9]]
    text_flag <- input$pca_text
    if(text_flag==TRUE){
      p <- plot_ly(
        x=PR$x[,rindex],
        y=PR$x[,cindex],
        type = "scatter",
        color=as.character(kmeans.result$cluster),
        mode="markers",
        colors = "Set1"
      ) %>% hide_colorbar() %>% 
        add_trace(
          x=PR$x[,rindex],
          y=PR$x[,cindex],
          type = 'scatter',
          mode = 'text', 
          text = names(kmeans.result$cluster), 
          textposition = 'top right'
        ) %>% layout(xaxis = list(title = xlabel), yaxis = list(title = ylabel),showlegend=FALSE)
    }else if(text_flag==FALSE){
      p <- plot_ly(
        x=PR$x[,rindex],
        y=PR$x[,cindex],
        type = "scatter",
        color=as.character(kmeans.result$cluster),
        mode="markers",
        text = names(kmeans.result$cluster),
        colors = "Set1"
      ) %>% hide_colorbar() %>% layout(xaxis = list(title = xlabel), yaxis = list(title = ylabel),showlegend=FALSE)
    }
  }
  
  
}

pca3dplot <- function(){
  li <- plotPCA()
  PR <- li[[1]]
  xlabel <- "PC1"
  ylabel <- "PC2"
  zlabel <- "PC3"
  cluster_flag <- li[[8]]
  if(cluster_flag==FALSE){
    p <- plot_ly(
      x=PR$x[,1],
      y=PR$x[,2],
      z=PR$x[,3],
      type="scatter3d",
      mode="markers",
      marker=list(size=5)
    ) %>% layout(scene=list(xaxis = list(title = xlabel), yaxis = list(title = ylabel),zaxis=list(title=zlabel)))
  }else if(cluster_flag==TRUE){
    kmeans.result <- li[[9]]
    text_flag <- input$pca_text
    if(text_flag==TRUE){
      p <- plot_ly(
        x=PR$x[,1],
        y=PR$x[,2],
        z=PR$x[,3],
        type = 'scatter3d',
        mode = 'text', 
        text = names(kmeans.result$cluster), 
        color = as.character(kmeans.result$cluster),
        textfont = list(size=10),
        textposition = 'top right'
      ) %>% layout(scene=list(xaxis = list(title = xlabel), yaxis = list(title = ylabel),zaxis=list(title=zlabel)),showlegend=FALSE)
    }else if(text_flag==FALSE){
      p <- plot_ly(
        x=PR$x[,1],
        y=PR$x[,2],
        z=PR$x[,3],
        type = "scatter3d",
        color=as.character(kmeans.result$cluster),
        mode="markers",
        marker=list(size=5),
        text = names(kmeans.result$cluster),
        colors = "Set1"
      ) %>% hide_colorbar() %>% layout(scene=list(xaxis = list(title = xlabel), yaxis = list(title = ylabel),zaxis=list(title=zlabel)),showlegend=FALSE)
    }
  }
}
