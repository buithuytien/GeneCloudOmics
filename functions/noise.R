SQCO <- function(MT){
  if(ncol(MT)==1){
    res <- matrix(0)
    return(res)
  }
  temp <- NULL
  for(i in 1:nrow(MT)){
    m <- sum(MT[i,])/length(MT[i,])
    v <- stats::var(MT[i,])
    if(m!=0)
      temp <- c(temp,v/(m*m))
  }
  res <- matrix(sum(temp)/length(temp))
  return (res)
}

noisePlot <- eventReactive(input$noise_plot, {
  noise.start.time <- Sys.time()
  type <- input$file_type
  if(type=='norm'){
    DS <- df_shiny()
  }else if(type=='raw'){
    DS <- df_raw_shiny()
  }
  numOfRep <- as.numeric(input$noise_numOfRep)
  numOfGeno <- ncol(DS)/numOfRep
  graph <- input$noise_graph_type
  names <- NULL
  for(i in 1:numOfGeno){
    id <- paste('noisetype',i,sep="")
    names <- c(names,input[[id]])
  }
  
  situation <- input$noise_situation
  if(situation=='a'){
    DS1 <- list()
    for(j in 1:numOfGeno){
      DS1[[j]] <- as.matrix(DS[,((j-1)*numOfRep+1):(j*numOfRep)]) 
    } 
    Noise <- NULL
    for(y in 1:numOfGeno){
      Noise <- c(Noise,SQCO(DS1[[y]]))
    }
    xform <- list(categoryorder = "array",
                  categoryarray = names)
    if(graph=="Bar chart"){
      p <- plot_ly(
        x = names,
        y = Noise,
        type = "bar"
      ) %>% layout(xaxis = xform)
    }else if(graph=="Line chart"){
      p <- plot_ly(
        x = names,
        y = Noise,
        type="scatter",
        mode="lines+markers"
      ) %>% layout(xaxis = xform,yaxis=list(range = c(0, max(Noise)+0.001)))
    }
  }else if(situation=='b'){
    DS_ave <- NULL
    for(j in 1:numOfGeno){
      part_DS <- as.matrix(DS[,((j-1)*numOfRep+1):(j*numOfRep)]) 
      DS_ave <- cbind(DS,data.frame(matrixStats::rowMeans2(part_DS)))
    }
    anchor <- input$noise_anchor_b
    names <- NULL
    for(i in 1:numOfGeno){
      id <- paste('noisetype',i,sep="")
      names <- c(names,input[[id]])
    }
    anchor_index <- match(anchor,names)
    Noise <- NULL
    for(i in 1:numOfGeno){
      if(i!=anchor_index){
        Noise <- c(Noise,SQCO(cbind(DS_ave[,anchor_index],DS_ave[,i])))
      }
    }
    names_wo_anchor <- NULL
    for(i in 1:numOfGeno){
      if(i!=anchor_index){
        id <- paste('noisetype',i,sep="")
        names_wo_anchor <- c(names_wo_anchor,input[[id]])
      }
    }
    xform <- list(categoryorder = "array",
                  categoryarray = names_wo_anchor)
    if(graph=="Bar chart"){
      p <- plot_ly(
        x = names_wo_anchor,
        y = Noise,
        type = "bar"
      ) %>% layout(xaxis = xform)
    }else if(graph=="Line chart"){
      p <- plot_ly(
        x = names_wo_anchor,
        y = Noise,
        type="scatter",
        mode="lines+markers"
      ) %>% layout(xaxis = xform,yaxis=list(range = c(0, max(Noise)+0.001)))
    }
  }else if(situation=='c'){
    anchor <- input$noise_anchor_c
    names <- colnames(DS)
    anchor_index <- match(anchor,names)
    Noise <- NULL
    for(i in 1:ncol(DS)){
      if(i!=anchor_index){
        Noise <- c(Noise,SQCO(cbind(DS[,anchor_index],DS[,i])))
      }
    }
    names_wo_anchor <- names[-anchor_index]
    xform <- list(categoryorder = "array",
                  categoryarray = names_wo_anchor)
    if(graph=="Bar chart"){
      p <- plot_ly(
        x = names_wo_anchor,
        y = Noise,
        type = "bar"
      ) %>% layout(xaxis = xform)
    }else if(graph=="Line chart"){
      p <- plot_ly(
        x = names_wo_anchor,
        y = Noise,
        type="scatter",
        mode="lines+markers"
      ) %>% layout(xaxis = xform,yaxis=list(range = c(0, max(Noise)+0.001)))
    }
  }
  noise.end.time <- Sys.time()
  print("noise time")
  print(noise.end.time - noise.start.time)
  return (p)
})
