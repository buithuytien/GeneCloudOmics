computeBin <- function(arr){  # Doane's rule
  n <- length(arr)
  gx <- moments::skewness(arr)
  sigmag <- sqrt(6*(n-2)/((n+1)*n+3))
  bin <- 1+log2(n)+log2(1+abs(gx)/sigmag)
  return (bin)
}

getBinCounts <- function(arr){
  vec <- entropy::discretize(arr,computeBin(arr),r=range(arr))
  return (vec)
}

entropyPlot <- reactive({
  entropy.start.time <- Sys.time()
  type <- input$file_type
  if(type=='norm'){
    DS <- df_shiny()
  }else if(type=='raw'){
    DS <- df_raw_shiny()
  }
  if(is.null(DS)==FALSE){
    tsflag <- input$tsflag
    graph <- input$entropy_graph_type
    names <- colnames(DS)
    xform <- list(categoryorder = "array",
                  categoryarray = names)
    entropy.vector <- NULL #entropy of each column
    for (i in 1:length(DS)){
      binCount <- getBinCounts(DS[,i])
      entropy <- entropy.empirical(binCount,unit="log2")
      entropy.vector <- c(entropy.vector,entropy)
    }
    if(tsflag==FALSE){
      if(graph=="Bar chart"){
        p <- plot_ly(
          x = names,
          y = entropy.vector,
          type = "bar"
        ) %>% layout(xaxis = xform)
      }else if(graph=="Line chart"){
        p <- plot_ly(
          x = names,
          y = entropy.vector,
          type="scatter",
          mode="lines+markers"
        ) %>% layout(xaxis = xform)
        #yaxis=list(range = c(0, max(ent)+0.002))
      }
    }else if(tsflag==TRUE){
      tp <- as.numeric(input$entropy_timepoints)
      numOfGeno <- ncol(DS)/tp
      names <- NULL
      for(i in 1:numOfGeno){
        id <- paste('entropytype',i,sep="")
        names <- c(names,input[[id]])
      }
      time_index <- c(1:tp)
      ent <- data.frame(time_index)
      for(j in 1:numOfGeno){
        part_ent <- entropy.vector[(tp*j-(tp-1)):(tp*j)]
        ent <- cbind(ent,part_ent)
      }
      if(graph=="Bar chart"){
        p <- plot_ly(x=ent[,1],y=ent[,2],name=names[1],type='bar')
        for(i in 1:(numOfGeno-1)){
          p <- add_trace(p,y=ent[,i+2],name=names[i+1],type='bar')
        }
        p <- layout(p,xaxis = list(title = 'Time'),yaxis=list(title='Entropy'))
      }else if(graph=="Line chart"){
        p <- plot_ly(x=ent[,1],y=ent[,2],name=names[1],type='scatter',mode='lines+markers')
        for(i in 1:(numOfGeno-1)){
          p <- add_trace(p,y=ent[,i+2],name=names[i+1],type='scatter',mode='lines+markers')
        }
        p <- layout(p,xaxis = list(title = 'Time'),yaxis=list(title='Entropy'))
      }
    }
    entropy.end.time <- Sys.time()
    print("entropy time")
    print(entropy.end.time - entropy.start.time)
    return (p)
  }
})
