plotHeatmap <- eventReactive(input$heatmap_plot, {#process and return data 
  heatmap.start.time <- Sys.time()
  
  type <- input$file_type
  value <- input$heatmap_value
  de_type <- input$heatmap_de_ind
  
  if(type=='norm'){
    DS <- df_shiny()
  }else if(type=='raw'){
    DS <- df_raw_shiny()
  }
  names <- NULL
  # for(i in 1:input$numOfGeno){
  #   id <- paste('type',i,sep="")
  #   names <- c(names,input[[id]])
  # }
  
  # numOfGeno <- input$numOfGeno
  # ref <- as.numeric(input$heatmap_anchor)
  clusterNum <- input$numOfCluster
  
  if(de_type == "ind"){
    fold <- as.numeric(input$fold)
    fold_ncol <- input$fold_ncol
    DS2 <- deWithoutStats(DS, FC=fold, n_col=fold_ncol)
    de_genes <- rownames(DS2)
    # print("from heatmap YT version")
    # print("de genes")
    # print(head(de_genes))
    # print(paste(length(de_genes),"genes"))
    
  } else if(de_type == "de"){
    res.df <- de_no_filt()
    if(is.null (res.df) ) return(NULL)
    p_val <- input$p_val
    fc <- input$fc
    rep_number <- input$n_rep
    res.df.filt <- de_filt(res.df, p_val,fc,rep_number) 
    de_genes <- res.df.filt$Gene
    # print("from line heatmap de result from DE analysis")
    # print("res.df.filt")
    # print(head(res.df.filt))
  }
  
  de_genes_exp <- DS[rownames(DS)%in%de_genes,]
  DS3 <- t(scale(t(de_genes_exp)))
  DS3 <- na.omit(DS3)
  # print("from line 1894 - heatmap de type")
  # print("DS3")
  # print(head(DS3))
  
  set.seed(110)
  a <- ComplexHeatmap::Heatmap(DS3, name="Normalized expression",
                               col = colorRamp2(c(min(DS3),0,max(DS3)), c("red","black", "green")),
                               row_names_gp = gpar(fontsize = 1), 
                               row_dend_gp = gpar(fontsize = 1),
                               row_title_gp = gpar(fontsize = 10),
                               cluster_columns = FALSE,
                               row_dend_width = unit(3, "cm"),
                               split = clusterNum, clustering_distance_rows = "pearson",                       
                               show_heatmap_legend = TRUE,
                               show_row_names = FALSE, show_column_names = T,
                               heatmap_legend_param = list(title = "Normalized expression") )
  set.seed(110)
  rcl.list <- row_order(a); 
  DS3.1 <- as.matrix(rownames(DS3))
  
  # Cluster <- NULL
  # for(i in 1:length(rcl.list)){
  #   for(j in 1:length(rcl.list[[i]])){
  #     pair <- c(i,DS3.1[rcl.list[[i]][j]])
  #     Cluster <- rbind(Cluster,pair)
  #   }
  # }
  # Cluster <- data.frame(Cluster,row.names = NULL)
  # colnames(Cluster) <- c("cluster","GeneID")
  
  rcl.list2 <- rcl.list
  for(i in 1:length(rcl.list)){
    rcl.list2[[i]] <- rownames(DS3)[rcl.list[[i]] ]
  }
  
  for(i in 1:length(rcl.list2)){
    genes <- rcl.list2[[i]]
    group_name <- rep(i, length(genes) )
    Cluster_i <- data.frame("GeneID"=genes, "Cluster"=i)
    if( i ==1 ){
      Cluster <- Cluster_i
    } else{
      Cluster <- rbind(Cluster, Cluster_i)
    }
  }
  
  print("line 2089, Cluster")
  print(head(Cluster))
  print(paste0("de_type = ",de_type))
  print("length of input DS")
  print(dim(DS3))
  
  # end heat map analysis
  heatmap.end.time <- Sys.time()
  print("heat map time")
  print(heatmap.end.time - heatmap.start.time)
  return (list(a,DS3,Cluster))
})

getCluster <- eventReactive(input$heatmap_plot, {
  set.seed(110) 
  ll <- plotHeatmap(); a <- ll[[1]]; DS3 <- ll[[2]]
  rcl.list <- row_order(a)
  DS3.1 <- as.matrix(rownames(DS3))
  
  Cluster <- NULL
  for(i in 1:length(rcl.list)){
    for(j in 1:length(rcl.list[[i]])){
      pair <- c(i,DS3.1[rcl.list[[i]][j]])
      Cluster <- rbind(Cluster,pair)
    }
  }
  Cluster <- data.frame(Cluster,row.names = NULL)
  colnames(Cluster) <- c("cluster","gene.id")
  
  return (Cluster[,c("gene.id", "cluster")])
})

mapPlot <- function(){
  myHeatmap <- plotHeatmap()[[1]]
  myHeatmap <- draw(myHeatmap)
}
