go_res <- eventReactive(input$submit_go, {
  go.start <- Sys.time()
  genes <- gene_list()
  if(is.null(genes)){
    showModal(modalDialog(
      title = "Error","Please check your input DE genes!"
    ))
    return (NULL)
  }
  go_method <- input$go_method
  if(go_method %in% c("clusterProfiler","GOstats")){
    dbs_name <- input$go_species
    orgDb <- DBS[[dbs_name]]
    keyType <- input$go_geneidtype
    ont <- input$subontology
    bg <- bg_list()
    if(go_method == "clusterProfiler"){
      temp <- goPrep(fg=genes,bg=bg,keyType=keyType,orgDb=orgDb) # list(fg_mapped, bg_mapped, fg_unmapped)
      if(is.null (temp) )
        return(NULL)
      res <- enrichgoApply(gene_list=genes, keyType, orgDb,ont=ont, pvalueCutoff=1, qvalueCutoff=1)
    } else if(go_method == "GOstats"){
      temp <- goPrep(fg=genes,bg=bg,keyType=keyType,orgDb=orgDb) # list(fg_mapped, bg_mapped, fg_unmapped)
      if(is.null (temp) )
        return(NULL)
      fg_mapped <- temp[[1]] ; bg_mapped <- temp[[2]]
      if(dbs_name == "org.Sc.sgd.db"){
        primary_id <-  "ORF"
      } else if(dbs_name == "org.At.tair.db"){
        primary_id <- "TAIR"
      } else{
        primary_id <- "ENTREZID"
      }
      res <- gostatsApply(fg_mapped, bg_mapped, keyType, orgDb, ont=ont,primary_id=primary_id, pvalueCutoff=1)
    }
  } else if(go_method == "enrichR"){
    if (!havingIP()){
      showModal(modalDialog(
        title = "Error","enrichR requires internet connection. Please check your internet connection and try again!"
      ))
      return(NULL)
    }
    dbs <- input$enrichR_dbs
    res <- enrichrApply(genes,dbs)
    if (nrow(res) == 0){
      showModal(modalDialog(
        title = "Warning","No term matched! Please ensure all input DE genes are in SYMBOL format!"
      ))
    }
  }
  
  if(is.null(res)){
    showModal(modalDialog(
      title = "Warning","NULL result given! Organism or identifier might not be correct! Or try another method"
    ))
  }
  go.end <- Sys.time()
  print("Go time")
  print(go.end - go.start)
  return(res)   # no filter by min counts in each go term
})

go_res_filt <- function(){
  res <- go_res()
  if(! is.null (res)){
    if(! is.data.frame(res) ) res <- as.data.frame(res)
    go_method <- input$go_method
    min_no <- input$go_min_no
    max_p <- input$go_max_p
    res_filt <- filter(res, Count>=min_no)   # filter by gene counts in go terms
    if(go_method == "clusterProfiler"){
      res_filt <- filter(res, p.adjust <= max_p)
    } else if(go_method == "GOstats"){
      res_filt <- filter(res, Pvalue <= max_p)
    }
    # res_filt
  } else{
    res_filt = NULL
  }
  return(res_filt)
}

output$go_table <- DT::renderDataTable({
  go_res_filt()
})

go_pie_res <- eventReactive(input$submit_go, {
  genes <- gene_list()
  if(is.null(genes))
    return (NULL)
  go_method <- input$go_method
  if(go_method == "enrichR")
    return (NULL)
  dbs_name <- input$go_species
  orgDb <- DBS[[dbs_name]]
  keyType <- input$go_geneidtype
  ont <- input$subontology
  bg <- bg_list()
  
  temp <- goPrep(fg=genes,bg=bg,keyType=keyType,orgDb=orgDb) # list(fg_mapped, bg_mapped, fg_unmapped)
  if(is.null (temp) )
    return(NULL)
  GOresult <- groupGO(genes, keyType= keyType, OrgDb = orgDb, ont = ont, level = 2,
                      readable = FALSE)
  res <- dplyr::filter(GOresult@result,Count!=0)
  res$Term <- paste(res$Description," (",res$ID,")",sep="")
  res <- dplyr::arrange(res,desc(Count))
  res2 <- res[,c("Term","Count")]
  res2 <- rbind(res2,c("Unclassified",length(temp[[3]])))
  return(res2)
})

goPie <- function(){
  res <- go_pie_res()
  print("from goPie line 2543 print res")
  print(res)
  if(! is.null(res) ){
    p <- plot_ly(res, labels = ~Term, values = ~Count, type = 'pie',
                 textposition = 'inside',
                 textinfo = 'label+percent',
                 insidetextfont = list(color = '#FFFFFF'),
                 #hoverinfo = 'text',
                 #text = ~paste('$', X1960, ' billions'),
                 marker = list(line = list(color = '#FFFFFF', width = 1)),
                 #The 'pull' attribute can also be used to create space between the sectors
                 showlegend = FALSE) %>%
      layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
    # title = 'Gene ontology',
    return (p)
  }
}


goList <- eventReactive(input$submit_go_graph, {
  go_method <- input$go_method
  go_terms <- input$go_term_slect
  res <- go_res_filt()
  # print("from goList")
  # print(go_method)
  # print(head(res))
  
  if(go_method=="clusterProfiler" & !is.null(res) ){
    res_filt <- res[res$Description%in%go_terms,]
    go_list <- goToList(res_filt)
    # print(head(go_list)) ##
    return(go_list)
  } else {
    return(NULL)
  }
})

go_graph_plot <- function() {
  go_list <-goList()
  show_gene_names <- input$show_gene_names
  if(! is.null (go_list)){
    temp <- graphGene(go_list); p <- temp[[1]]; q <- temp[[2]]
    if(show_gene_names){
      p <- p + geom_node_text(aes_(label=~name), repel=TRUE)
    }
  } else {
    return(NULL)
  }
  return(p)
}
