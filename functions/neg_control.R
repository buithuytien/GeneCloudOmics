neg_control <- reactive({
  if(is.null(input$spikes1))
    return(NULL)
  spikes <- read.csv(input$spikes1$datapath,header=F)
  spikes <- as.character(spikes[,1])
  # print(spikes[1:10])
  return(spikes)
})