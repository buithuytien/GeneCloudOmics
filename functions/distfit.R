plotDist <- reactive({
  dist.start <- Sys.time()
  dis <- input$distributions
  var <- input$dist.var
  DS <- distfit_df()
  fits <- list()
  distrs <- NULL
  numcol <- c(0,0,0,0,0,0)
  dist_zoom <- input$dist_zoom
  if(dist_zoom=='slider'){
    fit_range <- input$dist_range
  }else if(dist_zoom=='text input'){
    fit_range <- c(input$dist_range_min,input$dist_range_max)
  }
  if("Log-normal" %in% dis){
    fit_ln <- fitdist(DS[,var], "lnorm")
    fits <- c(fits,list(fit_ln))
    distrs <- c(distrs,"Log-normal")
    numcol[1]=1
  }
  if("Log-logistic" %in% dis){ 
    fit_ll <- fitdist(DS[,var], "llogis", start = list(shape = 10, scale = 10),lower=c(0,0))
    fits <- c(fits,list(fit_ll))
    distrs <- c(distrs,"Log-logistic")
    numcol[2]=1
  }
  if("Pareto" %in% dis){
    fit_P <- fitdist(DS[,var], "pareto", start = list(shape = 10, scale = 10),lower=c(0,0))
    fits <- c(fits,list(fit_P))
    distrs <- c(distrs,"Pareto")
    numcol[3]=1
  }
  if("Burr" %in% dis){
    fit_B <- fitdist(DS[,var], "burr", start = list(shape1 = 0.3, shape2 = 1, rate = 1),lower=c(0,0,0))
    fits <- c(fits,list(fit_B))
    distrs <- c(distrs,"Burr")
    numcol[4]=1
  }
  if("Weibull" %in% dis){
    fit_W <- fitdist(DS[,var], "weibull",lower=c(0,0))
    fits <- c(fits,list(fit_W))
    distrs <- c(distrs,"Weibull")
    numcol[5]=1
  }
  if("Gamma" %in% dis){
    fit_G <- fitdist(DS[,var], "gamma",lower=c(0, 0),start=list(scale=1,shape=1))
    fits <- c(fits,list(fit_G))
    distrs <- c(distrs,"Gamma")
    numcol[6]=1
  }
  dist.end <- Sys.time()
  print("Distribution fitting time")
  print(dist.end - dist.start)
  return (list(fits,distrs,numcol,var,fit_range))
  
})


distaic <- reactive({
  dist.start <- Sys.time()
  DS <- distfit_df()
  AIC.df <- as.data.frame(matrix(nrow=ncol(DS),ncol=6))
  rownames(AIC.df) <- colnames(DS)
  colnames(AIC.df) <- c("Log-normal","Log-logistic","Pareto", "Burr", "Weibull", "Gamma")
  for(i in 1:nrow(AIC.df)){
    fit_ln <- fitdist(DS[,i], "lnorm")
    fit_ll <- fitdist(DS[,i], "llogis", start = list(shape = 10, scale = 10),lower=c(0,0))
    fit_P <- fitdist(DS[,i], "pareto", start = list(shape = 10, scale = 10),lower=c(0,0))
    fit_B <- fitdist(DS[,i], "burr", start = list(shape1 = 0.3, shape2 = 1, rate = 1),lower=c(0,0,0))
    fit_W <- fitdist(DS[,i], "weibull",lower=c(0,0))
    fit_G <- fitdist(DS[,i], "gamma",lower=c(0, 0),start=list(scale=1,shape=1))
    fits <- list(fit_ln,fit_ll,fit_P,fit_B,fit_W,fit_G)
    AIC.df[i,] <- gofstat(fits)$aic
  }
  for(i in 1:nrow(AIC.df)){
    AIC.df$min.AIC[i]<-colnames(AIC.df)[which.min(AIC.df[i,1:6])]
  }
  dist.end <- Sys.time()
  print('distribution fitting time')
  print(dist.end - dist.start)
  return (AIC.df)
})


distplot <- function(){
  li <- plotDist()
  fits <- li[[1]]
  distrs <- li[[2]]
  numcol <- li[[3]]
  var <- li[[4]]
  fit_range <- li[[5]]
  line_types <- c(1,2,3,4,5,6) #par lty
  if(length(fits)!=0)
    cdfcomp(fits, xlogscale = TRUE, ylogscale = TRUE,
            ylab = "CDF", xlab = "Expression levels (log)", xlim = c(fit_range[1], fit_range[2]),
            legendtext = distrs, cex = 0.5 ,lwd=2, main = var,fitcol=rainbow(6)[which(numcol==1)],fitlty = line_types[which(numcol==1)])
} 

