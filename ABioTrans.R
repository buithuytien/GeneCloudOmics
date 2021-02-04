#Sys.setenv("plotly_username"=" your_plotly_username")
#Sys.setenv("plotly_api_key"="your_api_key")
## test repo



print("start loading")
start.load <- Sys.time()   ### time

if(length(find.package(package = 'shiny',quiet = T))>0){
  library(shiny)
}else{
  print("Package shiny not installed")
  install.packages("shiny")
  print("Package shiny installed")
  library(shiny)
}


if(length(find.package(package = 'shinythemes',quiet = T))>0){
  library(shinythemes)
}else{
  print("Package shinythemes not installed")
  install.packages("shinythemes")
  print("Package shinythemes installed")
  library(shinythemes)
}

if(length(find.package(package = 'rstudioapi',quiet = T))>0){
  library(rstudioapi)
}else{
  install.packages("rstudioapi")
  library(rstudioapi)
}


wd <- dirname(rstudioapi::getActiveDocumentContext()$path)  #set wd as the current folder
print(wd == getwd())
print(wd)
print(getwd())
if(! wd == getwd()){
  setwd(wd)
}

source(paste0("./www/utils.R"))
loadPkg()
end.load <- Sys.time()
print("loading time")
print(end.load-start.load)

##### UI from here ###########
source("ui.R")

####################################################
source("server.R")

shinyApp(ui,server)
