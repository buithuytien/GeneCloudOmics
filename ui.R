source(file.path("functions","home_tab.R"))
source(file.path("functions","preprocessing_tab.R"))
source(file.path("functions","scatter_tab.R"))
source(file.path("functions","distribution_tab.R"))
source(file.path("functions","correlation_tab.R"))
source(file.path("functions","PCA_tab.R"))
source(file.path("functions","DE_analysis_tab.R"))
source(file.path("functions","heatmap_tab.R"))
source(file.path("functions","noise_tab.R"))
source(file.path("functions","entropy_tab.R"))
source(file.path("functions","GO_analysis_tab.R"))

wd <- dirname(rstudioapi::getActiveDocumentContext()$path)  #set wd as the current folder
if(! wd == getwd()){
  setwd(wd)
}

# 
# ## sourcing util files
source(paste0("./www/utils.R"))

# 
loadPkg()

species.choices <<- c("Homo sapiens"='org.Hs.eg.db',"Mus musculus"='org.Mm.eg.db',"Rattus norvegicus"='org.Rn.eg.db',"Gallus gallus"='org.Gg.eg.db',"Danio rerio"='org.Dr.eg.db',"Drosophila melanogaster"='org.Dm.eg.db',"Caenorhabditis elegans"='org.Ce.eg.db',"Saccharomyces cereviasiae"='org.Sc.sgd.db',"Arabidopsis thaliana"='org.At.tair.db',"Escherichia coli (strain K12)"='org.EcK12.eg.db',"Escherichia coli (strain Sakai)"='org.EcSakai.eg.db',"Anopheles gambiae"='org.Ag.eg.db',"Bos taurus"='org.Bt.eg.db',"Canis familiaris"='org.Cf.eg.db',"Macaca mulatta"='org.Mmu.eg.db',"Plasmodium falciparum"='org.Pf.plasmo.db',"Pan troglodytes"='org.Pt.eg.db',"Sus scrofa"='org.Ss.eg.db',"Xenopus tropicalis"='org.Xl.eg.db')
DBS <<- list('org.Hs.eg.db'=org.Hs.eg.db,'org.Mm.eg.db'=org.Mm.eg.db,'org.Rn.eg.db'=org.Rn.eg.db,"org.Gg.eg.db"=org.Gg.eg.db,"org.Dr.eg.db"=org.Dr.eg.db,"org.Dm.eg.db"=org.Dm.eg.db,"org.Ce.eg.db"=org.Ce.eg.db,"org.Sc.sgd.db"=org.Sc.sgd.db,"org.At.tair.db"=org.At.tair.db,"org.EcK12.eg.db"=org.EcK12.eg.db,"org.EcSakai.eg.db"=org.EcSakai.eg.db,"org.Ag.eg.db"=org.Ag.eg.db,"org.Bt.eg.db"=org.Bt.eg.db,"org.Cf.eg.db"=org.Cf.eg.db,"org.Mmu.eg.db"=org.Mmu.eg.db,"org.Pf.plasmo.db"=org.Pf.plasmo.db,"org.Pt.eg.db"=org.Pt.eg.db,"org.Ss.eg.db"=org.Ss.eg.db,"org.Xl.eg.db"=org.Xl.eg.db)
enrichRdbs <- as.character(read.csv(paste0(wd,"/www/enrichRdbs.csv"))[,1])



ui <- navbarPage(id = "navbar",
                 theme = shinytheme("flatly"),
                 title = 'ABioTrans',
                 home_tab,
                 preprocessing_tab,
                 scatter_tab,
                 distribution_tab,
                 correlation_tab,
                 PCA_tab,
                 DE_analysis_tab,
                 heat_map_tab,
                 noise_tab, 
                 entropy_tab,
                 GO_analysis_tab
)
