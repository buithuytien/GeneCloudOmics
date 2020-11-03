GetSequences <- function (ProteinAccList, directorypath = NULL) 
{
  if (!has_internet()) {
    message("Please connect to the internet as the package requires internect connection.")
    return()
  }
  columns <- c("fragment,encodedon,comment(ALTERNATIVE PRODUCTS),comment(ERRONEOUS GENE MODEL PREDICTION),comment(ERRONEOUS INITIATION),comment(ERRONEOUS TERMINATION),comment(ERRONEOUS TRANSLATION),comment(FRAMESHIFT),comment(MASS SPECTROMETRY),comment(POLYMORPHISM),comment(RNA EDITING),comment(SEQUENCE CAUTION),length,mass,sequence,feature(ALTERNATIVE SEQUENCE),feature(NATURAL VARIANT),feature(NON ADJACENT RESIDUES),feature(NON STANDARD RESIDUE),feature(NON TERMINAL RESIDUE),feature(SEQUENCE CONFLICT),feature(SEQUENCE UNCERTAINTY),version(sequence)")
  baseUrl <- "http://www.uniprot.org/uniprot/"
  ProteinInfoParsed_total = data.frame()
  for (ProteinAcc in ProteinAccList) {
    Request <- tryCatch({
      GET(paste0(baseUrl, ProteinAcc, ".xml"))
    }, error = function(cond) {
      message("Internet connection problem occurs and the function will return the original error")
      message(cond)
    })
    ProteinName_url <- paste0("?query=accession:", ProteinAcc, 
                              "&format=tab&columns=", columns)
    RequestUrl <- paste0(baseUrl, ProteinName_url)
    RequestUrl <- URLencode(RequestUrl)
    if (length(Request) == 0) {
      message("Internet connection problem occurs")
      return()
    }
    if (Request$status_code == 200) {
      ProteinDataTable <- tryCatch(read.csv(RequestUrl, 
                                            header = TRUE, sep = "\t"), error = function(e) NULL)
      if (!is.null(ProteinDataTable)) {
        ProteinDataTable <- ProteinDataTable[1, ]
        ProteinInfoParsed <- as.data.frame(ProteinDataTable, 
                                           row.names = ProteinAcc)
        ProteinInfoParsed_total <- rbind(ProteinInfoParsed_total, 
                                         ProteinInfoParsed)
      }
    }
    else {
      HandleBadRequests(Request$status_code)
    }
  }
  if (!is.null(directorypath)) {
    write.csv(ProteinInfoParsed_total, paste0(directorypath, 
                                              "/", "Sequences Information.csv"))
  }
  return(ProteinInfoParsed_total)
}



#Get Sequence information from Uniprot
GETSeqFastaUniprot <- function(Accessions, FileName)
{
  OutNumber <- 0
  for (Acc in Accessions)
  {
    Request <- tryCatch(
      {
        GET(paste0("https://www.uniprot.org/uniprot/" , Acc , ".Fasta") , timeout(10))
      },error = function(cond)
      {
        message("Internet connection problem occurs and the function will return the original error")
        message(cond)
      }
    )
    if (Request$status_code == 200)
    {
      OutNumber <<- OutNumber + 1
      Fastadata <- read.csv(paste0("https://www.uniprot.org/uniprot/" , Acc , ".Fasta") , header = F , sep = "\t")
      Sequences <- paste0(as.character(unlist(Fastadata)) , collapse = "\n")
      write.table(x = Sequences , file = paste0(FileName ,".fasta") , quote = F , row.names = F , col.names = F, append = T)
    }
    
  }
  return(OutNumber)
}


PlotCharge <- function(SeqDataObjPath)
{
  seqdata <- select(SeqDataObjPath , "Sequence")
  result = aminoAcidProperties(seqdata , seq="Sequence")
  #Group data frame by charge
  ChargeGroup <- rep("Negative" , dim(result)[1])
  positive_index <- which(result$Sequence_aa_charge > 0)
  ChargeGroup[positive_index] <- "Positive"
  result <- cbind(result , ChargeGroup)
  
  #Get charge Sign ratios
  Chargecount <- table(sign(result$Sequence_aa_charge))
  if (length(Chargecount) == 1){
    Chargecount <- as.table(cbind(Chargecount , 0));
  }
  Chargeratio <- table(sign(result$Sequence_aa_charge))/dim(result)[1]*100
  if (length(Chargeratio) == 1){
    Chargeratio <- as.table(cbind(Chargeratio , 0));
  }
  #Construct Charge dataframe
  Chargedf = data.frame(x = c("Negative" , "Positive") , y = Chargeratio , z = Chargecount)
  #Charge plot
  Chargebarplot <- ggplot(result, aes(x=as.numeric(reorder(rownames(result) , result$Sequence_aa_charge)), y=result$Sequence_aa_charge, label=rownames(result))) +
    geom_bar(stat='identity', aes(fill=result$ChargeGroup)) + theme_classic() + theme(axis.title.x = element_blank()  , axis.ticks.x = element_blank() , plot.title = element_text(hjust = 0.5))+ylab("Protein charge")+
    guides(fill=guide_legend(title="Groups"))+ggtitle("Sequence Charge") + scale_x_continuous(limits = c(0,  dim(result)[1]), expand = c(0, 0))
  Chargebarplot
  
  Chargeframeplot <- ggplot(Chargedf , aes(x = Chargedf$x , y = Chargedf$z.Freq))+
    geom_bar(stat = "identity" , aes(fill = Chargedf$y.Freq)) + geom_text(aes(label = paste0(as.character(round(Chargedf$y.Freq),2), "%")) , size = 4, vjust = -1)+theme_bw()+
    theme(legend.position = "none" , plot.title = element_text(hjust = 0.5))+xlab("Groups") + ylab("Protein count") + ggtitle("Sequence Charge") + scale_y_continuous(limits = c(0,  dim(result)[1]), expand = c(0, 1))
  
  ChargeTotal <- ggarrange(Chargeframeplot, Chargebarplot , nrow = 1 , ncol = 2, align = "hv")
  plot(ChargeTotal)
}
PlotGravy <- function(SeqDataObjPath)
{
  seqdata <- select(SeqDataObjPath , "Sequence")
  result = aminoAcidProperties(seqdata , seq="Sequence")
  #Group data frame by GRAVY charge
  GravyGroup <- rep("Negative" , dim(result)[1])
  positive_index <- which(result$Sequence_aa_gravy > 0)
  GravyGroup[positive_index] <- "Positive"
  result <- cbind(result , GravyGroup)
  #Get GRAVY sign ratios
  GRAVYcount <- table(sign(result$Sequence_aa_gravy))
  if (length(GRAVYcount) == 1){
    GRAVYcount <- as.table(cbind(GRAVYcount , 0));
  }
  
  GRAVYratio <-  table(sign(result$Sequence_aa_gravy))/dim(result)[1]*100
  if (length(GRAVYratio) == 1){
    GRAVYratio <- as.table(cbind(GRAVYratio , 0));
  }
  #Construct GRAVY dataframe
  GRAVYdf = data.frame(x = c("Negative" , "Positive") , y = GRAVYratio , z = GRAVYcount)
  #GRAVY plot
  GRAVYbarplot <- ggplot(result, aes(x=as.numeric(reorder(rownames(result), result$Sequence_aa_gravy)), y=result$Sequence_aa_gravy, label=rownames(result))) +
    geom_bar(stat='identity', aes(fill=result$GravyGroup)) + theme_classic() + theme(axis.title.x = element_blank() , axis.ticks.x = element_blank() , plot.title = element_text(hjust = 0.5))+ylab("GRAVY index")+
    guides(fill=guide_legend(title="Groups"))+ggtitle("Sequence GRAVY index") + scale_x_continuous(limits = c(0,  dim(result)[1]), expand = c(0, 0))
  
  GRAVYframeplot <- ggplot(GRAVYdf , aes(x = GRAVYdf$x , y = GRAVYdf$z.Freq))+
    geom_bar(stat = "identity" , aes(fill = GRAVYdf$y.Freq)) + geom_text(aes(label = paste0(as.character(round(GRAVYdf$y.Freq),2), "%")) , size = 4, vjust = -0.6)+theme_bw()+
    theme(legend.position = "none" , plot.title = element_text(hjust = 0.5))+xlab("Groups") + ylab("Protein count") + ggtitle("GRAVY index") + scale_y_continuous(limits = c(0,  dim(result)[1]), expand = c(0, 1))
  
  GravyTotal <- ggarrange(GRAVYframeplot, GRAVYbarplot , nrow = 1 , ncol = 2, align = "hv")
  plot(GravyTotal)
}

PlotAcidity <- function(SeqDataObjPath)
{
  seqdata <- select(SeqDataObjPath , "Sequence")
  result = aminoAcidProperties(seqdata , seq="Sequence")
  #Plot acidic vs Basic Groups
  AcidicBasic <- c(result$Sequence_aa_acidic , result$Sequence_aa_basic)
  AcidicBasicGroup <- rep("Acidic" , length(result$Sequence_aa_acidic))
  AcidicBasicGroup <- c(AcidicBasicGroup , rep("Basic" , length(result$Sequence_aa_basic)))
  AcidicBasicframe <- data.frame(x = AcidicBasic , y = AcidicBasicGroup)
  
  
  p <- ggplot(AcidicBasicframe, aes(x=AcidicBasicframe$y, y=AcidicBasicframe$x, fill=AcidicBasicframe$y)) +
    geom_violin(alpha = 0.3)+ geom_boxplot(width=0.1) +
    guides(fill=guide_legend(title="Groups"))+ xlab("Groups") + ylab("Hydrophobicity")+
    theme_bw() + ggtitle("Hydrophobicity") + theme(plot.title = element_text(hjust = 0.5))
  plot(p)
}

ConstructPhylogeny <- function(ProteinDataObject)
{
  #Get sequences
  Seqlist <- as.character(ProteinDataObject$Sequence)
  
  #Apply multiple sequence alignment
  MSAresult <- msa(Seqlist , type = "protein")
  #Convert algniment results to tree
  MSAtree <- msaConvert(MSAresult, type="seqinr::alignment")
  # generate a distance matrix using seqinr package
  MSAdistance <- dist.alignment(MSAtree, "identity")
  #neighbor-joining tree estimation
  mTree <- nj(MSAdistance)
  mTree$tip.label <- as.character(rownames(ProteinDataObject))
  
  # pile up the functions to make a new tree
  NTree <- nj(dist.alignment(MSAtree, "identity"))
  NTree$tip.label <- mTree$tip.label
  #Plot new tree
  plot(NTree, "f", FALSE, cex = 0.7 , main="Phylogenetic Tree of proteins")
  #RegularTree <- plot(mTree)
}