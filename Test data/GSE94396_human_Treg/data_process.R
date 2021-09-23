setwd("C:\\Users\\BUITT\\Dropbox\\Biotrans\\ABioTrans\\Test data\\GSE94396_human_Treg")
exp_mat <- read.table("GSE94396_iTreg_raw_counts.txt", header = TRUE)
colnames(exp_mat)[1] <- "gene_id"

rnames <- exp_mat$gene_id
rnames <- sapply(rnames, function(x) strsplit(x, "\\.")[[1]][1] )
exp_mat$gene_id <- rnames

cnames <- c("G01_T01_A",	"G01_T01_B", "G01_T01_D",	
            "G03_T02_A",	"G03_T02_B", "G03_T02_D",	
            "G03_T03_A",	"G03_T03_B", "G03_T03_D",	
            "G03_T04_A",	"G03_T04_B", "G03_T04_D",	
            "G03_T05_A",	"G03_T05_B", "G03_T05_D",	
            "G03_T06_A",	"G03_T06_B", "G03_T06_D")
# map with human gtf features
human_gtf <- read.csv("D:\\Biotrans\\Sandro\\data\\gff\\human_gtf_extract.csv", row.names = 1)
exp_mat_combined <- merge(exp_mat, human_gtf, by = "gene_id")
rownames(exp_mat_combined) <- exp_mat_combined$gene_id
head(exp_mat_combined)


write.csv(exp_mat_combined[,cnames], "GSE94396_human_Treg_raw_count.csv")
# write.csv(exp_mat_combined[,"length"], "GSE94396_gene_length.csv")
