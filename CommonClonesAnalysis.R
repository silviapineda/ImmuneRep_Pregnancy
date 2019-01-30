rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies for pregnancy outcomes
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: Common Clones across Mother/Fetus
###         
###
### Author: Silvia Pineda
### Date: January, 2019
############################################################################################

library(pheatmap)
setwd("/Users/Pinedasans/ImmuneRep_Pregnancy/")
load("Data/BCR/BCR_data_summary.RData")

###################################################
## Common clones across mother and Fetus (BCR) ##
#################################################
isotype="IGHM"

############
## 1. Build the matrix with the clones by samples
data_qc_isotype<-data_qc[which(data_qc$isotype==isotype),]
clone_type<-t(as.data.frame(unclass(table(data_qc_isotype$V_J_lenghCDR3_Clone_igh,factor(data_qc_isotype$sample_label))))) # 31,547 (IGHA) #364755 (IGHD) #25591 (IGHG) #599338 (IGHM)
##Build present vs no present
clone_type_presence<-apply(clone_type,1,function(x) ifelse(x==0,0,1))
###Filter by clones that at least are share in 2 samples
clone_type_filter<-clone_type_presence[which(rowSums(clone_type_presence)>1),] #None (IGHA) 100 (IGHD) 1 (IGHG) 289 (IGHM)
tiff(paste0("Results/heatmap_common_clones_",isotype),width = 6000, height = 2000, res = 300)
pheatmap(t(clone_type_filter),cluster_rows = F,color = colorRampPalette(c("white", "red"))(50))
dev.off()

###################################################
## Common clones across mother and Fetus (TCR) ##
#################################################

load("Data/TCR/TCR_data_summary.RData")
############
## 1. Build the matrix with the clones by samples
clone_type<-t(as.data.frame(unclass(table(data_qc$V_J_length_CDR3_Clone_tcrb,factor(data_qc$sample_label))))) # 487669 (TCR) #1028849 (BCR)
##Build present vs no present
clone_type_presence<-apply(clone_type,1,function(x) ifelse(x==0,0,1))
###Filter by clones that at least are share in 2 samples
clone_type_filter<-clone_type_presence[which(rowSums(clone_type_presence)>1),] #2677 (TCR) #829 (BCR)
tiff("Results/heatmap_common_clones_TCR",width = 6000, height = 2000, res = 300)
pheatmap(t(clone_type_filter),cluster_rows = F,color = colorRampPalette(c("white", "red"))(50))
dev.off()

