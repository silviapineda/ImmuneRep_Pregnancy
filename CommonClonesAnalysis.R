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
library("RColorBrewer")

setwd("/Users/Pinedasans/ImmuneRep_Pregnancy/")
load("Data/BCR/BCR_data_summary.RData")

###################################################
## Common clones across mother and Fetus (BCR) ##
#################################################
isotype="IGHM"

############
## 1. Build the matrix with the clones by samples
data_qc_isotype<-data_qc[which(data_qc$isotype==isotype),]
clone_type<-t(as.data.frame(unclass(table(data_qc_isotype$V_J_lenghCDR3_Clone_igh,factor(data_qc_isotype$sample_label))))) # 31,547 (IGHA) #364755 (IGHD) #25591 (IGHG) #1 (IGHM)
##Build present vs no present
clone_type_presence<-apply(clone_type,1,function(x) ifelse(x==0,0,1))
###Filter by clones that at least are share in 2 samples
clone_type_filter<-clone_type_presence[which(rowSums(clone_type_presence)>1),] #None (IGHA) 100 (IGHD) 1 (IGHG) 289 (IGHM)
id<-match(rownames(clone_type_filter),colnames(clone_type))
write.csv(clone_type[,id], file = paste0("Results/common_clones_",isotype,".csv"))

brewer.pal(n = 3, name = "Accent")
COLOR=c("#BEAED4","#7FC97F")
fill2 = brewer.pal(10,"Set3")
#summary_data$sample<-ifelse(summary_data$sample=="Mother","Maternal","Fetal")
annotation.col <- data.frame(row.names = rownames(summary_data))
annotation.col$Type <- factor(summary_data$sample)
annotation.col$Pair <- factor(c("pair1","pair1","pair2","pair2","pair3","pair3","pair4","pair4","pair5","pair5","pair6","pair6",
                         "pair7","pair7","pair8","pair8","pair9","pair9","pair10","pair10"))

ann_colors = list(Type = c("Maternal" = COLOR [1] ,"Fetal" = COLOR [2]),
                  Pair = c("pair1" = fill2[1],"pair2" = fill2[2],"pair3" = fill2[3],"pair4" = fill2[4],"pair5" = fill2[5],
                           "pair6" = fill2[6],"pair7" = fill2[7],"pair8" = fill2[8],"pair9" = fill2[9],"pair10" = fill2[10]))


tiff(paste0("Results/heatmap_common_clones_",isotype),width = 2000, height = 4000, res = 300)
pheatmap(clone_type_filter,annotation = annotation.col,border_color=F, annotation_colors = ann_colors,color = colorRampPalette(c("white", "red"))(50))
dev.off()

###Shared among pairs
j<-1
k<-2
shared_pairs<-NULL
for (i in 1:10){
  shared_pairs[i]<- length(which(rowSums(clone_type_filter[,j:k])==2))
  j<-j+2
  k<-k+2
}
###Shared among mothers
shared_mother<- length(which(rowSums(clone_type_filter[,c(1,3,5,7,9,11,13,15,17,19)])==2))
shared_fetus<-length(which(rowSums(clone_type_filter[,c(2,4,6,8,10,12,14,16,18,20)])==2))

counts=matrix(data=c(7,8,8,13),nrow=2)
fisher.test(counts) #p-value = 0.7


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
id<-match(rownames(clone_type_filter),colnames(clone_type))
write.csv(clone_type[,id], file = "Results/common_clones_TCR.csv")

brewer.pal(n = 3, name = "Accent")
COLOR=c("#BEAED4","#7FC97F")
fill2 = brewer.pal(10,"Set3")
#summary_data$sample<-ifelse(summary_data$sample=="Mother","Maternal","Fetal")
annotation.col <- data.frame(row.names = rownames(summary_data))
annotation.col$Type <- factor(summary_data$sample)
annotation.col$Pair <- factor(c("pair1","pair1","pair2","pair2","pair3","pair3","pair4","pair4","pair5","pair5","pair6","pair6",
                                "pair7","pair7","pair8","pair8","pair9","pair9","pair10","pair10"))

ann_colors = list(Type = c("Maternal" = COLOR [1] ,"Fetal" = COLOR [2]),
                  Pair = c("pair1" = fill2[1],"pair2" = fill2[2],"pair3" = fill2[3],"pair4" = fill2[4],"pair5" = fill2[5],
                           "pair6" = fill2[6],"pair7" = fill2[7],"pair8" = fill2[8],"pair9" = fill2[9],"pair10" = fill2[10]))


tiff("Results/heatmap_common_clones_TCR",width = 2000, height = 4000, res = 300)
pheatmap(clone_type_filter,annotation = annotation.col,border_color=F, annotation_colors = ann_colors,color = colorRampPalette(c("white", "red"))(50))
dev.off()


###Shared among pairs
j<-1
k<-2
shared_pairs<-NULL
for (i in 1:10){
  shared_pairs[i]<- length(which(rowSums(clone_type_filter[,j:k])==2))
  j<-j+2
  k<-k+2
}

###Random shared
x<-1
shared_random<-NULL
for(j in c(1,3,5,7,9,11,13,15,17,19)){
  k<-2
  random_pairs<-NULL
  for (i in 1:10){
    random_pairs[i]<- length(which(rowSums(clone_type_filter[,c(j,k)])==2))
    k<-k+2
  }
  shared_random[x]<-sum(random_pairs)
  x<-x+1
}

hist(shared_random)
abline(v=152, col="red", lwd=2)

p=(sum(abs(shared_random) > abs(152)) + 1) / (length(shared_random) + 1) #0.36

###Shared among mothers
shared_mother<- length(which(rowSums(clone_type_filter[,c(1,3,5,7,9,11,13,15,17,19)])==2))
shared_fetus<-length(which(rowSums(clone_type_filter[,c(2,4,6,8,10,12,14,16,18,20)])==2))

