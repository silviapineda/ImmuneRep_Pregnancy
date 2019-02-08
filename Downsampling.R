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
### DESCRIP: Downsapmling
###         
###
### Author: Silvia Pineda
### Date: January, 2019
############################################################################################
library(ggplot2)
library("igraph")
library("ineq")
library("dplyr")
library("RColorBrewer")
library(lme4)
library(ggpubr)

setwd("/Users/Pinedasans/ImmuneRep_Pregnancy/")

################################################
####Prepare Downsampling data for clone calling 
#############################################

##BCR
load("Data/BCR/BCR_data_summary.RData")

###Minimum number of reads to do downsampling
min_reads<-min(summary_data$read_count) #291635
##Down-sampling to 291635 reads
sample_id<-unique(data_qc$sample_label)
for (b in 1:10){
  print(b)
  data_qc_downsampled<-NULL
  for (i in sample_id){
    print(i)
    data_qc_downsampled_reads<-data_qc[which(data_qc$sample_label==i),]
    id_down<-sample(1:dim(data_qc_downsampled_reads)[1],min_reads)
    data_qc_downsampled<-rbind(data_qc_downsampled,data_qc_downsampled_reads[id_down,])
  }
  save(data_qc_downsampled,file=paste0("Data/BCR/data_downsampled_",b,".Rdata"))
  
  
    
}

###Prepare the data for the python script using python3
for (i in 1:10){
  load(paste0("Data/BCR/data_downsampled_",i,".Rdata"))
  data_qc_downsampled$V_J_lenghCDR3 = paste(data_qc_downsampled$v_gene, data_qc_downsampled$j_gene, nchar(data_qc_downsampled$cdr3_seq),sep="_")
  data_qc_downsampled$unique_id<-seq(1,nrow(data_qc_downsampled))
  data_clonesInference<-data_qc_downsampled[,c("unique_id","sample_label","cdr3_seq","CDR3_length","v_gene","j_gene","V_J_lenghCDR3")]
  write.table(data_clonesInference,file=paste0("Data/BCR/data_clonesInference_down",i,".txt"),row.names = F,sep="\t")
}  


##TCR
load("Data/TCR/TCR_data_summary.RData")

###Minimum number of reads to do downsampling
min_reads<-min(summary_data$read_count) #219625
##Down-sampling to 219625 reads
sample_id<-unique(data_qc$sample_label)
for (b in 1:10){
  print(b)
  data_qc_downsampled<-NULL
  for (i in sample_id){
    print(i)
    data_qc_downsampled_reads<-data_qc[which(data_qc$sample_label==i),]
    id_down<-sample(1:dim(data_qc_downsampled_reads)[1],min_reads)
    data_qc_downsampled<-rbind(data_qc_downsampled,data_qc_downsampled_reads[id_down,])
  }
  data_qc_downsampled$V_J_lenghCDR3 = paste(data_qc_downsampled$v_gene, data_qc_downsampled$j_gene, nchar(data_qc_downsampled$cdr3_seq),sep="_")
  data_qc_downsampled$unique_id<-seq(1,nrow(data_qc_downsampled))
  data_clonesInference<-data_qc_downsampled[,c("unique_id","sample_label","cdr3_seq","CDR3_length","v_gene","j_gene","V_J_lenghCDR3")]
  write.table(data_clonesInference,file=paste0("Data/TCR/data_clonesInference_down",b,".txt"),row.names = F,sep="\t")
}



####################
###Read the data from the python script and prepare the summary matrix
####################
for (i in c(1:2,4:6)){
  ClonesInfered_downsampled<-read.csv(paste0("Data/BCR/ClonesInfered_down",i,".csv"))
  load(paste0("Data/BCR/data_downsampled_",i,".Rdata"))
  data_qc_downsampled$unique_id<-seq(1,nrow(data_qc_downsampled))
  id<-match(ClonesInfered_downsampled$unique_id,data_qc_downsampled$unique_id)
  ClonesInfered_downsampled$isotype<-data_qc_downsampled$isotype[id]
  ClonesInfered_downsampled$SHM_freq<-data_qc_downsampled$SHM_freq[id]
  ClonesInfered_downsampled$SHM<-data_qc_downsampled$SHM[id]
  assign(paste0("summary_data_down",i),summaryData(ClonesInfered_downsampled))
}


####Mean Clones
summary_clones_down<-matrix(NA,20,7)
k=1
for (i in c("clones_IGHA","clones_IGHD","clones_IGHE","clones_IGHG","clones_IGHM","clones_memory","clones_naive")){
  clones_down<-matrix(NA,20,5)
  for (j in c(1:2,4:5)){
    clones_down[,j]<-get(paste0("summary_data_down",j))[,i]
  }
  clones_down<-clones_down[,-3]
  summary_clones_down[,k]<-apply(clones_down,1,mean)
  k=k+1
}
colnames(summary_clones_down)<-c("clones_IGHA","clones_IGHD","clones_IGHE","clones_IGHG","clones_IGHM","clones_memory","clones_naive")
rownames(summary_clones_down)<-rownames(summary_data_down1)

####Mean Entropy
summary_entropy_down<-matrix(NA,20,7)
k=1
for (i in c("entropy_IGHA","entropy_IGHD","entropy_IGHE",
            "entropy_IGHG","entropy_IGHM","entropy_memory","entropy_naive")){
  entropy_down<-matrix(NA,20,5)
  for (j in c(1:2,4:5)){
    entropy_down[,j]<-get(paste0("summary_data_down",j))[,i]
  }
  entropy_down<-entropy_down[,-3]
  summary_entropy_down[,k]<-apply(entropy_down,1,mean)
  k=k+1
}
colnames(summary_entropy_down)<-c("entropy_IGHA","entropy_IGHD","entropy_IGHE","entropy_IGHG","entropy_IGHM","entropy_memory","entropy_naive")
rownames(summary_entropy_down)<-rownames(summary_data_down1)

####SHM
summary_SHM_down<-matrix(NA,20,5)
k=1
for (i in c("SHM_IGHA","SHM_IGHD","SHM_IGHE","SHM_IGHG","SHM_IGHM")){
  SHM_down<-matrix(NA,20,5)
  for (j in c(1:2,4:5)){
    SHM_down[,j]<-get(paste0("summary_data_down",j))[,i]
  }
  SHM_down<-SHM_down[,-3]
  summary_SHM_down[,k]<-apply(SHM_down,1,mean)
  k=k+1
}
colnames(summary_SHM_down)<-c("SHM_IGHA","SHM_IGHD","SHM_IGHE","SHM_IGHG","SHM_IGHM")
rownames(summary_SHM_down)<-rownames(summary_data_down1)

summary_data_down<-cbind(summary_clones_down,summary_entropy_down,summary_SHM_down)
summary_data_down<-data.frame(summary_data_down)
summary_data_down$sample<-summary_data_down1$sample
summary_data_down$pairs<-summary_data_down1$pairs
summary_data_down$NumCells<-summary_data_down1$NumCells
save(summary_data_down,file="Data/BCR/BCR_data_summary_down.RData")

#####
###Statistical Analysis on the downsampling
#####
load("Data/BCR/BCR_data_summary_down.RData")

###Function to plot summary plots with bar plots
####Summary plots
brewer.pal(n = 3, name = "Accent")
COLOR=c("#BEAED4","#7FC97F")

summary_data_down<-summary_data_down[order(summary_data_down$sample),]
for (i in c("clones_IGHA","clones_IGHD","clones_IGHE","clones_IGHG","clones_IGHM","clones_memory","clones_naive")){
  print(i)
  tiff(paste0("Results/boxplot_",i,".tiff"),res=300,w=1500,h=2000)
    g<-ggpaired(summary_data_down, x = "sample", y = i,
      color = "sample", line.color = "gray", line.size = 0.4,
      palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
      labs(title=i,x="Sample", y = "Number of clones") + theme(legend.position="none") +
      stat_compare_means(method = "wilcox.test",paired = TRUE)
    print(g)
  dev.off()
}

for (i in c("entropy_IGHA","entropy_IGHD","entropy_IGHE",
            "entropy_IGHG","entropy_IGHM","entropy_memory","entropy_naive")){
  print(i)
  tiff(paste0("Results/boxplot_",i,".tiff"),res=300,w=1500,h=2000)
  g<-ggpaired(summary_data_down, x = "sample", y = i,
               color = "sample", line.color = "gray", line.size = 0.4,
               palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
    labs(title=i,x="Sample", y = "Entropy") + theme(legend.position="none") +
    stat_compare_means(method = "wilcox.test",paired = TRUE)
  print(g)
  dev.off()
}

for (i in c("SHM_IGHA","SHM_IGHD","SHM_IGHE","SHM_IGHG","SHM_IGHM")){
  print(i)
  tiff(paste0("Results/boxplot_",i,".tiff"),res=300,w=1500,h=2000)
  g<-ggpaired(summary_data_down, x = "sample", y = i,
              color = "sample", line.color = "gray", line.size = 0.4,
              palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
    labs(title=i,x="Sample", y = "SHM") + theme(legend.position="none") +
    stat_compare_means(method = "wilcox.test",paired = TRUE)
  print(g)
  dev.off()
}


