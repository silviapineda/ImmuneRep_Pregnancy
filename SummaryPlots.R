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
### DESCRIP: Analysis of the VDJ data
###         
###
### Author: Silvia Pineda
### Date: December, 2018
############################################################################################
library("RColorBrewer")

setwd("/Users/Pinedasans/ImmuneRep_Pregnancy/")
load("Data/BCR_data_summary.RData")

###Function to plot summary plots with bar plots
plot <- function(statistic,color,sample,pair,label){
  barplot(statistic,col=color[sample],
          names.arg = pair,main=label,xlab = "Samples", ylab = label,las=2)
  legend("topright", legend=levels(sample),col=color[factor(sample)],pch=15, cex=0.8)
}  

####Summary plots
brewer.pal(n = 3, name = "Accent")
COLOR=c("#BEAED4","#7FC97F")

###Reads
for (i in c("UNMAPPED","IGHA","IGHD","IGHE","IGHG","IGHM","reads_memory","reads_naive")){
  print(i)
  tiff(paste0("Results/barplot_",i,".tiff"),res=300,w=4000,h=2000)
  plot(summary_data[,i],COLOR,summary_data$sample,summary_data$pairs,label=i)
  dev.off()
}
###Clones
for (i in c("clones_unmapped","clones_IGHA","clones_IGHD","clones_IGHE",
            "clones_IGHG","clones_IGHM","clones_memory","clones_naive")){
  print(i)
  tiff(paste0("Results/barplot_",i,".tiff"),res=300,w=4000,h=2000)
  plot(summary_data[,i],COLOR,summary_data$sample,summary_data$pairs,label=i)
  dev.off()
}

###Entropy
for (i in c("entropy_unmapped","entropy_IGHA","entropy_IGHD","entropy_IGHE",
            "entropy_IGHG","entropy_IGHM","entropy_memory","entropy_naive")){
  print(i)
  tiff(paste0("Results/barplot_",i,".tiff"),res=300,w=4000,h=2000)
  plot(summary_data[,i],COLOR,summary_data$sample,summary_data$pairs,label=i)
  dev.off()
}

###SHM
for (i in c("SHM_unmapped","SHM_IGHA","SHM_IGHD","SHM_IGHE",
            "SHM_IGHG","SHM_IGHM")){
  print(i)
  tiff(paste0("Results/barplot_",i,".tiff"),res=300,w=4000,h=2000)
  plot(summary_data[,i],COLOR,summary_data$sample,summary_data$pairs,label=i)
  dev.off()
}

###CDR3 length
for (i in c("CDR3_length_unmapped","CDR3_length_IGHA","CDR3_length_IGHD","CDR3_length_IGHE",
            "CDR3_length_IGHG","CDR3_length_IGHM")){
  print(i)
  tiff(paste0("Results/barplot_",i,".tiff"),res=300,w=4000,h=2000)
  plot(summary_data[,i],COLOR,summary_data$sample,summary_data$pairs,label=i)
  dev.off()
}
