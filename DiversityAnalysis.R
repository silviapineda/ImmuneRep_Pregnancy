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
### DESCRIP: Diversity Analysis
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
load("Data/BCR/BCR_data_summary.RData")

###Function to plot summary plots with bar plots
####Summary plots
brewer.pal(n = 3, name = "Accent")
COLOR=c("#BEAED4","#7FC97F")

summary_data<-summary_data[order(summary_data$pairs),]
for (i in c("clones_IGHA","clones_IGHD","clones_IGHE","clones_IGHG","clones_IGHM","clones_memory","clones_naive")){
  print(i)
  tiff(paste0("Results/boxplot_",i,".tiff"),res=300,w=1500,h=2000)
  g<-ggpaired(summary_data, x = "sample", y = i,
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
  g<-ggpaired(summary_data, x = "sample", y = i,
              color = "sample", line.color = "gray", line.size = 0.4,
              palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
    labs(title=i,x="Sample", y = "Entropy") + theme(legend.position="none") +
    stat_compare_means(method = "wilcox.test",paired = TRUE)
  print(g)
  dev.off()
}

for (i in c("clonality_IGHA","clonality_IGHD",
            "clonality_IGHG","clonality_IGHM","clonality_memory","clonality_naive")){
  print(i)
  tiff(paste0("Results/boxplot_",i,".tiff"),res=300,w=1500,h=2000)
  g<-ggpaired(summary_data, x = "sample", y = i,
              color = "sample", line.color = "gray", line.size = 0.4,
              palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
    labs(title=i,x="Sample", y = "clonality") + theme(legend.position="none") +
    stat_compare_means(method = "wilcox.test",paired = TRUE)
  print(g)
  dev.off()
}

for (i in c("SHM_IGHA","SHM_IGHD","SHM_IGHE","SHM_IGHG","SHM_IGHM")){
  print(i)
  tiff(paste0("Results/boxplot_",i,".tiff"),res=300,w=1500,h=2000)
  g<-ggpaired(summary_data, x = "sample", y = i,
              color = "sample", line.color = "gray", line.size = 0.4,
              palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
    labs(title=i,x="Sample", y = "SHM") + theme(legend.position="none") +
    stat_compare_means(method = "wilcox.test",paired = TRUE)
  print(g)
  dev.off()
}

########################################
######### TCR diversity analysis ######
######################################
load("Data/TCR/TCR_data_summary.RData")

###Function to plot summary plots with bar plots
####Summary plots
brewer.pal(n = 3, name = "Accent")
COLOR=c("#BEAED4","#7FC97F")

summary_data<-summary_data[order(summary_data$pairs),]
tiff("Results/boxplot_clones_TCR.tiff",res=300,w=1500,h=2000)
  g<-ggpaired(summary_data, x = "sample", y = "clones",
              color = "sample", line.color = "gray", line.size = 0.4,
              palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
    labs(title="TCR",x="Sample", y = "Number of clones") + theme(legend.position="none") +
    stat_compare_means(method = "wilcox.test",paired = TRUE)
  print(g)
  dev.off()

tiff("Results/boxplot_entropy_TCR.tiff",res=300,w=1500,h=2000)
  g<-ggpaired(summary_data, x = "sample", y = "entropy",
              color = "sample", line.color = "gray", line.size = 0.4,
              palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
    labs(title="TCR",x="Sample", y = "Entropy") + theme(legend.position="none") +
    stat_compare_means(method = "wilcox.test",paired = TRUE)
  print(g)
  dev.off()

tiff("Results/boxplot_clonality_TCR.tiff",res=300,w=1500,h=2000)
  g<-ggpaired(summary_data, x = "sample", y = "clonality",
              color = "sample", line.color = "gray", line.size = 0.4,
              palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
    labs(title="TCR",x="Sample", y = "clonality") + theme(legend.position="none") +
    stat_compare_means(method = "wilcox.test",paired = TRUE)
  print(g)
  dev.off()






