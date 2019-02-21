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
library(lmerTest)

setwd("/Users/Pinedasans/ImmuneRep_Pregnancy/")
load("Data/BCR/BCR_data_summary.RData")

###Function to plot summary plots with bar plots
####Summary plots
brewer.pal(n = 3, name = "Accent")
COLOR=c("#BEAED4","#7FC97F")
summary_data<-summary_data[order(rownames(summary_data)),]

##Clones
j<-1
p_adj<-NULL
for (i in c("clones_IGHA","clones_IGHD","clones_IGHE","clones_IGHG","clones_IGHM","clones_memory","clones_naive")){
  print(i)
  results<-coefficients(summary(lmer(summary_data[,i] ~ summary_data$sample + summary_data$NumCells + (1| summary_data$pairs))))
  p_adj[j]<-results[2,5]
  tiff(paste0("Results/boxplot_",i,".tiff"),res=300,w=1500,h=2000)
  g<-ggpaired(summary_data, x = "sample", y = i,
              color = "sample", line.color = "gray", line.size = 0.4,
              palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
    labs(title=i,x="Sample", y = "Number of clones") + theme(legend.position="none") +
    stat_compare_means(method = "wilcox.test",paired = TRUE)
  print(g)
  dev.off()
  j<-j+1
}
names(p_adj)<-c("clones_IGHA","clones_IGHD","clones_IGHE","clones_IGHG","clones_IGHM","clones_memory","clones_naive")
write.csv(p_adj,"Results/Diversity/p_adj_clones.csv")

plot(summary_data[,"clones_IGHM"], 
     summary_data[,"NumCells"],
     col = COLOR,pch=20,ylab = "Num cells",xlab = "Clones IGHM")
legend("topright",legend=c("Maternal","Fetal"), 
       col=COLOR,pch=20,cex=c(1.2),ncol=2)

##Entropy
j<-1
p_adj<-NULL
for (i in c("entropy_IGHA","entropy_IGHD","entropy_IGHE",
            "entropy_IGHG","entropy_IGHM","entropy_memory","entropy_naive")){
  print(i)
  results<-coefficients(summary(lmer(summary_data[,i] ~ summary_data$sample + summary_data$NumCells + (1| summary_data$pairs))))
  p_adj[j]<-results[2,5]
  tiff(paste0("Results/boxplot_",i,".tiff"),res=300,w=1500,h=2000)
  g<-ggpaired(summary_data, x = "sample", y = i,
              color = "sample", line.color = "gray", line.size = 0.4,
              palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
    labs(title=i,x="Sample", y = "Entropy") + theme(legend.position="none") +
    stat_compare_means(method = "wilcox.test",paired = TRUE)
  print(g)
  dev.off()
  j<-j+1
}
names(p_adj)<-c("entropy_IGHA","entropy_IGHD","entropy_IGHE","entropy_IGHG","entropy_IGHM","entropy_memory","entropy_naive")
write.csv(p_adj,"Results/Diversity/p_adj_entropy.csv")

tiff("Results/plot_IGHM_numcell.tiff",res=300,w=2000,h=2000)
summary(glm(summary_data[,"entropy_IGHM"]~summary_data[,"NumCells"]))
plot(summary_data[,"entropy_IGHM"], 
     summary_data[,"NumCells"],
     col = COLOR,cex=1.5,pch=20,ylab = "Num cells",xlab = "entropy IGHM (p=0.007)")
legend("topright",legend=c("Maternal","Fetal"), 
       col=COLOR,pch=20,cex=c(1.5),ncol=2)
dev.off()

##Clonality
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

##SHM
j<-1
p_adj<-NULL
for (i in c("SHM_IGHA","SHM_IGHD","SHM_IGHG","SHM_IGHM")){
  print(i)
  results<-coefficients(summary(lmer(summary_data[,i] ~ summary_data$sample + summary_data$NumCells + (1| summary_data$pairs))))
  p_adj[j]<-results[2,5]
  tiff(paste0("Results/boxplot_",i,".tiff"),res=300,w=1500,h=2000)
  g<-ggpaired(summary_data, x = "sample", y = i,
              color = "sample", line.color = "gray", line.size = 0.4,
              palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
    labs(title=i,x="Sample", y = "SHM") + theme(legend.position="none") +
    stat_compare_means(method = "wilcox.test",paired = TRUE)
  print(g)
  dev.off()
  j<-j+1
}
names(p_adj)<-c("SHM_IGHA","SHM_IGHD","SHM_IGHG","SHM_IGHM")
write.csv(p_adj,"Results/Diversity/p_adj_SHM.csv")

tiff("Results/plot_SHM_IGHM_numcell.tiff",res=300,w=2000,h=2000)
summary(glm(summary_data[,"SHM_IGHM"]~summary_data[,"NumCells"]))
plot(summary_data[,"SHM_IGHM"], 
     summary_data[,"NumCells"],
     col = COLOR,cex=1.5,pch=20,ylab = "Num cells",xlab = "SHM_IGHM (p=0.6)")
legend("topright",legend=c("Maternal","Fetal"), 
       col=COLOR,pch=20,cex=c(1.5),ncol=2)
dev.off()


##CDR3
j<-1
p_adj<-NULL
for (i in c("CDR3_length_IGHA","CDR3_length_IGHD","CDR3_length_IGHG","CDR3_length_IGHM","CDR3_length_memory","CDR3_length_naive")){
  print(i)
  results<-coefficients(summary(lmer(summary_data[,i] ~ summary_data$sample + summary_data$NumCells + (1| summary_data$pairs))))
  p_adj[j]<-results[2,5]
  tiff(paste0("Results/boxplot_",i,".tiff"),res=300,w=1500,h=2000)
  g<-ggpaired(summary_data, x = "sample", y = i,
              color = "sample", line.color = "gray", line.size = 0.4,
              palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
    labs(title=i,x="Sample", y = "CDR3 length") + theme(legend.position="none") +
    stat_compare_means(method = "wilcox.test",paired = TRUE)
  print(g)
  dev.off()
  j<-j+1
}
names(p_adj)<-c("CDR3_length_IGHA","CDR3_length_IGHD","CDR3_length_IGHG","CDR3_length_IGHM","CDR3_length_memory","CDR3_length_naive")
write.csv(p_adj,"Results/Diversity/p_adj_cdr3.csv")

tiff("Results/plot_CDR3_length_naive_numcell.tiff",res=300,w=2000,h=2000)
summary(glm(summary_data[,"CDR3_length_naive"]~summary_data[,"NumCells"]))
plot(summary_data[,"CDR3_length_naive"], 
     summary_data[,"NumCells"],
     col = COLOR,cex=1.5,pch=20,ylab = "Num cells",xlab = "CDR3_length_naive (p=0.7)")
legend("topright",legend=c("Maternal","Fetal"), 
       col=COLOR,pch=20,cex=c(1.5),ncol=2)
dev.off()



########################################
######### TCR diversity analysis ######
######################################
load("Data/TCR/TCR_data_summary.RData")
info_samples<-read.csv("Data/NormalSamplesScottBoydStanford9_13_18.csv")
id<-match(info_samples$Sample_ID,rownames(summary_data))
summary_data$pairs[id]<-info_samples$Pairs
summary_data$NumCells[id]<-info_samples$NumCells
summary_data$sample<-gsub("\\*", "", substr(rownames(summary_data), 1, 1))
summary_data$sample<-ifelse(summary_data$sample=="M","Maternal","Fetal")
summary_data$sample<-factor(summary_data$sample)
summary_data$sample<-relevel(summary_data$sample, ref = "Maternal")


###Function to plot summary plots with bar plots
####Summary plots
brewer.pal(n = 3, name = "Accent")
COLOR=c("#BEAED4","#7FC97F")

summary_data<-summary_data[order(rownames(summary_data)),]
tiff("Results/boxplot_clones_TCR.tiff",res=300,w=1500,h=2000)
  g<-ggpaired(summary_data, x = "sample", y = "clones",
              color = "sample", line.color = "gray", line.size = 0.4,
              palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
    labs(title="TCR",x="Sample", y = "Number of clones") + theme(legend.position="none") +
    stat_compare_means(method = "wilcox.test",paired = TRUE)
  print(g)
  dev.off()

summary(lmer(summary_data[,"clones"] ~ summary_data$sample + summary_data$NumCells + (1| summary_data$pairs)))

  
tiff("Results/boxplot_entropy_TCR.tiff",res=300,w=1500,h=2000)
  g<-ggpaired(summary_data, x = "sample", y = "entropy",
              color = "sample", line.color = "gray", line.size = 0.4,
              palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
    labs(title="TCR",x="Sample", y = "Entropy") + theme(legend.position="none") +
    stat_compare_means(method = "wilcox.test",paired = TRUE)
  print(g)
  dev.off()

  summary(lmer(summary_data[,"entropy"] ~ summary_data$sample + summary_data$NumCells + (1| summary_data$pairs)))
  
  
tiff("Results/boxplot_clonality_TCR.tiff",res=300,w=1500,h=2000)
  g<-ggpaired(summary_data, x = "sample", y = "clonality",
              color = "sample", line.color = "gray", line.size = 0.4,
              palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
    labs(title="TCR",x="Sample", y = "clonality") + theme(legend.position="none") +
    stat_compare_means(method = "wilcox.test",paired = TRUE)
  print(g)
  dev.off()

  
  tiff("Results/boxplot_mean_CDR3_length_TCR.tiff",res=300,w=1500,h=2000)
  g<-ggpaired(summary_data, x = "sample", y = "mean_CDR3_length",
              color = "sample", line.color = "gray", line.size = 0.4,
              palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
    labs(title="TCR",x="Sample", y = "mean_CDR3_length") + theme(legend.position="none") +
    stat_compare_means(method = "wilcox.test",paired = TRUE)
  print(g)
  dev.off()
  
  summary(lmer(summary_data[,"mean_CDR3_length"] ~ summary_data$sample + summary_data$NumCells + (1| summary_data$pairs))) #0.057
  



