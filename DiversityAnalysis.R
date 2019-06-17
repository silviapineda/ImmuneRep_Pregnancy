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
library(finalfit)
library(gridExtra)
library(grid)

setwd("/Users/Pinedasans/ImmuneRep_Pregnancy/")
load("Data/BCR_PTB_Term/BCR_PTB_Term_data_summary.RData")

###Function to plot summary plots with bar plots
####Summary plots
cols=brewer.pal(3,name = "Accent")

###Delete the twin for analysis
summary_data<-summary_data[-c(39:41),]
summary_data<-summary_data[order(rownames(summary_data)),]

##Clones
markers<-c("clones_IGHA","clones_IGHD","clones_IGHE","clones_IGHG","clones_IGHM","clones_mutated_IGHM","clones_unmutated_IGHM",
          "clones_mutated_IGHD","clones_unmutated_IGHD")

p_PTB<-NULL
p_Fetal<-NULL
for (i in 1:length(markers)){
  print(i)
  results<-coefficients(summary(glm(summary_data[,markers[i]] ~ summary_data$Outcome +summary_data$sample)))
  p_PTB[i]<-results[2,4]
  p_Fetal[i]<-results[3,4]
  summary_data$Marker<-summary_data[,markers[i]]
  tiff(paste0("Results/boxplot_Term_PTB",i,".tiff"),res=300,w=1500,h=2000)
  ggplot(summary_data) + geom_boxplot(aes(y = Marker, x = sample, fill = Outcome, color=Outcome),alpha = 0, position = position_dodge(width = .8)) +
    geom_point(aes(y=Marker, color=Outcome,x=sample,fill=Outcome), position = position_jitterdodge(dodge.width = 0.8)) +     
    scale_color_manual(values = c(cols[1], cols[2]), labels = c("Term", "PTB")) +
    scale_y_continuous(name=markers[i]) + stat_compare_means(aes(x=sample, y=Marker, color=Outcome,fill=Outcome,label.x.npc = "center"))

    
  
  
  #g<-ggpaired(summary_data, x = "Outcome", y = i,
  #            color = "Outcome", line.color = "gray", line.size = 0.4,
   #           palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
   # labs(title=i,x="Outcome", y = "Number of clones") + theme(legend.position="none") +
  #  stat_compare_means(method = "wilcox.test")
  #print(g)
  #dev.off()
  j<-j+1
}
names(p_adj)<-c("clones_IGHA","clones_IGHD","clones_IGHE","clones_IGHG","clones_IGHM","clones_mutated_IGHM","clones_unmutated_IGHM",
                "clones_mutated_IGHD","clones_unmutated_IGHD")
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
            "entropy_IGHG","entropy_IGHM","entropy_mutated_IGHM","entropy_unmutated_IGHM",
            "entropy_mutated_IGHD","entropy_unmutated_IGHD")){
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
names(p_adj)<-c("entropy_IGHA","entropy_IGHD","entropy_IGHE","entropy_IGHG","entropy_IGHM","entropy_mutated_IGHM","entropy_unmutated_IGHM",
                "entropy_mutated_IGHD","entropy_unmutated_IGHD")
write.csv(p_adj,"Results/Diversity/p_adj_entropy.csv")

tiff("Results/plot_entropy_IGHD_numcell.tiff",res=300,w=2000,h=2000)
ggplot(summary_data,aes(entropy_IGHD,NumCells,color=sample)) + 
  geom_point() + geom_smooth(method='lm')  +
  scale_color_manual(values=COLOR)
summary(glm(summary_data[,"entropy_IGHD"]~summary_data[,"NumCells"]))
dev.off()


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
ggplot(summary_data,aes(SHM_IGHM,NumCells,color=sample)) + 
  geom_point() + geom_smooth(method='lm')  +
  scale_color_manual(values=COLOR)
dev.off()


##CDR3
j<-1
p_adj<-NULL
for (i in c("CDR3_length_IGHA","CDR3_length_IGHD","CDR3_length_IGHG","CDR3_length_IGHM","CDR3_length_mutated_IGHM","CDR3_length_unmutated_IGHM",
            "CDR3_length_mutated_IGHD","CDR3_length_unmutated_IGHD")){
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
names(p_adj)<-c("CDR3_length_IGHA","CDR3_length_IGHD","CDR3_length_IGHG","CDR3_length_IGHM","CDR3_length_mutated_IGHM","CDR3_length_unmutated_IGHM",
                "CDR3_length_mutated_IGHD","CDR3_length_unmutated_IGHD")
write.csv(p_adj,"Results/Diversity/p_adj_cdr3.csv")

tiff("Results/plot_CDR3_length_naive_numcell.tiff",res=300,w=2000,h=2000)
summary(glm(summary_data[,"CDR3_length_unmutated_IGHD"]~summary_data[,"NumCells"]+summary_data[,"sample"]))
ggplot(summary_data,aes(CDR3_length_unmutated_IGHD,NumCells,color=sample)) + 
  geom_point() + geom_smooth(method='lm')  +
  scale_color_manual(values=COLOR)
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
  



