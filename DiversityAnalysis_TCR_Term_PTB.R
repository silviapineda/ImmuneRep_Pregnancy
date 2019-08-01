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
### DESCRIP: Diversity Analysis, TCR
###         
###
### Author: Silvia Pineda / adapted by Brian Le
### Date: July, 2019
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

load("Stanford_2/TCR/TCR_data_summary.RData")
directory <- "Stanford_2/TCR/Results/"

###Function to plot summary plots with bar plots
####Summary plots
cols=brewer.pal(3,name = "Set2")

###Delete the twin for analysis
summary_data<-summary_data[-c(39:41),]
summary_data<-summary_data[order(rownames(summary_data)),]

########################################
######### TCR diversity analysis ######
######################################
#load("Data/TCR/TCR_data_summary.RData")
#info_samples<-read.csv("Data/NormalSamplesScottBoydStanford9_13_18.csv")
#id<-match(info_samples$Sample_ID,rownames(summary_data))
#summary_data$pairs[id]<-info_samples$Pairs
#summary_data$NumCells[id]<-info_samples$NumCells
#summary_data$sample<-gsub("\\*", "", substr(rownames(summary_data), 1, 1))
#summary_data$sample<-ifelse(summary_data$sample=="M","Maternal","Fetal")
#summary_data$sample<-factor(summary_data$sample)
#summary_data$sample<-relevel(summary_data$sample, ref = "Maternal")


###Function to plot summary plots with bar plots
####Summary plots
brewer.pal(n = 3, name = "Accent")
COLOR=c("#BEAED4","#7FC97F")

summary_data<-summary_data[order(rownames(summary_data)),]
tiff(paste0(directory, "boxplot_clones_TCR.tiff"),res=300,w=1500,h=2000)
g<-ggpaired(summary_data, x = "sample", y = "clones",
            color = "sample", line.color = "gray", line.size = 0.4,
            palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
  labs(title="TCR",x="Sample", y = "Number of clones") + theme(legend.position="none") +
  stat_compare_means(method = "wilcox.test",paired = TRUE)
print(g)
dev.off()

summary(lmer(summary_data[,"clones"] ~ summary_data$sample + summary_data$NumCells + (1| summary_data$pairs)))

  
tiff(paste0(directory,"boxplot_entropy_TCR.tiff"),res=300,w=1500,h=2000)
g<-ggpaired(summary_data, x = "sample", y = "entropy",
            color = "sample", line.color = "gray", line.size = 0.4,
            palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
  labs(title="TCR",x="Sample", y = "Entropy") + theme(legend.position="none") +
  stat_compare_means(method = "wilcox.test",paired = TRUE)
print(g)
dev.off()

summary(lmer(summary_data[,"entropy"] ~ summary_data$sample + summary_data$NumCells + (1| summary_data$pairs)))
  
  
tiff(paste0(directory,"boxplot_clonality_TCR.tiff"),res=300,w=1500,h=2000)
g<-ggpaired(summary_data, x = "sample", y = "clonality",
            color = "sample", line.color = "gray", line.size = 0.4,
            palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
  labs(title="TCR",x="Sample", y = "clonality") + theme(legend.position="none") +
  stat_compare_means(method = "wilcox.test",paired = TRUE)
print(g)
dev.off()

  
tiff(paste0(directory,"boxplot_mean_CDR3_length_TCR.tiff"),res=300,w=1500,h=2000)
g<-ggpaired(summary_data, x = "sample", y = "mean_CDR3_length",
            color = "sample", line.color = "gray", line.size = 0.4,
            palette = c("#BEAED4","#7FC97F"))+theme(text = element_text(size=15)) +
  labs(title="TCR",x="Sample", y = "mean_CDR3_length") + theme(legend.position="none") +
  stat_compare_means(method = "wilcox.test",paired = TRUE)
print(g)
dev.off()

summary(lmer(summary_data[,"mean_CDR3_length"] ~ summary_data$sample + summary_data$NumCells + (1| summary_data$pairs))) #0.057




