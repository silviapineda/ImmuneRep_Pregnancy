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

####Summary plots
brewer.pal(n = 3, name = "Accent")
COLOR=c("#BEAED4","#7FC97F")

tiff("Results/barplot_reads.tiff",res=300,w=4000,h=2000)
barplot(summary_data$read_count,col=COLOR[summary_data$sample],
        names.arg = summary_data$pairs,main="Number of Reads",xlab = "Samples", ylab = "Reads",las=2)
legend("topright", legend=levels(summary_data$sample),col=COLOR[factor(summary_data$sample)],pch=15, cex=0.8)
dev.off()

tiff("Results/barplot_clones.tiff",res=300,w=4000,h=2000)
barplot(summary_data$clones,col=COLOR[summary_data$sample],
        names.arg = summary_data$pairs,main="Number of Clones",xlab = "Samples", ylab = "Clones",las=2)
legend("topright", legend=levels(summary_data$sample),col=COLOR[factor(summary_data$sample)],pch=15, cex=0.8)
dev.off()

tiff("Results/barplot_SHM.tiff",res=300,w=4000,h=2000)
barplot(summary_data$SHM,col=COLOR[summary_data$sample],
        names.arg = summary_data$pairs,main="SHM",xlab = "Samples", ylab = "SHM",las=2)
legend("topright", legend=levels(summary_data$sample),col=COLOR[factor(summary_data$sample)],pch=15, cex=0.8)
dev.off()

tiff("Results/barplot_cdr3length.tiff",res=300,w=4000,h=2000)
barplot(summary_data$mean_CDR3_length,col=COLOR[summary_data$sample],
        names.arg = summary_data$pairs,main="mean_CDR3_length",xlab = "Samples", ylab = "mean_CDR3_length",las=2)
legend("topright", legend=levels(summary_data$sample),col=COLOR[factor(summary_data$sample)],pch=15, cex=0.8)
dev.off()

tiff("Results/entropy.tiff",res=300,w=4000,h=2000)
barplot(summary_data$entropy,col=COLOR[summary_data$sample],
        names.arg = summary_data$pairs,main="entropy",xlab = "Samples", ylab = "entropy",las=2)
legend("topright", legend=levels(summary_data$sample),col=COLOR[factor(summary_data$sample)],pch=15, cex=0.8)
dev.off()

tiff("Results/simpson.tiff",res=300,w=4000,h=2000)
barplot(summary_data$simpson,col=COLOR[summary_data$sample],
        names.arg = summary_data$pairs,main="simpson",xlab = "Samples", ylab = "simpson",las=2)
legend("topright", legend=levels(summary_data$sample),col=COLOR[factor(summary_data$sample)],pch=15, cex=0.8)
dev.off()

tiff("Results/clonality.tiff",res=300,w=4000,h=2000)
barplot(summary_data$clonality,col=COLOR[summary_data$sample],
        names.arg = summary_data$pairs,main="clonality",xlab = "Samples", ylab = "clonality",las=2)
legend("topright", legend=levels(summary_data$sample),col=COLOR[factor(summary_data$sample)],pch=15, cex=0.8)
dev.off()
