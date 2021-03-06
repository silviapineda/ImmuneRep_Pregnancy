rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire. Analysis T cells antibodies for pregnancy outcomes
###
### CITATION: Modified from Silvia's SummaryPlots.R
###
### PROCESS: 
###           
### DESCRIP: Analysis of the VDJ data
###         
###
### Author: Brian Le
### Date: December, 2018
############################################################################################
library("RColorBrewer")

directory <- "TCR/"

load(paste(directory, 'TCR_data_summary.RData', sep = ""))

plot_directory <- "TCR/Results/"
dir.create(plot_directory)

####Summary plots
brewer.pal(n = 3, name = "Accent")
COLOR=c("#BEAED4","#7FC97F")

tiff(paste(plot_directory, "barplot_reads.tiff", sep = ""),res=300,w=4000,h=2000)
barplot(summary_data$read_count,col=COLOR[summary_data$sample],
        names.arg = summary_data$pairs,main="Number of Reads",xlab = "Samples", ylab = "Reads",las=2)
legend("topright", legend=levels(summary_data$sample),col=COLOR[factor(summary_data$sample)],pch=15, cex=0.8)
dev.off()

tiff(paste(plot_directory, "barplot_clones.tiff", sep = ""),res=300,w=4000,h=2000)
barplot(summary_data$clones,col=COLOR[summary_data$sample],
        names.arg = summary_data$pairs,main="clones_TCR",xlab = "Samples", ylab = "clones_TCR",las=2)
legend("topright", legend=levels(summary_data$sample),col=COLOR[factor(summary_data$sample)],pch=15, cex=0.8)
dev.off()

tiff(paste(plot_directory, "barplot_cdr3length.tiff", sep = ""),res=300,w=4000,h=2000)
barplot(summary_data$mean_CDR3_length,col=COLOR[summary_data$sample],
        names.arg = summary_data$pairs,main="mean_CDR3_length",xlab = "Samples", ylab = "mean_CDR3_length",las=2)
legend("topright", legend=levels(summary_data$sample),col=COLOR[factor(summary_data$sample)],pch=15, cex=0.8)
dev.off()

tiff(paste(plot_directory, "entropy.tiff", sep = ""),res=300,w=4000,h=2000)
barplot(summary_data$entropy,col=COLOR[summary_data$sample],
        names.arg = summary_data$pairs,main="entropy_TCR",xlab = "Samples", ylab = "entropy_TCR",las=2)
legend("topright", legend=levels(summary_data$sample),col=COLOR[factor(summary_data$sample)],pch=15, cex=0.8)
dev.off()

tiff(paste(plot_directory, "simpson.tiff", sep = ""),res=300,w=4000,h=2000)
barplot(summary_data$simpson,col=COLOR[summary_data$sample],
        names.arg = summary_data$pairs,main="simpson",xlab = "Samples", ylab = "simpson",las=2)
legend("topright", legend=levels(summary_data$sample),col=COLOR[factor(summary_data$sample)],pch=15, cex=0.8)
dev.off()

tiff(paste(plot_directory, "clonality.tiff", sep = ""),res=300,w=4000,h=2000)
barplot(summary_data$clonality,col=COLOR[summary_data$sample],
        names.arg = summary_data$pairs,main="clonality",xlab = "Samples", ylab = "clonality",las=2)
legend("topright", legend=levels(summary_data$sample),col=COLOR[factor(summary_data$sample)],pch=15, cex=0.8)
dev.off()
