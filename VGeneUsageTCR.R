###########################################################################################
### PROJECT: Immune Repertoire. Analysis T cells antibodies for pregnancy, v_gene usage
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: 
###         
###
### Author: Brian Le
### Date: January, 2019
############################################################################################

library(RColorBrewer)
library(ggplot2)
library(pheatmap)

directory <- "TCR/"

load(paste(directory, 'TCR_data_summary.RData', sep = ""))

plot_directory <- "TCR/Results/V_gene_usage/"
dir.create(plot_directory)

data_qc$v_gene_family <- gsub("\\-.*$", "", data_qc$v_gene)

##Isolate unique clones per sample to perform clone-level analysis
data_qc_clones<-data_qc[,c("sample_label", "V_J_length_CDR3_Clone_tcrb", "v_gene", "v_gene_family")]
mask<-!duplicated(data_qc_clones[,c("V_J_length_CDR3_Clone_tcrb", "sample_label")])
data_qc_clones<-data_qc_clones[mask,]

##Calculate v_gene counts and percents from reads for each sample
##percents calculated individually for each sample (mother or fetus)
v_gene_counts<-as.data.frame.matrix(table(data_qc_clones$v_gene, data_qc_clones$sample_label))

##filter out any sample which has fewer than 100 clones... not a problem here, but can be for BCR isotype analysis
v_gene_counts <- v_gene_counts[,colSums(v_gene_counts) >= 1000]
v_gene_percents<-data.frame(apply(v_gene_counts, 2, function(x) prop.table(x)))

#barplot(as.matrix(v_gene_percents),col=COLOR,
#        names.arg = summary_data$pairs,main="",xlab = "Samples", ylab = "%",las=2)

#barplot(t(as.matrix(v_gene_percents)),col=COLOR[summary_data$sample],
#        names.arg = rownames(v_gene_percents),main="",xlab = "V-Gene", ylab = "%",las=2)

##Look at v_gene percentage distribution for just mothers and for just babies
mothers<-summary_data[summary_data$sample == "Mother", "pairs"]
fetuses<-summary_data[summary_data$sample == "Fetus", "pairs"]


COLOR<-brewer.pal(n = length(mothers), name = "Paired")

tiff(paste(plot_directory, "barplot_mothers_genes.tiff", sep = ""),res=300,w=4000,h=2000)
barplot(t(as.matrix(v_gene_percents[,mothers])),col=COLOR,
        names.arg = rownames(v_gene_percents),main="V gene distribution - mothers",xlab = "", ylab = "% of clones",las=2)
legend("topright", legend=mothers,col=COLOR,pch=15, cex=0.8)
dev.off()

tiff(paste(plot_directory, "barplot_fetuses_genes.tiff", sep = ""),res=300,w=4000,h=2000)
barplot(t(as.matrix(v_gene_percents[,fetuses])),col=COLOR,
        names.arg = rownames(v_gene_percents),main="V gene distribution - fetuses",xlab = "", ylab = "% of clones",las=2)
legend("topright", legend=fetuses,col=COLOR,pch=15, cex=0.8)
dev.off()

tiff(paste(plot_directory, "barplot_mothers_and_fetuses_genes.tiff", sep = ""),res=300,w=4000,h=2000)
par(mfrow=c(2,1))
barplot(t(as.matrix(v_gene_percents[,mothers])),col=COLOR,
        main="V gene distribution - mothers",xaxt = 'n', xlab = NULL, ylab = "% of clones",las=2)
legend("topright", legend=mothers,col=COLOR,pch=15, cex=0.6)

barplot(t(as.matrix(v_gene_percents[,fetuses])),col=COLOR,
        names.arg = rownames(v_gene_percents),main="V gene distribution - fetuses",xlab = "", ylab = "% of clones",las=2)
legend("topright", legend=fetuses,col=COLOR,pch=15, cex=0.6)
dev.off()

##Calculate v_gene counts and percents from reads for POOLED mothers/fetuses
##No down-sampling, so samples with more reads have higher contribution
##Pooling these together results in effectively plotting MEANS which means we can calculate STANDARD DEVIATIONS
##which means we can also COMPARE MEANS AND STANDARD DEVIATIONS within groups

v_gene_percents_mothers <- v_gene_percents[,mothers]
v_gene_percents_fetuses <- v_gene_percents[,fetuses]

means_mothers <- colMeans(t(v_gene_percents_mothers))
stdev_mothers <- sapply(data.frame(t(v_gene_percents_mothers)), sd, na.rm = TRUE)

means_fetuses <- colMeans(t(v_gene_percents_fetuses))
stdev_fetuses <- sapply(data.frame(t(v_gene_percents_fetuses)), sd, na.rm = TRUE)

v_gene_percents_summary <- data.frame(means_mothers = means_mothers, stdev_mothers = stdev_mothers,
                                      means_fetuses = means_fetuses, stdev_fetuses = stdev_fetuses)
v_gene_percents_summary$gene <- rownames(v_gene_percents_summary)

v_gene_percents_summary_2 <- rbind(data.frame(gene = names(means_mothers), group = "Maternal", mean = means_mothers, std = stdev_mothers),
                                   data.frame(gene = names(means_fetuses), group = "Fetal", mean = means_fetuses, std = stdev_fetuses))

##unpaired t-test comparing gene distributions between mothers and fetuses
#n1 = ncol(v_gene_percents_mothers)
#n2 = ncol(v_gene_percents_fetuses)
#means_diff = (v_gene_percents_summary$means_mothers - v_gene_percents_summary$means_fetuses)
#SE = sqrt((v_gene_percents_summary$stdev_mothers)^2/n1 + (v_gene_percents_summary$stdev_fetuses)^2/n2)
#t = means_diff / SE
#p_t = pt(t, df = n1 + n2 - 2)
#p_t_sig = p_t < 0.05

for (gene in rownames(v_gene_percents)){
  
  if ((!is.na(v_gene_percents_summary[gene, "means_fetuses"])) &
      (!is.na(v_gene_percents_summary[gene, "means_fetuses"]))){
    v_gene_percents_summary[gene, "p_wilcoxon"] <- wilcox.test(unlist(v_gene_percents_mothers[gene,]),
                                                               unlist(v_gene_percents_fetuses[gene,]),
                                                               alternative = "two.sided")$p.value
  }
  else {
    v_gene_percents_summary[gene, "p_wilcoxon"] <- NA
  }
  #print(paste(gene, v_gene_percents_summary[gene, "p_wilcoxon"], sep = "  "))
}

v_gene_percents_summary$p_wilcoxon_adj <- p.adjust(v_gene_percents_summary$p_wilcoxon)
p_w_sig = v_gene_percents_summary$p_wilcoxon_adj < 0.05

sig_genes = v_gene_percents_summary[v_gene_percents_summary$p_wilcoxon_adj < 0.05, "gene"]
v_gene_percents_summary_2_sig = v_gene_percents_summary_2[v_gene_percents_summary_2$gene %in% sig_genes,]

label = ifelse(p_w_sig, "*", "")
label = c(label, rep("", length(p_w_sig)))

v_gene_percents_summary_2$sig <- label
v_gene_percents_summary_2$sig_ypos <- v_gene_percents_summary_2$mean + v_gene_percents_summary_2$std
##if value is NA, set it to 0
v_gene_percents_summary_2$sig_ypos <- ifelse(!is.na(v_gene_percents_summary_2$sig_ypos),
                                             v_gene_percents_summary_2$sig_ypos,
                                             0)
for (gene in unique(v_gene_percents_summary_2$gene)){
  
  ypos <- max(v_gene_percents_summary_2[(v_gene_percents_summary_2$gene == gene), "sig_ypos"])
  v_gene_percents_summary_2[(v_gene_percents_summary_2$gene == gene), "sig_ypos"] <- ypos
  
}

COLOR=c("#BEAED4","#7FC97F")

tiff(paste(plot_directory, "barplot_combined_genes_errorbars.tiff", sep = ""),res=300,w=4000,h=2000)
ggplot(v_gene_percents_summary_2, aes(x = gene, y = mean, fill = group)) +
  geom_bar(position = position_dodge(), stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Pooled v gene usage - TCR") +
  xlab("V Gene") +
  ylab("% of clones") +
  coord_cartesian(ylim = c(0,1.1*max(v_gene_percents_summary_2$mean + v_gene_percents_summary_2$std)), expand = FALSE) +
  scale_fill_manual(values = COLOR, name = "") +
  geom_text(aes(x = gene, y = sig_ypos), nudge_y = 0.02*max(v_gene_percents_summary_2$sig_ypos), 
            label = label, size = 5)
dev.off()

tiff(paste(plot_directory, "barplot_combined_genes_errorbars_sig.tiff", sep = ""),res=300,w=4000,h=2000)
ggplot(v_gene_percents_summary_2_sig, aes(x = gene, y = mean, fill = group)) +
  geom_bar(position = position_dodge(), stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Significantly different genes - TCR") + 
  xlab("V Gene") +
  ylab("% of clones") +
  coord_cartesian(ylim = c(0,1.1*max(v_gene_percents_summary_2_sig$mean + v_gene_percents_summary_2_sig$std)), expand = FALSE) +
  scale_fill_manual(values = COLOR, name = "")
dev.off()

#load in ordering on 7q34 chromosome for TRB genes, order barplot by ordering
filepath = "TRB_gene_order.csv"
v_gene_order = read.csv(filepath)
colnames(v_gene_order)[1] <- "gene"

#append gene_order from imgt to dataframe
v_gene_percents_summary_2 <- merge(v_gene_percents_summary_2, v_gene_order[,c("gene", "gene_order")],
                                   all.x = TRUE, all.y = FALSE)

##plotting by ordering the v genes by chromosome location order
tiff(paste(plot_directory, "barplot_combined_genes_chromosome_order.tiff", sep = ""),res=300,w=4000,h=2000)
ggplot(v_gene_percents_summary_2, aes(x = reorder(gene, gene_order), y = mean, fill = group)) +
  geom_bar(position = position_dodge(), stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Pooled v gene usage - TCR") +
  xlab("V Gene") +
  ylab("% of clones") +
  coord_cartesian(ylim = c(0,1.1*max(v_gene_percents_summary_2$mean + v_gene_percents_summary_2$std)), expand = FALSE) +
  scale_fill_manual(values = COLOR, name = "") +
  geom_text(aes(x = reorder(gene, gene_order), y = sig_ypos, label = sig), nudge_y = 0.02*max(v_gene_percents_summary_2$sig_ypos), 
            size = 5)
dev.off()

##OLD
data_qc$sample_label_level <- substr(data_qc$sample_label, 1, 1)

v_gene_counts_level<-as.data.frame.matrix(table(data_qc$v_gene, data_qc$sample_label_level))
v_gene_percents_level<-data.frame(apply(v_gene_counts_level, 2, function(x) prop.table(x)))
colnames(v_gene_counts_level) <- c("Fetal", "Maternal")
colnames(v_gene_percents_level) <- c("Fetal", "Maternal")

COLOR=c("#BEAED4","#7FC97F")

tiff(paste(plot_directory, "barplot_combined_genes.tiff", sep = ""),res=300,w=4000,h=2000)
barplot(t(as.matrix(v_gene_percents_level)),col=COLOR,beside=TRUE,
        names.arg = rownames(v_gene_percents_level),main="V gene distribution",xlab = "", ylab = "% of clones",las=2)
legend("topright", legend=colnames(v_gene_percents_level),col=COLOR,pch=15, cex=0.6)
dev.off()


##V_gene usage for mother/fetus pairs, plotted side-by-side as barplots
COLOR=c("#BEAED4","#7FC97F")

sample_labels <- rownames(summary_data)

for (i in 1:(length(sample_labels)/2)*2){
  labels = c(sample_labels[i-1], sample_labels[i])
  pair = paste(sample_labels[i-1], sample_labels[i], sep = "-")
  filename = paste("barplot_pair_", pair,".tiff", sep = "")
  
  tiff(paste(plot_directory, filename, sep = ""),res=300,w=4000,h=2000)
  barplot(t(as.matrix(v_gene_percents[,labels])),col=COLOR,beside=TRUE,
          names.arg = rownames(v_gene_percents),main=paste("V gene distribution - ", pair),xlab = "", ylab = "% of clones",las=2)
  legend("topright", legend=labels,col=COLOR,pch=15, cex=0.8)
  dev.off()
}



my_sample_col <- data.frame(pair = gsub("[FMa]", "", colnames(v_gene_percents)),
                            sample = ifelse(substr(colnames(v_gene_percents), 1, 1) == "F", "Fetal", "Maternal"))
rownames(my_sample_col) <- colnames(v_gene_percents)

tiff(paste(plot_directory, "heatmap.tiff", sep = ""),res=300,w=2000,h=2000)
pheatmap(v_gene_percents, annotation_col = my_sample_col,
         color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
         main = "Heatmap of Genes/Samples - TCR", fontsize = 8, fontsize_row = 6)
dev.off()

##heatmaps of just significant genes
tiff(paste(plot_directory, "heatmap_sig_genes.tiff", sep = ""),res=300,w=2000,h=2000)
pheatmap(v_gene_percents[sig_genes,], annotation_col = my_sample_col,
         color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
         main = "Heatmap of Significant Genes - TCR", fontsize = 8, fontsize_row = 6)
dev.off()


#heatmap, exclude genes with consistent low usage (mean: 0.005)
tiff(paste(plot_directory, "heatmap_usage_cutoff.tiff", sep = ""),res=300,w=2000,h=2000)
p <- pheatmap(v_gene_percents[(rowSums(v_gene_percents)/ncol(v_gene_percents)) > 0.005,], annotation_col = my_sample_col,
              color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
              main = paste("Heatmap of Genes/Samples - TCR"), fontsize = 8, fontsize_row = 6)
print(p)
dev.off()

##correlatin matrix heatmap
tiff(paste(plot_directory, "heatmap_correlation.tiff", sep = ""),res=300,w=2000,h=2000)
p <- pheatmap(cor(v_gene_percents), annotation_col = my_sample_col,
              color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
              main = paste("Sample Correlation - TCR"), fontsize = 8, fontsize_row = 6)
#corrplot.mixed(cor(v_gene_percents), lower = "number", upper = "shade", order = "hclust", hclust.method = "ward.D2", addrect = 2, 
#         tl.col = "black", tl.srt = 45)
print(p)
dev.off()

