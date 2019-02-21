###########################################################################################
### PROJECT: Immune Repertoire. Analysis T cells antibodies for pregnancy, j_segment usage
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: 
###         
###
### Author: Brian Le
### Date: February, 2019
############################################################################################

library(RColorBrewer)
library(ggplot2)
library(pheatmap)

directory <- "TCR/"

load(paste(directory, 'TCR_data_summary.RData', sep = ""))

plot_directory <- "TCR/Results/J_segment_usage/"
dir.create(plot_directory)

##Isolate unique clones per sample to perform clone-level analysis
data_qc_clones<-data_qc[,c("sample_label", "V_J_length_CDR3_Clone_tcrb", "j_segment")]
mask<-!duplicated(data_qc_clones[,c("V_J_length_CDR3_Clone_tcrb", "sample_label")])
data_qc_clones<-data_qc_clones[mask,]

##Calculate v_gene counts and percents from reads for each sample
##percents calculated individually for each sample (mother or fetus)

##c/p from VGeneUsageTCR.R, data_qc_clones$v_gene --> data_qc_clones$j_segment
v_gene_counts<-as.data.frame.matrix(table(data_qc_clones$j_segment, data_qc_clones$sample_label))

##filter out any sample which has fewer than 100 clones... not a problem here, but can be for BCR isotype analysis
v_gene_counts <- v_gene_counts[,colSums(v_gene_counts) >= 1000]
v_gene_percents<-data.frame(apply(v_gene_counts, 2, function(x) prop.table(x)))

##Look at v_gene percentage distribution for just mothers and for just babies
mothers<-summary_data[summary_data$sample == "Mother", "pairs"]
fetuses<-summary_data[summary_data$sample == "Fetus", "pairs"]

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

tiff(paste(plot_directory, "barplot_jSegment_usage_TCR.tiff", sep = ""),res=300,w=4000,h=2000)
ggplot(v_gene_percents_summary_2, aes(x = gene, y = mean, fill = group)) +
  geom_bar(position = position_dodge(), stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Pooled j segment usage - TCR") +
  xlab("J Segment") +
  ylab("% of clones") +
  coord_cartesian(ylim = c(0,1.1*max(v_gene_percents_summary_2$mean + v_gene_percents_summary_2$std)), expand = FALSE) +
  scale_fill_manual(values = COLOR, name = "") +
  geom_text(aes(x = gene, y = sig_ypos), nudge_y = 0.02*max(v_gene_percents_summary_2$sig_ypos), 
            label = label, size = 5)
dev.off()

tiff(paste(plot_directory, "barplot_jSegment_usage_sig_TCR.tiff", sep = ""),res=300,w=4000,h=2000)
ggplot(v_gene_percents_summary_2_sig, aes(x = gene, y = mean, fill = group)) +
  geom_bar(position = position_dodge(), stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = "Significantly different genes - TCR") + 
  xlab("J Segment") +
  ylab("% of clones") +
  coord_cartesian(ylim = c(0,1.1*max(v_gene_percents_summary_2_sig$mean + v_gene_percents_summary_2_sig$std)), expand = FALSE) +
  scale_fill_manual(values = COLOR, name = "")
dev.off()

my_sample_col <- data.frame(pair = gsub("[FMa]", "", colnames(v_gene_percents)),
                            sample = ifelse(substr(colnames(v_gene_percents), 1, 1) == "F", "Fetal", "Maternal"))
rownames(my_sample_col) <- colnames(v_gene_percents)

tiff(paste(plot_directory, "heatmap_jSegment.tiff", sep = ""),res=300,w=2000,h=2000)
pheatmap(v_gene_percents, annotation_col = my_sample_col,
         color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
         main = "J Segment Usage Heatmap - TCR", fontsize = 8, fontsize_row = 6)
dev.off()

##heatmaps of just significant genes
tiff(paste(plot_directory, "heatmap_jSegment_sig.tiff", sep = ""),res=300,w=2000,h=2000)
pheatmap(v_gene_percents[sig_genes,], annotation_col = my_sample_col,
         color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
         main = "J Segment Usage Heatmap (Significant) - TCR", fontsize = 8, fontsize_row = 6)
dev.off()


#heatmap, exclude genes with consistent low usage (mean: 0.005)
tiff(paste(plot_directory, "heatmap_jSegment_usage_cutoff.tiff", sep = ""),res=300,w=2000,h=2000)
p <- pheatmap(v_gene_percents[(rowSums(v_gene_percents)/ncol(v_gene_percents)) > 0.005,], annotation_col = my_sample_col,
              color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
              main = paste("J Segment Usage Heatmap - TCR"), fontsize = 8, fontsize_row = 6)
print(p)
dev.off()

##correlatin matrix heatmap
tiff(paste(plot_directory, "heatmap_jSegment_correlation.tiff", sep = ""),res=300,w=2000,h=2000)
p <- pheatmap(cor(v_gene_percents), annotation_col = my_sample_col,
              color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
              main = paste("J Segment Usage Sample Correlation - TCR"), fontsize = 8, fontsize_row = 6)
#corrplot.mixed(cor(v_gene_percents), lower = "number", upper = "shade", order = "hclust", hclust.method = "ward.D2", addrect = 2, 
#         tl.col = "black", tl.srt = 45)
print(p)
dev.off()
