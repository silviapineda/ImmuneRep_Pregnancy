###########################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies for pregnancy, v_gene usage
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

directory <- "BCR/"

load(paste(directory, 'BCR_data_summary.RData', sep = ""))

plot_directory <- "BCR/Results/V_gene_usage/"
dir.create(plot_directory)

##Extract *_genes from *_segments
data_qc$v_gene <- gsub("\\*.*$", "", data_qc$v_segment)
data_qc$j_gene <- gsub("\\*.*$", "", data_qc$j_segment)
data_qc$d_gene <- gsub("\\*.*$", "", data_qc$d_segment)

##Isolate unique clones per sample to perform clone-level analysis
data_qc_clones<-data_qc[,c("sample_label", "isotype", "V_J_lenghCDR3_Clone_igh", "v_gene")]
mask<-!duplicated(paste(data_qc_clones$V_J_lenghCDR3_Clone_igh, data_qc_clones$sample_label, sep = ''))
data_qc_clones<-data_qc_clones[mask,]

data_qc_clones[data_qc_clones$isotype == "", "isotype"] <- "noisotype"

summary_data$row <- rownames(summary_data)

mothers<-summary_data[summary_data$sample == "Mother", "row"]
fetuses<-summary_data[summary_data$sample == "Fetus", "row"]
isotypes <- unique(data_qc_clones$isotype)

##perform for each isotype
##NOTE: will FAIL for isotype = "IGHE", with only 45 reads and not covering every sample

for (isotype in isotypes){
  
  print(isotype)
  
  plot_directory_iso <- paste(plot_directory, isotype, "/", sep = "")
  dir.create(plot_directory_iso)
  
  data_qc_iso <- data_qc_clones[data_qc_clones$isotype == isotype,]
  
  ##Calculate v_gene counts and percents from reads for each sample
  ##percents calculated individually for each sample (mother or fetus)
  v_gene_counts<-as.data.frame.matrix(table(data_qc_iso$v_gene, data_qc_iso$sample_label))
  v_gene_percents<-data.frame(apply(v_gene_counts, 2, function(x) prop.table(x)))
  
  v_gene_percents_mothers <- v_gene_percents[,mothers]
  v_gene_percents_fetuses <- v_gene_percents[,fetuses]
  
  means_mothers <- colMeans(t(v_gene_percents_mothers))
  stdev_mothers <- sapply(data.frame(t(v_gene_percents_mothers)), sd, na.rm = TRUE)
  
  means_fetuses <- colMeans(t(v_gene_percents_fetuses))
  stdev_fetuses <- sapply(data.frame(t(v_gene_percents_fetuses)), sd, na.rm = TRUE)
  
  v_gene_percents_summary <- data.frame(means_mothers = means_mothers, stdev_mothers = stdev_mothers,
                                        means_fetuses = means_fetuses, stdev_fetuses = stdev_fetuses)
  v_gene_percents_summary$gene <- rownames(v_gene_percents_summary)
  
  v_gene_percents_summary_2 <- rbind(data.frame(gene = names(means_mothers), group = "Mothers", mean = means_mothers, std = stdev_mothers),
                                     data.frame(gene = names(means_fetuses), group = "Fetuses", mean = means_fetuses, std = stdev_fetuses))
  
  n1 = 8
  n2 = 8
  means_diff = (v_gene_percents_summary$means_mothers - v_gene_percents_summary$means_fetuses)
  SE = sqrt((v_gene_percents_summary$stdev_mothers)^2/n1 + (v_gene_percents_summary$stdev_fetuses)^2/n2)
  t = means_diff / SE
  p_t = pt(t, df = n1 + n2 - 2)
  p_t_sig = p_t < 0.05
  
  sig_genes = v_gene_percents_summary[p_t_sig, "gene"]
  v_gene_percents_summary_2_sig = v_gene_percents_summary_2[v_gene_percents_summary_2$gene %in% sig_genes,]
  
  label = ifelse(p_t_sig, "*", "")
  label = c(label, rep("", length(p_t_sig)))
  
  v_gene_percents_summary_2$sig <- label
  v_gene_percents_summary_2$sig_ypos <- v_gene_percents_summary_2$mean + v_gene_percents_summary_2$std
  for (gene in unique(v_gene_percents_summary_2$gene)){

    ypos <- max(v_gene_percents_summary_2[(v_gene_percents_summary_2$gene == gene), "sig_ypos"])
    v_gene_percents_summary_2[(v_gene_percents_summary_2$gene == gene), "sig_ypos"] <- ypos
    
    }
  
  tiff(paste(plot_directory_iso, "barplot_combined_genes_errorbars_", isotype, ".tiff", sep = ""),res=300,w=4000,h=2000)
  p <- ggplot(v_gene_percents_summary_2, aes(x = gene, y = mean, fill = group)) +
       geom_bar(position = position_dodge(), stat = "identity", color = "black") +
       geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
       theme(text = element_text(size = 15),
             axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
       labs(title = paste("Pooled v gene usage - ", isotype, sep = "")) +
       xlab("V Gene") +
       ylab("% of clones") +
       coord_cartesian(ylim = c(0,1.1*max(v_gene_percents_summary_2$mean + v_gene_percents_summary_2$std)), expand = FALSE) +
       scale_fill_discrete(name = "") +
       geom_text(aes(x = gene, y = sig_ypos), nudge_y = 0.02*max(v_gene_percents_summary_2$sig_ypos), 
                 label = label, size = 5)
  print(p)
  dev.off()
  
  tiff(paste(plot_directory_iso, "barplot_combined_genes_errorbars_sig_", isotype, ".tiff", sep = ""),res=300,w=4000,h=2000)
  p <- ggplot(v_gene_percents_summary_2_sig, aes(x = gene, y = mean, fill = group)) +
       geom_bar(position = position_dodge(), stat = "identity", color = "black") +
       geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
       theme(text = element_text(size = 15),
             axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
       labs(title = paste("Significantly different genes - ", isotype, sep = "")) + 
       xlab("V Gene") +
       ylab("% of clones") +
       coord_cartesian(ylim = c(0,1.1*max(v_gene_percents_summary_2$mean + v_gene_percents_summary_2$std)), expand = FALSE) +
       scale_fill_discrete(name = "")
  print(p)
  dev.off()
}

##create heatmaps for each isotype
##could be merged with above for loop

for (isotype in isotypes){

  print(isotype)
  
  plot_directory_iso <- paste(plot_directory, isotype, "/", sep = "")
  
  data_qc_iso <- data_qc_clones[data_qc_clones$isotype == isotype,]
  
  ##Calculate v_gene counts and percents from reads for each sample
  ##percents calculated individually for each sample (mother or fetus)
  v_gene_counts<-as.data.frame.matrix(table(data_qc_iso$v_gene, data_qc_iso$sample_label))
  v_gene_percents<-data.frame(apply(v_gene_counts, 2, function(x) prop.table(x)))
  
  my_sample_col <- data.frame(pair = rep(gsub("F", "", colnames(v_gene_percents)[1:10*2]),
                                         rep(2, c(10))))
  rownames(my_sample_col) <- colnames(v_gene_percents)
  
  tiff(paste(plot_directory_iso, "heatmap.tiff", sep = ""),res=300,w=2000,h=2000)
  p <- pheatmap(v_gene_percents, annotation_col = my_sample_col,
                color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
                main = paste("Heatmap of Genes/Samples - ", isotype), fontsize = 8, fontsize_row = 6)
  print(p)
  dev.off()
}
