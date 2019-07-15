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

data_qc$v_gene_family <- gsub("\\-.*$", "", data_qc$v_gene)

##Isolate unique clones per isotype per sample to perform clone-level analysis
data_qc_clones<-data_qc[,c("sample_label", "isotype", "V_J_lenghCDR3_Clone_igh", "v_gene", "v_gene_family", "IGHM_mutated", "IGHD_mutated", "d_gene")]
mask<-!duplicated(data_qc_clones[,c("V_J_lenghCDR3_Clone_igh", "sample_label", "isotype")])
data_qc_clones<-data_qc_clones[which(mask),]

data_qc_clones[which(data_qc_clones$isotype == ""), "isotype"] <- "noisotype"

summary_data$row <- rownames(summary_data)

mothers<-summary_data[summary_data$sample == "Mother", "row"]
fetuses<-summary_data[summary_data$sample == "Fetus", "row"]
isotypes <- unique(data_qc_clones$isotype)

##perform for each isotype
##NOTE: will FAIL for isotype = "IGHE", with only 45 reads and not covering every sample
#isotypes <- c('IGHM', 'IGHD')

for (isotype in isotypes){
  
  print(isotype)
  
  plot_directory_iso <- paste(plot_directory, isotype, "/", sep = "")
  dir.create(plot_directory_iso, showWarnings = FALSE)
  
  data_qc_iso <- data_qc_clones[which(data_qc_clones$isotype == isotype),]
  
  ##Calculate v_gene counts and percents from reads for each sample
  ##percents calculated individually for each sample (mother or fetus)
  v_gene_counts<-as.data.frame.matrix(table(data_qc_iso$v_gene, data_qc_iso$sample_label))
  
  ##filter out any sample which has fewer than 100 clones... not a problem here, but can be for BCR isotype analysis
  cutoff <- 1000
  v_gene_counts <- v_gene_counts[,colSums(v_gene_counts) >= cutoff]
  
  v_gene_percents<-data.frame(apply(v_gene_counts, 2, function(x) prop.table(x)))
  
  v_gene_percents_mothers <- v_gene_percents[,substr(colnames(v_gene_percents), 1, 1) == "M", drop = FALSE]
  v_gene_percents_fetuses <- v_gene_percents[,substr(colnames(v_gene_percents), 1, 1) == "F", drop = FALSE]

  means_mothers <- colMeans(t(v_gene_percents_mothers))
  stdev_mothers <- sapply(data.frame(t(v_gene_percents_mothers)), sd, na.rm = TRUE)

  means_fetuses <- colMeans(t(v_gene_percents_fetuses))
  stdev_fetuses <- sapply(data.frame(t(v_gene_percents_fetuses)), sd, na.rm = TRUE)

  v_gene_percents_summary <- data.frame(means_mothers = means_mothers, stdev_mothers = stdev_mothers,
                                        means_fetuses = means_fetuses, stdev_fetuses = stdev_fetuses)
  v_gene_percents_summary$gene <- rownames(v_gene_percents_summary)

  if(length(means_mothers) == 0 & length(means_fetuses) == 0){
         v_gene_percents_summary_2 <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("gene", "group", "mean", "std"))
  } else {
         v_gene_percents_summary_2 <- rbind(data.frame(gene = names(means_mothers), group = "Maternal", mean = means_mothers, std = stdev_mothers),
                                            data.frame(gene = names(means_fetuses), group = "Fetal", mean = means_fetuses, std = stdev_fetuses))
  }
  ##unpaired t-test comparing gene distributions between mothers and fetuses
  #n1 = ncol(v_gene_percents_mothers)
  #n2 = ncol(v_gene_percents_fetuses)
  #means_diff = (v_gene_percents_summary$means_mothers - v_gene_percents_summary$means_fetuses)
  #SE = sqrt((v_gene_percents_summary$stdev_mothers)^2/n1 + (v_gene_percents_summary$stdev_fetuses)^2/n2)
  #t = means_diff / SE
  #p_t = pt(t, df = n1 + n2 - 2)
  #p_t_sig = p_t < 0.05

  ##wilcoxon UNPAIRED rank sum test comparing gene distributions between mothers and fetuses
  ##unpaired because by excluding samples with <1000 clones, paired samples aren't always included

  for (gene in rownames(v_gene_percents)){

    if ((!is.na(v_gene_percents_summary[gene, "means_fetuses"])) &
        (!is.na(v_gene_percents_summary[gene, "means_mothers"]))){
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

  tiff(paste(plot_directory_iso, "barplot_combined_genes_errorbars_", isotype, ".tiff", sep = ""),res=300,w=4000,h=2000)
  p <- ggplot(v_gene_percents_summary_2, aes(x = gene, y = mean, fill = group)) +
       geom_bar(position = position_dodge(), stat = "identity", color = "black") +
       geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
       theme(text = element_text(size = 15),
             axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
       labs(title = paste("Pooled v gene usage - ", isotype, sep = "")) +
       xlab("V Gene") +
       ylab("% of clones") +
       coord_cartesian(ylim = c(0,1.1*max(0, v_gene_percents_summary_2$mean + v_gene_percents_summary_2$std,
                                          na.rm = TRUE)), #0 provided in max as a fail case for empty dataframes
                       expand = FALSE) +
    scale_fill_manual(values = COLOR, name = "") +
       geom_text(aes(x = gene, y = sig_ypos), nudge_y = 0.02*max(0, v_gene_percents_summary_2$sig_ypos, na.rm = TRUE),
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
       coord_cartesian(ylim = c(0,1.1*max(0, v_gene_percents_summary_2$mean + v_gene_percents_summary_2$std,
                                       na.rm = TRUE)),
                    expand = FALSE) +
    scale_fill_manual(values = COLOR, name = "")
  print(p)
  dev.off()

  filepath = "IGHV_gene_order.csv"
  v_gene_order = read.csv(filepath)
  colnames(v_gene_order)[1] <- "gene"

  v_gene_percents_summary_2 <- merge(v_gene_percents_summary_2, v_gene_order[,c("gene", "gene_order")],
                                     all.x = TRUE, all.y = FALSE)

  tiff(paste(plot_directory_iso, "barplot_combined_genes_chromosome_order_", isotype, ".tiff", sep = ""),res=300,w=4000,h=2000)
  p <- ggplot(v_gene_percents_summary_2, aes(x = reorder(gene, gene_order), y = mean, fill = group)) +
    geom_bar(position = position_dodge(), stat = "identity", color = "black") +
    geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
    theme(text = element_text(size = 15),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = paste("Pooled v gene usage - ", isotype, sep = "")) +
    xlab("V Gene") +
    ylab("% of clones") +
    coord_cartesian(ylim = c(0,1.1*max(0, v_gene_percents_summary_2$mean + v_gene_percents_summary_2$std,
                                       na.rm = TRUE)), #0 provided in max as a fail case for empty dataframes
                    expand = FALSE) +
    scale_fill_manual(values = COLOR, name = "") +
    geom_text(aes(x = reorder(gene, gene_order), y = sig_ypos, label = sig), nudge_y = 0.02*max(0, v_gene_percents_summary_2$sig_ypos, na.rm = TRUE),
              size = 5)
  print(p)
  dev.off()
  
  write.csv(v_gene_percents, file = paste0(plot_directory_iso, "v_gene_percents.csv"))
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
  
  cutoff <- 1000
  v_gene_counts <- v_gene_counts[,colSums(v_gene_counts) >= cutoff]
  v_gene_percents<-data.frame(apply(v_gene_counts, 2, function(x) prop.table(x)))
  
  
  my_sample_col <- data.frame(pair = gsub("[FMa]", "", colnames(v_gene_percents)),
                              sample = ifelse(substr(colnames(v_gene_percents), 1, 1) == "F", "Fetal", "Maternal"))
  rownames(my_sample_col) <- colnames(v_gene_percents)
  
  tiff(paste(plot_directory_iso, "heatmap.tiff", sep = ""),res=300,w=2000,h=2000)
  p <- pheatmap(v_gene_percents, annotation_col = my_sample_col,
                color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
                main = paste("Heatmap of Genes/Samples - ", isotype), fontsize = 8, fontsize_row = 6)
  print(p)
  dev.off()
  
  tiff(paste(plot_directory_iso, "heatmap_usage_cutoff.tiff", sep = ""),res=300,w=2000,h=2000)
  p <- pheatmap(v_gene_percents[(rowSums(v_gene_percents)/ncol(v_gene_percents)) > 0.005,], annotation_col = my_sample_col,
                color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
                main = paste("Heatmap of Genes/Samples - ", isotype), fontsize = 8, fontsize_row = 6)
  print(p)
  dev.off()
  
  tiff(paste(plot_directory_iso, "heatmap_correlation.tiff", sep = ""),res=300,w=2000,h=2000)
  p <- pheatmap(cor(v_gene_percents), annotation_col = my_sample_col,
                color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
                main = paste("Sample Correlation - ", isotype), fontsize = 8, fontsize_row = 6)
  #corrplot.mixed(cor(v_gene_percents), lower = "number", upper = "shade", order = "hclust", hclust.method = "ward.D2", addrect = 2, 
  #         tl.col = "black", tl.srt = 45)
  print(p)
  dev.off()
}

##repeat plots/heatmaps for IGHD and IGHM, splitting into unmutated (mutate_rate <= 0.01) and mutated (mutate_rate > 0.01)

for (isotype in c("IGHD", "IGHM")){
  
  for (mutate in c("unmutated", "mutated")){
    print(paste(isotype, mutate, sep =" "))
    
    plot_directory_iso <- paste(plot_directory, isotype, "/", mutate, "/", sep = "")
    dir.create(plot_directory_iso)
    
    data_qc_iso <- data_qc_clones[data_qc_clones$isotype == isotype,]
    if (isotype == "IGHD"){
      data_qc_iso <- data_qc_iso[data_qc_iso$IGHD_mutated == mutate,]
    }
    if (isotype == "IGHM"){
      data_qc_iso <- data_qc_iso[data_qc_iso$IGHM_mutated == mutate,]
    }
    
    ##Calculate v_gene counts and percents from reads for each sample
    ##percents calculated individually for each sample (mother or fetus)
    v_gene_counts<-as.data.frame.matrix(table(data_qc_iso$v_gene, data_qc_iso$sample_label))
    
    ##filter out any sample which has fewer than 100 clones... not a problem here, but can be for BCR isotype analysis
    cutoff <- 1000
    v_gene_counts <- v_gene_counts[,colSums(v_gene_counts) >= cutoff]
    
    v_gene_percents<-data.frame(apply(v_gene_counts, 2, function(x) prop.table(x)))
    
    v_gene_percents_mothers <- v_gene_percents[,substr(colnames(v_gene_percents), 1, 1) == "M", drop = FALSE]
    v_gene_percents_fetuses <- v_gene_percents[,substr(colnames(v_gene_percents), 1, 1) == "F", drop = FALSE]
    
    means_mothers <- colMeans(t(v_gene_percents_mothers))
    stdev_mothers <- sapply(data.frame(t(v_gene_percents_mothers)), sd, na.rm = TRUE)
    
    means_fetuses <- colMeans(t(v_gene_percents_fetuses))
    stdev_fetuses <- sapply(data.frame(t(v_gene_percents_fetuses)), sd, na.rm = TRUE)
    
    v_gene_percents_summary <- data.frame(means_mothers = means_mothers, stdev_mothers = stdev_mothers,
                                          means_fetuses = means_fetuses, stdev_fetuses = stdev_fetuses)
    v_gene_percents_summary$gene <- rownames(v_gene_percents_summary)
    
    if(length(means_mothers) == 0 & length(means_fetuses) == 0){
      v_gene_percents_summary_2 <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("gene", "group", "mean", "std"))
    } else {
      v_gene_percents_summary_2 <- rbind(data.frame(gene = names(means_mothers), group = "Maternal", mean = means_mothers, std = stdev_mothers),
                                         data.frame(gene = names(means_fetuses), group = "Fetal", mean = means_fetuses, std = stdev_fetuses))
    }
    ##unpaired t-test comparing gene distributions between mothers and fetuses
    #n1 = ncol(v_gene_percents_mothers)
    #n2 = ncol(v_gene_percents_fetuses)
    #means_diff = (v_gene_percents_summary$means_mothers - v_gene_percents_summary$means_fetuses)
    #SE = sqrt((v_gene_percents_summary$stdev_mothers)^2/n1 + (v_gene_percents_summary$stdev_fetuses)^2/n2)
    #t = means_diff / SE
    #p_t = pt(t, df = n1 + n2 - 2)
    #p_t_sig = p_t < 0.05
    
    ##wilcoxon UNPAIRED rank sum test comparing gene distributions between mothers and fetuses
    ##unpaired because by excluding samples with <1000 clones, paired samples aren't always included
    
    for (gene in rownames(v_gene_percents)){
      
      if ((!is.na(v_gene_percents_summary[gene, "means_fetuses"])) &
          (!is.na(v_gene_percents_summary[gene, "means_mothers"]))){
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
    
    tiff(paste(plot_directory_iso, "barplot_combined_genes_errorbars_sig_", isotype, ".tiff", sep = ""),res=300,w=4000,h=2000)
    p <- ggplot(v_gene_percents_summary_2_sig, aes(x = gene, y = mean, fill = group)) +
      geom_bar(position = position_dodge(), stat = "identity", color = "black") +
      geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
      theme(text = element_text(size = 15),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      labs(title = paste("Significantly different genes - ", isotype, " ", mutate, sep = "")) + 
      xlab("V Gene") +
      ylab("% of clones") +
      coord_cartesian(ylim = c(0,1.1*max(0, v_gene_percents_summary_2$mean + v_gene_percents_summary_2$std,
                                         na.rm = TRUE)),
                      expand = FALSE) +
      scale_fill_manual(values = COLOR, name = "")
    print(p)
    dev.off()
    
    filepath = "IGHV_gene_order.csv"
    v_gene_order = read.csv(filepath)
    colnames(v_gene_order)[1] <- "gene"
    
    v_gene_percents_summary_2 <- merge(v_gene_percents_summary_2, v_gene_order[,c("gene", "gene_order")],
                                       all.x = TRUE, all.y = FALSE)
    
    tiff(paste(plot_directory_iso, "barplot_combined_genes_chromosome_order_", isotype, ".tiff", sep = ""),res=300,w=4000,h=2000)
    p <- ggplot(v_gene_percents_summary_2, aes(x = reorder(gene, gene_order), y = mean, fill = group)) +
      geom_bar(position = position_dodge(), stat = "identity", color = "black") +
      geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
      theme(text = element_text(size = 15),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      labs(title = paste("Pooled v gene usage - ", isotype, " ", mutate, sep = "")) +
      xlab("V Gene") +
      ylab("% of clones") +
      coord_cartesian(ylim = c(0,1.1*max(0, v_gene_percents_summary_2$mean + v_gene_percents_summary_2$std,
                                         na.rm = TRUE)), #0 provided in max as a fail case for empty dataframes
                      expand = FALSE) +
      scale_fill_manual(values = COLOR, name = "") +
      geom_text(aes(x = reorder(gene, gene_order), y = sig_ypos, label = sig), nudge_y = 0.02*max(0, v_gene_percents_summary_2$sig_ypos, na.rm = TRUE), 
                size = 5)
    print(p)
    dev.off()
    
    my_sample_col <- data.frame(pair = gsub("[FMa]", "", colnames(v_gene_percents)),
                                sample = ifelse(substr(colnames(v_gene_percents), 1, 1) == "F", "Fetal", "Maternal"))
    rownames(my_sample_col) <- colnames(v_gene_percents)
    
    tiff(paste(plot_directory_iso, "heatmap.tiff", sep = ""),res=300,w=2000,h=2000)
    p <- pheatmap(v_gene_percents, annotation_col = my_sample_col,
                  color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
                  main = paste("Heatmap of Genes/Samples - ", isotype, " ", mutate, sep = ""), fontsize = 8, fontsize_row = 6)
    print(p)
    dev.off()
    
    ##heatmaps of just significant genes
    
    if(nrow(v_gene_percents[sig_genes,]) > 1){
      tiff(paste(plot_directory_iso, "heatmap_sig_genes.tiff", sep = ""),res=300,w=2000,h=2000)
      p <- pheatmap(v_gene_percents[sig_genes,], annotation_col = my_sample_col,
               color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
               main = paste0("Heatmap of Significant Genes - ", isotype, " ", mutate), fontsize = 8, fontsize_row = 6)
      print(p)
      dev.off()
    }
    
    tiff(paste(plot_directory_iso, "heatmap_usage_cutoff.tiff", sep = ""),res=300,w=2000,h=2000)
    p <- pheatmap(v_gene_percents[(rowSums(v_gene_percents)/ncol(v_gene_percents)) > 0.005,], annotation_col = my_sample_col,
                  color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
                  main = paste("Heatmap of Genes/Samples - ", isotype, " ", mutate, sep = ""), fontsize = 8, fontsize_row = 6)
    print(p)
    dev.off()
    
    tiff(paste(plot_directory_iso, "heatmap_correlation.tiff", sep = ""),res=300,w=2000,h=2000)
    p <- pheatmap(cor(v_gene_percents), annotation_col = my_sample_col,
                  color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
                  main = paste("Sample Correlation - ", isotype, " ", mutate, sep = ""), fontsize = 8, fontsize_row = 6)
    #corrplot.mixed(cor(v_gene_percents), lower = "number", upper = "shade", order = "hclust", hclust.method = "ward.D2", addrect = 2, 
    #         tl.col = "black", tl.srt = 45)
    print(p)
    dev.off()
  }
}


##correlation matrix

library("corrplot")
pheatmap(cor(v_gene_percents), color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)))

rcorr(as.matrix(v_gene_percents))
pheatmap(rcorr(as.matrix(v_gene_percents))$r)


corrplot(cor(v_gene_percents), type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

mat <- table(data_qc_clones$isotype, data_qc_clones$sample_label)
pheatmap(t(mat),
         cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = t(mat),
         color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "Greens")))(100)))


data_qc_clones$d_gene_status <- ifelse(data_qc_clones$d_gene == "", 0, 1)
table(data_qc_clones$isotype, data_qc_clones$d_gene_status)
