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

load('~/Box Sync/BCR_TCR-Maternal-Fetal/BCR data/BCR_PTB_Term_data_summary.RData')
summary_data$row <- rownames(summary_data)

directory <- "Stanford_2/BCR/"
plot_directory <- "Stanford_2/BCR/Results/V_gene_usage/"
dir.create(plot_directory)

##Isolate unique clones per isotype per sample to perform clone-level analysis
data_qc<-data_qc[,c("sample_label", "isotype", "V_J_lenghCDR3_Clone_igh", "v_gene", "IGHM_mutated", "IGHD_mutated")]
#rm(data_qc) #frees up memory just in case

#mask<-!duplicated(data_qc_clones[,c("V_J_lenghCDR3_Clone_igh", "sample_label", "isotype")])
#data_qc_clones<-data_qc_clones[mask,]

samples <- unique(data_qc$sample_label)
data_qc_clones <- data.frame()

for (sample in samples){
  print(sample)
  data_qc_temp <- data_qc[data_qc$sample_label == sample,]
  mask<-!duplicated(data_qc_temp[,c("V_J_lenghCDR3_Clone_igh", "sample_label", "isotype")])
  data_qc_temp<-data_qc_temp[mask,]
  
  data_qc_clones <- rbind(data_qc_clones, data_qc_temp)
}

#exclude samples M/F0786, per Renan
data_qc_clones <- data_qc_clones[data_qc_clones$sample_label != "M0786" &
                                 data_qc_clones$sample_label != "F0786",]

save(data_qc_clones, file = 'Stanford_2/BCR/BCR_unique_clones.RData')

isotypes <- unique(data_qc_clones$isotype)
sample_types <- unique(summary_data$sample)
sample_outcomes <- unique(summary_data$Outcome)

isotypes <- c('IGHM', 'IGHD') #having issues wiht last 3, A/E/G have no fetal samples with >1000

for (isotype in isotypes){
  
  print(isotype)
  
  plot_directory_iso <- paste(plot_directory, isotype, "/", sep = "")
  dir.create(plot_directory_iso, showWarnings = FALSE)
  
  data_qc_iso <- data_qc_clones[data_qc_clones$isotype == isotype,]
  
  ##Calculate v_gene counts and percents from reads for each sample
  ##percents calculated individually for each sample
  v_gene_counts<-as.data.frame.matrix(table(data_qc_iso$v_gene, data_qc_iso$sample))
  
  ###so few clones that most get washed out completely
  ##filter out any sample which has fewer than 1000 distinct clones
  #1000 is what was used before
  cutoff <- 1000
  v_gene_counts <- v_gene_counts[,colSums(v_gene_counts) >= cutoff]
  
  v_gene_percents<-data.frame(apply(v_gene_counts, 2, function(x) prop.table(x)))
  colnames(v_gene_percents) <- gsub(".twin.", "(twin)", colnames(v_gene_percents))
  
  genes <- rownames(v_gene_percents)
  
  for (type in sample_types){
    samples_preterm <- summary_data[summary_data$sample == type & 
                                      summary_data$Outcome == "PTB" &
                                      summary_data$row %in% colnames(v_gene_percents), "row"]
    samples_term <-summary_data[summary_data$sample == type & 
                                  summary_data$Outcome == "Term" &
                                  summary_data$row %in% colnames(v_gene_percents), "row"]
    
    v_gene_percents_preterm <- v_gene_percents[,as.character(samples_preterm)]
    v_gene_percents_term <- v_gene_percents[,as.character(samples_term)]
    
    means_preterm <- colMeans(t(v_gene_percents_preterm))
    stdev_preterm <- apply(v_gene_percents_preterm, 1, sd, na.rm = TRUE)
    
    means_term <- colMeans(t(v_gene_percents_term))
    stdev_term <- apply(v_gene_percents_term, 1, sd, na.rm = TRUE)
    
    v_gene_percents_summary <- data.frame(means_preterm = means_preterm, stdev_preterm = stdev_preterm,
                                          means_term = means_term, stdev_term = stdev_term,
                                          row.names = genes)
    
    v_gene_percents_summary$gene <- rownames(v_gene_percents_summary)
    
    v_gene_percents_summary_2 <- rbind(data.frame(gene = genes, group = "Preterm", mean = means_preterm, std = stdev_preterm),
                                       data.frame(gene = genes, group = "Term", mean = means_term, std = stdev_term))
    
    ##wilcoxon UNPAIRED rank sum test comparing gene distributions between term/preterm
    ##unpaired because not looking at paired samples
    
    for (gene in genes){
      
      if ((!is.na(v_gene_percents_summary[gene, "means_preterm"])) &
          (!is.na(v_gene_percents_summary[gene, "means_term"]))){
        v_gene_percents_summary[gene, "p_wilcoxon"] <- wilcox.test(unlist(v_gene_percents_preterm[gene,]),
                                                                   unlist(v_gene_percents_term[gene,]),
                                                                   #paired = TRUE,
                                                                   alternative = "two.sided")$p.value
      }
      else {
        v_gene_percents_summary[gene, "p_wilcoxon"] <- NA
      }
      #print(paste(gene, v_gene_percents_summary[gene, "p_wilcoxon"], sep = "  "))
    }
    
    v_gene_percents_summary$p_wilcoxon_adj <- p.adjust(v_gene_percents_summary$p_wilcoxon)
    p_w_sig = v_gene_percents_summary$p_wilcoxon_adj < 0.05
    
    sig_genes = v_gene_percents_summary[which(v_gene_percents_summary$p_wilcoxon_adj < 0.05), "gene"]
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
    
    tiff(paste(plot_directory_iso, type, "_", isotype, ".tiff", sep = ""),res=300,w=4000,h=2000)
    p <- ggplot(v_gene_percents_summary_2, aes(x = gene, y = mean, fill = group)) +
      geom_bar(position = position_dodge(), stat = "identity", color = "black") +
      geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
      theme(text = element_text(size = 15),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      labs(title = paste0("BCR - isotype ", isotype, " - ", type, " - pooled v gene usage by outcome")) +
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
    
    if (nrow(v_gene_percents_summary_2_sig) > 0){
      tiff(paste(plot_directory_iso, type, "_", isotype, "_sig.tiff", sep = ""),res=300,w=4000,h=2000)
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
    }

    #heatmaps
    v_gene_percents_combined <- cbind(v_gene_percents_preterm, v_gene_percents_term)
    
    col_labels <- data.frame(row.names = summary_data$row,
                             outcome = summary_data$Outcome,
                             type = summary_data$sample,
                             sample = gsub("[A-Za-z\\.]", "", gsub("_[0-9]", "", summary_data$row)))
    
    ann_colors <- list(
      outcome = c(Term = "#4DAF4A", PTB = "#984EA3"),
      type = c(Maternal = "#377EB8", Fetal = "#E41A1C")
    )
    
    tiff(paste(plot_directory_iso, type, "_", isotype, "heatmap.tiff", sep = ""),res=300,w=2000,h=2000)
    p <- pheatmap(v_gene_percents_combined, annotation_col = col_labels,
                  annotation_colors = ann_colors,
                  color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
                  main = paste("Heatmap of Genes/Samples - ", isotype, " ", mutate, sep = ""), fontsize = 8, fontsize_row = 6)
    print(p)
    dev.off()
    
    #heatmap of sig genes only
    if(nrow(v_gene_percents_combined[sig_genes,]) > 1){
      tiff(paste0(plot_directory_iso, type, "_", isotype, "heatmap_sig_genes.tiff"),res=300,w=2000,h=2000)
      p <- pheatmap(v_gene_percents_combined[sig_genes,], annotation_col = col_labels,
                    annotation_colors = ann_colors,
                    color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
                    main = paste0("Heatmap of Significant Genes - ", isotype, " ", mutate), fontsize = 8, fontsize_row = 6)
      print(p)
      dev.off()
    }
  }

  
  ###fetal vs. maternal blood
  
  for (outcome in sample_outcomes){
    samples_fetal <- summary_data[summary_data$Outcome == outcome & 
                                      summary_data$sample == "Fetal" &
                                      summary_data$row %in% colnames(v_gene_percents), "row"]
    
    samples_maternal <- summary_data[summary_data$Outcome == outcome &
                                     summary_data$sample == "Maternal" &
                                     summary_data$row %in% colnames(v_gene_percents), "row"]
    
    v_gene_percents_fetal <- v_gene_percents[,as.character(samples_fetal)]
    v_gene_percents_maternal <- v_gene_percents[,as.character(samples_maternal)]
    
    means_fetal <- colMeans(t(v_gene_percents_fetal))
    stdev_fetal <- apply(v_gene_percents_fetal, 1, sd, na.rm = TRUE)
    
    means_maternal <- colMeans(t(v_gene_percents_maternal))
    stdev_maternal <- apply(v_gene_percents_maternal, 1, sd, na.rm = TRUE)
    
    v_gene_percents_summary <- data.frame(means_fetal = means_fetal, stdev_fetal = stdev_fetal,
                                          means_maternal = means_maternal, stdev_maternal = stdev_maternal,
                                          row.names = genes)
    
    v_gene_percents_summary$gene <- rownames(v_gene_percents_summary)
    
    v_gene_percents_summary_2 <- rbind(data.frame(gene = genes, group = "fetal", mean = means_fetal, std = stdev_fetal),
                                       data.frame(gene = genes, group = "maternal", mean = means_maternal, std = stdev_maternal))
    
    ##wilcoxon UNPAIRED (filtering results in uneven amounts) rank sum test comparing
    #gene distributions between fetal and m_blood for the same samples
    
    for (gene in genes){
      
      if ((!is.na(v_gene_percents_summary[gene, "means_fetal"])) &
          (!is.na(v_gene_percents_summary[gene, "means_maternal"]))){
        v_gene_percents_summary[gene, "p_wilcoxon"] <- wilcox.test(unlist(v_gene_percents_fetal[gene,]),
                                                                   unlist(v_gene_percents_maternal[gene,]),
                                                                   #paired = TRUE,
                                                                   alternative = "two.sided")$p.value
      }
      else {
        v_gene_percents_summary[gene, "p_wilcoxon"] <- NA
      }
      #print(paste(gene, v_gene_percents_summary[gene, "p_wilcoxon"], sep = "  "))
    }
    
    v_gene_percents_summary$p_wilcoxon_adj <- p.adjust(v_gene_percents_summary$p_wilcoxon)
    p_w_sig = v_gene_percents_summary$p_wilcoxon_adj < 0.05
    
    sig_genes = v_gene_percents_summary[which(v_gene_percents_summary$p_wilcoxon_adj < 0.05), "gene"]
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
    #COLOR=c("#BEAED4","#7FC97F")
    #brewer.pal(n = 5, name = "Pastel1")[1:2] 
    COLOR = c("#FBB4AE", "#B3CDE3")
    
    
    tiff(paste(plot_directory_iso, outcome, "_", isotype, ".tiff", sep = ""),res=300,w=4000,h=2000)
    p <- ggplot(v_gene_percents_summary_2, aes(x = gene, y = mean, fill = group)) +
      geom_bar(position = position_dodge(), stat = "identity", color = "black") +
      geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
      theme(text = element_text(size = 15),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      labs(title = paste0("BCR - isotype ", isotype, " - ", outcome, " - pooled v gene usage by sample type")) +
      xlab("V Gene") +
      ylab("% of clones") +
      coord_cartesian(ylim = c(0,1.1*max(v_gene_percents_summary_2$mean + v_gene_percents_summary_2$std)), expand = FALSE) +
      scale_fill_manual(values = COLOR, name = "") +
      geom_text(aes(x = gene, y = sig_ypos, label = sig), nudge_y = 0.02*max(v_gene_percents_summary_2$sig_ypos), 
                size = 5)
    print(p)
    dev.off()
    
    if (nrow(v_gene_percents_summary_2_sig) > 0){
      tiff(paste(plot_directory_iso, outcome, "_", isotype, "_sig.tiff", sep = ""),res=300,w=4000,h=2000)
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
    }
    
    if (outcome == 'Term'){
      write.csv(v_gene_percents, file = paste0(plot_directory_iso, "Term_v_gene_percents.csv"))
    }
    
    #heatmaps
    v_gene_percents_combined <- cbind(v_gene_percents_fetal, v_gene_percents_maternal)
    
    col_labels <- data.frame(row.names = summary_data$row,
                             outcome = summary_data$Outcome,
                             type = summary_data$sample,
                             sample = gsub("[A-Za-z\\.]", "", gsub("_[0-9]", "", summary_data$row)))
    
    ann_colors <- list(
      outcome = c(Term = "#4DAF4A", PTB = "#984EA3"),
      type = c(Maternal = "#377EB8", Fetal = "#E41A1C")
    )
    
    tiff(paste0(plot_directory_iso, outcome, "_", isotype, "heatmap.tiff"),res=300,w=2000,h=2000)
    p <- pheatmap(v_gene_percents_combined, annotation_col = col_labels,
                  annotation_colors = ann_colors,
                  color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
                  main = paste("Heatmap of Genes/Samples - ", isotype, " ", mutate, sep = ""), fontsize = 8, fontsize_row = 6)
    print(p)
    dev.off()
    
    #heatmap of sig genes only
    if(nrow(v_gene_percents_combined[sig_genes,]) > 1){
      tiff(paste0(plot_directory_iso, outcome, "_", isotype, "heatmap_sig_genes.tiff"),res=300,w=2000,h=2000)
      p <- pheatmap(v_gene_percents_combined[sig_genes,], annotation_col = col_labels,
                    annotation_colors = ann_colors,
                    color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
                    main = paste0("Heatmap of Significant Genes - ", isotype, " ", mutate), fontsize = 8, fontsize_row = 6)
      print(p)
      dev.off()
    }
  }
}


###split mutaated and unmutated


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
    ##percents calculated individually for each sample
    v_gene_counts<-as.data.frame.matrix(table(data_qc_iso$v_gene, data_qc_iso$sample))
    
    ###so few clones that most get washed out completely
    ##filter out any sample which has fewer than 1000 distinct clones
    #1000 is what was used before
    cutoff <- 1000
    v_gene_counts <- v_gene_counts[,colSums(v_gene_counts) >= cutoff]
    
    v_gene_percents<-data.frame(apply(v_gene_counts, 2, function(x) prop.table(x)))
    colnames(v_gene_percents) <- gsub(".twin.", "(twin)", colnames(v_gene_percents))
    
    genes <- rownames(v_gene_percents)
    
    for (type in sample_types){
      samples_preterm <- summary_data[summary_data$sample == type & 
                                        summary_data$Outcome == "PTB" &
                                        summary_data$row %in% colnames(v_gene_percents), "row"]
      samples_term <-summary_data[summary_data$sample == type & 
                                    summary_data$Outcome == "Term" &
                                    summary_data$row %in% colnames(v_gene_percents), "row"]
      
      v_gene_percents_preterm <- v_gene_percents[,as.character(samples_preterm)]
      v_gene_percents_term <- v_gene_percents[,as.character(samples_term)]
      
      means_preterm <- colMeans(t(v_gene_percents_preterm))
      stdev_preterm <- apply(v_gene_percents_preterm, 1, sd, na.rm = TRUE)
      
      means_term <- colMeans(t(v_gene_percents_term))
      stdev_term <- apply(v_gene_percents_term, 1, sd, na.rm = TRUE)
      
      v_gene_percents_summary <- data.frame(means_preterm = means_preterm, stdev_preterm = stdev_preterm,
                                            means_term = means_term, stdev_term = stdev_term,
                                            row.names = genes)
      
      v_gene_percents_summary$gene <- rownames(v_gene_percents_summary)
      
      v_gene_percents_summary_2 <- rbind(data.frame(gene = genes, group = "Preterm", mean = means_preterm, std = stdev_preterm),
                                         data.frame(gene = genes, group = "Term", mean = means_term, std = stdev_term))
      
      ##wilcoxon UNPAIRED rank sum test comparing gene distributions between term/preterm
      ##unpaired because not looking at paired samples
      
      for (gene in genes){
        
        if ((!is.na(v_gene_percents_summary[gene, "means_preterm"])) &
            (!is.na(v_gene_percents_summary[gene, "means_term"]))){
          v_gene_percents_summary[gene, "p_wilcoxon"] <- wilcox.test(unlist(v_gene_percents_preterm[gene,]),
                                                                     unlist(v_gene_percents_term[gene,]),
                                                                     #paired = TRUE,
                                                                     alternative = "two.sided")$p.value
        }
        else {
          v_gene_percents_summary[gene, "p_wilcoxon"] <- NA
        }
        #print(paste(gene, v_gene_percents_summary[gene, "p_wilcoxon"], sep = "  "))
      }
      
      v_gene_percents_summary$p_wilcoxon_adj <- p.adjust(v_gene_percents_summary$p_wilcoxon)
      p_w_sig = v_gene_percents_summary$p_wilcoxon_adj < 0.05
      
      sig_genes = v_gene_percents_summary[which(v_gene_percents_summary$p_wilcoxon_adj < 0.05), "gene"]
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
      
      tiff(paste(plot_directory_iso, type, "_", isotype, mutate, ".tiff", sep = ""),res=300,w=4000,h=2000)
      p <- ggplot(v_gene_percents_summary_2, aes(x = gene, y = mean, fill = group)) +
        geom_bar(position = position_dodge(), stat = "identity", color = "black") +
        geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
        theme(text = element_text(size = 15),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        labs(title = paste0("BCR - isotype ", isotype, " - ", mutate, " - ", type, " - pooled v gene usage by outcome")) +
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
      
      if (nrow(v_gene_percents_summary_2_sig) > 0){
        tiff(paste(plot_directory_iso, type, "_", isotype, mutate, "_sig.tiff", sep = ""),res=300,w=4000,h=2000)
        p <- ggplot(v_gene_percents_summary_2_sig, aes(x = gene, y = mean, fill = group)) +
          geom_bar(position = position_dodge(), stat = "identity", color = "black") +
          geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
          theme(text = element_text(size = 15),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
          labs(title = paste("Significantly different genes - ", isotype, " - ", mutate, sep = "")) + 
          xlab("V Gene") +
          ylab("% of clones") +
          coord_cartesian(ylim = c(0,1.1*max(0, v_gene_percents_summary_2$mean + v_gene_percents_summary_2$std,
                                             na.rm = TRUE)),
                          expand = FALSE) +
          scale_fill_manual(values = COLOR, name = "")
        print(p)
        dev.off()
      }
      
      #heatmaps
      v_gene_percents_combined <- cbind(v_gene_percents_preterm, v_gene_percents_term)
      
      col_labels <- data.frame(row.names = summary_data$row,
                               outcome = summary_data$Outcome,
                               type = summary_data$sample,
                               sample = gsub("[A-Za-z\\.]", "", gsub("_[0-9]", "", summary_data$row)))
      
      ann_colors <- list(
        outcome = c(Term = "#4DAF4A", PTB = "#984EA3"),
        type = c(Maternal = "#377EB8", Fetal = "#E41A1C")
      )
      
      tiff(paste(plot_directory_iso, type, "_", isotype, mutate, "heatmap.tiff", sep = ""),res=300,w=2000,h=2000)
      p <- pheatmap(v_gene_percents_combined, annotation_col = col_labels,
                    annotation_colors = ann_colors,
                    color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
                    main = paste("Heatmap of Genes/Samples - ", isotype, " ", mutate, sep = ""), fontsize = 8, fontsize_row = 6)
      print(p)
      dev.off()
      
      #heatmap of sig genes only
      if(nrow(v_gene_percents_combined[sig_genes,]) > 1){
        tiff(paste0(plot_directory_iso, type, "_", isotype, mutate, "heatmap_sig_genes.tiff"),res=300,w=2000,h=2000)
        p <- pheatmap(v_gene_percents_combined[sig_genes,], annotation_col = col_labels,
                      annotation_colors = ann_colors,
                      color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
                      main = paste0("Heatmap of Significant Genes - ", isotype, " ", mutate), fontsize = 8, fontsize_row = 6)
        print(p)
        dev.off()
      }
    }
    ###fetal vs. maternal blood
    
    for (outcome in sample_outcomes){
      samples_fetal <- summary_data[summary_data$Outcome == outcome & 
                                      summary_data$sample == "Fetal" &
                                      summary_data$row %in% colnames(v_gene_percents), "row"]
      
      samples_maternal <- summary_data[summary_data$Outcome == outcome &
                                         summary_data$sample == "Maternal" &
                                         summary_data$row %in% colnames(v_gene_percents), "row"]
      
      v_gene_percents_fetal <- v_gene_percents[,as.character(samples_fetal)]
      v_gene_percents_maternal <- v_gene_percents[,as.character(samples_maternal)]
      
      means_fetal <- colMeans(t(v_gene_percents_fetal))
      stdev_fetal <- apply(v_gene_percents_fetal, 1, sd, na.rm = TRUE)
      
      means_maternal <- colMeans(t(v_gene_percents_maternal))
      stdev_maternal <- apply(v_gene_percents_maternal, 1, sd, na.rm = TRUE)
      
      v_gene_percents_summary <- data.frame(means_fetal = means_fetal, stdev_fetal = stdev_fetal,
                                            means_maternal = means_maternal, stdev_maternal = stdev_maternal,
                                            row.names = genes)
      
      v_gene_percents_summary$gene <- rownames(v_gene_percents_summary)
      
      v_gene_percents_summary_2 <- rbind(data.frame(gene = genes, group = "fetal", mean = means_fetal, std = stdev_fetal),
                                         data.frame(gene = genes, group = "maternal", mean = means_maternal, std = stdev_maternal))
      
      ##wilcoxon UNPAIRED (filtering results in uneven amounts) rank sum test comparing
      #gene distributions between fetal and m_blood for the same samples
      
      for (gene in genes){
        
        if ((!is.na(v_gene_percents_summary[gene, "means_fetal"])) &
            (!is.na(v_gene_percents_summary[gene, "means_maternal"]))){
          v_gene_percents_summary[gene, "p_wilcoxon"] <- wilcox.test(unlist(v_gene_percents_fetal[gene,]),
                                                                     unlist(v_gene_percents_maternal[gene,]),
                                                                     #paired = TRUE,
                                                                     alternative = "two.sided")$p.value
        }
        else {
          v_gene_percents_summary[gene, "p_wilcoxon"] <- NA
        }
        #print(paste(gene, v_gene_percents_summary[gene, "p_wilcoxon"], sep = "  "))
      }
      
      v_gene_percents_summary$p_wilcoxon_adj <- p.adjust(v_gene_percents_summary$p_wilcoxon)
      p_w_sig = v_gene_percents_summary$p_wilcoxon_adj < 0.05
      
      sig_genes = v_gene_percents_summary[which(v_gene_percents_summary$p_wilcoxon_adj < 0.05), "gene"]
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
      #COLOR=c("#BEAED4","#7FC97F")
      #brewer.pal(n = 5, name = "Pastel1")[1:2] 
      COLOR = c("#FBB4AE", "#B3CDE3")
      
      
      tiff(paste(plot_directory_iso, outcome, "_", isotype, mutate, ".tiff", sep = ""),res=300,w=4000,h=2000)
      p <- ggplot(v_gene_percents_summary_2, aes(x = gene, y = mean, fill = group)) +
        geom_bar(position = position_dodge(), stat = "identity", color = "black") +
        geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
        theme(text = element_text(size = 15),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        labs(title = paste0("BCR - isotype ", isotype, " - ", mutate, " - ", outcome, " - pooled v gene usage by sample type")) +
        xlab("V Gene") +
        ylab("% of clones") +
        coord_cartesian(ylim = c(0,1.1*max(v_gene_percents_summary_2$mean + v_gene_percents_summary_2$std)), expand = FALSE) +
        scale_fill_manual(values = COLOR, name = "") +
        geom_text(aes(x = gene, y = sig_ypos, label = sig), nudge_y = 0.02*max(v_gene_percents_summary_2$sig_ypos), 
                  size = 5)
      print(p)
      dev.off()
      
      if (nrow(v_gene_percents_summary_2_sig) > 0){
        tiff(paste(plot_directory_iso, outcome, "_", isotype, mutate, "_sig.tiff", sep = ""),res=300,w=4000,h=2000)
        p <- ggplot(v_gene_percents_summary_2_sig, aes(x = gene, y = mean, fill = group)) +
          geom_bar(position = position_dodge(), stat = "identity", color = "black") +
          geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
          theme(text = element_text(size = 15),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
          labs(title = paste("Significantly different genes - ", isotype, " - ", mutate, sep = "")) + 
          xlab("V Gene") +
          ylab("% of clones") +
          coord_cartesian(ylim = c(0,1.1*max(0, v_gene_percents_summary_2$mean + v_gene_percents_summary_2$std,
                                             na.rm = TRUE)),
                          expand = FALSE) +
          scale_fill_manual(values = COLOR, name = "")
        print(p)
        dev.off()
      }
      
      #heatmaps
      v_gene_percents_combined <- cbind(v_gene_percents_fetal, v_gene_percents_maternal)
      
      col_labels <- data.frame(row.names = summary_data$row,
                               outcome = summary_data$Outcome,
                               type = summary_data$sample,
                               sample = gsub("[A-Za-z\\.]", "", gsub("_[0-9]", "", summary_data$row)))
      
      ann_colors <- list(
        outcome = c(Term = "#4DAF4A", PTB = "#984EA3"),
        type = c(Maternal = "#377EB8", Fetal = "#E41A1C")
        )
      
      tiff(paste0(plot_directory_iso, outcome, "_", isotype, mutate, "heatmap.tiff"),res=300,w=2000,h=2000)
      p <- pheatmap(v_gene_percents_combined, annotation_col = col_labels,
                    annotation_colors = ann_colors,
                    color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
                    main = paste("Heatmap of Genes/Samples - ", isotype, " ", mutate, sep = ""), fontsize = 8, fontsize_row = 6)
      print(p)
      dev.off()
      
      #heatmap of sig genes only
      if(nrow(v_gene_percents_combined[sig_genes,]) > 1){
        tiff(paste0(plot_directory_iso, outcome, "_", isotype, mutate, "heatmap_sig_genes.tiff"),res=300,w=2000,h=2000)
        p <- pheatmap(v_gene_percents_combined[sig_genes,], annotation_col = col_labels,
                      annotation_colors = ann_colors,
                      color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
                      main = paste0("Heatmap of Significant Genes - ", isotype, " ", mutate), fontsize = 8, fontsize_row = 6)
        print(p)
        dev.off()
      }
      
      if (outcome == 'Term'){
        write.csv(v_gene_percents, file = paste0(plot_directory_iso, mutate, "_Term_v_gene_percents.csv"))
      }
    }
  }
}
