###########################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies for pregnancy, v_gene usage
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: after doing clone analysis, look at v gene usage among unique clones
###         
###
### Author: Brian Le
### Date: April, 2019
############################################################################################

library(RColorBrewer)
library(ggplot2)
library(pheatmap)

directory <- "iRepertoire/"

summary_file_path <- paste(directory, 'iRepertoire_data/iRepertoire_summary_data.csv', sep = "")
summary_data <- read.csv(summary_file_path)

#input filename of file that processed clones without downsampling
data_qc <- read.csv("iRepertoire_clones_BCR_processed.csv")

plot_directory <- "iRepertoire/Results/V_gene_usage_clones/"
dir.create(plot_directory)

##Extract *_genes from *_segments
data_qc$v_gene <- gsub("\\*.*$", "", data_qc$V)
data_qc$j_gene <- gsub("\\*.*$", "", data_qc$J)

get_isotype <- function(class){
  
  class <- gsub("hTRB", "", class)
  class <- gsub("hIGH", "", class)
  class <- gsub("\\*.*", "", class)
  
  return(class)
}

#this leaves G1, G3, and G4 as separate isotype subclasses
data_qc$isotype <- with(data_qc, get_isotype(C))
data_qc$V_J_lengthCDR3_Clone <- paste(data_qc$V_J_lengthCDR3, data_qc$numberClone)

#data_qc$v_gene_family <- gsub("\\-.*$", "", data_qc$v_gene)

##Isolate unique clones per isotype per sample to perform clone-level analysis
data_qc_clones<-data_qc[,c("sample", "isotype", "V_J_lengthCDR3_Clone", "v_gene")]
mask<-!duplicated(data_qc_clones[,c("V_J_lengthCDR3_Clone", "sample", "isotype")])
data_qc_clones<-data_qc_clones[mask,]

#drop the few missing isotype clones
data_qc_clones<-data_qc_clones[data_qc_clones$isotype != "-",]

isotypes <- unique(data_qc_clones$isotype)

##perform for each isotype

for (isotype in isotypes){
  
  print(isotype)
  
  plot_directory_iso <- paste(plot_directory, isotype, "/", sep = "")
  dir.create(plot_directory_iso)
  
  data_qc_iso <- data_qc_clones[data_qc_clones$isotype == isotype,]
  
  ##Calculate v_gene counts and percents from reads for each sample
  ##percents calculated individually for each sample
  v_gene_counts<-as.data.frame.matrix(table(data_qc_iso$v_gene, data_qc_iso$sample))
  
  ###filter out any sample which has fewer than 50 clones
  cutoff <- 50
  v_gene_counts <- v_gene_counts[,colSums(v_gene_counts) >= cutoff]
  
  v_gene_percents<-data.frame(apply(v_gene_counts, 2, function(x) prop.table(x)))
  colnames(v_gene_percents) <- gsub("X", "", colnames(v_gene_percents))
  
  genes <- rownames(v_gene_percents)
  samples_preterm <- summary_data[summary_data$outcome == "ptl" &
                                    summary_data$Sample.id %in% colnames(v_gene_percents), "Sample.id"]
  samples_term <-summary_data[summary_data$outcome == "normal" &
                                summary_data$Sample.id %in% colnames(v_gene_percents), "Sample.id"]
  
  v_gene_percents_preterm <- v_gene_percents[,as.character(samples_preterm)]
  v_gene_percents_term <- v_gene_percents[,as.character(samples_term)]
  
  means_preterm <- colMeans(t(v_gene_percents_preterm))
  stdev_preterm <- sapply(transpose(v_gene_percents_preterm), sd, na.rm = TRUE)
  
  means_term <- colMeans(t(v_gene_percents_term))
  stdev_term <- sapply(transpose(v_gene_percents_term), sd, na.rm = TRUE)
  
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
  
  tiff(paste(plot_directory_iso, "outcome_barplot_genes_", isotype, ".tiff", sep = ""),res=300,w=4000,h=2000)
  p <- ggplot(v_gene_percents_summary_2, aes(x = gene, y = mean, fill = group)) +
       geom_bar(position = position_dodge(), stat = "identity", color = "black") +
       geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
       theme(text = element_text(size = 15),
             axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
       labs(title = paste("Pooled v gene usage by outcome - isotype ", isotype, sep = "")) +
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
  
  tiff(paste(plot_directory_iso, "outcome_barplot_genes_sig_", isotype, ".tiff", sep = ""),res=300,w=4000,h=2000)
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
  
  ###decidua vs. maternal blood
  
  sample_outcomes <- unique(summary_data$outcome)
  
  for (outcome in sample_outcomes){
    samples_decidua <- summary_data[summary_data$outcome == outcome & 
                                      summary_data$sample_type == "decidua" &
                                      summary_data$Sample.id %in% colnames(v_gene_percents), "Sample.id"]
    
    samples_mblood <- summary_data[summary_data$outcome == outcome &
                                     summary_data$sample_type == "maternal_blood" &
                                     summary_data$Sample.id %in% colnames(v_gene_percents), "Sample.id"]
    
    v_gene_percents_decidua <- v_gene_percents[,as.character(samples_decidua)]
    v_gene_percents_mblood <- v_gene_percents[,as.character(samples_mblood)]
    
    means_decidua <- colMeans(t(v_gene_percents_decidua))
    stdev_decidua <- sapply(transpose(v_gene_percents_decidua), sd, na.rm = TRUE)
    
    means_mblood <- colMeans(t(v_gene_percents_mblood))
    stdev_mblood <- sapply(transpose(v_gene_percents_mblood), sd, na.rm = TRUE)
    
    v_gene_percents_summary <- data.frame(means_decidua = means_decidua, stdev_decidua = stdev_decidua,
                                          means_mblood = means_mblood, stdev_mblood = stdev_mblood,
                                          row.names = genes)
    
    v_gene_percents_summary$gene <- rownames(v_gene_percents_summary)
    
    v_gene_percents_summary_2 <- rbind(data.frame(gene = genes, group = "decidua", mean = means_decidua, std = stdev_decidua),
                                       data.frame(gene = genes, group = "maternal_blood", mean = means_mblood, std = stdev_mblood))
    
    ##wilcoxon UNPAIRED (filtering results in uneven amounts) rank sum test comparing gene distributions between decidua and m_blood for the same samples
    
    for (gene in genes){
      
      if ((!is.na(v_gene_percents_summary[gene, "means_decidua"])) &
          (!is.na(v_gene_percents_summary[gene, "means_mblood"]))){
        v_gene_percents_summary[gene, "p_wilcoxon"] <- wilcox.test(unlist(v_gene_percents_decidua[gene,]),
                                                                   unlist(v_gene_percents_mblood[gene,]),
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
    ##plotting by ordering the v genes by chromosome location order
    tiff(paste(plot_directory_iso, "sample-type_barplot_genes_", isotype, ".tiff", sep = ""),res=300,w=4000,h=2000)
    p <- ggplot(v_gene_percents_summary_2, aes(x = gene, y = mean, fill = group)) +
      geom_bar(position = position_dodge(), stat = "identity", color = "black") +
      geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
      theme(text = element_text(size = 15),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      labs(title = paste("Pooled v gene usage by sample type - isotype ", isotype, sep = "")) +
      xlab("V Gene") +
      ylab("% of clones") +
      coord_cartesian(ylim = c(0,1.1*max(v_gene_percents_summary_2$mean + v_gene_percents_summary_2$std)), expand = FALSE) +
      scale_fill_manual(values = COLOR, name = "") +
      geom_text(aes(x = gene, y = sig_ypos, label = sig), nudge_y = 0.02*max(v_gene_percents_summary_2$sig_ypos), 
                size = 5)
    print(p)
    dev.off()
  }
}