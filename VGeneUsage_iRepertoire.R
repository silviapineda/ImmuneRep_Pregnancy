###########################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies for pregnancy
###          
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: looking at v gene usage from iRepertoire _1_usage.csv file (normalized)
###          
###
### Author: Brian Le
### Date: February, 2019
############################################################################################

library(reshape)
library(data.table)
library(ggplot2)

directory <- "iRepertoire/"

load(file=paste(directory, 'iRepertoire_data/iRepertoire_V_1_usage_data.RData', sep = ""))

summary_file_path <- paste(directory, 'iRepertoire_data/iRepertoire_summary_data.csv', sep = "")
summary_data <- read.csv(summary_file_path)

##do some clean up
colnames(data)[1:2] <- c("v_gene", "usage")
data$v_gene <- gsub("h", "", data$v_gene)

v_gene_data <- merge(data, summary_data[,c("Sample.id", "sample_name", "sample_type",
                                           "sample_population_type", "outcome")],
                     by.x = "sample", by.y = "Sample.id",
                     all.x = TRUE, all.y = FALSE
                     )

##let's do a loop, comparing term to preterm for CD4, CD8, CD19, TCR, each by m_blood and decidua

plot_directory <- "iRepertoire/Results/V_gene_usage/"
dir.create(plot_directory)

sample_pops <- c('CD4', 'CD8', 'CD19', 'TCR')
sample_types <- c('decidua', 'maternal_blood')

samples_term <- unique(summary_data[summary_data$outcome == "normal", c("sample_name", "outcome")])$sample_name
samples_preterm <- unique(summary_data[summary_data$outcome == "ptl", c("sample_name", "outcome")])$sample_name

for (pop in sample_pops){
  for (type in sample_types){
    
    if (type == 'decidua'){
      samples_preterm_selected <- samples_preterm[substr(samples_preterm, 1, 1) == "D"]
      samples_term_selected <- samples_term[substr(samples_term, 1, 1) == "D"]
    } else {
      samples_preterm_selected <- samples_preterm[substr(samples_preterm, 1, 1) == "M"]
      samples_term_selected <- samples_term[substr(samples_term, 1, 1) == "M"]
    }
    
    ##reduce data to that of sample_population and sample_type
    v_gene_data_subset <- v_gene_data[v_gene_data$sample_population_type == pop &
                                        v_gene_data$sample_type == type,]
    
    ##use cast() from reshape package to create table of percentages
    v_gene_percents <- cast(v_gene_data_subset, v_gene ~ sample_name, value = "usage")
    rownames(v_gene_percents) <- v_gene_percents$v_gene
    v_gene_percents <- v_gene_percents[,2:ncol(v_gene_percents)]
    
    genes <- rownames(v_gene_percents)
    
    v_gene_percents_preterm <- v_gene_percents[,as.character(samples_preterm_selected)]
    v_gene_percents_term <- v_gene_percents[,as.character(samples_term_selected)]
    
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
    
    ##wilcoxon UNPAIRED rank sum test comparing gene distributions between decidua and m_blood
    #paired makes no sense in this context; 10 independent samples being pooled into groups of 4 and 6
    ##unpaired to match BCR analysis where not all pairs of samples included due to 1000 clone minimum
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
    ##plotting by ordering the v genes by chromosome location order
    tiff(paste(plot_directory, "pooled_", type, "_", pop, ".tiff", sep = ""),res=300,w=4000,h=2000)
    p <- ggplot(v_gene_percents_summary_2, aes(x = gene, y = mean, fill = group)) +
      geom_bar(position = position_dodge(), stat = "identity", color = "black") +
      geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
      theme(text = element_text(size = 15),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      labs(title = paste("Pooled v gene usage - ", type, " - ", pop, sep = "")) +
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

###let's do another loop!
##comparing m_blood to decidua for CD4, CD8, CD19, TCR, each by term and preterm

plot_directory <- "iRepertoire/Results/V_gene_usage/"
dir.create(plot_directory)

sample_pops <- c('CD4', 'CD8', 'CD19', 'TCR')
sample_types <- c('decidua', 'maternal_blood')
sample_outcomes <- c('normal', 'ptl')

for (pop in sample_pops){
  for (outcome in sample_outcomes){
    
    sample_numbers <- summary_data[which(summary_data$sample_population_type == pop &
                                         summary_data$outcome == outcome), "Sample.id"]
    
    samples_decidua <- summary_data[which(summary_data$sample_population_type == pop &
                                            summary_data$outcome == outcome &
                                            summary_data$sample_type == "decidua"), "Sample.id"]
    
    samples_mblood <- summary_data[which(summary_data$sample_population_type == pop &
                                            summary_data$outcome == outcome &
                                            summary_data$sample_type == "maternal_blood"), "Sample.id"]
    
    ##reduce data to the chosen samples
    v_gene_data_subset <- v_gene_data[v_gene_data$sample %in% sample_numbers,]
    
    ##use cast() from reshape package to create table of percentages
    v_gene_percents <- cast(v_gene_data_subset, v_gene ~ sample, value = "usage")
    rownames(v_gene_percents) <- v_gene_percents$v_gene
    v_gene_percents <- v_gene_percents[,2:ncol(v_gene_percents)]
    
    genes <- rownames(v_gene_percents)
    
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
    
    ##wilcoxon PAIRED rank sum test comparing gene distributions between decidua and m_blood for the same samples

    for (gene in genes){
      
      if ((!is.na(v_gene_percents_summary[gene, "means_decidua"])) &
          (!is.na(v_gene_percents_summary[gene, "means_mblood"]))){
        v_gene_percents_summary[gene, "p_wilcoxon"] <- wilcox.test(unlist(v_gene_percents_decidua[gene,]),
                                                                   unlist(v_gene_percents_mblood[gene,]),
                                                                   paired = TRUE,
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
    tiff(paste(plot_directory, "pooled_", outcome, "_", pop, ".tiff", sep = ""),res=300,w=4000,h=2000)
    p <- ggplot(v_gene_percents_summary_2, aes(x = gene, y = mean, fill = group)) +
      geom_bar(position = position_dodge(), stat = "identity", color = "black") +
      geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
      theme(text = element_text(size = 15),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      labs(title = paste("Pooled v gene usage - ", outcome, " - ", pop, sep = "")) +
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
