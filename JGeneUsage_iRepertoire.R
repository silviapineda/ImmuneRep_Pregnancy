###########################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies for pregnancy
###          
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: looking at j gene usage from iRepertoire _1_usage.csv file (normalized)
###          
###
### Author: Brian Le
### Date: February, 2019
############################################################################################

library(reshape)
library(data.table)
library(ggplot2)

directory <- "iRepertoire/"

load(file=paste(directory, 'iRepertoire_J_1_usage_data.RData', sep = ""))

summary_file_path <- paste(directory, 'iRepertoire_summary_data.csv', sep = "")
summary_data <- read.csv(summary_file_path)

##do some clean up
colnames(data)[1:2] <- c("j_gene", "usage")
data$j_gene <- gsub("h", "", data$j_gene)

j_gene_data <- merge(data, summary_data[,c("Sample.id", "sample_name", "sample_type",
                                           "sample_population_type", "outcome")],
                     by.x = "sample", by.y = "Sample.id",
                     all.x = TRUE, all.y = FALSE
)

##let's do a loop, comparing term to preterm for CD4, CD8, CD19, TCR, each by m_blood and decidua

plot_directory <- "iRepertoire/Results/j_gene_usage/"
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
    j_gene_data_subset <- j_gene_data[j_gene_data$sample_population_type == pop &
                                        j_gene_data$sample_type == type,]
    
    ##use cast() from reshape package to create table of percentages
    j_gene_percents <- cast(j_gene_data_subset, j_gene ~ sample_name, value = "usage")
    rownames(j_gene_percents) <- j_gene_percents$j_gene
    j_gene_percents <- j_gene_percents[,2:ncol(j_gene_percents)]
    
    genes <- rownames(j_gene_percents)
    
    j_gene_percents_preterm <- j_gene_percents[,as.character(samples_preterm_selected)]
    j_gene_percents_term <- j_gene_percents[,as.character(samples_term_selected)]
    
    means_preterm <- colMeans(t(j_gene_percents_preterm))
    stdej_preterm <- sapply(transpose(j_gene_percents_preterm), sd, na.rm = TRUE)
    
    means_term <- colMeans(t(j_gene_percents_term))
    stdej_term <- sapply(transpose(j_gene_percents_term), sd, na.rm = TRUE)
    
    j_gene_percents_summary <- data.frame(means_preterm = means_preterm, stdej_preterm = stdej_preterm,
                                          means_term = means_term, stdej_term = stdej_term,
                                          row.names = genes)
    
    j_gene_percents_summary$gene <- rownames(j_gene_percents_summary)
    
    j_gene_percents_summary_2 <- rbind(data.frame(gene = genes, group = "Preterm", mean = means_preterm, std = stdej_preterm),
                                       data.frame(gene = genes, group = "Term", mean = means_term, std = stdej_term))
    
    ##wilcoxon UNPAIRED rank sum test comparing gene distributions between decidua and m_blood
    #paired makes no sense in this context; 10 independent samples being pooled into groups of 4 and 6
    ##unpaired to match BCR analysis where not all pairs of samples included due to 1000 clone minimum
    for (gene in genes){
      
      if ((!is.na(j_gene_percents_summary[gene, "means_preterm"])) &
          (!is.na(j_gene_percents_summary[gene, "means_term"]))){
        j_gene_percents_summary[gene, "p_wilcoxon"] <- wilcox.test(unlist(j_gene_percents_preterm[gene,]),
                                                                   unlist(j_gene_percents_term[gene,]),
                                                                   #paired = TRUE,
                                                                   alternative = "two.sided")$p.value
      }
      else {
        j_gene_percents_summary[gene, "p_wilcoxon"] <- NA
      }
      #print(paste(gene, j_gene_percents_summary[gene, "p_wilcoxon"], sep = "  "))
    }
    
    j_gene_percents_summary$p_wilcoxon_adj <- p.adjust(j_gene_percents_summary$p_wilcoxon)
    p_w_sig = j_gene_percents_summary$p_wilcoxon_adj < 0.05
    
    sig_genes = j_gene_percents_summary[j_gene_percents_summary$p_wilcoxon_adj < 0.05, "gene"]
    j_gene_percents_summary_2_sig = j_gene_percents_summary_2[j_gene_percents_summary_2$gene %in% sig_genes,]
    
    label = ifelse(p_w_sig, "*", "")
    label = c(label, rep("", length(p_w_sig)))
    
    j_gene_percents_summary_2$sig <- label
    j_gene_percents_summary_2$sig_ypos <- j_gene_percents_summary_2$mean + j_gene_percents_summary_2$std
    ##if value is NA, set it to 0
    j_gene_percents_summary_2$sig_ypos <- ifelse(!is.na(j_gene_percents_summary_2$sig_ypos),
                                                 j_gene_percents_summary_2$sig_ypos,
                                                 0)
    for (gene in unique(j_gene_percents_summary_2$gene)){
      
      ypos <- max(j_gene_percents_summary_2[(j_gene_percents_summary_2$gene == gene), "sig_ypos"])
      j_gene_percents_summary_2[(j_gene_percents_summary_2$gene == gene), "sig_ypos"] <- ypos
      
    }
    
    COLOR=c("#BEAED4","#7FC97F")
    
    
    ##plotting by ordering the v genes by chromosome location order
    tiff(paste(plot_directory, "pooled_", type, "_", pop, ".tiff", sep = ""),res=300,w=4000,h=2000)
    p <- ggplot(j_gene_percents_summary_2, aes(x = gene, y = mean, fill = group)) +
      geom_bar(position = position_dodge(), stat = "identity", color = "black") +
      geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width= 0.2, position = position_dodge(0.9)) +
      theme(text = element_text(size = 15),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      labs(title = paste("Pooled v gene usage - ", type, " - ", pop, sep = "")) +
      xlab("V Gene") +
      ylab("% of clones") +
      coord_cartesian(ylim = c(0,1.1*max(j_gene_percents_summary_2$mean + j_gene_percents_summary_2$std)), expand = FALSE) +
      scale_fill_manual(values = COLOR, name = "") +
      geom_text(aes(x = gene, y = sig_ypos, label = sig), nudge_y = 0.02*max(j_gene_percents_summary_2$sig_ypos), 
                size = 5)
    print(p)
    dev.off()
    
    
  }
}
