###########################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies for pregnancy
###          
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: plots of basic features
###          
###
### Author: Brian Le
### Date: February, 2019
############################################################################################

library(reshape)
library(data.table)
library(ggplot2)

directory <- "iRepertoire/"

summary_file_path <- paste(directory, 'iRepertoire_data/iRepertoire_summary_data.csv', sep = "")
summary_data <- read.csv(summary_file_path)

####Summary plots
sample_pops <- c('CD4', 'CD8', 'CD19', 'TCR')
sample_types <- c('decidua', 'maternal_blood')

COLOR=c("#BEAED4","#7FC97F")

summary_data$order <- ifelse(summary_data$sample_population_type == "CD19", 1,
                             summary_data$order <- ifelse(summary_data$sample_population_type == "TCR", 2,
                                                          summary_data$order <- ifelse(summary_data$sample_population_type == "CD4", 3,
                                                                                       summary_data$order <- ifelse(summary_data$sample_population_type == "CD8", 4, 5))))
summary_data[summary_data$outcome == "ptl", "order"] <- summary_data[summary_data$outcome == "ptl", "order"] + 0.5

tiff(paste(plot_directory, "Reads_all.tiff", sep = ""),res=300,w=4000,h=2000)
p <- ggplot(summary_data, aes(x = reorder(Sample.id, order), fill = outcome)) +
  geom_bar(aes(y = Reads), stat = 'identity') +
  labs(title = paste("Number of reads", sep = "")) +
  xlab("Sample") +
  ylab("Reads") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(p)
dev.off()

tiff(paste(plot_directory, "Unique_CDR3_all.tiff", sep = ""),res=300,w=4000,h=2000)
p <- ggplot(summary_data, aes(x = reorder(Sample.id, order), fill = outcome)) +
  geom_bar(aes(y = Unique.CDR3), stat = 'identity') +
  labs(title = paste("Unique.CDR3 per sample", sep = "")) +
  xlab("Sample") +
  ylab("Unique.CDR3") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(p)
dev.off()

for (pop in sample_pops){
  plot_directory <- "iRepertoire/Results/summary_plots/"
  dir.create(plot_directory)
  
  for (type in sample_types){
    
    summary_data_subset <- summary_data[summary_data$sample_population_type == pop &
                                          summary_data$sample_type == type,]
    summary_data_subset$order <- ifelse(summary_data_subset$outcome == "normal", 1, 2)
    
    tiff(paste(plot_directory, "Entropy_", type, "_", pop, ".tiff", sep = ""),res=300,w=4000,h=2000)
    p <- ggplot(summary_data_subset, aes(x = reorder(sample_name, order), fill = outcome)) +
      geom_bar(aes(y = Entropy), stat = 'identity') +
      labs(title = paste("Shannon entropy - ", type, " - ", pop, sep = "")) +
      xlab("Sample") +
      ylab("Entropy")
    print(p)
    dev.off()

    tiff(paste(plot_directory, "D50_", type, "_", pop, ".tiff", sep = ""),res=300,w=4000,h=2000)
    p <- ggplot(summary_data_subset, aes(x = reorder(sample_name, order), fill = outcome)) +
      geom_bar(aes(y = D50), stat = 'identity') +
      labs(title = paste("D50 - ", type, " - ", pop, sep = "")) +
      xlab("Sample") +
      ylab("Diversity Index (D50)")
    print(p)
    dev.off()
    
    tiff(paste(plot_directory, "Diversity_index_", type, "_", pop, ".tiff", sep = ""),res=300,w=4000,h=2000)
    p <- ggplot(summary_data_subset, aes(x = reorder(sample_name, order), fill = outcome)) +
      geom_bar(aes(y = Diversity.Index), stat = 'identity') +
      labs(title = paste("Diversity Index - ", type, " - ", pop, sep = "")) +
      xlab("Sample") +
      ylab("Diversity Index (D50)")
    print(p)
    dev.off()
    
    tiff(paste(plot_directory, "Reads_", type, "_", pop, ".tiff", sep = ""),res=300,w=4000,h=2000)
    p <- ggplot(summary_data_subset, aes(x = reorder(sample_name, order), fill = outcome)) +
      geom_bar(aes(y = Reads), stat = 'identity') +
      labs(title = paste("Reads - ", type, " - ", pop, sep = "")) +
      xlab("Sample") +
      ylab("Reads")
    print(p)
    dev.off()
    
    tiff(paste(plot_directory, "Clones_", type, "_", pop, ".tiff", sep = ""),res=300,w=4000,h=2000)
    p <- ggplot(summary_data_subset, aes(x = reorder(sample_name, order), fill = outcome)) +
      geom_bar(aes(y = Unique.CDR3), stat = 'identity') +
      labs(title = paste("Clones - ", type, " - ", pop, sep = "")) +
      xlab("Sample") +
      ylab("Clones")
    print(p)
    dev.off()
  }
}

#boxplots of entropy and diversity index

for (pop in sample_pops){

  plot_directory <- paste0(directory, "Results/summary_plots/", pop, "/")
  dir.create(plot_directory, showWarnings = FALSE)
  
  plot_data <- summary_data[summary_data$sample_population_type == pop,]
  
  #decidua: comparing CDR3s term vs ptb
  tiff(paste0(plot_directory, "entropy_", pop, "_decidua.tiff"),res=300,w=1400,h=2000)
  p <- ggboxplot(plot_data[plot_data$sample_type == "decidua",], x = "outcome", y = "Entropy",
                 color = "outcome", add = "jitter", title = paste0(pop, " - Decidua - Entropy"),
                 xlab = "Outcome", ylab = "Entropy") + stat_compare_means()
  print(p)
  dev.off()
  
  #maternal blood: CDR3s, term vs ptb
  tiff(paste0(plot_directory, "entropy_", pop, "_maternalblood.tiff"),res=300,w=1400,h=2000)
  p <- ggboxplot(plot_data[plot_data$sample_type == "maternal_blood",], x = "outcome", y = "Entropy",
                 color = "outcome", add = "jitter", title = paste0(pop, " - Maternal Blood - Entropy"),
                 xlab = "Outcome", ylab = "Entropy") + stat_compare_means()
  print(p)
  dev.off()
  
  #preterm: paired CDR3s, maternal blood vs decidua
  tiff(paste0(plot_directory, "entropy_", pop, "_paired_preterm.tiff"),res=300,w=1400,h=2000)
  p <- ggpaired(plot_data[plot_data$outcome == "ptl",], x = "sample_type", y = "Entropy", 
                id = "sample_number", line.color = "gray",
                color = "sample_type", add = "jitter", title = paste0(pop, " - Preterm - Entropy"),
                xlab = "Sample Type", ylab = "Entropy") +
    stat_compare_means(paired = TRUE)
  print(p)
  dev.off()
  
  #preterm: paired CDR3s, maternal blood vs decidua
  tiff(paste0(plot_directory, "entropy_", pop, "_paired_term.tiff"),res=300,w=1400,h=2000)
  p <- ggpaired(plot_data[plot_data$outcome == "normal",], x = "sample_type", y = "Entropy",
                id = "sample_number", line.color = "gray",
                color = "sample_type", add = "jitter", title = paste0(pop, " - Term - Entropy"),
                xlab = "Sample Type", ylab = "Entropy") +
    stat_compare_means(paired = TRUE)
  print(p)
  dev.off()
}

boxplot(Reads ~ sample_number, summary_data,
        main = "Reads per Pregnancy",
        xlab = "Sample Number")
boxplot(Reads ~ sample_population_type, summary_data,
        main = "Reads per Population",
        xlab = "Sample Population")
boxplot(Reads ~ sample_type * outcome, summary_data,
        main = "Reads per Sample Type / Outcome")

hist(summary_data$Entropy,
     main = "Histogram of Entropy Distribution",
     xlab = "Entropy")
