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
plot_directory <- "iRepertoire/Results/summary_plots/"
dir.create(plot_directory)

summary_file_path <- paste(directory, 'iRepertoire_summary_data.csv', sep = "")
summary_data <- read.csv(summary_file_path)

####Summary plots
sample_pops <- c('CD4', 'CD8', 'CD19', 'TCR')
sample_types <- c('decidua', 'maternal_blood')

COLOR=c("#BEAED4","#7FC97F")

for (pop in sample_pops){
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
      xlab("Reads") +
      ylab("Entropy")
    print(p)
    dev.off()
    
    tiff(paste(plot_directory, "Clones_", type, "_", pop, ".tiff", sep = ""),res=300,w=4000,h=2000)
    p <- ggplot(summary_data_subset, aes(x = reorder(sample_name, order), fill = outcome)) +
      geom_bar(aes(y = Unique.CDR3), stat = 'identity') +
      labs(title = paste("Shannon entropy - ", type, " - ", pop, sep = "")) +
      xlab("Sample") +
      ylab("Entropy")
    print(p)
    dev.off()
  }
}
