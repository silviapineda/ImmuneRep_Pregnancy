###########################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies for pregnancy
###          iRepertoire data (prepared files)
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: extracting summaries prepared by iRepertoire
###          files ending in _0 - non-normalized (all reads); ending in _1 - normalized (unique clones only)
###
### Author: Brian Le
### Date: February, 2019
############################################################################################

directory <- "iRepertoire/"

summary_path_1 <- "iRepertoire/x5714-JH-0003-UCSFSperDUOA_summary.csv"
summary_path_2 <- "iRepertoire/x5714-JH-0004-UCSFSperDUOB_summary.csv"
sample_outcome_path <- "iRepertoire/sample_outcome.csv"

##combine iRepertoire's summary data
summary_data <- rbind(read.csv(summary_path_1),
                      read.csv(summary_path_2))

summary_data$sample_name <- gsub("_.*", "", summary_data$Sample)
summary_data$sample_number <- gsub("[A-Z]", "", summary_data$sample_name)
summary_data$sample_type <- ifelse(substr(summary_data$sample_name, 1, 1)=="M","maternal_blood","decidua")
summary_data$sample_population <- gsub("^[A-Z]*[0-9]*_", "", summary_data$Sample)

#formatting sample_population to be consistent between CD19/CD8/CD4[non Treg or Treg-]/Treg

summary_data$sample_population_type <- ifelse(grepl("CD19", summary_data$sample_population), "CD19",
                                              ifelse(grepl("CD8", summary_data$sample_population), "CD8",
                                                     ifelse(grepl("CD4", summary_data$sample_population), "CD4", "TCR"
                                                     )))

##read in outcomes (normal vs preterm labor)
sample_outcome <- read.csv(sample_outcome_path)
summary_data <- merge(summary_data, sample_outcome,
                      by = "sample_number",
                      all.x = TRUE, all.y = FALSE)

samples <- as.character(summary_data$Sample.id)
sample_locations <- paste("iRepertoire", summary_data$Project, summary_data$Sample.id, sep = "/")

summaries <- c("0_CDR3Length", "1_CDR3Length",
               "0_Naddition", "1_Naddition",
               "CDR3_list_1", "CDR3_list_2", "CDRs",
               "J_0_trim", "J_0_usage", "J_1_trim", "J_1_usage",
               "V_0_trim", "V_0_usage", "V_1_trim", "V_1_usage")

for (summary in summaries){
  
  cat(summary, "\n")
  
  data <- c()
  for(i in sample_locations) {
    t <- read.csv(paste(i, "_", summary, ".csv", sep = ""), header = FALSE)
    t$sample <- gsub(".*/", "", i)
    data <- rbind(data, t)
  }
  
  write.csv(data, file = paste(directory, 'iRepertoire_', summary, '_allsamples.csv', sep = ""))
  save(data, file=paste(directory, 'iRepertoire_', summary, '_data.RData', sep = ""))
  
}

write.csv(summary_data, file = paste(directory, 'iRepertoire_summary_data.csv', sep = ""))
