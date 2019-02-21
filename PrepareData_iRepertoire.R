###########################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies for pregnancy
###          iRepertoire data
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

data <- c()
for(i in sample_locations) {
  cat(i, "\n")
  t <- read.csv(paste(i, "_pep.csv", sep = ""))
  t$sample <- gsub(".*/", "", i)
  data <- rbind(data, t)
}


save(data, summary_data, file=paste(directory, 'iRepertoire_data.RData', sep = ""))
