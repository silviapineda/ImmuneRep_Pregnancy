###########################################################################################
### PROJECT: Immune Repertoire. Shared clone analysis of pregnancy outcome iRep data
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: after doing clone analysis (without downsampling)
###          for TCR / CD4 / CD8 / CD19 (listed as BCR)
###
### Author: Brian Le
### Date: April, 2019
############################################################################################
library(RColorBrewer)
library(pheatmap)

directory <- "iRepertoire/"

summary_file_path <- paste(directory, 'iRepertoire_data/iRepertoire_summary_data.csv', sep = "")
summary_data <- read.csv(summary_file_path)

get_isotype <- function(class){
  
  class <- gsub("hTRB", "", class)
  class <- gsub("hIGH", "", class)
  class <- gsub("\\*.*", "", class)
  class <- gsub("[0-9]", "", class) #combines C1/C2 in TCR, G1/G3/G4 in BCR
  
  return(class)
}

#start with BCR
pop <- "CD4"

data_qc <- read.csv(paste0("iRepertoire_clones_",pop,"_processed.csv"))
data_qc$isotype <- with(data_qc, get_isotype(C))
data_qc$V_J_lengthCDR3_Clone <- paste(data_qc$V_J_lengthCDR3, data_qc$numberClone)
#data_all <- rbind(data_all, data_qc[,c("sample", "copy", "isotype", "V_J_lengthCDR3_Clone")])
# 
# length(unique(data_qc$V_J_lengthCDR3_Clone)) #69274 clones
# nrow(unique(data_qc[, c("isotype", "V_J_lengthCDR3_Clone")])) #134256 clones by isotype
# nrow(unique(data_qc[, c("sample", "isotype", "V_J_lengthCDR3_Clone")])) #246631 clones by isotype and sample

clones_unique <- unique(data_qc[, c("sample", "isotype", "V_J_lengthCDR3_Clone")])

#clones_by_sample <- table(clones_unique$sample, clones_unique$V_J_lengthCDR3_Clone)
#table(clones_by_sample) #436 each shared across 5 samples, ignoring isotype :o

##look across all 20 samples, then stratify from there, per isotype
clones_subsample <- clones_unique[which(clones_unique$isotype == "C"),]
#print(length(clones_subsample$V_J_lengthCDR3_Clone)) #total clones
print(length(unique(clones_subsample$V_J_lengthCDR3_Clone))) #total unique clones
clones_by_sample <- as.data.frame.matrix(table(clones_subsample$V_J_lengthCDR3_Clone, clones_subsample$sample))
clones_by_sample <- clones_by_sample[rowSums(clones_by_sample) > 1,] 
#filter to only clones present across multiple samples
print(nrow(clones_by_sample)) #shared by 2 or more

samples_term <- summary_data[summary_data$outcome == "normal", "Sample.id"]
samples_preterm <- summary_data[summary_data$outcome == "ptl", "Sample.id"]
clones_by_sample_subset <- clones_by_sample[, colnames(clones_by_sample) %in% samples_term]
clones_by_sample_subset <- clones_by_sample_subset[rowSums(clones_by_sample_subset) > 1,]
print(paste("term", nrow(clones_by_sample_subset)))
clones_by_sample_subset <- clones_by_sample[, colnames(clones_by_sample) %in% samples_preterm]
clones_by_sample_subset <- clones_by_sample_subset[rowSums(clones_by_sample_subset) > 1,]
print(paste("ptb", nrow(clones_by_sample_subset)))

#across all 20 samples, looking at pairwise data...
pairwise_data <- summary_data[sort(summary_data$sample_name)]
summary_data$pair_name <- gsub("[A-Z]", "", summary_data$sample_name)
sample_name <- unique(summary_data$pair_name)
sample_decidua <- c()
sample_maternalblood <- c()
for (sample in sample_name){
  sample_decidua <- c(sample_decidua,
                      as.character(summary_data[(summary_data$pair_name == sample) &
                                   (summary_data$sample_type == 'decidua') &
                                   (summary_data$sample_population_type == pop), "Sample.id"])) #if pop is BCR, need to fill in CD19
  sample_maternalblood <- c(sample_maternalblood,
                      as.character(summary_data[(summary_data$pair_name == sample) &
                                     (summary_data$sample_type == 'maternal_blood') &
                                     (summary_data$sample_population_type == pop), "Sample.id"]))
}

for (i in 1:length(sample_name)){
  print(i)
  #if column doesn't exist because sample has no reads of this type, populate with 0s
  if(!(sample_decidua[i] %in% colnames(clones_by_sample))){
    clones_by_sample[,sample_decidua[i]] <- 0
  }
  if(!(sample_maternalblood[i] %in% colnames(clones_by_sample))){
    clones_by_sample[,sample_maternalblood[i]] <- 0
  }
  clones_by_sample[, sample_name[i]] <- clones_by_sample[, sample_decidua[i]] + clones_by_sample[, sample_maternalblood[i]]
}
clones_by_sample_pairs <- clones_by_sample[, sample_name]
clones_by_sample_pairs[clones_by_sample_pairs == 1] <- 0
clones_by_sample_pairs <- clones_by_sample_pairs[rowSums(clones_by_sample_pairs) > 1,]
print(nrow(clones_by_sample_pairs))


#start by looking at maternal_blood BCR samples, isotype IgM
samples <- summary_data[summary_data$sample_type == "decidua", "Sample.id"]
clones_subsample <- clones_unique[clones_unique$sample %in% samples,]
clones_subsample <- clones_subsample[clones_subsample$isotype == "C",]
length(unique(clones_subsample$V_J_lengthCDR3_Clone)) #54,664 total clones for maternal blood, IgM
clones_by_sample <- as.data.frame.matrix(table(clones_subsample$V_J_lengthCDR3_Clone, clones_subsample$sample))
clones_by_sample <- clones_by_sample[rowSums(clones_by_sample) > 1,] 
#clones_by_sample$total <- rowSums(clones_by_sample)
#filter to only clones present across multiple samples
nrow(clones_by_sample)
#23,124 rows/clones

###

samples_term <- summary_data[summary_data$outcome == "normal", "Sample.id"]
samples_preterm <- summary_data[summary_data$outcome == "ptl", "Sample.id"]
clones_by_sample_subset <- clones_by_sample[, colnames(clones_by_sample) %in% samples_term]
clones_by_sample_subset <- clones_by_sample_subset[rowSums(clones_by_sample_subset) > 1,]
print(paste("term", nrow(clones_by_sample_subset)))
clones_by_sample_subset <- clones_by_sample[, colnames(clones_by_sample) %in% samples_preterm]
clones_by_sample_subset <- clones_by_sample_subset[rowSums(clones_by_sample_subset) > 1,]
print(paste("ptb", nrow(clones_by_sample_subset)))


########
#let's do term vs. preterm!
term_samples <- summary_data[summary_data$outcome == "normal", "Sample.id"]
preterm_samples <- summary_data[summary_data$outcome == "ptl", "Sample.id"]

clones_by_sample_term <- clones_by_sample[, colnames(clones_by_sample) %in% term_samples]
clones_by_sample_term$total_term <- rowSums(clones_by_sample_term)
clones_by_sample_preterm <- clones_by_sample[, colnames(clones_by_sample) %in% preterm_samples]
clones_by_sample_preterm$total_preterm <- rowSums(clones_by_sample_preterm)

clones_by_sample$total_term <- clones_by_sample_term$total_term
clones_by_sample$total_preterm <- clones_by_sample_preterm$total_preterm
clones_by_sample$total <- clones_by_sample$total_term + clones_by_sample$total_preterm

#sample size too small, not at all a useful measurement... think of something else
clones_by_sample$phyper <- phyper(clones_by_sample$total_term, 6, 4, clones_by_sample$total, lower.tail = TRUE)
clones_by_sample$q <- qvalue(clones_by_sample$phyper)$qvalues
#phyper(5, 6, 4, 6, lower.tail = TRUE) #5 term (or fewer, inclusive) selected, 6 term total, 4 preterm total, 6 balls drawn

table(clones_by_sample$total)
table(clones_by_sample$total, clones_by_sample$total_preterm)

dhyper(4, 4, 6, 5)*75

#CD4 population
#2 samples: observed, 217-542-71, expected: 276-442-110
#3 samples: observed 33-187-100-6, expected: 54-163-98-10
#4 samples: observed 4-52-91-7-0, expected: 11-59-66-18-1
#5 samples: observed 0-11-52-12-0, expectd: 2-18-36-18-2

##for each pair of samples, how many shared clones do they have?
samples <- colnames(clones_by_sample[1:10])

pairwise_clones <- NULL

for (i in 1:(length(samples)-1)){
  for (j in (i+1):length(samples)){
    print(paste(i, j))
    
    clones_selected <- clones_by_sample[, c(samples[i], samples[j])]
    clones_selected$total <- rowSums(clones_selected)
    
    pairwise_clones <- rbind(pairwise_clones,
                             data.frame(sample1 = samples[i],
                                        sample2 = samples[j],
                                        total_shared = sum(clones_selected$total > 1)))
  }
}

pairwise_clones <- merge(pairwise_clones, summary_data[,c("Sample.id", "outcome")],
                         by.x = "sample1",
                         by.y = "Sample.id",
                         all.x = TRUE, all.y = FALSE)

colnames(pairwise_clones)[4] <- "outcome1"

pairwise_clones <- merge(pairwise_clones, summary_data[,c("Sample.id", "outcome")],
                         by.x = "sample2",
                         by.y = "Sample.id",
                         all.x = TRUE, all.y = FALSE)
colnames(pairwise_clones)[5] <- "outcome2"

#maybe plot this as a heatmap with samples on the axes and number of shared clones as bubbles?

pop <- 'BCR'
isotype <- 'M'

data_qc <- read.csv(paste0("iRepertoire_clones_",pop,"_processed.csv"))
data_qc$isotype <- with(data_qc, get_isotype(C))
data_qc$V_J_lengthCDR3_Clone <- paste(data_qc$V_J_lengthCDR3, data_qc$numberClone)
clones_unique <- unique(data_qc[, c("sample", "isotype", "V_J_lengthCDR3_Clone")])
##look across all 20 samples, then stratify from there, per isotype
clones_subsample <- clones_unique[which(clones_unique$isotype == isotype),]
clones_by_sample <- as.data.frame.matrix(table(clones_subsample$V_J_lengthCDR3_Clone, clones_subsample$sample))
clones_by_sample <- clones_by_sample[rowSums(clones_by_sample) > 1,] 
#filter to only clones present across multiple samples


pairwise_clones <- matrix(data = 0, nrow = ncol(clones_by_sample), ncol = ncol(clones_by_sample))
samples <- colnames(clones_by_sample)

for (i in 1:length(samples)){
  for (j in 1:length(samples)){
    print(paste(i, j))
    
    clones_selected <- clones_by_sample[, c(samples[i], samples[j])]
    clones_selected$total <- rowSums(clones_selected)
    
    pairwise_clones[i,j] <- sum(clones_selected$total > 1)
  }
}

pairwise_clones <- data.frame(pairwise_clones)

colnames(pairwise_clones) <- samples
rownames(pairwise_clones) <- samples

col_labels <- summary_data[summary_data$Sample.id %in% colnames(pairwise_clones), c("Sample.id", "outcome", "sample_type")]
rownames(col_labels) <- col_labels$Sample.id
col_labels <- subset(col_labels, select = c("outcome", "sample_type"))

#brewer.pal(n = 4, name = "Set1")[3:4] #muted green/purple, [2:1] for sample_type
ann_colors = list(
  outcome = c(normal = "#4DAF4A", ptl = "#984EA3"),
  sample_type = c(maternal_blood = "#E41A1C", decidua = "#377EB8")
)

pheatmap(pairwise_clones,
         annotation_col = col_labels,
         #annotation_colors = ann_colors,
         color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
         #treeheight_row = FALSE,
         show_rownames = FALSE,
         main = "test")


# clones_sample_shared <- clones_by_sample[,c("total_term", "total_preterm", "total")]
# 
# test <- clones_sample_shared[clones_sample_shared$total == 2,]

# clones_subsample <- clones_unique[clones_unique$sample %in% samples,]
# clones_subsample <- clones_subsample[clones_subsample$isotype == "A",]
# clones_by_sample <- data.frame(table(clones_subsample$sample, clones_subsample$V_J_lengthCDR3_Clone))
# 
# 
# #take a step back, ignoring isotype, look at the sample level
# clones_by_sample <- data.frame(table(clones_unique$sample, clones_unique$V_J_lengthCDR3_Clone))


#20 samples, 5 isotypes, let's see how much is shared and how much isn't

##testing area
#BCR: heatmap of all samples split into isotypes, see if they cluster by shared clones
#jk runs out of memory

# pop <- 'BCR'
# 
# data_qc <- read.csv(paste0("iRepertoire_clones_",pop,"_processed.csv"))
# data_qc$isotype <- with(data_qc, get_isotype(C))
# data_qc <- data_qc[data_qc$isotype != "-",]
# data_qc$V_J_lengthCDR3_Clone <- paste(data_qc$V_J_lengthCDR3, data_qc$numberClone)
# 
# sample_types <- c("maternal_blood", "decidua")
# isotypes <- unique(data_qc$isotype)
# 
# data_qc$sampleisotype <- paste0(data_qc$sample, data_qc$isotype)
# 
# clones_unique <- unique(data_qc[, c("sampleisotype", "V_J_lengthCDR3_Clone")])
# 
# clones_by_sample <- as.data.frame.matrix(table(clones_unique$V_J_lengthCDR3_Clone, clones_unique$sampleisotype))
# clones_by_sample <- clones_by_sample[rowSums(clones_by_sample) > 1,] 
# #filter to only clones present across multiple samples
# 
# 
# col_labels <- data.frame(sample = colnames(clones_by_sample), isotype = gsub("[0-9]", "", colnames(clones_by_sample)))
# rownames(col_labels) <- col_labels$sample
# col_labels <- subset(col_labels, select = c("isotype"))
# 
# pheatmap(clones_by_sample,
#          annotation_col = col_labels,
#          #annotation_colors = ann_colors,
#          color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
#          #treeheight_row = FALSE,
#          show_rownames = FALSE,
#          main = paste(pop, "-", sample_type, "- isotype", isotype, "- shared clones"))

#####PLOTS

###creating heatmaps of shared clones (1/red if present in sample, 0/yellow if not in sample) and clustering

pops <- c("BCR", "CD4", "CD8", "TCR")

for (pop in pops){
  print(pop)
  data_qc <- read.csv(paste0("iRepertoire_clones_",pop,"_processed.csv"))
  data_qc$isotype <- with(data_qc, get_isotype(C))
  data_qc <- data_qc[data_qc$isotype != "-",]
  data_qc$V_J_lengthCDR3_Clone <- paste(data_qc$V_J_lengthCDR3, data_qc$numberClone)

  sample_types <- c("maternal_blood", "decidua")
  isotypes <- unique(data_qc$isotype)
  
  clones_unique <- unique(data_qc[, c("sample", "isotype", "V_J_lengthCDR3_Clone")])
  
  
  
  for (sample_type in sample_types){
    print(sample_type)
    #start by looking at maternal_blood BCR samples, isotype IgM
    samples <- summary_data[summary_data$sample_type == sample_type, "Sample.id"]
    
    for (isotype in isotypes){
      print(isotype)
      plot_directory <- paste0(directory, "Results/shared_clones/", pop, "/")
      dir.create(plot_directory, showWarnings = FALSE)
      clones_subsample <- clones_unique[clones_unique$sample %in% samples,]
      clones_subsample <- clones_subsample[clones_subsample$isotype == isotype,]

      clones_by_sample <- as.data.frame.matrix(table(clones_subsample$V_J_lengthCDR3_Clone, clones_subsample$sample))
      clones_by_sample <- clones_by_sample[rowSums(clones_by_sample) > 1,] 
      #filter to only clones present across multiple samples

      if(nrow(clones_by_sample) > 0){ #if shared clones exist, plot heatmap
        col_labels <- summary_data[summary_data$Sample.id %in% samples &
                                     summary_data$Sample.id %in% colnames(clones_by_sample), c("Sample.id", "outcome")]
        rownames(col_labels) <- col_labels$Sample.id
        col_labels <- subset(col_labels, select = "outcome")
        
        #brewer.pal(n = 4, name = "Set1")[3:4] #muted green/purple
        ann_colors = list(
          outcome = c(normal = "#4DAF4A", ptl = "#984EA3")
        )
        
        tiff(paste0(plot_directory, sample_type, "_isotype", isotype, "_sharedclones_heatmap.tiff"),res=300,w=2000,h=2000)
        p <- pheatmap(clones_by_sample,
                 annotation_col = col_labels,
                 annotation_colors = ann_colors,
                 color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
                 #treeheight_row = FALSE,
                 show_rownames = FALSE,
                 main = paste(pop, "-", sample_type, "- isotype", isotype, "- shared clones"))
        print(p)
        dev.off()
      }
    }
  }
}

###keeping maternal_blood and decidua together as an additional pheatmap annotation

for (pop in pops){
  print(pop)
  data_qc <- read.csv(paste0("iRepertoire_clones_",pop,"_processed.csv"))
  data_qc$isotype <- with(data_qc, get_isotype(C))
  data_qc <- data_qc[data_qc$isotype != "-",]
  data_qc$V_J_lengthCDR3_Clone <- paste(data_qc$V_J_lengthCDR3, data_qc$numberClone)
  
  isotypes <- unique(data_qc$isotype)
  
  clones_unique <- unique(data_qc[, c("sample", "isotype", "V_J_lengthCDR3_Clone")])
  
  # for (sample_type in sample_types){
  #   print(sample_type)
  #   #start by looking at maternal_blood BCR samples, isotype IgM
  #   samples <- summary_data[summary_data$sample_type == sample_type, "Sample.id"]
    
    for (isotype in isotypes){
      print(isotype)
      plot_directory <- paste0(directory, "Results/shared_clones/", pop, "/")
      dir.create(plot_directory, showWarnings = FALSE)

      clones_subsample <- clones_unique[clones_unique$isotype == isotype,]
      
      clones_by_sample <- as.data.frame.matrix(table(clones_subsample$V_J_lengthCDR3_Clone, clones_subsample$sample))
      clones_by_sample <- clones_by_sample[rowSums(clones_by_sample) > 1,] 
      #filter to only clones present across multiple samples
      
      if(nrow(clones_by_sample) > 0){ #if shared clones exist, plot heatmap
        col_labels <- summary_data[summary_data$Sample.id %in% colnames(clones_by_sample), c("Sample.id", "outcome", "sample_type", "pair_name")]
        rownames(col_labels) <- col_labels$Sample.id
        col_labels <- subset(col_labels, select = c("outcome", "sample_type", "pair_name"))
        
        #brewer.pal(n = 4, name = "Set1")[3:4] #muted green/purple, [2:1] for sample_type
        ann_colors = list(
          outcome = c(normal = "#4DAF4A", ptl = "#984EA3"),
          sample_type = c(maternal_blood = "#E41A1C", decidua = "#377EB8")
        )
        
        tiff(paste0(plot_directory, "combined_isotype", isotype, "_sharedclones_heatmap.tiff"),res=300,w=2000,h=2000)
        p <- pheatmap(clones_by_sample,
                      annotation_col = col_labels,
                      annotation_colors = ann_colors,
                      color = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
                      #treeheight_row = FALSE,
                      show_rownames = FALSE,
                      main = paste(pop, "- isotype", isotype, "- shared clones"))
        print(p)
        dev.off()
      }
    }
  # }
}

