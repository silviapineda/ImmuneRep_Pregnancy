rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire. Analysis T cells antibodies for pregnancy
###
### CITATION: Modified from Silvia's PrepareData_BCR.R
###
### PROCESS: 
###           
### DESCRIP: 
###         
###
### Author: Brian Le
### Date: December, 2018
############################################################################################

#Run once to combine data files together and save as R workspace

##Read all the files and save into and Rdata all together
directory <- "TCR/"

files <- list.files(directory)

data <- c()
for(i in files) {
  cat(i, "\n")
  t <- read.delim(paste(directory, i, sep = ""))
  data <- rbind(data, t)
}
save(data, file=paste(directory, 'TCR_data.RData', sep = ""))

##############################
###Some Quality Control
#############################
directory <- "TCR/"

load(paste(directory, 'TCR_data.RData', sep = ""))
##75% of v-scores are below 124, discarding all below 200 results in no leftover data
##Instead of filtering by v-score, filter by num_mismatches / length

mismatch_threshold <- 0.99
mask <- ((nchar(as.character(data$v_sequence)) -
         lengths(regmatches(data$v_sequence,
                            gregexpr("[A-Z]",
                                     data$v_sequence))))/
           nchar(as.character(data$v_sequence))) > mismatch_threshold

data_qc<-data[mask,]
##Discard the non-functional sequences
data_qc<-data_qc[which(data_qc$productive=="t"),]

##Extract the gene from the segment with the allele
#v_segment is variable in length: 8, 9, 10, 11, or 15 characters
#j_segment is empty or 8 characters ending with *01 (discard that segment)
#d_segment is empty or 10 characters ending with *0x
#in all cases, keep gene abbreviation appearing prior to *
data_qc$v_gene <- gsub("\\*.*$", "", data_qc$v_segment)
data_qc$j_gene <- gsub("\\*.*$", "", data_qc$j_segment)
data_qc$d_gene <- gsub("\\*.*$", "", data_qc$d_segment)

###Extract the CDR3 region
data_qc$cdr3_seq <- gsub(" ","", data_qc$cdr3_seq_nt_q)
###Extract the isotype
#isosubtype column doesn't exist in this data
#data_qc$isotype <- substr(data_qc$isosubtype, 1, 4)

##Count the somatic hypermutations
#SHM does not occur in t-cells
#data_qc$Vlength<-nchar(as.character(data_qc$v_sequence))
#v_sequence = as.character(data_qc$v_sequence)
#data_qc$SHM<-sapply(regmatches(v_sequence, gregexpr("[A-Z]", v_sequence, perl=TRUE)), length)
#data_qc$SHM_freq<-data_qc$SHM/data_qc$Vlength

##count the CDR3 length
data_qc$CDR3_length<-nchar(as.character(data_qc$cdr3_seq))

##count the read length
data_qc$read_length<-nchar(as.character(data_qc$trimmed_sequence))

##Read counts per sample (no isosubtypes in TCR)
read_count <- table(data_qc$sample_label)
reads <- cbind(read_count)

##Count number of clones per sample 
data_qc$V_J_length_CDR3_Clone_tcrb = paste(data_qc$v_gene, data_qc$j_gene, nchar(data_qc$cdr3_seq), data_qc$tcrb_clone_id,sep="_")
read_count_tcrbClones<- unique(data_qc[,c("sample_label","V_J_length_CDR3_Clone_tcrb")])
clones_tcrb<-data.matrix(table(read_count_tcrbClones$sample_label))
colnames(clones_tcrb)<-c("clones")

reads_clones_tcrb<-cbind(reads,clones_tcrb)

###length of CDR3
cdr3_length<-aggregate(data_qc$CDR3_length,by=list(data_qc$sample_label), FUN=mean)
id<-match(rownames(reads_clones_tcrb),cdr3_length$Group.1)
reads_clones_tcrb_cdr3length<-cbind(reads_clones_tcrb,cdr3_length$x[id])
colnames(reads_clones_tcrb_cdr3length)[3]<-c("mean_CDR3_length")

summary_data<-data.frame(reads_clones_tcrb_cdr3length)

###Introduce a variable to design mother and fetal
#info_samples<-read.csv("Data/NormalSamplesScottBoydStanford9_13_18.csv")
summary_data$pairs<-row.names(summary_data)
summary_data$sample<-gsub("\\*", "", substr(summary_data$pairs, 1, 1))
summary_data$sample<-ifelse(summary_data$sample=="M","Mother","Fetus")
summary_data$sample<-factor(summary_data$sample)
summary_data$sample<-relevel(summary_data$sample, ref = "Mother")

## Diversity measures
sample<-rownames(summary_data)
entropy<-rep(NA,length(sample))
simpson<-rep(NA,length(sample))
for (i in 1:length(sample)){
  print(i)
  data_sample_unique<-data_qc[which(data_qc$sample_label==sample[i]),]
  clones_sample<-data_sample_unique[,"V_J_length_CDR3_Clone_tcrb"]
  #To write file to run with Recon
  #write.delim(data.frame(table(table(clones_sample))),file=paste("clones_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  fi<-as.numeric(table(clones_sample))/length(clones_sample)
  hi<-fi*log2(fi)
  entropy[i]=-sum(hi) #entropy(table(clones_specimen)) returns the same result but by default is the natural logarithm (log)
  simpson[i]=sum(fi*fi) 
}
entropy_norm<-entropy/max(entropy,na.rm = T)
clonality<-(1-entropy_norm)
names(clonality)<-sample
diversity<-cbind(clonality,entropy,simpson)

summary_data<-cbind(summary_data,diversity)

save(data_qc,summary_data,file=paste(directory, "TCR_data_summary.RData", sep = ""))

