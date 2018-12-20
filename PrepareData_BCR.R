rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies for pregnancy
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: 
###         
###
### Author: Silvia Pineda
### Date: December, 2018
############################################################################################

setwd("/Users/Pinedasans/ImmuneRep_Pregnancy/")

##Read all the files and save into and Rdata all together
files <- list.files("Data/BCR/")

data <- c()
for(i in files) {
  cat(i, "\n")
  t <- read.delim(paste("Data/BCR/", i, sep = ""))
  data <- rbind(data, t)
}
save(data, file="Data/BCR_data.RData")

##############################
###Some Quality Control
#############################
load("Data/BCR_data.RData")
##Discard all the V_score < 200
data_qc<-data[which(data$v_score>=200),]
##Discard the non-functional sequences
data_qc<-data_qc[which(data_qc$productive=="t"),]

##Extract the gene from the segment with the allele
data_qc$v_gene <- gsub("\\*", "", substr(data_qc$v_segment, 1, 8))
data_qc$j_gene <- gsub("\\*", "", substr(data_qc$j_segment, 1, 5))
data_qc$d_gene <- gsub("\\*", "", substr(data_qc$d_segment, 1, 8))

###Extract the CDR3 region
data_qc$cdr3_seq <- gsub(" ","", data_qc$cdr3_seq_nt_q)
###Extract the isotype
data_qc$isotype <- substr(data_qc$isosubtype, 1, 4)

##Count the somatic hypermutations
data_qc$Vlength<-nchar(as.character(data_qc$v_sequence))
v_sequence = as.character(data_qc$v_sequence)
data_qc$SHM<-sapply(regmatches(v_sequence, gregexpr("[A-Z]", v_sequence, perl=TRUE)), length)
data_qc$SHM_freq<-data_qc$SHM/data_qc$Vlength

##count the CDR3 length
data_qc$CDR3_length<-nchar(as.character(data_qc$cdr3_seq))

##count the read length
data_qc$read_length<-nchar(as.character(data_qc$trimmed_sequence))

##Read counts and clones per sample and data point
read_count <- table(data_qc$sample_label)
read_count_isotype <- table(data_qc$sample_label, data_qc$isotype)
colnames(read_count_isotype)[1] = "UNMAPPED"
reads <- cbind(read_count,read_count_isotype)

##Count number of clones per sample 
data_qc$V_J_lenghCDR3_Clone_igh = paste(data_qc$v_gene, data_qc$j_gene, nchar(data_qc$cdr3_seq),data_qc$igh_clone_id,sep="_")
read_count_ighClones<- unique(data_qc[,c("sample_label","V_J_lenghCDR3_Clone_igh")])
clones_igh<-data.matrix(table(read_count_ighClones$sample_label))
colnames(clones_igh)<-c("clones")

reads_clones_igh<-cbind(reads,clones_igh)

#####SHM
sum_SHM_freq<-aggregate(data_qc$SHM_freq, by=list(data_qc$sample_label), FUN=sum)
id<-match(rownames(reads_clones_igh),sum_SHM_freq$Group.1)
colnames(sum_SHM_freq)[2]<-c("SHM")
reads_clones_igh_SHM<-cbind(reads_clones_igh,sum_SHM_freq$SHM[id])
colnames(reads_clones_igh_SHM)[9]<-c("SHM")

###length of CDR3
cdr3_length<-aggregate(data_qc$CDR3_length,by=list(data_qc$sample_label), FUN=mean)
id<-match(rownames(reads_clones_igh_SHM),cdr3_length$Group.1)
reads_clones_igh_SHM_cdr3length<-cbind(reads_clones_igh_SHM,cdr3_length$x[id])
colnames(reads_clones_igh_SHM_cdr3length)[10]<-c("mean_CDR3_length")

summary_data<-reads_clones_igh_SHM_cdr3length
save(data_qc,summary_data,file="Data/BCR_data_summary.RData")

