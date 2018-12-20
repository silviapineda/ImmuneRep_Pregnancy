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

summary_data<-data.frame(reads_clones_igh_SHM_cdr3length)

###Introduce a variable to design mother and fetal
info_samples<-read.csv("Data/NormalSamplesScottBoydStanford9_13_18.csv")
summary_data$pairs<-info_samples$Pairs
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
  clones_sample<-data_sample_unique[,"V_J_lenghCDR3_Clone_igh"]
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

save(data_qc,summary_data,file="Data/BCR_data_summary.RData")

