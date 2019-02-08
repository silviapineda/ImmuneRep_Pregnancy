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
load("Data/BCR/BCR_data.RData")
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

###Delete the unmapped isotypes
data_qc<-data_qc[which(data_qc$isotype!=""),]

##Read counts and clones per sample and data point
read_count <- table(data_qc$sample_label)
read_count_isotype <- table(data_qc$sample_label, data_qc$isotype)
#colnames(read_count_isotype)[1] = "UNMAPPED"
reads <- cbind(read_count,read_count_isotype)

##Count number of clones per sample 
data_qc$V_J_lenghCDR3_Clone_igh = paste(data_qc$v_gene, data_qc$j_gene, nchar(data_qc$cdr3_seq),data_qc$igh_clone_id,sep="_")
read_count_ighClones<- unique(data_qc[,c("sample_label","V_J_lenghCDR3_Clone_igh")])
clones_count<-data.matrix(table(read_count_ighClones$sample_label))
colnames(clones_count)<-c("clones_count")
read_count_ighClones_isotype<- unique(data_qc[,c("sample_label","V_J_lenghCDR3_Clone_igh","isotype")])
clones_igh<-data.matrix(table(read_count_ighClones_isotype$sample_label,read_count_ighClones_isotype$isotype))
colnames(clones_igh)<-c("clones_IGHA","clones_IGHD","clones_IGHE","clones_IGHG","clones_IGHM")
reads_clones_igh<-cbind(reads,clones_count,clones_igh)

#####SHM
sum_SHM_freq<-aggregate(data_qc$SHM_freq, by=list(data_qc$sample_label,data_qc$isotype), FUN=sum)
#sum_SHM_unmapped<-sum_SHM_freq[which(sum_SHM_freq$Group.2==""),]
sum_SHM_IGHA<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="IGHA"),]
sum_SHM_IGHD<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="IGHD"),]
sum_SHM_IGHE<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="IGHE"),]
sum_SHM_IGHG<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="IGHG"),]
sum_SHM_IGHM<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="IGHM"),]
sum_SHM_IGHM<-sum_SHM_freq[which(sum_SHM_freq$Group.2=="IGHM"),]

#id_unmapped<-match(rownames(reads_clones_igh),sum_SHM_unmapped$Group.1)
id_IGHA<-match(rownames(reads_clones_igh),sum_SHM_IGHA$Group.1)
id_IGHD<-match(rownames(reads_clones_igh),sum_SHM_IGHD$Group.1)
id_IGHE<-match(rownames(reads_clones_igh),sum_SHM_IGHE$Group.1)
id_IGHG<-match(rownames(reads_clones_igh),sum_SHM_IGHG$Group.1)
id_IGHM<-match(rownames(reads_clones_igh),sum_SHM_IGHM$Group.1)
reads_clones_igh_SHM<-cbind(reads_clones_igh,sum_SHM_IGHA$x[id_IGHA],sum_SHM_IGHD$x[id_IGHD],
                            sum_SHM_IGHE$x[id_IGHE],sum_SHM_IGHG$x[id_IGHG],sum_SHM_IGHM$x[id_IGHM])
colnames(reads_clones_igh_SHM)[13:17]<-c("SHM_IGHA","SHM_IGHD","SHM_IGHE","SHM_IGHG","SHM_IGHM")

###length of CDR3
cdr3_length<-aggregate(data_qc$CDR3_length,by=list(data_qc$sample_label,data_qc$isotype), FUN=mean)
#cdr3_length_unmapped<-cdr3_length[which(cdr3_length$Group.2==""),]
cdr3_length_IGHA<-cdr3_length[which(cdr3_length$Group.2=="IGHA"),]
cdr3_length_IGHD<-cdr3_length[which(cdr3_length$Group.2=="IGHD"),]
cdr3_length_IGHE<-cdr3_length[which(cdr3_length$Group.2=="IGHE"),]
cdr3_length_IGHG<-cdr3_length[which(cdr3_length$Group.2=="IGHG"),]
cdr3_length_IGHM<-cdr3_length[which(cdr3_length$Group.2=="IGHM"),]

#id_unmapped<-match(rownames(reads_clones_igh_SHM),cdr3_length_unmapped$Group.1)
id_IGHA<-match(rownames(reads_clones_igh_SHM),cdr3_length_IGHA$Group.1)
id_IGHD<-match(rownames(reads_clones_igh_SHM),cdr3_length_IGHD$Group.1)
id_IGHE<-match(rownames(reads_clones_igh_SHM),cdr3_length_IGHE$Group.1)
id_IGHG<-match(rownames(reads_clones_igh_SHM),cdr3_length_IGHG$Group.1)
id_IGHM<-match(rownames(reads_clones_igh_SHM),cdr3_length_IGHM$Group.1)

reads_clones_igh_SHM_cdr3length<-cbind(reads_clones_igh_SHM,cdr3_length_IGHA$x[id_IGHA],cdr3_length_IGHD$x[id_IGHD],
                                       cdr3_length_IGHE$x[id_IGHE],cdr3_length_IGHG$x[id_IGHG],cdr3_length_IGHM$x[id_IGHM])
colnames(reads_clones_igh_SHM_cdr3length)[18:22]<-c("CDR3_length_IGHA","CDR3_length_IGHD","CDR3_length_IGHE","CDR3_length_IGHG","CDR3_length_IGHM")
#

### Add naive and memory B-cells clones from the IGHM isotype
## SHM<=4 is navie SHM>4 is memory
data_qc$IGHM_naive_memory<-ifelse(data_qc$isotype=="IGHM" & data_qc$SHM<=4,"naive",
                                  ifelse(data_qc$isotype=="IGHM" & data_qc$SHM>4,"memory",NA))

###Extract reads
read_count_naive_memory <- data.matrix(table(data_qc$sample_label, data_qc$IGHM_naive_memory))
colnames(read_count_naive_memory) = c("reads_memory","reads_naive")
reads_clones_igh_SHM_cdr3length_naive_memory<-cbind(reads_clones_igh_SHM_cdr3length,read_count_naive_memory)

###Extract clones
clones_naive_memory<- unique(data_qc[,c("sample_label","V_J_lenghCDR3_Clone_igh","IGHM_naive_memory")])
clones_naive_memory_matrix<-data.matrix(table(clones_naive_memory$sample_label,clones_naive_memory$IGHM_naive_memory))
colnames(clones_naive_memory_matrix)<-c("clones_memory","clones_naive")
reads_clones_igh_SHM_cdr3length_naive_memory<-cbind(reads_clones_igh_SHM_cdr3length_naive_memory,
                                                    clones_naive_memory_matrix[,1])
reads_clones_igh_SHM_cdr3length_naive_memory<-cbind(reads_clones_igh_SHM_cdr3length_naive_memory,
                                                    clones_naive_memory_matrix[,2])
colnames(reads_clones_igh_SHM_cdr3length_naive_memory)[25:26]<-c("clones_memory","clones_naive")

summary_data<-data.frame(reads_clones_igh_SHM_cdr3length_naive_memory)

###Introduce a variable to design mother and fetal
info_samples<-read.csv("Data/NormalSamplesScottBoydStanford9_13_18.csv")
id<-match(info_samples$Sample_ID,rownames(summary_data))
summary_data$pairs[id]<-info_samples$Pairs
summary_data$NumCells[id]<-info_samples$NumCells
summary_data$sample<-gsub("\\*", "", substr(summary_data$pairs, 1, 1))
summary_data$sample<-ifelse(summary_data$sample=="M","Mother","Fetus")
summary_data$sample<-factor(summary_data$sample)
summary_data$sample<-relevel(summary_data$sample, ref = "Mother")

## Diversity measures
sample<-rownames(summary_data)
entropy_IGHA<-NULL
entropy_IGHD<-NULL
entropy_IGHE<-NULL
entropy_IGHG<-NULL
entropy_IGHM<-NULL
entropy_naive<-NULL
entropy_memory<-NULL
for (i in 1:length(sample)){
  print(i)
  data_sample_unique<-data_qc[which(data_qc$sample_label==sample[i]),]
  clones_sample<-data_sample_unique[,"V_J_lenghCDR3_Clone_igh"]
  clones_specimen_IGHA<-data_sample_unique[which(data_sample_unique$isotype=="IGHA"),"V_J_lenghCDR3_Clone_igh"]
  clones_specimen_IGHD<-data_sample_unique[which(data_sample_unique$isotype=="IGHD"),"V_J_lenghCDR3_Clone_igh"]
  clones_specimen_IGHE<-data_sample_unique[which(data_sample_unique$isotype=="IGHE"),"V_J_lenghCDR3_Clone_igh"]
  clones_specimen_IGHG<-data_sample_unique[which(data_sample_unique$isotype=="IGHG"),"V_J_lenghCDR3_Clone_igh"]
  clones_specimen_IGHM<-data_sample_unique[which(data_sample_unique$isotype=="IGHM"),"V_J_lenghCDR3_Clone_igh"]
  clones_specimen_naive<-data_sample_unique[which(data_sample_unique$IGHM_naive_memory=="naive"),"V_J_lenghCDR3_Clone_igh"]
  clones_specimen_memory<-data_sample_unique[which(data_sample_unique$IGHM_naive_memory=="memory"),"V_J_lenghCDR3_Clone_igh"]
  #To write file to run with Recon
  #write.delim(data.frame(table(table(clones_specimen_IGHA))),file=paste("clones_IGHA_",specimen_unique[i],".txt",sep=""),sep="\t",col.names=F)
  #write.delim(data.frame(table(table(clones_specimen_IGHD))),file=paste("clones_IGHD_",specimen_unique[i],".txt",sep=""),sep="\t",col.names=F)
  #write.delim(data.frame(table(table(clones_specimen_IGHE))),file=paste("clones_IGHE_",specimen_unique[i],".txt",sep=""),sep="\t",col.names=F)
  #write.delim(data.frame(table(table(clones_specimen_IGHG))),file=paste("clones_IGHG_",specimen_unique[i],".txt",sep=""),sep="\t",col.names=F)
  #write.delim(data.frame(table(table(clones_specimen_IGHM))),file=paste("clones_IGHM_",specimen_unique[i],".txt",sep=""),sep="\t",col.names=F)
  
  fi_IGHA<-as.numeric(table(clones_specimen_IGHA))/length(clones_specimen_IGHA)
  fi_IGHD<-as.numeric(table(clones_specimen_IGHD))/length(clones_specimen_IGHD)
  fi_IGHE<-as.numeric(table(clones_specimen_IGHE))/length(clones_specimen_IGHE)
  fi_IGHG<-as.numeric(table(clones_specimen_IGHG))/length(clones_specimen_IGHG)
  fi_IGHM<-as.numeric(table(clones_specimen_IGHM))/length(clones_specimen_IGHM)
  fi_naive<-as.numeric(table(clones_specimen_naive))/length(clones_specimen_naive)
  fi_memory<-as.numeric(table(clones_specimen_memory))/length(clones_specimen_memory)
  
  hi_IGHA<-fi_IGHA*log2(fi_IGHA)
  hi_IGHD<-fi_IGHD*log2(fi_IGHD)
  hi_IGHE<-fi_IGHE*log2(fi_IGHE)
  hi_IGHG<-fi_IGHG*log2(fi_IGHG)
  hi_IGHM<-fi_IGHM*log2(fi_IGHM)
  hi_naive<-fi_naive*log2(fi_naive)
  hi_memory<-fi_memory*log2(fi_memory)
  
  entropy_IGHA[i]=-sum(hi_IGHA)
  entropy_IGHD[i]=-sum(hi_IGHD)
  entropy_IGHE[i]=-sum(hi_IGHE)
  entropy_IGHG[i]=-sum(hi_IGHG)
  entropy_IGHM[i]=-sum(hi_IGHM)
  entropy_naive[i]=-sum(hi_naive)
  entropy_memory[i]=-sum(hi_memory)
  
}

entropy_norm_IGHA<-entropy_IGHA/max(entropy_IGHA,na.rm = T)
clonality_IGHA<-(1-entropy_norm_IGHA)
entropy_norm_IGHD<-entropy_IGHD/max(entropy_IGHD,na.rm = T)
clonality_IGHD<-(1-entropy_norm_IGHD)
entropy_norm_IGHG<-entropy_IGHG/max(entropy_IGHG,na.rm = T)
clonality_IGHG<-(1-entropy_norm_IGHG)
entropy_norm_IGHM<-entropy_IGHM/max(entropy_IGHM,na.rm = T)
clonality_IGHM<-(1-entropy_norm_IGHM)
entropy_norm_memory<-entropy_memory/max(entropy_memory,na.rm = T)
clonality_memory<-(1-entropy_norm_memory)
entropy_norm_naive<-entropy_naive/max(entropy_naive,na.rm = T)
clonality_naive<-(1-entropy_norm_naive)

diversity<-cbind(entropy_unmapped,entropy_IGHA,entropy_IGHD,entropy_IGHE,entropy_IGHG,entropy_IGHM,entropy_naive,entropy_memory,
                 clonality_IGHA,clonality_IGHD,clonality_IGHG,clonality_IGHM,clonality_memory,clonality_naive)
rownames(diversity)<-sample
summary_data<-cbind(summary_data,diversity)

save(data_qc,summary_data,file="Data/BCR_data_summary.RData")
