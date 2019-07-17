#rm(list = ls(all = TRUE))
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
files <- list.files("Data/BCR_preT_and_cntrls//")

data <- c()
for(i in files) {
  cat(i, "\n")
  t <- read.delim(paste("Data/BCR_preT_and_cntrls/", i, sep = ""))
  data <- rbind(data, t)
}
save(data, file="Data/BCR_PTB_Term_data.RData")

##############################
###Some Quality Control
#############################
load("Data/BCR_PTB_Term_data.RData")
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

####CDR3 missing filter out
data_qc<-data_qc[which(data_qc$CDR3_length>4),]

##count the read length
data_qc$read_length<-nchar(as.character(data_qc$trimmed_sequence))

###Delete the unmapped isotypes
data_qc<-data_qc[which(data_qc$isotype!=""),]
save(data_qc, file="Data/BCR_PTB_Term_data_qc.RData")

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

### Add unmutated and mutated B-cells clones from the IGHM isotype
## SHM<=1% is navie SHM>1% is memory
data_qc$mutated_rate<-data_qc$SHM/data_qc$Vlength
data_qc$IGHM_mutated<-ifelse(data_qc$isotype=="IGHM" & data_qc$mutated_rate<=0.01,"unmutated",
                                  ifelse(data_qc$isotype=="IGHM" & data_qc$mutated_rate>0.01,"mutated",NA))
data_qc$IGHD_mutated<-ifelse(data_qc$isotype=="IGHD" & data_qc$mutated_rate<=0.01,"unmutated",
                             ifelse(data_qc$isotype=="IGHD" & data_qc$mutated_rate>0.01,"mutated",NA))

###Extract reads
read_count_mutated_IGHM <- data.matrix(table(data_qc$sample_label, data_qc$IGHM_mutated))
colnames(read_count_mutated_IGHM) = c("mutated_IGHM","unmutaed_IGHM")
read_count_mutated_IGHD <- data.matrix(table(data_qc$sample_label, data_qc$IGHD_mutated))
colnames(read_count_mutated_IGHD) = c("mutated_IGHD","unmutaed_IGHD")
reads_clones_igh_SHM_cdr3length_naive_memory<-cbind(reads_clones_igh_SHM_cdr3length,read_count_mutated_IGHM,read_count_mutated_IGHD)

###Extract clones
clones_mutated_IGHM<- unique(data_qc[,c("sample_label","V_J_lenghCDR3_Clone_igh","IGHM_mutated")])
clones_mutated_IGHM_matrix<-data.matrix(table(clones_mutated_IGHM$sample_label,clones_mutated_IGHM$IGHM_mutated))
colnames(clones_mutated_IGHM_matrix)<-c("clones_mutated_IGHM","clones_unmutated_IGHM")
clones_mutated_IGHD<- unique(data_qc[,c("sample_label","V_J_lenghCDR3_Clone_igh","IGHD_mutated")])
clones_mutated_IGHD_matrix<-data.matrix(table(clones_mutated_IGHD$sample_label,clones_mutated_IGHD$IGHD_mutated))
colnames(clones_mutated_IGHD_matrix)<-c("clones_mutated_IGHD","clones_unmutated_IGHD")
reads_clones_igh_SHM_cdr3length_naive_memory<-cbind(reads_clones_igh_SHM_cdr3length_naive_memory,
                                                    clones_mutated_IGHM_matrix[,1:2],clones_mutated_IGHD_matrix[,1:2])


####cdr3 length for naive and memory
cdr3_length_IGHM<-aggregate(data_qc$CDR3_length,by=list(data_qc$sample_label,data_qc$IGHM_mutated), FUN=mean)
cdr3_length_unmutated_IGHM<-cdr3_length_IGHM[which(cdr3_length_IGHM$Group.2=="unmutated"),]
cdr3_length_mutated_IGHM<-cdr3_length_IGHM[which(cdr3_length_IGHM$Group.2=="mutated"),]
cdr3_length_IGHD<-aggregate(data_qc$CDR3_length,by=list(data_qc$sample_label,data_qc$IGHD_mutated), FUN=mean)
cdr3_length_unmutated_IGHD<-cdr3_length_IGHD[which(cdr3_length_IGHD$Group.2=="unmutated"),]
cdr3_length_mutated_IGHD<-cdr3_length_IGHD[which(cdr3_length_IGHD$Group.2=="mutated"),]

id_unmutated_IGHM<-match(rownames(reads_clones_igh_SHM_cdr3length_naive_memory),cdr3_length_unmutated_IGHM$Group.1)
id_mutated_IGHM<-match(rownames(reads_clones_igh_SHM_cdr3length_naive_memory),cdr3_length_mutated_IGHM$Group.1)
id_unmutated_IGHD<-match(rownames(reads_clones_igh_SHM_cdr3length_naive_memory),cdr3_length_unmutated_IGHD$Group.1)
id_mutated_IGHD<-match(rownames(reads_clones_igh_SHM_cdr3length_naive_memory),cdr3_length_mutated_IGHD$Group.1)

reads_clones_igh_SHM_cdr3length_naive_memory_cdr3length<-cbind(reads_clones_igh_SHM_cdr3length_naive_memory,cdr3_length_unmutated_IGHM$x[id_unmutated_IGHM],
                                                               cdr3_length_mutated_IGHM$x[id_mutated_IGHM],cdr3_length_unmutated_IGHD$x[id_unmutated_IGHD],
                                                               cdr3_length_mutated_IGHD$x[id_mutated_IGHD])
colnames(reads_clones_igh_SHM_cdr3length_naive_memory_cdr3length)[31:34]<-c("CDR3_length_unmutated_IGHM","CDR3_length_mutated_IGHM",
                                                                            "CDR3_length_unmutated_IGHD","CDR3_length_mutated_IGHD")

summary_data<-data.frame(reads_clones_igh_SHM_cdr3length_naive_memory_cdr3length)

###Introduce a variable to design mother and fetal
#info_samples<-read.csv("Data/NormalSamplesScottBoydStanford9_13_18.csv") ##First batch
info_samples<-read.csv("Data/BCR_PTB_Term/BCR_PTB_Term.csv") ##Second batch
id<-match(info_samples$Sample_ID,rownames(summary_data))
summary_data$pairs[id]<-info_samples$Pairs
summary_data$NumCells[id]<-info_samples$NumCells
summary_data$Outcome[id]<-info_samples$X
summary_data$sample<-gsub("\\*", "", substr(rownames(summary_data), 1, 1))
summary_data$sample<-ifelse(summary_data$sample=="M","Maternal","Fetal")
summary_data$sample<-factor(summary_data$sample)
summary_data$sample<-relevel(summary_data$sample, ref = "Maternal")
summary_data$Outcome<-factor(summary_data$Outcome)
summary_data$Outcome<-relevel(summary_data$Outcome, ref = "Term")

## Diversity measures
sample<-rownames(summary_data)
entropy_IGHA<-NULL
entropy_IGHD<-NULL
entropy_IGHE<-NULL
entropy_IGHG<-NULL
entropy_IGHM<-NULL
entropy_unmutated_IGHM<-NULL
entropy_mutated_IGHM<-NULL
entropy_unmutated_IGHD<-NULL
entropy_mutated_IGHD<-NULL
for (i in 1:length(sample)){
  print(i)
  data_sample_unique<-data_qc[which(data_qc$sample_label==sample[i]),]
  clones_sample<-data_sample_unique[,"V_J_lenghCDR3_Clone_igh"]
  clones_specimen_IGHA<-data_sample_unique[which(data_sample_unique$isotype=="IGHA"),"V_J_lenghCDR3_Clone_igh"]
  clones_specimen_IGHD<-data_sample_unique[which(data_sample_unique$isotype=="IGHD"),"V_J_lenghCDR3_Clone_igh"]
  clones_specimen_IGHE<-data_sample_unique[which(data_sample_unique$isotype=="IGHE"),"V_J_lenghCDR3_Clone_igh"]
  clones_specimen_IGHG<-data_sample_unique[which(data_sample_unique$isotype=="IGHG"),"V_J_lenghCDR3_Clone_igh"]
  clones_specimen_IGHM<-data_sample_unique[which(data_sample_unique$isotype=="IGHM"),"V_J_lenghCDR3_Clone_igh"]
  clones_specimen_unmutated_IGHM<-data_sample_unique[which(data_sample_unique$IGHM_mutated=="unmutated"),"V_J_lenghCDR3_Clone_igh"]
  clones_specimen_mutated_IGHM<-data_sample_unique[which(data_sample_unique$IGHM_mutated=="mutated"),"V_J_lenghCDR3_Clone_igh"]
  clones_specimen_unmutated_IGHD<-data_sample_unique[which(data_sample_unique$IGHD_mutated=="unmutated"),"V_J_lenghCDR3_Clone_igh"]
  clones_specimen_mutated_IGHD<-data_sample_unique[which(data_sample_unique$IGHD_mutated=="mutated"),"V_J_lenghCDR3_Clone_igh"]
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
  fi_unmutated_IGHM<-as.numeric(table(clones_specimen_unmutated_IGHM))/length(clones_specimen_unmutated_IGHM)
  fi_mutated_IGHM<-as.numeric(table(clones_specimen_mutated_IGHM))/length(clones_specimen_mutated_IGHM)
  fi_unmutated_IGHD<-as.numeric(table(clones_specimen_unmutated_IGHD))/length(clones_specimen_unmutated_IGHD)
  fi_mutated_IGHD<-as.numeric(table(clones_specimen_mutated_IGHD))/length(clones_specimen_mutated_IGHD)
  
  hi_IGHA<-fi_IGHA*log2(fi_IGHA)
  hi_IGHD<-fi_IGHD*log2(fi_IGHD)
  hi_IGHE<-fi_IGHE*log2(fi_IGHE)
  hi_IGHG<-fi_IGHG*log2(fi_IGHG)
  hi_IGHM<-fi_IGHM*log2(fi_IGHM)
  hi_unmutated_IGHM<-fi_unmutated_IGHM*log2(fi_unmutated_IGHM)
  hi_mutated_IGHM<-fi_mutated_IGHM*log2(fi_mutated_IGHM)
  hi_unmutated_IGHD<-fi_unmutated_IGHD*log2(fi_unmutated_IGHD)
  hi_mutated_IGHD<-fi_mutated_IGHD*log2(fi_mutated_IGHD)
  
  entropy_IGHA[i]=-sum(hi_IGHA)
  entropy_IGHD[i]=-sum(hi_IGHD)
  entropy_IGHE[i]=-sum(hi_IGHE)
  entropy_IGHG[i]=-sum(hi_IGHG)
  entropy_IGHM[i]=-sum(hi_IGHM)
  entropy_unmutated_IGHM[i]=-sum(hi_unmutated_IGHM)
  entropy_mutated_IGHM[i]=-sum(hi_mutated_IGHM)
  entropy_unmutated_IGHD[i]=-sum(hi_unmutated_IGHD)
  entropy_mutated_IGHD[i]=-sum(hi_mutated_IGHD)  
}

# entropy_norm_IGHA<-entropy_IGHA/max(entropy_IGHA,na.rm = T)
# clonality_IGHA<-(1-entropy_norm_IGHA)
# entropy_norm_IGHD<-entropy_IGHD/max(entropy_IGHD,na.rm = T)
# clonality_IGHD<-(1-entropy_norm_IGHD)
# entropy_norm_IGHG<-entropy_IGHG/max(entropy_IGHG,na.rm = T)
# clonality_IGHG<-(1-entropy_norm_IGHG)
# entropy_norm_IGHM<-entropy_IGHM/max(entropy_IGHM,na.rm = T)
# clonality_IGHM<-(1-entropy_norm_IGHM)
# entropy_norm_memory<-entropy_memory/max(entropy_memory,na.rm = T)
# clonality_memory<-(1-entropy_norm_memory)
# entropy_norm_naive<-entropy_naive/max(entropy_naive,na.rm = T)
# clonality_naive<-(1-entropy_norm_naive)
# 
diversity<-cbind(entropy_IGHA,entropy_IGHD,entropy_IGHE,entropy_IGHG,entropy_IGHM,entropy_unmutated_IGHM,entropy_mutated_IGHM,
                 entropy_unmutated_IGHD,entropy_mutated_IGHD)
rownames(diversity)<-sample
summary_data<-cbind(summary_data,diversity)

save(data_qc,summary_data,file="Data/BCR_PTB_Term_data_summary.RData")


