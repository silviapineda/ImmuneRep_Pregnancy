rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies for pregnancy outcomes
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: Network analysis
###         
###
### Author: Silvia Pineda
### Date: January, 2019
############################################################################################
library(ggplot2)
library("igraph")
library("ineq")
library("dplyr")
library("RColorBrewer")

setwd("/Users/Pinedasans/ImmuneRep_Pregnancy/")
load("Data/BCR/BCR_data_summary.RData")

##1.Obtain the vertex and edges
Obtain_vertex_edges<-function(data,isotype){
 
  if(isotype=="memory" | isotype=="naive"){
    data<-data[which(data$IGHM_naive_memory==isotype),]
  } else{
    data<-data[which(data$isotype==isotype),]
  }
  
  data$CloneId_CDR3<-paste0(data[,c("V_J_lenghCDR3_Clone_igh")],data[,c("cdr3_seq_aa_q")])
  sample<-unique(data$sample_label)
  for (i in sample){
      print(i)
      data_sample<-data[which(data$sample_label==i),]
      df_sample<-data_sample[,c("CloneId_CDR3","V_J_lenghCDR3_Clone_igh")]
      groups <- group_by(df_sample,CloneId_CDR3,V_J_lenghCDR3_Clone_igh)
      assign(paste0("edges",i),unique(data.frame(groups)))
      df_vertex<-data.frame(table(data_sample$CloneId_CDR3))
      assign(paste0("vertex",i),df_vertex[which(df_vertex$Freq!=0),])
      write.table(get(paste0("edges",i)),paste0("Results/edges_",isotype,"_",i,".txt"),sep="\t",row.names = F)
      write.table(get(paste0("vertex",i)),paste0("Results/vertex_",isotype,"_",i,".txt"),sep="\t",row.names = F)
  }
}

network_IGHM<-Obtain_vertex_edges(data_qc,"IGHM")
network_IGHD<-Obtain_vertex_edges(data_qc,"IGHD")
network_IGHG<-Obtain_vertex_edges(data_qc,"IGHG")
network_IGHA<-Obtain_vertex_edges(data_qc,"IGHA")
network_memory<-Obtain_vertex_edges(data_qc,"memory")
network_naive<-Obtain_vertex_edges(data_qc,"naive")



##2.Apply the nucleotides-assembly-1.0.jar made by Mikel using the Network.sh 

##3.Plot the network
sample<-unique(data_qc$sample_label)

isotype="naive"

##Mother
sample_M<-sample[which(summary_data$sample=="Mother")] 
for(i in sample_M) {
  print(i)
    edges <- read.delim(paste("Results/Network/",isotype,"/edges_",isotype,"_",i,".txt.outcome.txt",sep = ""))
    vertex <- read.delim(paste("Results/Network/",isotype,"/vertex_",isotype,"_",i,".txt",sep = ""))
    net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
    V(net)$size <- V(net)$Freq/100
    V(net)$color <- c("#BEAED4")
    net <- simplify(net, remove.multiple = F, remove.loops = T) 
    E(net)$arrow.mode <- 0
    E(net)$width <- 0.4
    E(net)$color <- c("black")
    tiff(paste("Results/Network/",isotype,"/network_",isotype,"_",i,".tiff",sep=""),res=300,h=3000,w=3000)
    plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
    dev.off()
}

##Fetus
sample_F<-sample[which(summary_data$sample=="Fetus")] 
for(i in sample_F) {
  print(i)
  edges <- read.delim(paste("Results/Network/",isotype,"/edges_",isotype,"_",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("Results/Network/",isotype,"/vertex_",isotype,"_",i,".txt",sep = ""))
  net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
  V(net)$size <- V(net)$Freq/100
  V(net)$color <- c("#7FC97F")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  E(net)$width <- 0.4
  E(net)$color <- c("black")
  tiff(paste("Results/Network/",isotype,"/network_",isotype,"_",i,".tiff",sep=""),res=300,h=3000,w=3000)
  plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
  dev.off()
}


##4. Obtain the vertex and cluster Gini index
Obtain_gini_index<-function(data,isotype,summary_data){
  if(isotype=="memory" | isotype=="naive"){
    data<-data[which(data$IGHM_naive_memory==isotype),]
  } else{
    data<-data[which(data$isotype==isotype),]
  }
  
  data$CloneId_CDR3<-paste0(data[,c("V_J_lenghCDR3_Clone_igh")],data[,c("cdr3_seq_aa_q")])
  sample<-unique(data$sample_label)
  
  vertex_max<-NULL
  vertex_gini<-NULL
  cluster_max<-NULL
  cluster_gini<-NULL
  num_reads_max_cluster<-NULL
  clusters<-NULL
  j<-1
  for (i in sample){
    assign(paste0("edges",i),read.delim(paste0("Results/Network/",isotype,"/edges_",isotype,"_",i,".txt")))
    assign(paste0("vertex",i),read.delim(paste0("Results/Network/",isotype,"/vertex_",isotype,"_",i,".txt")))
    vertex_max[j]<-max(get(paste0("vertex",i))$Freq)
    vertex_gini[j]<-Gini(get(paste0("vertex",i))$Freq)
    cluster_max[j]<-max(table(get(paste0("edges",i))$V_J_lenghCDR3_Clone_igh))
    clusters[j]<-sum(table(table(get(paste0("edges",i))$V_J_lenghCDR3_Clone_igh)))
    num_reads_max_cluster[j]<-tail(table(table(get(paste0("edges",i))$V_J_lenghCDR3_Clone_igh)),1)
    cluster_gini[j]<-Gini(table(get(paste0("edges",i))$V_J_lenghCDR3_Clone_igh))
    j=j+1
  }
  
  #clonal_expansion<-(num_reads_max_cluster/summary_data[,isotype])*100
  results<-cbind(cluster_gini,vertex_gini,vertex_max,cluster_max,num_reads_max_cluster,clusters)
  
  rownames(results)<-sample
  
  return(results)
}

###Plot Gini boxplot
isotype = "naive"
assign(paste0("cluster_gini_",isotype),data.frame(Obtain_gini_index(data_qc,isotype,summary_data)))

brewer.pal(n = 3, name = "Accent")
COLOR=c("#BEAED4","#7FC97F")

tiff(paste0("Results/network_vertex_cluster_gini_",isotype,".tiff"),h=2000,w=2000,res=300)
par(fig=c(0,0.8,0,0.8))
plot(get(paste0("cluster_gini_",isotype))[,"cluster_gini"], 
         get(paste0("cluster_gini_",isotype))[,"vertex_gini"],
     col = COLOR,pch=20,ylab = "Gini (Vextex)",xlab = "Gini (Cluster)")
legend("bottomright",legend=c("Maternal","Fetal"), 
       col=COLOR,pch=20,cex=c(1.2),ncol=2)

par(fig=c(0,0.8,0.55,1), new=TRUE)
summary(lm(get(paste0("cluster_gini_",isotype))[,"cluster_gini"]~summary_data$sample))
boxplot(get(paste0("cluster_gini_",isotype))[,"cluster_gini"]~summary_data$sample,
        col=COLOR, horizontal=TRUE, axes=FALSE)

par(fig=c(0.65,1,0,0.8),new=TRUE)
summary(lm( get(paste0("cluster_gini_",isotype))[,"vertex_gini"]~summary_data$sample))
boxplot(get(paste0("cluster_gini_",isotype))[,"vertex_gini"]~summary_data$sample,
        col=COLOR,axes=FALSE)
dev.off()

####################################
#### Network Analysis for TCR ######
###################################
load("Data/TCR/TCR_data_summary.RData")

##1.Obtain the vertex and edges
data_qc$CloneId_CDR3<-paste0(data_qc[,c("V_J_length_CDR3_Clone_tcrb")],data_qc[,c("cdr3_seq_aa_q")])
sample<-unique(data_qc$sample_label)
for (i in sample){
  print(i)
  data_sample<-data_qc[which(data_qc$sample_label==i),]
  df_sample<-data_sample[,c("CloneId_CDR3","V_J_length_CDR3_Clone_tcrb")]
  groups <- group_by(df_sample,CloneId_CDR3,V_J_length_CDR3_Clone_tcrb)
  assign(paste0("edges",i),unique(data.frame(groups)))
  df_vertex<-data.frame(table(data_sample$CloneId_CDR3))
  assign(paste0("vertex",i),df_vertex[which(df_vertex$Freq!=0),])
  write.table(get(paste0("edges",i)),paste0("Results/edges_",i,".txt"),sep="\t",row.names = F)
  write.table(get(paste0("vertex",i)),paste0("Results/vertex_",i,".txt"),sep="\t",row.names = F)
}


##2.Apply the nucleotides-assembly-1.0.jar made by Mikel using the Network.sh 

##3.Plot the network
sample<-unique(data_qc$sample_label)

##Mother
sample_M<-sample[which(summary_data$sample=="Mother")] 
for(i in sample_M) {
  print(i)
  edges <- read.delim(paste("Results/edges_",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("Results/vertex_",i,".txt",sep = ""))
  net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
  V(net)$size <- V(net)$Freq/100
  V(net)$color <- c("#BEAED4")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  E(net)$width <- 0.4
  E(net)$color <- c("black")
  tiff(paste("Results/network_",i,".tiff",sep=""),res=300,h=3000,w=3000)
  plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
  dev.off()
}

##Fetus
sample_F<-sample[which(summary_data$sample=="Fetus")] 
for(i in sample_F) {
  print(i)
  edges <- read.delim(paste("Results/edges_",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("Results/vertex_",i,".txt",sep = ""))
  net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
  V(net)$size <- V(net)$Freq/100
  V(net)$color <- c("#7FC97F")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  E(net)$width <- 0.4
  E(net)$color <- c("black")
  tiff(paste("Results/network_",i,".tiff",sep=""),res=300,h=3000,w=3000)
  plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
  dev.off()
}

###Obtain Gini Cluster and Vertex
data_qc$CloneId_CDR3<-paste0(data_qc[,c("V_J_length_CDR3_Clone_tcrb")],data_qc[,c("cdr3_seq_aa_q")])
sample<-unique(data_qc$sample_label)

vertex_max<-NULL
vertex_gini<-NULL
cluster_max<-NULL
cluster_gini<-NULL
num_reads_max_cluster<-NULL
clusters<-NULL
j<-1
for (i in sample){
  assign(paste0("edges",i),read.delim(paste0("Results/Network/TCR/edges_",i,".txt")))
  assign(paste0("vertex",i),read.delim(paste0("Results/Network/TCR/vertex_",i,".txt")))
  vertex_max[j]<-max(get(paste0("vertex",i))$Freq)
  vertex_gini[j]<-Gini(get(paste0("vertex",i))$Freq)
  cluster_max[j]<-max(table(get(paste0("edges",i))$V_J_length_CDR3_Clone_tcrb))
  clusters[j]<-sum(table(table(get(paste0("edges",i))$V_J_length_CDR3_Clone_tcrb)))
  num_reads_max_cluster[j]<-tail(table(table(get(paste0("edges",i))$V_J_length_CDR3_Clone_tcrb)),1)
  cluster_gini[j]<-Gini(table(get(paste0("edges",i))$V_J_length_CDR3_Clone_tcrb))
  j=j+1
}
clonal_expansion<-(num_reads_max_cluster/summary_data[,"read_count"])*100
cluster_gini_TCR<-data.frame(cbind(cluster_gini,vertex_gini,vertex_max,cluster_max,num_reads_max_cluster,clusters,clonal_expansion))

rownames(cluster_gini_TCR)<-sample

###Plot Gini boxplot
brewer.pal(n = 3, name = "Accent")
COLOR=c("#BEAED4","#7FC97F")

tiff("Results/network_vertex_cluster_gini_TCR.tiff",h=2000,w=2000,res=300)
par(fig=c(0,0.8,0,0.8))
plot(cluster_gini_TCR[,"cluster_gini"], 
     cluster_gini_TCR[,"vertex_gini"],
     col = COLOR,pch=20,ylab = "Gini (Vextex)",xlab = "Gini (Cluster)")
legend("bottomright",legend=c("Maternal","Fetal"), 
       col=COLOR,pch=20,cex=c(1.2),ncol=2)

par(fig=c(0,0.8,0.55,1), new=TRUE)
summary(lm(cluster_gini_TCR[,"cluster_gini"]~summary_data$sample))
boxplot(cluster_gini_TCR[,"cluster_gini"]~summary_data$sample,
        col=COLOR, horizontal=TRUE, axes=FALSE)

par(fig=c(0.65,1,0,0.8),new=TRUE)
summary(lm( cluster_gini_TCR[,"vertex_gini"]~summary_data$sample))
boxplot(cluster_gini_TCR[,"vertex_gini"]~summary_data$sample,
        col=COLOR,axes=FALSE)
dev.off()

