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

setwd("/Users/Pinedasans/ImmuneRep_Pregnancy/")
load("Data/BCR_data_summary.RData")

##1.Obtain the vertex and edges
Obtain_vertex_edges<-function(data,isotype){
  data<-data[which(data$isotype==isotype),]
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

##2.Apply the nucleotides-assembly-1.0.jar made by Mikel using the Network.sh 

##3.Plot the network
sample<-unique(data_qc$sample_label)

isotype="IGHM"

##Mother
sample_M<-sample[which(summary_data$sample=="Mother")] 
for(i in sample_M) {
  print(i)
    edges <- read.delim(paste("Results/edges_",isotype,"_",i,".txt.outcome.txt",sep = ""))
    vertex <- read.delim(paste("Results/vertex_",isotype,"_",i,".txt",sep = ""))
    net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
    V(net)$size <- V(net)$Freq/100
    V(net)$color <- c("#BEAED4")
    net <- simplify(net, remove.multiple = F, remove.loops = T) 
    E(net)$arrow.mode <- 0
    E(net)$width <- 0.4
    E(net)$color <- c("black")
    tiff(paste("Results/network_",isotype,"_",i,".tiff",sep=""),res=300,h=3000,w=3000)
    plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
    dev.off()
}

##Fetus
sample_F<-sample[which(summary_data$sample=="Fetus")] 
for(i in sample_F) {
  print(i)
  edges <- read.delim(paste("Results/edges_",isotype,"_",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("Results/vertex_",isotype,"_",i,".txt",sep = ""))
  net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
  V(net)$size <- V(net)$Freq/100
  V(net)$color <- c("#7FC97F")
  net <- simplify(net, remove.multiple = F, remove.loops = T) 
  E(net)$arrow.mode <- 0
  E(net)$width <- 0.4
  E(net)$color <- c("black")
  tiff(paste("Results/network_",isotype,"_",i,".tiff",sep=""),res=300,h=3000,w=3000)
  plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
  dev.off()
}
