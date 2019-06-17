rm(list = ls(all = TRUE))
x <-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies for pregnancy outcomes
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: Network analysis of all (not downsampled) BCR iRepertoire data
###         
###
### Author: Silvia Pineda / Brian Le
### Date: January, 2019
############################################################################################
library(ggplot2)
library("igraph")
library("ineq")
library("dplyr")
library("RColorBrewer")

#setwd("/Users/Pinedasans/ImmuneRep_Pregnancy/")
directory <- "iRepertoire/"
load(file=paste(directory, 'iRepertoire_data.RData', sep = ""))

summary_file_path <- paste(directory, 'iRepertoire_summary_data.csv', sep = "")
summary_data <- read.csv(summary_file_path)

#this data is the downsampled data that assigned clone numbers to the CD19 data only
#this is also NO DOWNSAMPLING
data_BCR_clones <- read.csv("iRepertoire_clones_BCR_processed.csv")
data_BCR_clones$V_J_lengthCDR3_clone <- paste(data_BCR_clones$V_J_lengthCDR3, data_BCR_clones$numberClone)

get_isotype <- function(class){
  
  class <- gsub("hTRB", "", class)
  class <- gsub("hIGH", "", class)
  class <- gsub("\\*.*", "", class)
  
  return(class)
}

#this leaves G1, G3, and G4 as separate isotype subclasses
data_BCR_clones$isotype <- with(data_BCR_clones, get_isotype(C))

##1.Obtain the vertex and edges
#Obtain_vertex_edges<-function(data,isotype){
# 
#  data<-data[which(data$isotype==isotype),]
#  
#  data$CloneID_CDR3pep<-paste0(data[,c("V_J_lengthCDR3_clone")],data[,c("CDR3.pep.")])
#  sample<-unique(data$sample_label)
#  for (i in sample){
#      print(i)
#      data_sample<-data[which(data$sample_label==i),]
#      df_sample<-data_sample[,c("CloneID_CDR3pep","V_J_lengthCDR3_clone")]
#      groups <- group_by(df_sample,CloneID_CDR3pep,V_J_lengthCDR3_clone)
#      assign(paste0("edges",i),unique(data.frame(groups))) #unique to draw out edges between unique pairs
#      #df_vertex<-data.frame(table(data_sample$CloneID_CDR3pep)) #frequency of each vertex... need to include copy here
#      df_vertex<-aggregate(copy ~ CloneID_CDR3pep, data_sample, sum)
#      assign(paste0("vertex",i),df_vertex[which(df_vertex$Freq!=0),])
#      write.table(get(paste0("edges",i)),paste0("Results/edges_",isotype,"_",i,".txt"),sep="\t",row.names = F)
#      write.table(get(paste0("vertex",i)),paste0("Results/vertex_",isotype,"_",i,".txt"),sep="\t",row.names = F)
#  }
#}

#trying to update it to combine steps #1 and #2
Obtain_vertex_edges<-function(data,isotype,directorypath){
  
  dir.create(directorypath, showWarnings = FALSE)
  
  data<-data[which(data$isotype==isotype),]
  
  data$CloneID_CDR3pep<-paste0(data[,c("V_J_lengthCDR3_clone")],data[,c("CDR3.pep.")])
  sample<-unique(data$sample)
  
  print(sample)
  for (sample_number in sample){
    print(sample_number)
    data_sample<-data[which(data$sample==sample_number),c("CloneID_CDR3pep","V_J_lengthCDR3_clone","copy")]
    #df_sample<-data_sample[,c("CloneID_CDR3pep","V_J_lengthCDR3_clone")]
    #groups <- group_by(df_sample,CloneID_CDR3pep,V_J_lengthCDR3_clone)
    data_edges <- unique(data_sample[,c("CloneID_CDR3pep","V_J_lengthCDR3_clone")])
    #unique to draw out edges between unique pairs
    #assign(paste0("edges",i),data_edges)
    
    print("edges start")
    
    group_V_J_CDR3 <- unique(data_edges$V_J_lengthCDR3_clone)
    edge_all <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("vertex_A", "vertex_B", "clone"))
    
    for (edge_connector in group_V_J_CDR3){
      
      edge_group <- data_edges[data_edges$V_J_lengthCDR3_clone == edge_connector, "CloneID_CDR3pep"]
      edge_group_length <- length(edge_group)
      if (edge_group_length > 1){
        edge_list_to <- c()
        edge_list_from <- c()
        for (i in 1:edge_group_length){
          for (j in 1:i){
            if (i != j){
              edge_list_to <- c(edge_list_to, edge_group[i])
              edge_list_from <- c(edge_list_from, edge_group[j])
            }
          }
        }
        #test <- data.frame(vertex_A = edge_list_to, vertex_B = edge_list_from, clone = edge_connector)
        edge_all <- rbind(edge_all, data.frame(vertex_A = edge_list_to, vertex_B = edge_list_from, clone = edge_connector))
        
      }
    }
    print("edges finish")
    
    #assign(paste0("edges",i),edge_all)
    
    print("vertex start")
    #df_vertex<-data.frame(table(data_sample$CloneID_CDR3pep)) #frequency of each vertex... need to include copy here
    df_vertex<-aggregate(copy ~ CloneID_CDR3pep, data_sample, sum)
    #assign(paste0("vertex",i),df_vertex[which(df_vertex$Freq!=0),])
    df_vertex<-df_vertex[which(df_vertex$copy!=0),]
    
    print("vertex finish; writing results")
    write.csv(edge_all,paste0(directorypath,"edges_",isotype,"_",sample_number,".csv"),row.names = F)
    write.csv(df_vertex,paste0(directorypath,"vertex_",isotype,"_",sample_number,".csv"),row.names = F)
  }
}

#####

isotypes <- unique(data_BCR_clones$isotype)
isotypes <- c("M", "G1", "G3", "G4", "D", "E", "A2")

for (isotype in isotypes){
  print(isotype)
  Obtain_vertex_edges(data_BCR_clones, isotype, directorypath = paste0(directory, "network_analysis/"))
}

#network_IGHM<-Obtain_vertex_edges(data_qc,"IGHM",type_mutated = 0)
#network_IGHD<-Obtain_vertex_edges(data_qc,"IGHD")
#network_IGHG<-Obtain_vertex_edges(data_qc,"IGHG")
#network_IGHA<-Obtain_vertex_edges(data_qc,"IGHA")
#network_IGHD_mutated<-Obtain_vertex_edges(data_qc,"IGHD","mutated")
#network_IGHD_unmutated<-Obtain_vertex_edges(data_qc,"IGHD","unmutated")
#network_IGHM_mutated<-Obtain_vertex_edges(data_qc,"IGHM","mutated")
#network_IGHM_unmutated<-Obtain_vertex_edges(data_qc,"IGHM","unmutated")

##2.Apply the nucleotides-assembly-1.0.jar made by Mikel using the Network.sh 
#need to write something here that will connect all the vertices by their shared property (edge)

#testing with edgesM0974 generated through step 1 analysis with BCR

#since the only thing linking the clones (vertices) together is their base VJCDR3, it necessarily follows that the groups
#will be disjoint, and the groups will only be connected to every other clone in the same group...
#so instead, let's iterate through each unique VJCDR3, and then for the edges list, populate it with each corresponding pair
#within that set of unique VJCDR3

# group_V_J_CDR3 <- unique(edgesM0974$V_J_lengthCDR3_clone)
# edge_all <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("vertex_A", "vertex_B", "clone"))
# 
# for (edge_connector in group_V_J_CDR3){
# 
#   edge_group <- edgesM0974[edgesM0974$V_J_lengthCDR3_clone == edge_connector, "CloneID_CDR3pep"]
#   edge_group_length <- length(edge_group)
#   if (edge_group_length > 1){
#     edge_list_to <- c()
#     edge_list_from <- c()
#     for (i in 1:edge_group_length){
#       for (j in 1:i){
#         if (i != j){
#           edge_list_to <- c(edge_list_to, edge_group[i])
#           edge_list_from <- c(edge_list_from, edge_group[j])
#         }
#       }
#     }
#     #test <- data.frame(vertex_A = edge_list_to, vertex_B = edge_list_from, clone = edge_connector)
#     edge_all <- rbind(edge_all, data.frame(vertex_A = edge_list_to, vertex_B = edge_list_from, clone = edge_connector))
#   }
# }
# 
# #qa check
# nrow(edge_all[edge_all$clone == "IGHV3-9_IGHJ4_39_1493",]) #should be 3003, not 6006

##3.Plot the network


##term and preterm
samples_BCR_term <- summary_data[which(summary_data$sample_population_type == "CD19" &
                                   summary_data$outcome == "normal"), "Sample.id"]
samples_BCR_preterm <- summary_data[which(summary_data$sample_population_type == "CD19" &
                                      summary_data$outcome == "ptl"), "Sample.id"]


#isotypes <- c("M", "D", "E", "A2")
#isotypes <- c("D", "E", "A2")

#term
dir.create(paste0(directory, "network_analysis/preterm/"))
dir.create(paste0(directory, "network_analysis/term/"))
for (isotype in isotypes){
  print(isotype)
    for(sample_number in samples_BCR_term) {
      if (file.exists(paste0(directory, "network_analysis/","edges_",isotype,"_",sample_number,".csv"))){
        print(sample_number)
        edges <- read.csv(paste0(directory, "network_analysis/","edges_",isotype,"_",sample_number,".csv"))
        vertex <- read.csv(paste0(directory, "network_analysis/","vertex_",isotype,"_",sample_number,".csv"))
        net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
        V(net)$size <- V(net)$copy/100
        V(net)$color <- c("#BEAED4")
        net <- simplify(net, remove.multiple = F, remove.loops = T)
        E(net)$arrow.mode <- 0
        E(net)$width <- 0.4
        E(net)$color <- c("black")
        tiff(paste(directory,"network_analysis/term/",isotype,"_",sample_number,".tiff",sep=""),res=300,h=3000,w=3000)
        plot(net,vertex.label=NA,layout=layout_with_graphopt(net,niter=800))
        dev.off()
    }
  }

  for(sample_number in samples_BCR_preterm) {
    if (file.exists(paste0(directory, "network_analysis/","edges_",isotype,"_",sample_number,".csv"))){
      print(sample_number)
      edges <- read.csv(paste0(directory, "network_analysis/","edges_",isotype,"_",sample_number,".csv"))
      vertex <- read.csv(paste0(directory, "network_analysis/","vertex_",isotype,"_",sample_number,".csv"))
      net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
      V(net)$size <- V(net)$copy/100
      V(net)$color <- c("#BEAED4")
      net <- simplify(net, remove.multiple = F, remove.loops = T) 
      E(net)$arrow.mode <- 0
      E(net)$width <- 0.4
      E(net)$color <- c("black")
      tiff(paste(directory,"network_analysis/preterm/",isotype,"_",sample_number,".tiff",sep=""),res=300,h=3000,w=3000)
      plot(net,vertex.label=NA,layout=layout_with_graphopt(net,niter=800))
      dev.off()
    }
  }
}

#brian-testing-area
# edge <- read.csv("iRepertoire/network_analysis/edges_A2_95313.csv")
# vertex <- read.csv("iRepertoire/network_analysis/vertex_A2_95313.csv")
# net<-graph_from_data_frame(d=edge[,c("vertex_A", "vertex_B")],vertices = vertex,directed=F)



##4. Obtain the vertex and cluster Gini index
Obtain_gini_index<-function(isotype,samples){
  
  #sample<-rownames(summary_data)
  
  vertex_max<-NULL
  vertex_gini<-NULL
  cluster_max<-NULL
  cluster_gini<-NULL
  num_reads_max_cluster<-NULL
  clusters<-NULL
  samples_processed<-NULL
  j<-1
  for (sample_number in samples){
    print(sample_number)
    if (file.exists(paste0(directory, "network_analysis/","edges_",isotype,"_",sample_number,".csv"))){
      # assign(paste0("edges",i),read.delim(paste0("Results/Network/",isotype,"/edges_",isotype2,"_",i,".txt")))
      # assign(paste0("vertex",i),read.delim(paste0("Results/Network/",isotype,"/vertex_",isotype2,"_",i,".txt")))
      edges <- read.csv(paste0(directory, "network_analysis/","edges_",isotype,"_",sample_number,".csv"))
      
      #skip sample is edges file exists but is empty
      if (nrow(edges) == 0){
        next
      }
      
      vertex <- read.csv(paste0(directory, "network_analysis/","vertex_",isotype,"_",sample_number,".csv"))
      vertex_max[j]<-max(vertex$copy)
      vertex_gini[j]<-Gini(vertex$copy)
      # cluster_max[j]<-max(table(get(paste0("edges",i))$V_J_lengthCDR3_clone))
      edges_table <- table(c(edges$vertex_A, edges$vertex_B))
      cluster_max[j]<-max(edges_table)
      # clusters[j]<-sum(table(table(get(paste0("edges",i))$V_J_lengthCDR3_clone)))
      clusters[j]<-sum(table(edges_table))
      num_reads_max_cluster[j]<-tail(table(edges_table),1)
      cluster_gini[j]<-Gini(table(edges_table))
      j=j+1
      samples_processed <- c(samples_processed, sample_number)
    }
  }
  
  #clonal_expansion<-(num_reads_max_cluster/summary_data[,isotype])*100
  results<-cbind(cluster_gini,vertex_gini,vertex_max,cluster_max,num_reads_max_cluster,clusters)
  
  rownames(results)<-samples_processed
  
  return(results)
}


#Gini calculation/plotting for BCR data
sample_types <- unique(summary_data$sample_type)

for (sample_type in sample_types){
  print(sample_type)
  
  dir.create(paste0("iRepertoire/network_analysis/", sample_type, "/"))
  
  samples_term <- summary_data[which(summary_data$sample_type == sample_type &
                                       summary_data$outcome == "normal" &
                                       summary_data$sample_population_type == "CD19"), "Sample.id"]
  samples_preterm <- summary_data[which(summary_data$sample_type == sample_type &
                                          summary_data$outcome == "ptl" &
                                          summary_data$sample_population_type == "CD19"), "Sample.id"]
  
  for (isotype in isotypes){
  
    print(isotype)
    ###Plot Gini boxplot
    #isotype = "G4"
    cluster_gini_term <- data.frame(Obtain_gini_index(isotype, samples_term))
    cluster_gini_preterm <- data.frame(Obtain_gini_index(isotype, samples_preterm))
    cluster_gini_term$outcome <- "normal"
    cluster_gini_preterm$outcome <- "ptl"
    summary_data_gini <- rbind(cluster_gini_term, cluster_gini_preterm)
    summary_data_gini$outcome <- factor(summary_data_gini$outcome)
    summary_data_gini$sample <- rownames(summary_data_gini)
  
    #brewer.pal(n = 3, name = "Accent")[2:1]
    palette(brewer.pal(n = 3, name = "Accent")[2:1]) #set palette colors for plotting by factor
    
    tiff(paste0("iRepertoire/network_analysis/", sample_type, "/network_vertex_cluster_gini_",isotype,".tiff"),h=2000,w=2000,res=300)
    par(fig=c(0,0.8,0,0.8))
    plot(summary_data_gini[,"cluster_gini"], 
         summary_data_gini[,"vertex_gini"],
         col = summary_data_gini$outcome,pch=20,ylab = "Gini (Vertex)",xlab = "Gini (Cluster)")
    legend("bottomright",legend=c("ptl","term"), 
           col=palette()[1:2],pch=20,cex=c(1.2),ncol=2)
    
    par(fig=c(0,0.8,0.55,1), new=TRUE)
    summary(lm(summary_data_gini[,"cluster_gini"]~summary_data_gini$outcome))
    boxplot(summary_data_gini[,"cluster_gini"]~summary_data_gini$outcome,
            col=palette()[1:2], horizontal=TRUE, axes=FALSE)
    
    par(fig=c(0.65,1,0,0.8),new=TRUE)
    summary(lm(summary_data_gini[,"vertex_gini"]~summary_data_gini$outcome))
    boxplot(summary_data_gini[,"vertex_gini"]~summary_data_gini$outcome,
            col=palette()[1:2],axes=FALSE)
    dev.off()
  }
}
# ####################################
# #### Network Analysis for TCR ######
# ###################################
# load("Data/TCR/TCR_data_summary.RData")
# 
# ##1.Obtain the vertex and edges
# data_qc$CloneID_CDR3pep<-paste0(data_qc[,c("V_J_length_CDR3_Clone_tcrb")],data_qc[,c("cdr3_seq_aa_q")])
# sample<-unique(data_qc$sample_label)
# for (i in sample){
#   print(i)
#   data_sample<-data_qc[which(data_qc$sample_label==i),]
#   df_sample<-data_sample[,c("CloneID_CDR3pep","V_J_length_CDR3_Clone_tcrb")]
#   groups <- group_by(df_sample,CloneID_CDR3pep,V_J_length_CDR3_Clone_tcrb)
#   assign(paste0("edges",i),unique(data.frame(groups)))
#   df_vertex<-data.frame(table(data_sample$CloneID_CDR3pep))
#   assign(paste0("vertex",i),df_vertex[which(df_vertex$Freq!=0),])
#   write.table(get(paste0("edges",i)),paste0("Results/edges_",i,".txt"),sep="\t",row.names = F)
#   write.table(get(paste0("vertex",i)),paste0("Results/vertex_",i,".txt"),sep="\t",row.names = F)
# }
# 
# 
# ##2.Apply the nucleotides-assembly-1.0.jar made by Mikel using the Network.sh 
# 
# ##3.Plot the network
# sample<-unique(data_qc$sample_label)
# 
# ##Mother
# sample_M<-sample[which(summary_data$sample=="Mother")] 
# for(i in sample_M) {
#   print(i)
#   edges <- read.delim(paste("Results/edges_",i,".txt.outcome.txt",sep = ""))
#   vertex <- read.delim(paste("Results/vertex_",i,".txt",sep = ""))
#   net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
#   V(net)$size <- V(net)$Freq/100
#   V(net)$color <- c("#BEAED4")
#   net <- simplify(net, remove.multiple = F, remove.loops = T) 
#   E(net)$arrow.mode <- 0
#   E(net)$width <- 0.4
#   E(net)$color <- c("black")
#   tiff(paste("Results/network_",i,".tiff",sep=""),res=300,h=3000,w=3000)
#   plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
#   dev.off()
# }
# 
# ##Fetus
# sample_F<-sample[which(summary_data$sample=="Fetus")] 
# for(i in sample_F) {
#   print(i)
#   edges <- read.delim(paste("Results/edges_",i,".txt.outcome.txt",sep = ""))
#   vertex <- read.delim(paste("Results/vertex_",i,".txt",sep = ""))
#   net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
#   V(net)$size <- V(net)$Freq/100
#   V(net)$color <- c("#7FC97F")
#   net <- simplify(net, remove.multiple = F, remove.loops = T) 
#   E(net)$arrow.mode <- 0
#   E(net)$width <- 0.4
#   E(net)$color <- c("black")
#   tiff(paste("Results/network_",i,".tiff",sep=""),res=300,h=3000,w=3000)
#   plot(net,vertex.label=NA,layout=layout_with_graphopt(net))
#   dev.off()
# }
# 
# ###Obtain Gini Cluster and Vertex
# data_qc$CloneID_CDR3pep<-paste0(data_qc[,c("V_J_length_CDR3_Clone_tcrb")],data_qc[,c("cdr3_seq_aa_q")])
# sample<-unique(data_qc$sample_label)
# 
# vertex_max<-NULL
# vertex_gini<-NULL
# cluster_max<-NULL
# cluster_gini<-NULL
# num_reads_max_cluster<-NULL
# clusters<-NULL
# j<-1
# for (i in sample){
#   assign(paste0("edges",i),read.delim(paste0("Results/Network/TCR/edges_",i,".txt")))
#   assign(paste0("vertex",i),read.delim(paste0("Results/Network/TCR/vertex_",i,".txt")))
#   vertex_max[j]<-max(get(paste0("vertex",i))$Freq)
#   vertex_gini[j]<-Gini(get(paste0("vertex",i))$Freq)
#   cluster_max[j]<-max(table(get(paste0("edges",i))$V_J_length_CDR3_Clone_tcrb))
#   clusters[j]<-sum(table(table(get(paste0("edges",i))$V_J_length_CDR3_Clone_tcrb)))
#   num_reads_max_cluster[j]<-tail(table(table(get(paste0("edges",i))$V_J_length_CDR3_Clone_tcrb)),1)
#   cluster_gini[j]<-Gini(table(get(paste0("edges",i))$V_J_length_CDR3_Clone_tcrb))
#   j=j+1
# }
# clonal_expansion<-(num_reads_max_cluster/summary_data[,"read_count"])*100
# cluster_gini_TCR<-data.frame(cbind(cluster_gini,vertex_gini,vertex_max,cluster_max,num_reads_max_cluster,clusters,clonal_expansion))
# rownames(cluster_gini_TCR)<-sample
# summary_data_gini_TCR<-cbind(summary_data,cluster_gini_TCR)
# write.csv(summary_data_gini_TCR,file="Results/Network/summary_data_gini_TCR.csv")
# 
# 
# ###Plot Gini boxplot
# brewer.pal(n = 3, name = "Accent")
# COLOR=c("#BEAED4","#7FC97F")
# 
# tiff("Results/network_vertex_cluster_gini_TCR.tiff",h=2000,w=2000,res=300)
# par(fig=c(0,0.8,0,0.8))
# plot(cluster_gini_TCR[,"cluster_gini"], 
#      cluster_gini_TCR[,"vertex_gini"],
#      col = COLOR,pch=20,ylab = "Gini (Vextex)",xlab = "Gini (Cluster)")
# legend("bottomright",legend=c("Maternal","Fetal"), 
#        col=COLOR,pch=20,cex=c(1.2),ncol=2)
# 
# par(fig=c(0,0.8,0.55,1), new=TRUE)
# summary(lm(cluster_gini_TCR[,"cluster_gini"]~summary_data$sample))
# boxplot(cluster_gini_TCR[,"cluster_gini"]~summary_data$sample,
#         col=COLOR, horizontal=TRUE, axes=FALSE)
# 
# par(fig=c(0.65,1,0,0.8),new=TRUE)
# summary(lm( cluster_gini_TCR[,"vertex_gini"]~summary_data$sample))
# boxplot(cluster_gini_TCR[,"vertex_gini"]~summary_data$sample,
#         col=COLOR,axes=FALSE)
# dev.off()

