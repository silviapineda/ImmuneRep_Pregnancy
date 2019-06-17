###########################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies for pregnancy outcomes
###
### CITATION: Based off of Silvia Pineda's original network analysis code
###
### PROCESS: 
###           
### DESCRIP: iRepertoire data - network analysis of data by type
###          for downsampled data (see ClonesInferenceBrianBCRDownsampling.py)
###
### Author: Brian Le
### Date: April, 2019
############################################################################################
library(ggplot2)
library("igraph")
library("ineq")
library("dplyr")
library("RColorBrewer")

directory <- "iRepertoire/network_analysis/downsampled/"

summary_file_path <- 'iRepertoire/iRepertoire_summary_data.csv'
summary_data <- read.csv(summary_file_path)

#load in downsampled BCR data
data_BCR_clones <- read.csv("iRepertoire_clones_BCR_processed_downsampled.csv")
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
Obtain_vertex_edges<-function(data,isotype,directorypath){
  
  dir.create(directorypath, showWarnings = FALSE)
  
  data<-data[which(data$isotype==isotype),]
  
  data$CloneID_CDR3pep<-paste0(data[,c("V_J_lengthCDR3_clone")],data[,c("CDR3.pep.")])
  sample<-unique(data$sample)
  
  print(sample)
  for (sample_number in sample){
    print(sample_number)
    data_sample<-data[which(data$sample==sample_number),c("CloneID_CDR3pep","V_J_lengthCDR3_clone")]
    #data_sample<-data[which(data$sample==sample_number),c("CloneID_CDR3pep","V_J_lengthCDR3_clone","copy")]
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
    df_vertex<-data.frame(table(data_sample$CloneID_CDR3pep)) #frequency of each vertex... need to include copy here
    #df_vertex<-aggregate(copy ~ CloneID_CDR3pep, data_sample, sum)
    #assign(paste0("vertex",i),df_vertex[which(df_vertex$Freq!=0),])
    #df_vertex<-df_vertex[which(df_vertex$copy!=0),]
    df_vertex<-df_vertex[which(df_vertex$Freq!=0),]
    
    print("vertex finish; writing results")
    #write.csv(edge_all,paste0(directorypath,"edges_",isotype,"_",sample_number,".csv"),row.names = F)
    #write.csv(df_vertex,paste0(directorypath,"vertex_",isotype,"_",sample_number,".csv"),row.names = F)
  }
}

#isotypes <- unique(data_BCR_clones$isotype)
isotypes <- c("M", "G1", "G3", "G4", "D", "E", "A2")

for (isotype in isotypes){
  print(isotype)
  Obtain_vertex_edges(data_BCR_clones, isotype, directorypath = paste0(directory, "edges_and_vertices/"))
}

##2. Obtain the vertex and cluster Gini index
Obtain_gini_index<-function(isotype,samples){
  
  #samples are numerical sample numbers (e.g. Sample.id from summary_data)
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
    if (file.exists(paste0(directory, "edges_and_vertices/","edges_",isotype,"_",sample_number,".csv"))){
      # assign(paste0("edges",i),read.delim(paste0("Results/Network/",isotype,"/edges_",isotype2,"_",i,".txt")))
      # assign(paste0("vertex",i),read.delim(paste0("Results/Network/",isotype,"/vertex_",isotype2,"_",i,".txt")))
      edges <- read.csv(paste0(directory, "edges_and_vertices/","edges_",isotype,"_",sample_number,".csv"))
      
      #skip sample is edges file exists but is empty
      if (nrow(edges) == 0){
        next
      }
      
      vertex <- read.csv(paste0(directory, "edges_and_vertices/","vertex_",isotype,"_",sample_number,".csv"))
      vertex_max[j]<-max(vertex$Freq)
      vertex_gini[j]<-Gini(vertex$Freq)
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


#Gini calculation/plotting for BCR data, split into decidua and maternal_blood
isotypes <- c("M", "G1", "G3", "G4", "D", "E", "A2")

for (sample_type in sample_types){
  print(sample_type)
  
  samples_term <- summary_data[which(summary_data$sample_type == sample_type &
                                         summary_data$outcome == "normal"), "Sample.id"]
  samples_preterm <- summary_data[which(summary_data$sample_type == sample_type &
                                         summary_data$outcome == "ptl"), "Sample.id"]
  
  for (isotype in isotypes){
    
    print(isotype)
    ###Plot Gini boxplot
    cluster_gini_term <- data.frame(Obtain_gini_index(isotype, samples_term))
    cluster_gini_preterm <- data.frame(Obtain_gini_index(isotype, samples_preterm))
    cluster_gini_term$outcome <- "normal"
    cluster_gini_preterm$outcome <- "ptl"
    summary_data_gini <- rbind(cluster_gini_term, cluster_gini_preterm)
    summary_data_gini$outcome <- factor(summary_data_gini$outcome)
    summary_data_gini$sample <- rownames(summary_data_gini)
    
    #brewer.pal(n = 3, name = "Accent")[2:1]
    palette(brewer.pal(n = 3, name = "Accent")[2:1]) #set palette colors for plotting by factor
    
    tiff(paste0(directory, sample_type, "/network_vertex_cluster_gini_",isotype,".tiff"),h=2000,w=2000,res=300)
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

##3.Plot the network
#takes a long time to run!!!
#separate into groups: decidua and maternal_blood, term and preterm

sample_types <- unique(summary_data$sample_type)
outcomes <- unique(summary_data$outcome)

#isotypes <- c("M", "D", "E", "A2")
#isotypes <- c("D", "E", "A2")
isotypes <- c("M", "G1", "G3", "G4", "D", "E", "A2")

for (sample_type in sample_types){
  print(sample_type)
  
  for (outcome in outcomes){
    save_path <- paste0(directory, "network_analysis/", sample_type, "/", outcome, "/")
    dir.create(save_path)
    
    sample_numbers <- summary_data[which(summary_data$sample_type == sample_type &
                                           summary_data$outcome == outcome), "Sample.id"]
    data_sample <- data_BCR_clones[data_BCR_clones$sample %in% sample_numbers,]    
    
    for (isotype in isotypes){
      print(isotype)
      for(sample_number in sample_numbers) {
        if (file.exists(paste0(directory, "edges_and_vertices/","edges_",isotype,"_",sample_number,".csv"))){
          print(sample_number)
          edges <- read.csv(paste0(directory, "edges_and_vertices/","edges_",isotype,"_",sample_number,".csv"))
          vertex <- read.csv(paste0(directory, "edges_and_vertices/","vertex_",isotype,"_",sample_number,".csv"))
          net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
          V(net)$size <- V(net)$Freq/100
          V(net)$color <- c("#BEAED4")
          net <- simplify(net, remove.multiple = F, remove.loops = T)
          E(net)$arrow.mode <- 0
          E(net)$width <- 0.4
          E(net)$color <- c("black")
          tiff(paste(save_path,isotype,"_",sample_number,".tiff",sep=""),res=300,h=3000,w=3000)
          plot(net,vertex.label=NA,layout=layout_with_graphopt(net,niter=800))
          dev.off()
        }
      }
    }
  }
}