rm(list=ls())
setwd("~/Data/multilayerNetwork/data/")

# Select only top-scoring edges - inspect distribution first
allEdges <- read.delim("STRING_11.0_human_cleaned_edges.act_inh_directional.txt",header=F)
# edgeScores <- as.numeric(as.character(allEdges[,3]))
# hist(edgeScores,100)
# edgeScores.sorted <- sort(abs(edgeScores),index.return=T)
# # topNegativeEdges <- edgeScores.sorted$ix[1:round(propToKeep*nrow(allEdges))]
# # topPositiveEdges <- edgeScores.sorted$ix[(nrow(allEdges)-round(propToKeep*nrow(allEdges))+1):nrow(allEdges)]
# keptNodes <- c()
# pSeq <- seq(0,1,by=0.05)
# for (propToKeep in pSeq) {
#   topEdges <- rev(edgeScores.sorted$ix)[1:round(propToKeep*nrow(allEdges))]
#   keptNodes <- c(keptNodes,length(unique(c(allEdges[topEdges,1],allEdges[topEdges,2]))))
# }
# plot(pSeq,keptNodes,type="l")
# propToKeep <- 0.05
# topEdges <- rev(edgeScores.sorted$ix)[1:round(propToKeep*nrow(allEdges))]
filteredEdges <- allEdges
# The distribution of scores is very skewed and concentrated at certain values. 
# This is not good for modeling a weighted protein signaling network
# An alternative is to augment the FUNCOUP network using edge sign information
# from STRING
# First map gene symbols
allNodes <- unique(c(as.character(allEdges[,1]),as.character(allEdges[,2])))
write.table(allNodes,"STRING.allENSPIDs.txt",sep="\t",
            row.names=F,col.names=F,quote=F)
IDMapper <- read.delim("ENSP_geneSym.txt",sep=",",header=T)
mappableProteins <- intersect(allNodes,as.character(IDMapper$Protein.ID))
# Induce a subgraph from mappable proteins
library(igraph)
g.all <- graph.edgelist(as.matrix(allEdges[,1:2]))
E(g.all)$weight <- as.numeric(as.character(allEdges[,3]))
g.mapped <- subgraph(g.all,mappableProteins)
mapped.Edgelist <- get.edgelist(g.mapped)
mapped.Edgelist <- cbind(mapped.Edgelist,weight=E(g.mapped)$weight)
write.table(mapped.Edgelist,"geneSym.mapped.STRING.v11.act_inh.edges.txt",
            sep="\t",row.names=F,col.names=F,quote=F)
# Use Python to map node names since there may be one-to-many mappings
# from ENSG to ENSP IDs.




