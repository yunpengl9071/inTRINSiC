rm(list=ls())

setwd("~/Data/multilayerNetwork/data/piq")

# Compare DHS and FANTOM networks.
# We see little overlap between the two, and need
# to know why this is the case.
# One possible way to look into this is to see if the degree
# of overlap nonlinearly varies with different percent cutoffs
# when selecting motif-region pairs from PIQ results.
# Full set of motif-region pairs is way too large. For now we
# look at just a couple of different cutoffs.

#celllines <- c("Daoy","H4","A172","M059J")
celllines <- c("Daoy","H4","A172")
cutoffs <- c("0.01","0.05","0.1")
networksByCutoff.list <- lapply(cutoffs,function(cutoff){
  network.list <- lapply(celllines,function(cellline){
    fn <- paste0("diffThresNetworks/",cellline,".expressedCGNet_",cutoff,".txt")
    edgelist <- read.delim(fn,header=F,check.names=F)
    cat(".")
    edgeset <- apply(edgelist,1,function(edge){
      paste0(edge[1],"_",edge[2])
    })
  })
})

networksByCutoff.cleanedOnly.list <- lapply(cutoffs,function(cutoff){
  network.list <- lapply(celllines,function(cellline){
    fn <- paste0("diffThresNetworks/",cellline,".cleanedCGNet_",cutoff,".txt")
    edgelist <- read.delim(fn,header=F,check.names=F)
    cat(".")
    edgeset <- apply(edgelist,1,function(edge){
      paste0(edge[1],"_",edge[2])
    })
  })
})

networkUnionsByCutoff.list <- lapply(networksByCutoff.list,function(x){
  unique(Reduce("union",x[1:3]))
})
networkUnionsByCutoff.cleanedOnly.list <- lapply(networksByCutoff.cleanedOnly.list,function(x){
  unique(Reduce("union",x[1:3]))
})


FANTOM5Nets <- read.delim("../data/FANTOM5_brain_0.1_expressed.txt",
                          header=F,check.names=F)
#FANTOM5Nets <- read.delim("FANTOM5_0.1_only_edges.txt",header=F,check.names=F)
FANTOM5Nets.edgeset <- unique(sapply(1:nrow(FANTOM5Nets),function(x){
  paste0(FANTOM5Nets[x,1],"_",FANTOM5Nets[x,2])
}))

FANTOM5CoverageInDHS <- lapply(networkUnionsByCutoff.list,function(n){
  length(intersect(n,FANTOM5Nets.edgeset))
})

plot(sapply(networkUnionsByCutoff.list,length),unlist(FANTOM5CoverageInDHS))

