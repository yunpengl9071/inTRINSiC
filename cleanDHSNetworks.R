rm(list=ls())

setwd("~/Data/multilayerNetwork/")

CGNet.files <- list.files("piq/")
CGNet.files <- CGNet.files[grep("CGNet.txt",CGNet.files,fixed=T)]

JASPAR.IDs <- read.delim("piq/JasparIDs2018.txt",header=F,check.names=F)

# Now format ID-TFName pairs to match that in the raw network file
formattedIDs <- apply(JASPAR.IDs,1,function(x){
  if (length(grep("(",as.character(x[2]),fixed=T)) & 
      length(grep(":",as.character(x[2]),fixed=T))) {
    regname <- gsub("(.*?)::(.*?)\\((.*?)\\.(.*?)\\)","\\1\\2\\3\\4",as.character(x[2]))
  } else if (length(grep(":",as.character(x[2]),fixed=T))) {
    regname <- gsub("(.*?)::(.*?)","\\1\\2",as.character(x[2]))
  } else if (length(grep("(",as.character(x[2]),fixed=T))) {
    regname <- gsub("(.*?)\\((.*?)\\.(.*?)\\)","\\1\\2\\3",as.character(x[2]))
  }  else {
    regname <- as.character(x[2])
  }
  regname <- gsub("-","",regname,fixed=T)
  regid <- gsub(".","",as.character(x[1]),fixed=T)
  paste0(regid," ",regname)
})
JASPARNames <- as.character(JASPAR.IDs[,2])
names(JASPARNames) <- formattedIDs
write.table(JASPARNames,"piq/JASPARNames.txt",sep='\t',
            quote=F,row.names = names(JASPARNames),col.names=F)
JASPARNamesToTFs <- lapply(JASPARNames,function(x){
  TFs <- strsplit(x,"(",fixed = T)[[1]][1]
  TFSplit <- strsplit(TFs,"::",fixed=T)[[1]]
  toupper(TFSplit)
})
allTFs <- unique(unlist(JASPARNamesToTFs))
# Now map TF symbols to gene symbols, splitting TF-TF complexes into single TF-gene edges

CGNet.list <- lapply(CGNet.files,function(x){
  read.delim(paste0("piq/",x),header=F,check.names=F)
})
allGenes <- Reduce("union",lapply(CGNet.list,function(x){
  unique(as.character(x[,2]))
}))

sampleNames <- sapply(CGNet.files,function(x){
  strsplit(x,".",fixed=T)[[1]][1]
})

cleanedDHSNets <- lapply(sampleNames,function(x){
  fn <- paste0("piq/",x,".expressedCGNet.txt")
  exprsCGNet <- read.delim(fn,header=F,check.names=F)
  allEdges <- apply(exprsCGNet,1,function(y){
    paste0(y[1],"_",y[2])
  })
  cat(".")
  allEdges
})
cleanedDHSNets.union <- Reduce("union",cleanedDHSNets)
cleanedDHSNets.inters <- Reduce("intersect",cleanedDHSNets)

# x <- CGNet.files[1]
# CGNet <- read.delim(paste0("piq/",x),header=F,check.names=F)
# allEdges <- lapply(1:nrow(CGNet),function(x){
#   TFs <- JASPARNamesToTFs[[as.character(CGNet[x,1])]]
#   gene <- as.character(CGNet[x,2])
#   sapply(TFs,function(y){
#     paste0(y,"_",gene)
#   })
# })


