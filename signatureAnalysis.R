setwd('~/Data/multilayerNetwork/')
source("src/coreFunctions.R")

subtypes <- c("Classical","Neural","Proneural","Mesenchymal")
lambda <- "0_0"
CV.fold <- 5

expr <- read.csv("data/hgMatrix.txt",sep="\t")
labels <- read.table("data/class_labels.txt",sep="\t",header=F)
labels <- as.character(labels[,1])
names(labels) <- colnames(expr)
commonSamples <- colnames(expr)
commonLabels <- labels[commonSamples]
commonExpr <- expr[,commonSamples]
colnamevec <- c(names(commonLabels)[commonLabels=="Classical"],
                names(commonLabels)[commonLabels=="Neural"],
                names(commonLabels)[commonLabels=="Proneural"],
                names(commonLabels)[commonLabels=="Mesenchymal"])
# colnamevec <- gsub(".","-",colnamevec,fixed=T)
# colnamevec <- gsub("-01\\b","",colnamevec)
sampleByLabels <- lapply(names(table(commonLabels)),function(x){
  commonLabels[commonLabels==x]
})
names(sampleByLabels) <- names(table(commonLabels))

lambda <- "0_0"
capacities.list <- lapply(subtypes,function(x){
  fn <- paste0(,x,"_TFGeneCapacities_ri_nomi_DHSFANTOM__lambda",lambda,".txt")
  temp <- data.matrix(read.delim(fn,header=T,row.names=1,check.names=F))
  constid <- which(rownames(temp)=="const")
  temp[-constid,]
})
names(capacities.list) <- subtypes

names(capacities.list) <- subtypes
capacities.log2.list <- lapply(capacities.list,log2)
# A F value is significant only when it induces gene expression changes
# of 2 fold or greater
abslog2FCutoff <- 1

unionNet <- read.table("data/piq/TFs_for_genes.txt",sep="\t",header=F)
allTFs <- unique(as.character(unionNet[,2]))
allGenes <- unique(as.character(unionNet[,1]))
unionNet.byGene <- split(unionNet,f=unionNet[,1])
unionNet.byTF <- split(unionNet,f=unionNet[,2])
regulators <- lapply(unionNet.byGene,function(x){
  as.character(x[,2])
})

source("src/sigAnalysisUtilityFunctions.R")

# The following are precomputed F value / correlation matrices and their
# corresponding signature subtype matrices without any filtering.

subtypeSpecificMats.bgAdj <- lapply(subtypes,function(x){
  data.matrix(read.delim(paste0(x,"_subtypeSpecificMats.bgAdj.txt"),
                         header=T,row.names=1,check.names=F))
})
sigSubtype <- data.matrix(read.delim("sigSubtype.mat.txt",
                                     header=T,row.names=1,check.names=F))

outDegrees <- sapply(unionNet.byTF,nrow)
outDegrees <- outDegrees[allTFs]
sigOutDegrees <- (outDegrees >= 1)
percentSigCutoff <- 0.005
sigSubtype.filtered <- sigSubtype[sigOutDegrees,]
domSigSubtype.filtered <- sapply(rownames(sigSubtype.filtered),function(x){
  temp <- table(sigSubtype.filtered[x,])
  if (length(temp)>1) {
    temp <- temp[names(temp)!="0"]
    if(max(temp)/outDegrees[x] >= percentSigCutoff) {
      names(temp)[temp==max(temp)][1]
    } else {
      "0"
    }
  } else {
    "0"
  }
})
domSigSubtype.filtered.multipleTgts <- sapply(names(domSigSubtype.filtered[domSigSubtype.filtered!="0"]),function(x){
  tgts <- tgtGenes <- which(sigSubtype[x,]==as.numeric(domSigSubtype.filtered[x]))
  length(tgts) > 5
})
GSE57872.all <- read.delim("v2.1/GSE57872_GBM_data_matrix.txt",header=T,row.names=1)
GSE57872.all <- data.matrix(GSE57872.all)

sigTFs <- domSigSubtype.filtered[intersect(rownames(GSE57872.all),
                                           names(domSigSubtype.filtered[domSigSubtype.filtered!="0"]
                                                 [domSigSubtype.filtered.multipleTgts]))]

geneRegs.TFs.coreg <- lapply(subtypes,function(x){
  fn <- paste0(x,"_TF_TF_cormat_DFNet_rectGate_031318.txt")
  data.matrix(read.delim(fn,header=T,row.names=1,check.names=F))
})

sigTFCorSubtype <- data.matrix(read.delim("sigTFCorSubtype.mat.txt",
                                          header=T,row.names=1,check.names=F))


# Examine if signature TFs and their interaction partners
# bear consistent mutational signatures
# Extend analysis to first neighbours in STRING interaction map

mutMat <- data.matrix(read.delim("unified_nonsilent_mut_0.05_GBM.txt",
                                 header=T,row.names=1,check.names=T))
mutMat <- mutMat[,intersect(colnames(commonExpr),colnames(mutMat))]
mutMat.bySubtype <- lapply(subtypes,function(x){
  subtypeSamples <- names(commonLabels)[commonLabels==x]
  mutMat[,intersect(subtypeSamples,colnames(mutMat))]
})
names(mutMat.bySubtype) <- subtypes

sigTF.mutSig <- list()
for (sigTF in names(sigTFs)) {
  if (sigTF %in% rownames(mutMat)) {
    sigTF.mutSig[[sigTF]] <- sapply(mutMat.bySubtype,function(x){
      sum(x[sigTF,]!=0)
    })
  }
}

nPartners <- 5
sigTF.partners <- lapply(names(sigTFs),function(sigTF){
  ss <- as.numeric(sigTFs[sigTF])
  corVec <- geneRegs.TFs.coreg[[ss]][sigTF,]
  topPartners <- names(sort(corVec,decreasing = T))[1:nPartners]
  topPartners
})
names(sigTF.partners) <- names(sigTFs)

sigTF.partners.mutSig <- list()
for (sigTF in names(sigTFs)) {
  availMut <- intersect(rownames(mutMat),sigTF.partners[[sigTF]])
  if (length(availMut) > 0) {
    temp <- lapply(availMut,function(x){
      sapply(mutMat.bySubtype,function(y){
        sum(y[x,]!=0)
      })
    })
    sigTF.partners.mutSig[[sigTF]] <- data.matrix(do.call(rbind.data.frame,temp))
  }
}

STRING.network.df <- read.delim("../STRING/geneSym.mapped.STRING.edges.txt",
                                header=F,check.names=F)
ENSP.ID.mapper <- read.table("../ENSP_geneSym.txt",sep=",",header=T,check.names=F)
STRING.network.geneSym.1 <- as.character(STRING.network.df[,1])
STRING.network.geneSym.2 <- as.character(STRING.network.df[,2])
STRING.network.geneSym.1 <- match(STRING.network.geneSym.1,
                                  as.character(ENSP.ID.mapper$`Protein ID`))
STRING.network.geneSym.1 <- as.character(ENSP.ID.mapper$`Associated Gene Name`[STRING.network.geneSym.1])
STRING.network.geneSym.2 <- match(STRING.network.geneSym.2,
                                  as.character(ENSP.ID.mapper$`Protein ID`))
STRING.network.geneSym.2 <- as.character(ENSP.ID.mapper$`Associated Gene Name`[STRING.network.geneSym.2])
STRING.network.geneSym <- cbind(STRING.network.geneSym.1,STRING.network.geneSym.2)
STRING.network <- STRING.network.geneSym[complete.cases(STRING.network.geneSym),]
STRING.network <- data.frame(STRING.network)
library(igraph)
STRING.g <- graph.data.frame(STRING.network)
STRING.adj <- data.matrix(get.adjacency(STRING.g))

sigTF.neighbours <- lapply(names(sigTFs),function(sigTF){
  if (sigTF %in% rownames(STRING.adj)) {
    union(names(which(STRING.adj[sigTF,]!=0)),
          names(which(STRING.adj[,sigTF]!=0)))
  } else {
    c()
  }
})
names(sigTF.neighbours) <- names(sigTFs)

sigTF.neighbours.mutSig <- list()
for (sigTF in names(sigTFs)) {
  availMut <- intersect(rownames(mutMat),sigTF.neighbours[[sigTF]])
  if (length(availMut) > 0) {
    temp <- lapply(availMut,function(x){
      sapply(mutMat.bySubtype,function(y){
        sum(y[x,]!=0)
      })
    })
    temp <- data.matrix(do.call(rbind.data.frame,temp))
    colnames(temp) <- subtypes
    rownames(temp) <- availMut
    sigTF.neighbours.mutSig[[sigTF]] <- temp
  }
}
