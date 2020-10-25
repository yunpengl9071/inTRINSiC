rm(list=ls())

setwd('/ahg/regevdata/projects/txnRegModeling/regression/perturbation/')
source("coreFunctions.R")
#.libPaths('/ahg/regevdata/projects/txnRegModeling/regression/Rlibs')

args = commandArgs(trailingOnly = T)

subtypes <- c("Classical","Neural","Proneural","Mesenchymal")
lambda <- "0.1_0"

cat("Loading data...\n")
noniter <- 0
curriter <- 0
expr <- read.csv("hgMatrix.txt",sep="\t")
labels <- read.table("class_labels.txt",sep="\t",header=F)
labels <- as.character(labels[,1])
names(labels) <- colnames(expr)
commonSamples <- colnames(expr)
commonLabels.all <- labels[commonSamples]
commonLabels <- commonLabels.all
commonExpr <- expr[,commonSamples]
colnamevec <- c(names(commonLabels)[commonLabels=="Classical"],
                names(commonLabels)[commonLabels=="Neural"],
                names(commonLabels)[commonLabels=="Proneural"],
                names(commonLabels)[commonLabels=="Mesenchymal"])
subtypeColors <- rep("red",length(commonLabels))
subtypeColors[which(commonLabels=="Neural")] <- "blue"
subtypeColors[which(commonLabels=="Proneural")] <- "green"
subtypeColors[which(commonLabels=="Mesenchymal")] <- "orange"


unionNet <- read.table("TFs_for_genes.txt",sep="\t",header=F)
allTFs <- unique(as.character(unionNet[,2]))
allGenes <- unique(as.character(unionNet[,1]))
unionNet.byGene <- split(unionNet,f=unionNet[,1])
unionNet.byTF <- split(unionNet,f=unionNet[,2])
regulators <- lapply(unionNet.byGene,function(x){
  as.character(x[,2])
})
targets <- lapply(unionNet.byTF,function(x){
  as.character(x[,1])
})
# Build a TF-TF regulatory network.
library(igraph)

unionNet.TFs <- sapply(1:nrow(unionNet),function(i){
  (as.character(unionNet[i,1]) %in% allTFs) &
    (as.character(unionNet[i,2]) %in% allTFs)
})
unionNet.TFs <- unionNet[unionNet.TFs,]
# unionNet.TFs <- data.frame(unionNet.TFs,stringsAsFactors = F)
unionNet.TFs.ch <- lapply(unionNet.TFs,as.character)
unionNet.TFs.ch <- cbind(unionNet.TFs.ch[[1]],unionNet.TFs.ch[[2]])
g.TFNet <- graph.edgelist((unionNet.TFs.ch),directed = T)


# Load the output matrices from the regression pipeline.

lambda <- "0.1_0"
capacities.list <- lapply(subtypes,function(x){
  fn <- paste0(x,"_TFGeneCapacities_ri_nomi__lambda",lambda,".txt")
  temp <- data.matrix(read.delim(fn,header=T,row.names=1,check.names=F))
  #constid <- which(rownames(temp)=="const")
  #temp[-constid,]
})
names(capacities.list) <- subtypes
capacities.log2.list <- lapply(capacities.list,log2)

# Map protein network to available nodes
# PrNet <- read.csv("../FUNCOUP_signed.txt",sep="\t",header=F)
library(igraph)
# g <- graph.data.frame(PrNet)
# adjPrNet <- get.adjacency(g)
# Extend protein & mRNA list to include all TFs in the random walk model,
# regardless of whether they interact in the FunCoup network
# Use full FC network to fill up TF-TF interactions
print("Preparing adjacency matrix...\n")
fc.full <- read.csv("FUNCOUP_STRING_signed_act_inh_inters.txt",sep="\t")
fc.full <- data.frame(gene1=fc.full[,1],gene2=fc.full[,2],weight=fc.full[,3])
quantileCutoff <- 1
fc.full <- fc.full[fc.full$weight>quantile(fc.full$weight,1-quantileCutoff),]
geneIDs <- unique(c(as.character(fc.full[,1]),as.character(fc.full[,2])))
idmapper <- read.delim("hg19.IDMapper.txt",header=F)
matched <- match(geneIDs,idmapper[,1])
notinDB <- geneIDs[which(is.na(matched))]
inDB <- geneIDs[which(!is.na(matched))]
matchedSyms <- as.character(idmapper[matched[!is.na(matched)],2])
lens <- sapply(matchedSyms,nchar)
woSyms <- inDB[lens==0]
unmapped <- c(notinDB,woSyms)
symbolMap <- cbind(ensembl_gene_id=inDB[lens!=0],hgnc_gene_symbol=matchedSyms[lens!=0])
TFMatches <- symbolMap[match(allTFs,symbolMap[,2]),1]
TFsInFunCoup <- TFMatches[!is.na(TFMatches)]
TFsNotInFunCoup <- allTFs[is.na(TFMatches)]
allNodes <- union(unique(c(as.character(fc.full[,1]),as.character(fc.full[,2]))),
                  TFsInFunCoup)
allNodesWithExprs <- intersect(symbolMap[match(allNodes,symbolMap[,1]),2],
                               rownames(expr))
allIDsWithExprs <- symbolMap[match(allNodesWithExprs,symbolMap[,2]),1]
fullGraph <- graph.data.frame(fc.full)
trimmedGraph <- induced.subgraph(fullGraph,allIDsWithExprs)
adj <- get.adjacency(trimmedGraph,attr="weight")
adjPrNet <- data.matrix(adj)
rownames(adjPrNet) <- colnames(adjPrNet) <- c(symbolMap[match(rownames(adj),symbolMap[,1]),2])
protNetGenes <- rownames(adjPrNet)

model.pred <- function(tfvec, mivec, fvecmi, fvectf, const) {
  if (length(fvecmi)==0 | sum(fvecmi)==length(fvecmi)) {
    temp <- log2(const) + sum(log2((1 + fvectf * (tfvec^2)) / (1 + tfvec^2)))
    temp
  } else {
    temp <- log2(const) + sum(0.9 * log2((1 + fvectf * (tfvec)^2) / (1 + tfvec^2))) + sum(0.1 * log2((1 + fvecmi * (mivec)^2) / (1 + mivec^2)))
    temp
  }
}

exprEstSubtype <- function(i){
  type <- subtypes[i]
  idx <- names(commonLabels)[commonLabels==type]
  tfmat <- expr[allTFs,idx]
  capmat <- capacities.list[[type]]
  exprEstimates <- lapply(colnames(capmat),function(x){
    regs <- regulators[[x]]
    tfsubmat <- tfmat[regs,,drop=F]
    capvec <- capmat[regs,x]
    const <- capmat["const",x]
    # cat(x)
    # cat('\n')
    exprsEst <- as.numeric(apply(tfsubmat,2,model.pred,
                                 mivec=c(),fvecmi=c(),fvectf=capvec,const=const))
  })
  exprEstimates <- do.call(rbind.data.frame,exprEstimates)
  exprEstimates <- data.matrix(exprEstimates)
  colnames(exprEstimates) <- idx
  rownames(exprEstimates) <- colnames(capmat)
  cat(".")
  exprEstimates
}

computeDeltaExprs <- function(target,subtypeExprsEst,currTFExprs,capmat) {
  TFs <- regulators[[target]]
  TFSubmat <- currTFExprs[TFs,colnames(subtypeExprsEst),drop=F]
  capvec <- capmat[TFs,target]
  const <- capmat["const",target]
  exprsEstNew <- as.numeric(apply(TFSubmat,2,model.pred,
                               mivec=c(),fvecmi=c(),fvectf=capvec,const=const))
  deltaExprsEst <- exprsEstNew - subtypeExprsEst
  deltaExprsEst
}

# Use a hierarchical model for perturbations - use recursion.
# Store current state of the expression matrix in a global variable.
# Also keep track of whether a TF has already been perturbed.

print("Estimating expression...\n")
exprEstSubtypes <- lapply(c(1:length(subtypes)),exprEstSubtype)
names(exprEstSubtypes) <- subtypes

hPerturb <- function(perturbedTF) {
  targetTFs <- intersect(targets[[perturbedTF]],allTFs)
  #cat(perturbedTF)
  cat(".")
  #print(sum(perturbedTFList))
  # Check if targets contain TFs
  if (length(targetTFs)){
    # Do not perturb those target TFs that have already been perturbed and have
    # perturbed the TF (equivalent to being connected in the TF subnetwork) in
    # question (i.e. found a feedback loop!)
    # Also try to be more efficient: a TF shouldn't need be perturbed unless all of its
    # regulators that are on the path from the topmost TF (knockdown target) have already
    # been perturbed - this will avoid perturbing a TF more than once.
    doNotPerturb <- targetTFs[is.finite(igraph::distances(g.TFNet,targetTFs,perturbedTF,mode="in"))]
    doNotPerturb <- intersect(doNotPerturb, allTFs[perturbedTFList==1])
    # doNotPerturb <- union(doNotPerturb,targetTFs[sapply(regulatorsOnPath[targetTFs],function(regs){
    #   length(regs) == 0 | sum(perturbedTFList[regs]) == length(regs)
    # })])
    doNotPerturb <- union(doNotPerturb,targetTFs[sapply(regulatorsOnPath[targetTFs],function(regs){
      length(regs) != 0 & sum(perturbedTFList[regs]) < length(regs)
    })])
    targetsToPerturb <- setdiff(targets[[perturbedTF]],doNotPerturb)
    targetTFsToPerturb <- intersect(targetsToPerturb,targetTFs)
    temp <- lapply(targetsToPerturb,function(target){
      targetExprsEst <- exprEstSubtypes[[s]][target,,drop=F]
      #cat(".")
      xxx <- subtypeExprs.current
      xxx[target,] <- computeDeltaExprs(target,
                                                         subtypeExprsEst = targetExprsEst,
                                                      currTFExprs = subtypeExprs.current[allTFs,],
                                                      capmat = capacities.list[[s]]) +
        subtypeExprs.current[target,]
      subtypeExprs.current <<- xxx
      rm(xxx)
    })
    # If all TFs within the target gene set have already been perturbed, also exit.
    if (length(targetTFsToPerturb)==0) {
      return()
    } else {
      # Otherwise look recursively through all of the unperturbed target TFs.
      perturbedTFList[targetTFsToPerturb] <<- 1
      for (targetTFToPerturb in targetTFsToPerturb) {
        print(sum(perturbedTFList))
        #cat(" ")
        #cat(paste0(perturbedTF," ---> ",targetTFToPerturb))
        hPerturb(targetTFToPerturb)
      }
    }
  } else { 
    # If none of the targets are TFs then perturb and return.
    temp <- lapply(targets[[perturbedTF]],function(target){
      targetExprsEst <- exprEstSubtypes[[s]][target,,drop=F]
      xxx <- subtypeExprs.current
      xxx[target,] <- computeDeltaExprs(target,
                                        subtypeExprsEst = targetExprsEst,
                                        currTFExprs = subtypeExprs.current[allTFs,],
                                        capmat = capacities.list[[s]]) +
        subtypeExprs.current[target,]
      subtypeExprs.current <<- xxx
    })
    return()
  }
}

#######################################################################################
drugTargets <- read.delim("druggableGenome",header=T,check.names=F)
drugTargets <- as.character(drugTargets$hgnc_names)
#druggableTFs <- intersect(allTFs,drugTargets)
#druggableTFs <- c("AKT1","AKT2","AKT3","PTEN","YY1","PIK3CA")
druggableTFs <- unique(c(allTFs,"PTEN","PIK3CA","AKT1","YY1"))
batchSize <- as.numeric(args[1])
batch <- as.numeric(args[2])
startIdx <- (batch-1)*batchSize+1
endIdx <- min(c(length(druggableTFs),batchSize*batch))
#startIdx <- 5
#endIdx <- 5
cat("Perturbing TFs...\n")
#perturbedTF <- druggableTFs[as.numeric(args[1])]
perturbedExprs.list <- list()
toPerturb <- druggableTFs[startIdx:endIdx]

for (perturbedTF in druggableTFs[startIdx:endIdx]) {
  cat(paste0("Now perturbing ",perturbedTF,"\n"))
  fn <- paste0("perturbedMatrix_TF_",perturbedTF)
  #if (T) {
  if (!file.exists(fn)) {
    # Build a list of regulators that are on the path from the topmost TF (KD target)
    # to the TF in question.
    if (perturbedTF %in% allTFs) {
    regulatorsOnPath <- lapply(allTFs,function(x){
      xxx <- sapply(names(V(g.TFNet)),function(TF){
        are.connected(g.TFNet,perturbedTF,TF) & are.connected(g.TFNet,TF,x)
      })
      cat(".")
      names(V(g.TFNet))[xxx]
    })
    names(regulatorsOnPath) <- allTFs}
    cat("Start multilayer perturbation...\n")
    finalMatrix <- list()
    log2FoldDecrease <- 5
    for(i in 1:length(subtypes)) {
      s <- subtypes[i]
      subtypeExprs.current <- data.matrix(expr[,(commonLabels)==s])
      subtypeExprs.current[perturbedTF,] <- subtypeExprs.current[perturbedTF,] - log2FoldDecrease
      perturbedTFList <- rep(0,length(allTFs))
      names(perturbedTFList) <- allTFs
      perturbedTFList[perturbedTF] <- 1
      hPerturb(perturbedTF)
      finalMatrix[[i]] <- subtypeExprs.current #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
    }
    finalMatrix <- data.matrix(do.call(cbind.data.frame,finalMatrix))
    # Change this to the specific subtypes/cell types used
    colnames(finalMatrix) <- c(colnames(expr[,(commonLabels)=="Classical"]),
                               colnames(expr[,(commonLabels)=="Neural"]),
                               colnames(expr[,(commonLabels)=="Proneural"]),
                               colnames(expr[,(commonLabels)=="Mesenchymal"]))
    finalMatrix[finalMatrix > 15] <- 15
    write.table(finalMatrix,fn,sep="\t",
                quote=F)
    perturbedExprs.list[[perturbedTF]] <- finalMatrix
  } else {
    perturbedExprs.list[[perturbedTF]] <- data.matrix(read.delim(fn,header=T,row.names=1,
                                                                 check.names = F))
  }
}
gc()

cat("Running exponential ranking...\n")
protNetGenesID <- match(protNetGenes,rownames(perturbedExprs.list[[1]]))
prScores.perturbed.list <- lapply(toPerturb,function(x){
  cat("\n")
  prScores.perturbed <- apply(perturbedExprs.list[[x]],2,function(exprVec){
    currmat <- updateTransitionProbExpRank(2^exprVec[protNetGenesID],adjPrNet)
    prvec <- 2^exprVec[protNetGenesID]/sum(2^exprVec[protNetGenesID])
    res <- expRank(prvec,currmat,mu=max(currmat)-min(currmat))
    cat(".")
    res[[1]]
  })
})
names(prScores.perturbed.list) <- toPerturb
fn.suffix.ppl <- "_lambda_0_0.1_prScores.perturbed_act_inh.txt"
count <- 1
for (prScores.perturbed in prScores.perturbed.list) {
  fn <- paste0("analyses/",toPerturb[[count]],fn.suffix.ppl)
  write.table(prScores.perturbed,fn,
              quote=F)
  count <- count + 1
}

fn.suffix.ppal <- "_lambda_0_0.1_prScores.perturbed_act_inh.adjAdj0.1.txt"
prScores.perturbed.list.adjAdj <- lapply(toPerturb[1:length(toPerturb)],function(x){
  prScores.perturbed.adj <- lapply(1:ncol(perturbedExprs.list[[x]]),function(i){
    exprVec <- as.numeric(perturbedExprs.list[[x]][protNetGenesID,i])
    currmat <- updateTransitionProbExpRank(2^exprVec,adjPrNet)
    scores <- t(currmat) %*% as.numeric(prScores.perturbed.list[[x]][,i])
    cat(".")
    # if (i %% 50 == 0) {
    #   cat(".")
    # }
    scores
  })
  #fn <- paste0(x,fn.suffix.ppal)
  temp <- do.call(cbind.data.frame,prScores.perturbed.adj)
  temp <- data.matrix(temp)
  cat("\n")
  colnames(temp) <- colnames(perturbedExprs.list[[x]])
  temp
})

names(prScores.perturbed.list.adjAdj) <- toPerturb
temp <- lapply(toPerturb[1:length(toPerturb)],function(x){
  fn <- paste0(x,fn.suffix.ppal)
  (write.table(prScores.perturbed.list.adjAdj[[x]],
               fn,sep="\t",quote=F))
})

