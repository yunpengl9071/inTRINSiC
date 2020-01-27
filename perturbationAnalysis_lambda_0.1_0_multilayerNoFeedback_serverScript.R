srm(list=ls())

setwd('~/Data/multilayerNetwork/')
source("src/coreFunctions.R")

args = commandArgs(trailingOnly = T)

subtypes <- c("Classical","Neural","Proneural","Mesenchymal")
lambda <- "0.1_0"

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


unionNet <- read.table("data/piq/TFs_for_genes.txt",sep="\t",header=F)
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
fc.full <- read.csv("data/FUNCOUP_STRING_signed_inters.txt",sep="\t")
fc.full <- data.frame(gene1=fc.full[,1],gene2=fc.full[,2],weight=fc.full[,3])
quantileCutoff <- 1
fc.full <- fc.full[fc.full$weight>quantile(fc.full$weight,1-quantileCutoff),]
geneIDs <- unique(c(as.character(fc.full[,1]),as.character(fc.full[,2])))
idmapper <- read.delim("data/hg19.IDMapper.txt",header=F)
matched <- match(geneIDs,idmapper[,1])
notinDB <- geneIDs[which(is.na(matched))]
inDB <- geneIDs[which(!is.na(matched))]
matchedSyms <- as.character(idmapper[matched[!is.na(matched)],2])
lens <- sapply(matchedSyms,nchar)
woSyms <- inDB[lens==0]
unmapped <- c(notinDB,woSyms)
symbolMap <- cbind(ensembl_gene_id=inDB[lens!=0],hgnc_gene_symbol=matchedSyms[lens!=0])
#TFList <- allTFs
TFMatches <- symbolMap[match(allTFs,symbolMap[,2]),1]
TFsInFunCoup <- TFMatches[!is.na(TFMatches)]
# But these 2 have expression values and are in the ENCODE-Net
TFsNotInFunCoup <- allTFs[is.na(TFMatches)]
# allNodes <- union(unique(c(as.character(PrNet[,1]),as.character(PrNet[,2]))),
#                   TFsInFunCoup)
allNodes <- union(unique(c(as.character(fc.full[,1]),as.character(fc.full[,2]))),
                  TFsInFunCoup)
allNodesWithExprs <- intersect(symbolMap[match(allNodes,symbolMap[,1]),2],
                               rownames(expr))
allIDsWithExprs <- symbolMap[match(allNodesWithExprs,symbolMap[,2]),1]
fullGraph <- graph.data.frame(fc.full)
trimmedGraph <- induced.subgraph(fullGraph,allIDsWithExprs)
adj <- get.adjacency(trimmedGraph,attr="weight")
adjPrNet <- matrix(0,{nrow(adj)+length(TFsNotInFunCoup)},{ncol(adj)+length(TFsNotInFunCoup)})
rownames(adjPrNet) <- c(symbolMap[match(rownames(adj),symbolMap[,1]),2],TFsNotInFunCoup)
colnames(adjPrNet) <- rownames(adjPrNet)
adjPrNet[rownames(adjPrNet)[1:dim(adj)[1]],rownames(adjPrNet)[1:dim(adj)[1]]] <- as.matrix(adj)
adjPrNet[{length(rownames(adj))+1}:dim(adjPrNet)[1],
         {length(colnames(adj))+1}:dim(adjPrNet)[2]] <- matrix(diag(length(TFsNotInFunCoup)),
                                                               length(TFsNotInFunCoup),
                                                               length(TFsNotInFunCoup))
protNetGenes <- rownames(adjPrNet)
# write.table(adjPrNet,"adjPrNet_DHSFANTOM_top0.1FUNCOUP.txt",
#             sep="\t",quote=F)

# adjPrNet <- read.csv("adjPrNet_signed_allmapped.txt",sep="\t",
#                      check.names=F)

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

# What's a good way to simulate expression changes?
# Should we start from the perturbed gene and work downstream,
# tracking all subsequently TFs one-by-one?

# Try a hierarchical model for perturbations - use recursion.
# Store current state of the expression matrix in a global variable.
# Also keep track of whether a TF has already been perturbed.

exprEstSubtypes <- lapply(c(1:length(subtypes)),exprEstSubtype)
names(exprEstSubtypes) <- subtypes

hPerturb <- function(perturbedTF) {
  targetTFs <- intersect(targets[[perturbedTF]],allTFs)
  #cat(perturbedTF)
  cat(" ")
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
      length(regs) != 0 & sum(perturbedTFList[regs]) == length(regs)
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
        cat(paste0(perturbedTF," ---> ",targetTFToPerturb))
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

# perturbedTF <- "HNF4A"
# # Build a list of regulators that are on the path from the topmost TF (KD target)
# # to the TF in question.
# regulatorsOnPath <- lapply(allTFs,function(x){
#   #regs <- regulators[[x]]
#   xxx <- sapply(names(V(g.TFNet)),function(TF){
#     are.connected(g.TFNet,perturbedTF,TF) & are.connected(g.TFNet,TF,x)
#   })
#   cat(".")
#   names(V(g.TFNet))[xxx]
# })
# names(regulatorsOnPath) <- allTFs
# 
# s <- "Classical"
# log2FoldDecrease <- 5
# subtypeExprs.current <- data.matrix(expr[,(commonLabels)==s])
# subtypeExprs.current[perturbedTF,] <- subtypeExprs.current[perturbedTF,] - log2FoldDecrease
# perturbedTFList <- rep(0,length(allTFs))
# names(perturbedTFList) <- allTFs
# perturbedTFList[perturbedTF] <- 1
# 
# system.time(hPerturb(perturbedTF))

#######################################################################################
drugTargets <- read.delim("data/druggableGenome",header=T,check.names=F)
drugTargets <- as.character(drugTargets$hgnc_names)
druggableTFs <- intersect(allTFs,drugTargets)
batchSize <- as.numeric(args[1])
batch <- as.numeric(args[2])
startIdx <- (batch-1)*batchSize+1
endIdx <- min(c(length(druggableTFs),batchSize*batch))

#perturbedTF <- druggableTFs[as.numeric(args[1])]
perturbedExprs.list <- list()

for (perturbedTF in druggableTFs[startIdx:endIdx]) {
  fn <- paste0("perturbedMatrix_TF_",perturbedTF)
  if (!file.exists(fn)) {
    # Build a list of regulators that are on the path from the topmost TF (KD target)
    # to the TF in question.
    regulatorsOnPath <- lapply(allTFs,function(x){
      xxx <- sapply(names(V(g.TFNet)),function(TF){
        are.connected(g.TFNet,perturbedTF,TF) & are.connected(g.TFNet,TF,x)
      })
      cat(".")
      names(V(g.TFNet))[xxx]
    })
    names(regulatorsOnPath) <- allTFs
    # library(foreach)
    # library(doParallel)
    # cl <- makeCluster(4)
    # registerDoParallel(cl)
    # subtypeExprs.current <<- data.matrix(expr[,(commonLabels)=="Classical"])
    # finalMatrix <- foreach(i=1:length(subtypes), .combine=cbind) %dopar% {
    #   s <- subtypes[i]
    #   log2FoldDecrease <- 5
    #   subtypeExprs.current <<- data.matrix(expr[,(commonLabels)==s])
    #   subtypeExprs.current[perturbedTF,] <<- subtypeExprs.current[perturbedTF,] - log2FoldDecrease
    #   perturbedTFList <- rep(0,length(allTFs))
    #   names(perturbedTFList) <- allTFs
    #   perturbedTFList[perturbedTF] <- 1
    #   hPerturb(perturbedTF)
    #   subtypeExprs.current #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
    # }
    # #stop cluster
    # stopCluster(cl)
    
    #Non-parallel version
    finalMatrix <- list()
    #subtypeExprs.current <<- data.matrix(expr[,(commonLabels)=="Classical"])
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
    colnames(finalMatrix) <- c(colnames(expr[,(commonLabels)=="Classical"]),
                               colnames(expr[,(commonLabels)=="Neural"]),
                               colnames(expr[,(commonLabels)=="Proneural"]),
                               colnames(expr[,(commonLabels)=="Mesenchymal"]))
    
    write.table(finalMatrix,fn,sep="\t",
                quote=F)
    perturbedExprs.list[[perturbedTF]] <- finalMatrix
  } else {
    perturbedExprs.list[[perturbedTF]] <- data.matrix(read.delim(fn,header=T,row.names=1,
                                                                 check.names = F))
  }
}
#perturbedTF <- "ZKSCAN3"

prScores.perturbed.list <- lapply(toPerturb,function(x){
  cat("\n")
  prScores.perturbed <- apply(perturbedExprs.list[[x]],2,function(exprVec){
    currmat <- updateTransitionProbExpRank(2^exprVec[protNetGenes],adjPrNet)
    prvec <- 2^exprVec[protNetGenes]/sum(2^exprVec[protNetGenes])
    res <- expRank(prvec,currmat,mu=max(currmat)-min(currmat))
    cat(".")
    res[[1]]
  })
})
names(prScores.perturbed.list) <- toPerturb
fn.suffix.ppl <- "_lambda_0_0.1_prScores.perturbed.txt"
count <- 1
for (prScores.perturbed in prScores.perturbed.list) {
  fn <- paste0("analyses/",toPerturb[[count]],fn.suffix.ppl)
  write.table(prScores.perturbed,fn,
              quote=F)
  count <- count + 1
}

fn.suffix.ppal <- "_lambda_0_0.1_prScores.perturbed.adjAdj0.1.txt"
prScores.perturbed.list.adjAdj <- lapply(toPerturb[1:length(toPerturb)],function(x){
  prScores.perturbed.adj <- lapply(1:ncol(perturbedExprs.list[[x]]),function(i){
    exprVec <- as.numeric(perturbedExprs.list[[x]][protNetGenes,i])
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

#######################################################################################

# drugTargets <- read.delim("data/druggableGenome",header=T,check.names=F)
# drugTargets <- as.character(drugTargets$hgnc_names)
# # Druggable targets have no overlap with significant TFs
# # Try top partners of significant TFs instead
# 
# #iteration <- 1
# 
# 
# 
# geneRegs.TFs.coreg <- lapply(subtypes,function(x){
#   fn <- paste0("analyses/",x,"_TF_TF_cormat_lambda_0_0_rectGate_100918.txt")
#   data.matrix(read.delim(fn,header=T,row.names=1,check.names=F))
# })
# 
# # Estimate change in expression in response to perturbing druggable
# # transcription factors
# #temp <- topSigCoreg.druggable
# druggableTFs <- intersect(allTFs,drugTargets)
# #domSigSubtype.filtered <- read.delim("domSigSubtype.filtered.txt",header=T,check.names=F)
# domSigSubtype.filtered <- read.delim("analyses/domSigSubtype.filtered_percentCutoff0.005.multipleTarget.withSCExprsTop3000.txt",
#                                      header=F,row.names=1,check.names=F)
# temp <- domSigSubtype.filtered
# domSigSubtype.filtered <- as.numeric(as.character(domSigSubtype.filtered[,1]))
# names(domSigSubtype.filtered) <- rownames(temp)
# druggableSigTFs <- intersect(names(domSigSubtype.filtered),druggableTFs)
# topSigCoreg <- lapply(names(domSigSubtype.filtered),function(x){
#   subtype <- domSigSubtype.filtered[x]
#   cormat <- geneRegs.TFs.coreg[[subtype]]
#   topCor <- names(rev(sort(cormat[x,])))[1:5]
#   topCor
# })
# # tsc <- lapply(c("FOXM1"),function(x){
# #   subtype <- 2
# #   cormat <- geneRegs.TFs.coreg[[subtype]]
# #   topCor <- names(rev(sort(cormat[x,])))[1:10]
# #   topCor
# # })
# topSigCoreg.druggable <- intersect(Reduce("union",topSigCoreg),druggableTFs)
# topSigCoreg.druggable <- union(topSigCoreg.druggable,druggableSigTFs)
# #topSigCoreg.druggable <- c(topSigCoreg.druggable,"ESR2","PPARG")
# #additionalTargets <- c("FOS","HIF1A","RXRA","HNF4A")
# #topSigCoreg.druggable <- additionalTargets
# # Are the druggable TFs coregulation partners with any signature TFs?
# # sigTFs.topCoregs <- lapply(names(domSigSubtype),function(x){
# #   TF.idx <- match(x,rownames(geneRegs.TFs.coreg[[1]]))
# #   geneRegs.topCoreg[[as.numeric(domSigSubtype[x])]][[TF.idx]]
# # })
# # sigTFs.topCoregs.druggableTFs <- lapply(sigTFs.topCoregs,intersect,druggableTFs)
# # dump("sigTFs.topCoregs.druggableTFs",
# #      file="sigTFs.topCoregs.druggableTFs.R")
# #source("sigTFs.topCoregs.druggableTFs.R")
# 
# toPerturb <- topSigCoreg.druggable
# toPerturb <- c("PIK3CA","AKT1","AKT2","AKT3","PTEN")
# toPerturb <- c("EIF3B","RAD51","ACTB","SF1","YY1")
# toPerturb <- c("RPS6","RPL7A","RPS11")
# toPerturb <- druggableTFs
# 
# exprEstSubtypes <- lapply(c(1:length(subtypes)),exprEstSubtype)
# names(exprEstSubtypes) <- subtypes
# foldDecrease <- 5
# estExprs.perturb.list <- lapply(toPerturb,estExprsPerturb,foldDecrease=foldDecrease)
# names(estExprs.perturb.list) <- toPerturb
# 
# estDeltaExprs.perturb.list <- lapply(names(estExprs.perturb.list),function(x){
#   deltaExprEstSubtypes <- lapply(c(1:length(subtypes)),function(i){
#     temp <- estExprs.perturb.list[[x]][[i]] - exprEstSubtypes[[i]]
#     temp[x,] <- rep((0 - foldDecrease),ncol(temp))
#     temp
#   })
# })
# 
# perturbedExprs.list <- lapply(estDeltaExprs.perturb.list,function(x){
#   deltaExprs <- do.call(cbind.data.frame,x)
#   deltaExprs <- data.matrix(deltaExprs)
#   exprs.perturbed <- expr[,colnames(deltaExprs)]
#   exprs.perturbed[rownames(deltaExprs),] <- exprs.perturbed[rownames(deltaExprs),] +
#     deltaExprs
#   exprs.perturbed
# })
# names(perturbedExprs.list) <- toPerturb
# 
# prScores.perturbed.list <- lapply(toPerturb,function(x){
#   cat("\n")
#   prScores.perturbed <- apply(perturbedExprs.list[[x]],2,function(exprVec){
#     currmat <- updateTransitionProbExpRank(2^exprVec[protNetGenes],adjPrNet)
#     prvec <- 2^exprVec[protNetGenes]/sum(2^exprVec[protNetGenes])
#     res <- expRank(prvec,currmat,mu=max(currmat)-min(currmat))
#     cat(".")
#     res[[1]]
#   })
# })
# names(prScores.perturbed.list) <- toPerturb
# fn.suffix.ppl <- "_lambda_0_0_prScores.perturbed.txt"
# count <- 1
# for (prScores.perturbed in prScores.perturbed.list) {
#   fn <- paste0("analyses/",toPerturb[[count]],fn.suffix.ppl)
#   write.table(prScores.perturbed,fn,
#               quote=F)
#   count <- count + 1
# }
# 
# fn.suffix.ppal <- "_lambda_0_0_prScores.perturbed.adjAdj0.1.txt"
# prScores.perturbed.list.adjAdj <- lapply(toPerturb[1:length(toPerturb)],function(x){
#   prScores.perturbed.adj <- lapply(1:ncol(perturbedExprs.list[[x]]),function(i){
#     exprVec <- as.numeric(perturbedExprs.list[[x]][protNetGenes,i])
#     currmat <- updateTransitionProbExpRank(2^exprVec,adjPrNet)
#     scores <- t(currmat) %*% as.numeric(prScores.perturbed.list[[x]][,i])
#     cat(".")
#     # if (i %% 50 == 0) {
#     #   cat(".")
#     # }
#     scores
#   })
#   fn <- paste0(x,fn.suffix.ppal)
#   write.table(prScores.perturbed.adj,fn,sep="\t",quote=F)
#   cat("\n")
#   temp <- do.call(cbind.data.frame,prScores.perturbed.adj)
#   temp <- data.matrix(temp)
#   colnames(temp) <- colnames(perturbedExprs.list[[x]])
#   temp
# })
# 
# names(prScores.perturbed.list.adjAdj) <- toPerturb
# temp <- lapply(toPerturb[1:length(toPerturb)],function(x){
#   fn <- paste0("analyses/",x,fn.suffix.ppal)
#   (write.table(prScores.perturbed.list.adjAdj[[x]],
#                fn,sep="\t",quote=F))
# })
# 
# prScores.perturbed.list.adjAdj <- lapply(toPerturb[1:length(toPerturb)],function(x){
#   fn <- paste0("analyses/",x,fn.suffix.ppal)
#   data.matrix(read.delim(fn,header=T,row.names=1,check.names=F))
# })
# 
# names(prScores.perturbed.list.adjAdj) <- toPerturb
# 
# # Also get protein scores for the unperturbed matrix
# prScores.unperturbed.list <- lapply(colnames(commonExpr),function(x){
#   cat(".")
#   exprVec <- as.numeric(commonExpr[protNetGenes,x])
#   currmat <- updateTransitionProbExpRank(2^exprVec,adjPrNet)
#   prvec <- 2^exprVec/sum(2^exprVec)
#   res <- expRank(prvec,currmat,mu=max(currmat)-min(currmat))
#   res[[1]]
# })
# prScores.unperturbed.mat <- do.call(cbind.data.frame,prScores.unperturbed.list)
# colnames(prScores.unperturbed.mat) <- colnames(commonExpr)
# prScores.unperturbed.mat <- prScores.unperturbed.mat[,c(names(commonLabels)[commonLabels=="Classical"],
#                                                         names(commonLabels)[commonLabels=="Neural"],
#                                                         names(commonLabels)[commonLabels=="Proneural"],
#                                                         names(commonLabels)[commonLabels=="Mesenchymal"])]
# prScores.unperturbed.adj.mat <- lapply(colnames(commonExpr),function(x){
#   cat(".")
#   exprVec <- as.numeric(commonExpr[protNetGenes,x])
#   currmat <- updateTransitionProbExpRank(2^exprVec,adjPrNet)
#   scores <- t(currmat) %*% as.numeric(prScores.unperturbed.mat[,x])
#   scores
# })
# prScores.unperturbed.adj.mat <- do.call(cbind.data.frame,prScores.unperturbed.adj.mat)
# colnames(prScores.unperturbed.adj.mat) <- colnames(commonExpr)
# prScores.unperturbed.adj.mat <- data.matrix(prScores.unperturbed.adj.mat)
# fn.puam <- "data/prScores.unperturbed.adj.mat_FUNCOUP_STRING_inters_0.1.txt"
# # write.table(prScores.unperturbed.adj.mat,fn.puam,
# #             sep="\t",quote=F)
# prScores.unperturbed.adj.mat <- data.matrix(read.delim(fn.puam,
#                                                        row.names=1,header=T,check.names = F))
# 
# 
# colID_order <- lapply(subtypes,function(s){
#   idx <- names(commonLabels)[commonLabels==s]
# })
# colID_order <- unlist(colID_order)
# prScores.perturbed.list.adjAdj <- lapply(prScores.perturbed.list.adjAdj,function(x){
#   colnames(x) <- colID_order
#   x
# })
# names(prScores.perturbed.list.adjAdj) <- toPerturb
# 
# prScores.unperturbed.adj.mat <- prScores.unperturbed.adj.mat[,colID_order]
# 
# compareGeneActivity <- function(gene.fav,gene.unfav,unptb,ptb,subtype) {
#   sampleIds <- names(commonLabels)[commonLabels==subtype]
#   unptb.activity <- (colMeans(unptb[gene.fav,sampleIds,drop=F])-colMeans(unptb[gene.unfav,sampleIds,drop=F]))
#   #ptb.activity <- as.numeric(ptb[match(gene,rownames(adjPrNet)),sampleIds])
#   ptb.activity <- (colMeans(ptb[gene.fav,sampleIds,drop=F])-colMeans(ptb[gene.unfav,sampleIds,drop=F]))
#   # boxplot(log2(unptb.activity+1),log2(ptb.activity+1),main=paste0(gene,"-",subtype))
#   boxplot((unptb.activity),(ptb.activity),main=paste0("readout","-",subtype))
#   # list(unptb=log2(unptb.activity+1),
#   #      ptb=log2(ptb.activity+1))
#   list(unptb=unptb.activity,
#        ptb=ptb.activity)
# }
# 
# readout <- c("CASP3","CASP9")
# readout.unfav <- c("BCL2")
# allActivityChanges <- lapply(names(prScores.perturbed.list.adjAdj),function(x){
#   res <- lapply(subtypes,function(subtype){
#     compareGeneActivity(readout,readout.unfav,
#                         prScores.unperturbed.adj.mat,
#                         prScores.perturbed.list.adjAdj[[x]],
#                         subtype)
#   })
#   subtypeDiff <- sapply(res,function(r){
#     mean((r$ptb)) - mean((r$unptb))
#   })
# })
# names(allActivityChanges) <- names(prScores.perturbed.list.adjAdj)
# 
# # Using CASP3/9/BCL2 combined readout, the following look interesting
# # MLX - proneural
# # ARNT - M (suboptimal)
# # ESRRG - C (interactor of FOSB)
# # RXRB - N
# # AR - M (suboptimal signature)
# # ESR1 - C (suboptimal signature)
# # STAT5B - C (subopt. sig.)
# # NR3C1 - M (subopt. sig.)
# # ATF1 - M (subopt. sig.)
# 
# readout <- "MLX"
# subtype <- "Classical"
# 
# temp <- compareGeneActivity(readout,readout.unfav,
#                             prScores.unperturbed.adj.mat,
#                             prScores.perturbed.list.adjAdj[[perturbedGene]],subtype)
# perturbedGene <- "ATF1"
# 
# reslist <- lapply(subtypes,function(x){
#   temp <- compareGeneActivity(readout,readout.unfav,
#                               prScores.unperturbed.adj.mat,
#                               prScores.perturbed.list.adjAdj[[perturbedGene]],x)
#   df <- cbind(temp[[1]],temp[[2]])
#   #fn <- paste0("","_",x,"_activity_",perturbedGene,"_unptb_ptb.txt")
#   write.table(df,paste0("analyses/CASP3CASP9BCL2_",x,"_activity_",perturbedGene,"_unptb_ptb.txt"),
#               sep="\t",row.names=F,col.names=F,quote=F)
# })
# 
# compareGeneExprs <- function(gene,unptb,ptb,subtype) {
#   sampleIds <- names(commonLabels)[commonLabels==subtype]
#   unptb.activity <- as.numeric(unptb[gene,sampleIds])
#   ptb.activity <- as.numeric(ptb[gene,sampleIds])
#   boxplot(unptb.activity,ptb.activity,main=paste0(gene,"-",subtype))
# }
# 
# 
# 
# 
# 
# perturbedTF.idx <- 1
# ptb.unptb.mat <- data.matrix(cbind(prScores.perturbed.list[[perturbedTF.idx]],
#                                    prScores.unperturbed.mat))
# unptb.colnames <- sapply(colnames(prScores.unperturbed.mat),function(x){
#   paste0(x,"_U")
# })
# colnames(ptb.unptb.mat) <- c(colnames(prScores.perturbed.list[[1]]),
#                              unptb.colnames)
# library(Rtsne)
# set.seed(38947329)
# # distMat2 <- log2(distMat + 0.00001)
# # diag(distMat2) <- 0
# topGenes <- rownames(ptb.unptb.mat)[rev(order(apply(ptb.unptb.mat,1,mad)))][1:200]
# tsne.clustMat.log2.scale <- Rtsne(t(t(scale(t(log2(ptb.unptb.mat))))),is_distance = F,
#                        perplexity=15,pca=F,theta=0,check_duplicates = F)
# colVec <- c(rep("red",sum(commonLabels=="Classical")),
#             rep("blue",sum(commonLabels=="Neural")),
#             rep("green",sum(commonLabels=="Proneural")),
#             rep("orange",sum(commonLabels=="Mesenchymal")))
# colVec <- rep(colVec,2)
# shapeVec <- c(rep(16,ncol(commonExpr)),
#               rep(17,ncol(commonExpr)))
# plot(tsne.clustMat.log2.scale$Y,pch=shapeVec,main="",col=colVec)
# #text(tsne.clustMat$Y,colnames(perturbMat.all.qn),cex=0.5)
# 
# compareGeneActivity <- function(gene,unptb,ptb,subtype) {
#   sampleIds <- names(commonLabels)[commonLabels==subtype]
#   unptb.activity <- as.numeric(unptb[gene,sampleIds])
#   ptb.activity <- as.numeric(ptb[match(gene,rownames(adjPrNet)),sampleIds])
#   boxplot(log2(unptb.activity+1),log2(ptb.activity+1),main=paste0(gene,"-",subtype))
#   list(unptb=log2(unptb.activity+1),
#        ptb=log2(ptb.activity+1))
# }
# 
# 
# # ARNT (all except Neural) --> FOXM1 seems interesting 
# # NR1H4 (Classical) --> CXCL12
# # RXRB (Classical SOX9 interactor) --> APEX1
# # PPARG (C/N) --> CASP3
# # ESR2 (P) --> CASP3
# # ATF1 (M) --> CASP9
# readout <- "CASP3"
# allActivityChanges <- lapply(names(prScores.perturbed.list.adjAdj),function(x){
#   res <- lapply(subtypes,function(subtype){
#     compareGeneActivity(readout,
#                         prScores.unperturbed.adj.mat,
#                         prScores.perturbed.list.adjAdj[[x]],
#                         subtype)
#   })
#   subtypeDiff <- sapply(res,function(r){
#     mean((r$ptb)) - mean((r$unptb))
#   })
# })
# 
# # RXRA and HNF4A seem interesting
# # RXRA is significant coregulator for 2 Neural and 1 Mesenchymal signature TFs
# temp <- compareGeneActivity("SHC1",
#                                       prScores.unperturbed.adj.mat,
#                                       prScores.perturbed.list.adjAdj[[2]],"Classical")
# reslist <- lapply(subtypes,function(x){
#   temp <- compareGeneActivity("CASP3",
#                       prScores.unperturbed.adj.mat,
#                       prScores.perturbed.list.adjAdj[["ESR2"]],x)
#   df <- cbind(temp[[1]],temp[[2]])
#   write.table(df,paste0("CASP3_",x,"_activity_ATF1_unptb_ptb.txt"),
#               sep="\t",row.names=F,col.names=F,quote=F)
# })
# 
# temp <- as.numeric(PRScores.all)
# temp[temp>=0.0002241] <- 0.0002241
# hist(temp,50,col="gray")
# 
# 
# 
# 
# compareGeneExprs <- function(gene,unptb,ptb,subtype) {
#   sampleIds <- names(commonLabels)[commonLabels==subtype]
#   unptb.activity <- as.numeric(unptb[gene,sampleIds])
#   ptb.activity <- as.numeric(ptb[gene,sampleIds])
#   boxplot(unptb.activity,ptb.activity,main=paste0(gene,"-",subtype))
# }
# 
# # shapeVec <- c(rep(16,length(toPerturb)*length(table(commonLabels))),
# #               rep(17,length(commonLabels)))
# 
# # Proof of concept perturbations: known signaling relationships
# 
# toPerturb <- c("EGFR") # Assay for stuff like AKT
# toPerturb <- c("MDM2","EGFR")
# toPerturb <- c("MYC")
# foldDecrease <- 5
# estExprs.perturb.list <- lapply(toPerturb,function(x){
#   expr.perturb <- expr
#   expr.perturb[x,] <- expr.perturb[x,] - foldDecrease
#   #temp <- expr.perturb[x,]
#   exprEstSubtypes <- lapply(1:length(subtypes), function(i) {
#     cat(".")
#     type <- subtypes[i]
#     idx <- names(commonLabels)[commonLabels==type]
#     tfmat <- expr.perturb[allTFs,idx]
#     capmat <- capacities.list[[type]]
#     exprEstimates <- lapply(colnames(capmat),function(x){
#       regs <- regulators[[x]]
#       tfsubmat <- tfmat[regs,,drop=F]
#       capvec <- capmat[regs,x]
#       const <- capmat["const",x]
#       # cat(x)
#       # cat('\n')
#       exprsEst <- as.numeric(apply(tfsubmat,2,model.pred,
#                                    mivec=c(),fvecmi=c(),fvectf=capvec,const=const))
#     })
#     exprEstimates <- do.call(rbind.data.frame,exprEstimates)
#     exprEstimates <- data.matrix(exprEstimates)
#     colnames(exprEstimates) <- idx
#     rownames(exprEstimates) <- colnames(capmat)
#     # if (x %in% rownames(exprEstimates)) {
#     #   exprEstimates[x,] <- exprEstimates
#     # }
#     exprEstimates
#   })
#   cat("\n")
#   exprEstSubtypes
# })
# names(estExprs.perturb.list) <- toPerturb
# # 
# # # Then apply change to original expression matrix (to avoid systematic bias)
# estDeltaExprs.perturb.list <- lapply(names(estExprs.perturb.list),function(x){
#   deltaExprEstSubtypes <- lapply(c(1:length(subtypes)),function(i){
#     temp <- estExprs.perturb.list[[x]][[i]] - exprEstSubtypes[[i]]
#     temp[x,] <- rep((0 - foldDecrease),ncol(temp))
#     temp
#   })
# })
# perturbedExprs.list <- lapply(estDeltaExprs.perturb.list,function(x){
#   deltaExprs <- do.call(cbind.data.frame,x)
#   deltaExprs <- data.matrix(deltaExprs)
#   exprs.perturbed <- expr[,colnames(deltaExprs)]
#   exprs.perturbed[rownames(deltaExprs),] <- exprs.perturbed[rownames(deltaExprs),] +
#     deltaExprs
#   exprs.perturbed
# })
# names(perturbedExprs.list) <- toPerturb
# 
# prScores.perturbed.list <- lapply(toPerturb,function(x){
#   cat("\n")
#   prScores.perturbed <- apply(perturbedExprs.list[[x]],2,function(exprVec){
#     currmat <- updateTransitionProbExpRank(2^exprVec[protNetGenes],adjPrNet)
#     prvec <- 2^exprVec[protNetGenes]/sum(2^exprVec[protNetGenes])
#     res <- expRank(prvec,currmat,mu=max(currmat)-min(currmat))
#     cat(".")
#     res[[1]]
#   })
# })
# 
# names(prScores.perturbed.list) <- toPerturb
# 
# count <- 1
# for (prScores.perturbed in prScores.perturbed.list) {
#   fn <- paste0(toPerturb[[count]],"_prScores.perturbed0.1.txt")
#   write.table(prScores.perturbed,fn,
#               quote=F)
#   count <- count + 1
# }
# # # dump("prScores.perturbed.list",file = "prScores.perturbed.list.R")
# prScores.perturbed.list <- lapply(toPerturb,function(x){
#   fn <- paste0(x,"_prScores.perturbed0.1.txt")
#   data.matrix(read.table(fn,sep=" ",row.names=1,header=T,check.names = F))
# })
# names(prScores.perturbed.list) <- toPerturb
# 
# prScores.perturbed.list.adjAdj <- lapply(toPerturb,function(x){
#   prScores.perturbed.adj <- lapply(1:ncol(perturbedExprs.list[[x]]),function(i){
#     exprVec <- as.numeric(perturbedExprs.list[[x]][protNetGenes,i])
#     currmat <- updateTransitionProbExpRank(2^exprVec,adjPrNet)
#     scores <- t(currmat) %*% as.numeric(prScores.perturbed.list[[x]][,i])
#     if (i %% 50 == 0) {
#       cat(".")
#     }
#     scores
#   })
#   cat("\n")
#   temp <- do.call(cbind.data.frame,prScores.perturbed.adj)
#   temp <- data.matrix(temp)
#   colnames(temp) <- colnames(perturbedExprs.list[[x]])
#   temp
# })
# names(prScores.perturbed.list.adjAdj) <- toPerturb
# for (g in toPerturb) {
#   fn <- paste0(g,"_prScores.perturbed.adjAdj0.1.txt")
#   write.table(prScores.perturbed.list.adjAdj[[g]],
#               fn,sep="\t",quote=F)
# }
# 
# prScores.perturbed.list.adjAdj <- lapply(toPerturb,function(x){
#   fn <- paste0(x,"_prScores.perturbed.adjAdj0.1.txt")
#   data.matrix(read.delim(fn,row.names=1,header=T,check.names=F))
# })
# names(prScores.perturbed.list.adjAdj) <- toPerturb
# 
# # capPr2miRNA.list <- lapply(cl,function(z){
# #   data.matrix(read.csv(paste0(z,"_TFmiRNACapacities_ri_nomi_",paste0("_lambda",lambda),
# #                               "_iter_",as.character(iteration-1),".txt"),
# #                        sep="\t",check.names=F))
# # })
# # names(capPr2miRNA.list) <- cl
# # capmi2mRNA.list <- lapply(cl,function(z){
# #   data.matrix(read.csv(paste0(z,"_mimRNACapacities_ri_nomi_",paste0("_lambda",lambda),
# #                               "_iter_",as.character(iteration-1),".txt"),
# #                        sep="\t",check.names=F))
# # })
# # names(capmi2mRNA.list) <- cl
# 
# # model.pred <- function(tfvec, mivec, fvecmi, fvectf, const) {
# #   if (length(fvecmi)==0 | sum(fvecmi)==length(fvecmi)) {
# #     temp <- log2(const) + sum(log2((1 + fvectf * (tfvec^2)) / (1 + tfvec^2)))
# #     temp
# #   } else {
# #     temp <- log2(const) + sum(0.9 * log2((1 + fvectf * (tfvec)^2) / (1 + tfvec^2))) + sum(0.1 * log2((1 + fvecmi * (mivec)^2) / (1 + mivec^2)))
# #     temp
# #   }
# # }
# # 
# # exprEstSubtypes <- lapply(1:length(subtypes), function(i) {
# #   type <- subtypes[i]
# #   idx <- names(commonLabels)[commonLabels==type]
# #   tfmat <- expr[allTFs,idx]
# #   capmat <- capacities.list[[type]]
# #   exprEstimates <- lapply(colnames(capmat),function(x){
# #     regs <- regulators[[x]]
# #     tfsubmat <- tfmat[regs,,drop=F]
# #     capvec <- capmat[regs,x]
# #     const <- capmat["const",x]
# #     # cat(x)
# #     # cat('\n')
# #     exprsEst <- as.numeric(apply(tfsubmat,2,model.pred,
# #                                  mivec=c(),fvecmi=c(),fvectf=capvec,const=const))
# #   })
# #   exprEstimates <- do.call(rbind.data.frame,exprEstimates)
# #   exprEstimates <- data.matrix(exprEstimates)
# #   colnames(exprEstimates) <- idx
# #   rownames(exprEstimates) <- colnames(capmat)
# #   cat(".")
# #   exprEstimates
# # })
# estExprsPerturb <- function(x, foldDecrease){
#   expr.perturb <- expr
#   expr.perturb[x,] <- expr.perturb[x,] - foldDecrease
#   #temp <- expr.perturb[x,]
#   exprEstSubtypes <- lapply(1:length(subtypes), function(i) {
#     type <- subtypes[i]
#     idx <- names(commonLabels)[commonLabels==type]
#     tfmat <- expr.perturb[allTFs,idx]
#     capmat <- capacities.list[[type]]
#     exprEstimates <- lapply(colnames(capmat),function(x){
#       regs <- regulators[[x]]
#       tfsubmat <- tfmat[regs,,drop=F]
#       capvec <- capmat[regs,x]
#       const <- capmat["const",x]
#       # cat(x)
#       # cat('\n')
#       exprsEst <- as.numeric(apply(tfsubmat,2,model.pred,
#                                    mivec=c(),fvecmi=c(),fvectf=capvec,const=const))
#     })
#     exprEstimates <- do.call(rbind.data.frame,exprEstimates)
#     exprEstimates <- data.matrix(exprEstimates)
#     colnames(exprEstimates) <- idx
#     rownames(exprEstimates) <- colnames(capmat)
#     cat(".")
#     # if (x %in% rownames(exprEstimates)) {
#     #   exprEstimates[x,] <- exprEstimates
#     # }
#     exprEstimates
#   })
#   cat("\n")
#   exprEstSubtypes
# }
