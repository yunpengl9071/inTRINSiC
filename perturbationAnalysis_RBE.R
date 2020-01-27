rm(list=ls())
setwd('~/Data/multilayerNetwork/')
source("src/coreFunctions.R")

subtypes <- c("Classical","Neural","Proneural","Mesenchymal")
lambda <- "0.1_0"

noniter <- 0
curriter <- 0
expr <- read.csv("data/combinedExprsRBE.txt",sep="\t")
labels <- read.table("data/combinedLabels.txt",sep="\t",header=F)
labels <- as.character(labels[,2])
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

labels <- labels[1:613]
commonLabels <- commonLabels[1:613]
commonExpr <- commonExpr[,1:613]
expr <- expr[,1:613]
ccle.class.labels <- commonLabels

lambda <- "0.1_0"
capacities.list <- lapply(subtypes,function(x){
  fn <- paste0("data/regrOutput/withCCLE/",x,"_TFGeneCapacities_ri_nomi__lambda",lambda,".txt")
  temp <- data.matrix(read.delim(fn,header=T,row.names=1,check.names=F))
  #constid <- which(rownames(temp)=="const")
  #temp[-constid,]
})
names(capacities.list) <- subtypes
capacities.log2.list <- lapply(capacities.list,log2)

options(stringsAsFactors = F)
unionNet <- read.table("data/piq/TFs_for_genes.txt",sep="\t",header=F)
allTFs <- setdiff(rownames(capacities.list[[1]]),c('const'))
allGenes <- colnames(capacities.list[[1]])
unionNet.TFs <- sapply(1:nrow(unionNet),function(i){
  (as.character(unionNet[i,1]) %in% allGenes) &
    (as.character(unionNet[i,2]) %in% allTFs)
})
unionNet <- unionNet[unionNet.TFs,]

#allTFs <- unique(as.character(unionNet[,2]))
#allGenes <- unique(as.character(unionNet[,1]))
unionNet.byGene <- split(unionNet,f=unionNet[,1],drop = T)
unionNet.byTF <- split(unionNet,f=unionNet[,2],drop = T)
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
unionNet.TFs.ch <- cbind(unionNet.TFs.ch[[2]],unionNet.TFs.ch[[1]])
g.TFNet <- graph.edgelist((unionNet.TFs.ch),directed = T)

celllines.info <- read.csv("data/DepMap-2018q3-celllines.csv",header=T)
celllines.cns <- celllines.info[grep("CENTRAL_NERVOUS_SYSTEM",celllines.info$CCLE_Name),]

library(igraph)
# g <- graph.data.frame(PrNet)
# adjPrNet <- get.adjacency(g)
# Extend protein & mRNA list to include all TFs in the random walk model,
# regardless of whether they interact in the FunCoup network
# Use full FC network to fill up TF-TF interactions
print("Preparing adjacency matrix...\n")
fc.full <- read.csv("data/FUNCOUP_STRING_signed_act_inh_inters.txt",sep="\t")
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
#adjPrNet <- matrix(0,{nrow(adj)+length(TFsNotInFunCoup)},{ncol(adj)+length(TFsNotInFunCoup)})
#rownames(adjPrNet) <- c(symbolMap[match(rownames(adj),symbolMap[,1]),2],TFsNotInFunCoup)
#colnames(adjPrNet) <- rownames(adjPrNet)
#adjPrNet[rownames(adjPrNet)[1:dim(adj)[1]],rownames(adjPrNet)[1:dim(adj)[1]]] <- as.matrix(adj)
#adjPrNet[{length(rownames(adj))+1}:dim(adjPrNet)[1],
#         {length(colnames(adj))+1}:dim(adjPrNet)[2]] <- matrix(diag(length(TFsNotInFunCoup)),
#                                                               length(TFsNotInFunCoup),
#                                                               length(TFsNotInFunCoup))
adjPrNet <- data.matrix(adj)
rownames(adjPrNet) <- colnames(adjPrNet) <- c(symbolMap[match(rownames(adj),symbolMap[,1]),2])
protNetGenes <- rownames(adjPrNet)


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
fn.puam <- "data/prScores.combinedRBE.unperturbed.adj.mat_FUNCOUP_STRING_inters_0.1.txt"
# write.table(prScores.unperturbed.adj.mat,fn.puam,
#             sep="\t",quote=F)
prScores.unperturbed.adj.mat <- data.matrix(read.delim(fn.puam,header=T,row.names=1,
                                                       check.names = F))

ccle.DepMap.scores <- read.csv('data/gene_effect.csv',
                               header=T, row.names = 1, check.names = F)
colnames(ccle.DepMap.scores) <- sapply(colnames(ccle.DepMap.scores),function(x){
  strsplit(x," ",fixed=T)[[1]][1]
})
ccle.withDepMap <- match(celllines.cns$Broad_ID,
                         rownames(ccle.DepMap.scores))
ccle.withDepMap.idx <- which(complete.cases(ccle.withDepMap))
ccle.withDepMap.names <- as.character(celllines.cns$CCLE_Name)[ccle.withDepMap.idx]
ccle.class.labels.withDepMap <- ccle.class.labels[ccle.withDepMap.names]
ccle.withDepMap <- ccle.withDepMap[complete.cases(ccle.withDepMap)]
ccle.DepMap.cns.scores <- ccle.DepMap.scores[ccle.withDepMap,]
rownames(ccle.DepMap.cns.scores) <- ccle.withDepMap.names
ccle.DepMap.cns.scores <- t(ccle.DepMap.cns.scores)
ccle.DepMap.cns.TF.scores <- ccle.DepMap.cns.scores[intersect(rownames(ccle.DepMap.cns.scores),
                                                              allTFs),]

# Key assumption is that certain TFs show differential essentiality scores across
# CCLE brain tumor cell lines of different subtypes.

# Identify perturbed expRanking scores that are correlated with essentiality scores
# for each TF perturbed and select ones which are constantly correlated.

# May need to combine scores from individual genes into 'metagenes' to denoise.

# Use all cell lines in initial exploratory analysis, but need to perform
# some sort of training-validation split afterwards taking facility of experimental
# validation into consideration.

fn.suffix <- "_lambda_0_0.1_RBE_CCLE_prScores.perturbed_act_inh.adjAdj0.1.txt"
availTFs <- intersect(rownames(ccle.DepMap.cns.scores),
                      allTFs)
#prScores.unperturbed.adj.scaled.mat <- t(scale(t(prScores.unperturbed.adj.mat)))
ccle.prScores.corWithEssScores <- lapply(availTFs,function(TF){
  fn <- paste0("data/perturbation/",TF,fn.suffix)
  TF.prScores <- data.matrix(read.delim(fn,sep="\t",header=T,row.names=1,check.names=F))
  combined.deltaScores <- cbind(TF.prScores,prScores.unperturbed.adj.mat[,colnames(TF.prScores)])
  # Probably should not scale the scores (???). Will introduce bias when predicting scores on new samples (???).
  # 091219: log2 transform positive and negative PRScore values separately
  combined.deltaScores[combined.deltaScores>0] <- log2(combined.deltaScores[combined.deltaScores>0]+1)
  combined.deltaScores[combined.deltaScores<0] <- (-log2(-combined.deltaScores[combined.deltaScores<0]+1))
  #combined.deltaScores <- t(scale(t(combined.deltaScores)))
  TF.deltaPRScores <- combined.deltaScores[,1:ncol(TF.prScores)] - 
    combined.deltaScores[,(ncol(TF.prScores)+1):ncol(combined.deltaScores)]
  fn.out <- paste0("data/perturbation/",TF,"_lambda_0_0.1_RBE_CCLE_deltaPRScores.perturbed_act_inh.adjAdj0.1.txt")
  write.table(TF.deltaPRScores,fn.out,sep="\t",quote=F)
  availLines <- intersect(ccle.withDepMap.names,colnames(TF.prScores))
  TF.prScores <- TF.prScores[,availLines]
  TF.DepMepScores <- ccle.DepMap.cns.TF.scores[TF,availLines]
  cat(".")
  # apply(TF.prScores,1,function(x){
  #   cor(x,TF.DepMepScores)
  # })
})
# ccle.prScores.corWithEssScores <- lapply(ccle.prScores.corWithEssScores,function(x){
#   sort(x[complete.cases(x)],decreasing = F)
# })
# names(ccle.prScores.corWithEssScores) <- availTF
TF <- "CTCF"

fn <- paste0("data/perturbation/",TF,fn.suffix)
TF.prScores <- data.matrix(read.delim(fn,sep="\t",header=T,row.names=1,check.names=F))
availLines <- intersect(ccle.withDepMap.names,colnames(TF.prScores))
write.table(ccle.DepMap.cns.TF.scores[,availLines],'data/CCLE_DepMap_cns_scores_avail.txt',
            sep='\t',quote=F)
# ccle.prScores.corWithEssScores <- lapply(ccle.prScores.corWithEssScores,function(x){
#   temp <- rep(0,nrow(TF.prScores))
#   names(temp) <- rownames(TF.prScores)
#   temp[names(x)] <- x
#   temp
#   #sort(x[complete.cases(x)],decreasing = F)
# })
# names(ccle.prScores.corWithEssScores) <- availTFs
# ccle.prScores.corWithEssScores.mat <- data.matrix(do.call(cbind.data.frame,
#                                                           ccle.prScores.corWithEssScores))
# head(rev(sort(rowMeans(ccle.prScores.corWithEssScores.mat))))
# 
# consistentTopScores <- Reduce('intersect',lapply(ccle.prScores.corWithEssScores,function(x){
#   names(x)[1:500]
# }))
# 
# myheatmapfun <- colorRampPalette((c(rgb(15,61,111,maxColorValue=255),
#                                     rgb(125,182,212,maxColorValue=255),
#                                     rgb(255,255,255,maxColorValue=255),
#                                     rgb(234,145,118,maxColorValue=255),
#                                     rgb(102,3,32,maxColorValue=255))))
# library(pheatmap)
# heatmap.data <- ccle.prScores.corWithEssScores.mat[rowSums(ccle.prScores.corWithEssScores.mat)!=0,]
# xxx <- pheatmap(heatmap.data,
#                 breaks=seq(-1,1,length.out=76),
#                 clustering_method = "ward.D",
#                 #clustering_distance_cols = "euclidean",
#                 color = myheatmapfun(75),dendrogram=F,
#                 labels_row = "",labels_col = "")

# Looks like a simple regression method does not capture good features that explain
# essentiality scores well. Try using multitask elastic net regression with CV in
# scikit-learn instead.
# Prep data for scikit-learn

testSamples <- c('GB1_CENTRAL_NERVOUS_SYSTEM',
                 'F5_CENTRAL_NERVOUS_SYSTEM',
                 'T98G_CENTRAL_NERVOUS_SYSTEM',
                 'DAOY_CENTRAL_NERVOUS_SYSTEM',
                 'LN18_CENTRAL_NERVOUS_SYSTEM',
                 'D283MED_CENTRAL_NERVOUS_SYSTEM',
                 'HS683_CENTRAL_NERVOUS_SYSTEM',
                 'U87MG_CENTRAL_NERVOUS_SYSTEM')
testSamples <- intersect(availLines,testSamples)

fn.suffix <- "_lambda_0_0.1_RBE_CCLE_deltaPRScores.perturbed_act_inh.adjAdj0.1.txt"
ccle.deltaPRScores.list <- lapply(availTFs,function(TF){
  fn <- paste0("data/perturbation/",TF,fn.suffix)
  TF.deltaPRScores <- data.matrix(read.delim(fn,sep="\t",header=T,row.names=1,check.names=F))
  TF.deltaPRScores <- TF.deltaPRScores[rowSums(TF.deltaPRScores)!=0,]
  TF.deltaPRScores <- TF.deltaPRScores[complete.cases(TF.deltaPRScores),]
  TF.availCovariates <- rownames(TF.deltaPRScores)
  #TF.deltaPRScores <- t(scale(t(TF.deltaPRScores[,availLines]))) # Should not scale the differences
  cat(".")
  #write.table(TF.deltaPRScores,paste0("data/perturbation/avail_unscaled_",TF,fn.suffix),sep='\t',quote=F)
  #write.table(TF.availCovariates,paste0("data/perturbation/",TF,"availCovariates.txt"),
  #            sep="\t",row.names=F,col.names=F,quote=F)
  TF.deltaPRScores
})
# ccle.zeroSDs.list <- lapply(ccle.deltaPRScores.list,function(x){
#   cat(".")
#   rownames(x)[which(apply(x,1,sd)==0)]
# })
# allZeroSDs <- Reduce("union",ccle.zeroSDs.list)
names(ccle.deltaPRScores.list) <- availTFs

write.table(ccle.DepMap.cns.TF.scores[,setdiff(availLines,testSamples)],
            'data/CCLE_DepMap_cns_scores_avail_TV.txt',
            sep='\t',quote=F)

library(pheatmap)
library(viridis)
heatmap.data <- ccle.DepMap.cns.TF.scores[,setdiff(availLines,testSamples)]
pheatmap(heatmap.data,
         breaks=seq(-0.5,0.5,length.out=76),
         clustering_method = "ward.D",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         #clustering_distance_cols = "euclidean",
         color = magma(75),dendrogram=F,treeheight_row = 0,
         treeheight_col = 0,labels_row = F,labels_col = F)
#labels_row = "",labels_col = "")

temp <- lapply(availTFs,function(TF){
  fn <- paste0("data/perturbation/avail_unscaled_delta_",TF,fn.suffix)
  fn.test <- paste0("data/perturbation/avail_unscaled_delta_test_",TF,fn.suffix)
  TF.deltaPRScores <- ccle.deltaPRScores.list[[TF]]
  #TF.deltaPRScores <- t(scale(t(TF.deltaPRScores[setdiff(rownames(TF.deltaPRScores),allZeroSDs),])))
  #TF.deltaPRScores <- TF.deltaPRScores[setdiff(rownames(TF.deltaPRScores),allZeroSDs),]
  TF.deltaPRScores.TV <- TF.deltaPRScores[,setdiff(availLines,testSamples)]
  TF.deltaPRScores.test <- TF.deltaPRScores[,testSamples]
  cat(".")
  write.table(TF.deltaPRScores.TV,fn,sep='\t',quote=F)
  write.table(TF.deltaPRScores.test,fn.test,sep='\t',quote=F)
})

fn.suffix <- "_lambda_0_0.1_RBE_CCLE_deltaPRScores.perturbed_act_inh.adjAdj0.1.txt"
prScores.test.list <- lapply(availTFs,function(TF){
  fn.test <- paste0("data/perturbation/avail_unscaled_delta_test_",TF,fn.suffix)
  cat(".")
  data.matrix(read.delim(fn.test,header=T,row.names=1,check.names=F))
})
names(prScores.test.list) <- availTFs

# elnetRes <- read.delim('data/perturbation/elasticNetCVRes.txt',header=F,row.names=1,
#                        check.names=F)
# elnetRes <- read.delim('data/perturbation/deltaScores_elasticNetCVRes.txt',header=F,row.names=1,
#                        check.names=F)
elnetRes <- read.delim('data/perturbation/deltaScores_elasticNetCVRes_unscaled_NZCovariates.txt',header=F,row.names=1,
                       check.names=F)
elnetRes <- data.matrix(elnetRes)
intercept <- elnetRes[,1]
elnetRes <- elnetRes[,2:ncol(elnetRes)]
rownames(elnetRes) <- rownames(ccle.DepMap.cns.TF.scores)
colnames(elnetRes) <- rownames(prScores.test.list[[1]])
names(intercept) <- rownames(ccle.DepMap.cns.TF.scores)
model.TFs <- rownames(ccle.DepMap.cns.TF.scores)
hist(elnetRes[elnetRes!=0],100)
elnetRes.sig <- elnetRes[apply((abs(elnetRes)),1,max)>0.01,apply(abs(elnetRes),2,max)>0.01]
dim(elnetRes.sig)
library(pheatmap)
library(viridis)
heatmap.data <- elnetRes.sig
xxx <- pheatmap(heatmap.data,
                breaks=seq(-0.02,0.02,length.out=76),
                clustering_method = "median",
                #clustering_distance_cols = "euclidean",
                color = viridis(75),dendrogram=F)
                #labels_row = "",labels_col = "")

elnetRes.nonZeros <- elnetRes[rowSums(elnetRes)!=0,]
elnetRes.cormat <- cor(t(elnetRes.nonZeros))
heatmap.data <- elnetRes.cormat
xxx <- pheatmap(heatmap.data,
                breaks=seq(-0.5,0.5,length.out=76),
                clustering_method = "median",
                #clustering_distance_cols = "euclidean",
                color = magma(75),dendrogram=F,fontsize = 1)
#labels_row = "",labels_col = "")


# Now predict essentiality scores using the trained model
# on the test cell lines
geneNames <- colnames(elnetRes)
ccle.essentiality.test <- ccle.DepMap.cns.TF.scores[,testSamples]
ccle.essentiality.test.pred <- lapply(model.TFs,function(x){
  cat(".")
  params <- c(intercept[x],elnetRes[x,])
  pred <- sapply(testSamples,function(ts){
    covs <- prScores.test.list[[x]][,ts]
    params[1] + sum(params[2:length(params)] * covs)
  })
})
ccle.essentiality.test.pred <- data.matrix(do.call(rbind.data.frame,
                                                   ccle.essentiality.test.pred))
rownames(ccle.essentiality.test.pred) <- model.TFs
colnames(ccle.essentiality.test.pred) <- testSamples

# Also perform prediction across all CCLE lines for final
# presentation
ccle.essentiality.allSamples <- ccle.DepMap.cns.TF.scores[,availLines]
ccle.essentiality.allSamples.pred <- lapply(model.TFs,function(x){
  cat(".")
  params <- c(intercept[x],elnetRes[x,])
  pred <- sapply(colnames(ccle.essentiality.allSamples),function(ts){
    covs <- ccle.deltaPRScores.list[[x]][,ts]
    params[1] + sum(params[2:length(params)] * covs)
  })
})
ccle.essentiality.allSamples.pred <- data.matrix(do.call(rbind.data.frame,
                                                   ccle.essentiality.allSamples.pred))
rownames(ccle.essentiality.allSamples.pred) <- model.TFs
colnames(ccle.essentiality.allSamples.pred) <- availLines


ccle.essentiality.test.cor <- sapply(model.TFs,function(TF){
  cor(ccle.essentiality.test[TF,],ccle.essentiality.test.pred[TF,])
})
ccle.essentiality.allSamples.cor <- sapply(model.TFs,function(TF){
  cor(ccle.essentiality.allSamples[TF,],ccle.essentiality.allSamples.pred[TF,])
})

# Adjusted with mean of absolute values for fairness of comparison
ccle.essentiality.test.rmse <- sapply(model.TFs,function(TF){
  (mean((ccle.essentiality.test[TF,]-ccle.essentiality.test.pred[TF,])^2))^0.5/
    mean(abs(ccle.essentiality.test[TF,]))
})

# ccle.essentiality.test.sMAPE <- sapply(model.TFs,function(TF){
#   sMAPEs <- apply(c(1:ncol(ccle.essentiality.test.pred)),function(x){
#     
#   })
#   (mean((ccle.essentiality.test[TF,]-ccle.essentiality.test.pred[TF,])^2))^0.5
# })

drugTargets <- read.delim("data/druggableGenome",header=T,check.names=F)
drugTargets <- as.character(drugTargets$hgnc_names)

ccle.essentiality.test.cor <- ccle.essentiality.test.cor[rownames(elnetRes)[rowSums(abs(elnetRes))>0.01]]
head(rev(sort(ccle.essentiality.test.cor)),50)
ccle.essentiality.allSamples.cor <- ccle.essentiality.allSamples.cor[rownames(elnetRes)
                                                                     [rowSums(abs(elnetRes))>0.01]]
head(rev(sort(ccle.essentiality.allSamples.cor)),50)

ccle.essentiality.test.rmse <- ccle.essentiality.test.rmse[rownames(elnetRes)[rowSums(abs(elnetRes))>0.01]]
head((sort(ccle.essentiality.test.rmse)),50)

intersect(drugTargets,
          names(ccle.essentiality.test.cor[which(ccle.essentiality.test.cor>0.5)]))

intersect(drugTargets,
          names(ccle.essentiality.test.rmse[which(ccle.essentiality.test.rmse<0.54)]))

#names(ccle.essentiality.test.cor[which(ccle.essentiality.test.cor>0.5)])
head(names(ccle.essentiality.test.cor)
     [complete.cases(ccle.essentiality.test.cor)][rev(order(ccle.essentiality.test.cor
                                                            [complete.cases(ccle.essentiality.test.cor)]))],15)
head(names(ccle.essentiality.test.rmse)
     [complete.cases(ccle.essentiality.test.rmse)][(order(ccle.essentiality.test.rmse
                                                          [complete.cases(ccle.essentiality.test.rmse)]))],15)


plot(ccle.essentiality.test["RXRA",],
     ccle.essentiality.test.pred["RXRA",],pch=16)

#TFoI <- "STAT1"
plotTFRes <- function(TFoI) {
  df <- data.frame(X=ccle.essentiality.test.pred[TFoI,],
                   Y=ccle.essentiality.test[TFoI,])
  testSamples.labels <- ccle.class.labels[testSamples]
  names(testSamples.labels) <- gsub("_CENTRAL_NERVOUS_SYSTEM","",
                                    names(testSamples.labels),fixed=T)
  rownames(df) <- names(testSamples.labels)
  
  fit <- lm(Y~X, data=df)
  lineFunc <- function(x){
    y <- fit$coefficients[1] + fit$coefficients[2] * x
  }
  library(gplots)
  rangeX <- range(df$X)[2] - range(df$X)[1]
  rangeY <- range(df$Y)[2] - range(df$Y)[1]
  ggplot(data = df, aes(x = df$X,y = df$Y,color = as.factor(testSamples.labels))) + 
    geom_point(size = 5) + geom_text(size = 5, label = names(testSamples.labels),
                                     nudge_x = 0.00,nudge_y=rangeY*0.05) +
    scale_color_manual(values=c('brown3','goldenrod1','cornflowerblue','mediumseagreen')) +
    theme(panel.background = element_blank()) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA,size = 2),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20)) + stat_function(fun=lineFunc,size=2,color='gray') +
    # ylim(min(c(min(df$Y),min(df$X)))-0.01,max(c(max(df$Y),max(df$X)))+0.01) +
    # xlim(min(c(min(df$Y),min(df$X)))-0.01,max(c(max(df$Y),max(df$X)))+0.01)
  ylim(min(c(min(df$Y)))-0.05*rangeY,max(c(max(df$Y)))+0.05*rangeY) +
    xlim(min(c(min(df$X)))-0.05*rangeX,max(c(max(df$X)))+0.05*rangeX)
}

plotTFRes.TrnVlnTst <- function(TFoI) {
  df <- data.frame(X=ccle.essentiality.allSamples.pred[TFoI,availLines],
                   Y=ccle.essentiality.allSamples[TFoI,availLines])
  samples.labels <- ccle.class.labels[availLines]
  names(samples.labels) <- gsub("_CENTRAL_NERVOUS_SYSTEM","",
                                    names(samples.labels),fixed=T)
  rownames(df) <- names(samples.labels)
  
  fit <- lm(Y~X, data=df)
  lineFunc <- function(x){
    y <- fit$coefficients[1] + fit$coefficients[2] * x
  }
  library(gplots)
  rangeX <- range(df$X)[2] - range(df$X)[1]
  rangeY <- range(df$Y)[2] - range(df$Y)[1]
  ggplot(data = df, aes(x = df$X,y = df$Y,color = as.factor(samples.labels))) + 
    geom_point(size = 5) +
    scale_color_manual(values=c('brown3','goldenrod1','cornflowerblue','mediumseagreen')) +
    theme(panel.background = element_blank()) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA,size = 2)
          #axis.text.x = element_text(size=20),
          #axis.text.y = element_text(size=20)
          ) + stat_function(fun=lineFunc,size=2,color='gray') +
    # ylim(min(c(min(df$Y),min(df$X)))-0.01,max(c(max(df$Y),max(df$X)))+0.01) +
    # xlim(min(c(min(df$Y),min(df$X)))-0.01,max(c(max(df$Y),max(df$X)))+0.01)
    ylim(min(c(min(df$Y)))-0.05*rangeY,max(c(max(df$Y)))+0.05*rangeY) +
    xlim(min(c(min(df$X)))-0.05*rangeX,max(c(max(df$X)))+0.05*rangeX)
}
head(names(ccle.essentiality.test.cor)
     [complete.cases(ccle.essentiality.test.cor)][rev(order(ccle.essentiality.test.cor
                                                            [complete.cases(ccle.essentiality.test.cor)]))],15)
head(names(ccle.essentiality.test.rmse)
     [complete.cases(ccle.essentiality.test.rmse)][(order(ccle.essentiality.test.rmse
                                                          [complete.cases(ccle.essentiality.test.rmse)]))],15)

TFoI <- "NFE2"
#plotTFRes(TFoI)
plotTFRes.TrnVlnTst(TFoI)
# 0.9177996
cor(ccle.essentiality.allSamples[TFoI,],
    ccle.essentiality.allSamples.pred[TFoI,])

TFoI <- "MYBL2"
#plotTFRes(TFoI)
plotTFRes.TrnVlnTst(TFoI)
# 0.7551809
cor(ccle.essentiality.test[TFoI,],
    ccle.essentiality.test.pred[TFoI,])
cor(ccle.essentiality.allSamples[TFoI,],
    ccle.essentiality.allSamples.pred[TFoI,])


#geom_text(size = 2, color = 'black')
yyy <- sapply(1:1000000,function(x){
  y <- sample(xxxl,21,replace = F)
  sum(y=='red') < 9 & sum(y=='green') < 7 & sum(y=='yellow') < 5 & sum(y=='blue') < 3 & sum(y=='brown') < 1
})
# ccle.prScores.bySample <- lapply(availLines,function(line){
#   allTFData <- lapply(ccle.prScores.list,function(ccle.prScores){
#     ccle.prScores[,line]
#   })
#   allTFData <- data.matrix(do.call(cbind.data.frame,allTFData))
#   rownames(allTFData) <- rownames(TF.prScores)
#   colnames(allTFData) <- availTFs
#   # # Select top 500 variable genes as predictors
#   # topMADs <- rev(order(sapply(allTFData)))
#   fn <- paste0('data/perturbation/CCLE_RBE_',line,'.prScores.allTFs.txt')
#   cat(".")
#   write.table(allTFData,fn,sep='\t',quote=F)
# })

# These TCGA tumor sample scores were computed by concatenating
# the TCGA and CCLE expression data after RemoveBatchEffects in the limma package,
# sending it down the regression pipeline, and running the perturbation algorithm
# only on the TCGA subset. In this way, we can utilize the trained weights (from 
# CCLE data) for predicting gene essentiality scores using protein activity scores.

# Should subtract unperturbed data

fn.suffix <- "_lambda_0_0.1_RBE_TCGA_prScores.perturbed_act_inh.adjAdj0.1.txt"
xxx <- list.files("data/perturbation/")
xxx <- xxx[grep(fn.suffix,xxx,fixed=T)]
names(xxx) <- sapply(strsplit(xxx,"_",fixed=T),function(x){
  x[1]
})
missingTFs <- setdiff(availTFs,names(xxx))
write.table(missingTFs,'data/perturbation/missingTFs.txt',
            sep='\t',row.names=F,col.names=F,quote=F)
availTFs <- intersect(availTFs,names(xxx))

# Process data the same way as done on CCLE training/test data
TCGA.prScores.ess.pred <- lapply(availTFs,function(TF){
  fn <- paste0("data/perturbation/",xxx[TF])
  #availLines <- intersect(ccle.withDepMap.names,colnames(TF.prScores))
  #TF.prScores <- TF.prScores[,availLines]
  #TF.DepMepScores <- ccle.DepMap.cns.TF.scores[TF,availLines]
  cat(".")
  TF.prScores <- data.matrix(read.delim(fn,sep="\t",header=T,row.names=1,check.names=F))
  combined.deltaScores <- cbind(TF.prScores,prScores.unperturbed.adj.mat[,colnames(TF.prScores)])
  combined.deltaScores[combined.deltaScores>0] <- log2(combined.deltaScores[combined.deltaScores>0]+1)
  combined.deltaScores[combined.deltaScores<0] <- (-log2(-combined.deltaScores[combined.deltaScores<0]+1))
  TF.deltaPRScores <- combined.deltaScores[,1:ncol(TF.prScores)] - 
    combined.deltaScores[,(ncol(TF.prScores)+1):ncol(combined.deltaScores)]
  #zeroSDs <- apply(TF.deltaPRScores,1,sd)
  #zeroSDs <- which(zeroSDs==0)
  #TF.deltaPRScores <- t(scale(t(TF.deltaPRScores)))
  #TF.deltaPRScores[zeroSDs,] <- 0
  params <- c(intercept[TF],elnetRes[TF,])
  pred <- sapply(colnames(TF.prScores),function(s){
    covs <- combined.deltaScores[rownames(prScores.test.list[[1]]),s]
    params[1] + sum(params[2:length(params)] * covs)
  })
})
TCGA.prScores.ess.pred.mat <- do.call(rbind.data.frame,TCGA.prScores.ess.pred)
TCGA.prScores.ess.pred.mat <- data.matrix(TCGA.prScores.ess.pred.mat)
rownames(TCGA.prScores.ess.pred.mat) <- availTFs
sampleByLabels <- c(names(labels)[1:544][labels[1:544]=="Classical"],
                    names(labels)[1:544][labels[1:544]=="Neural"],
                    names(labels)[1:544][labels[1:544]=="Proneural"],
                    names(labels)[1:544][labels[1:544]=="Mesenchymal"])
orderedLabels <- c(rep("Classical",sum(labels[1:544]=="Classical")),
                   rep("Neural",sum(labels[1:544]=="Neural")),
                   rep("Proneural",sum(labels[1:544]=="Proneural")),
                   rep("Mesenchymal",sum(labels[1:544]=="Mesenchymal")))
names(orderedLabels) <- sampleByLabels
colnames(TCGA.prScores.ess.pred.mat) <- sampleByLabels
write.table(TCGA.prScores.ess.pred.mat,"data/perturbation/TCGA.prScores.ess.pred.mat.doubleExpTransf.unscaled.txt",
            sep='\t',quote=F)

ccle.essentiality.allSamples.cor.valid <- ccle.essentiality.allSamples.cor[complete.cases(ccle.essentiality.allSamples.cor)]
hist(rev(sort(ccle.essentiality.allSamples.cor.valid)),
     col="gray",
     xlab="Pearson correlation coefficient",
     main="Correlation between predicted\nand experimental CCLE essentiality scores")

topTFs <- intersect(names(rev(sort(ccle.essentiality.allSamples.cor.valid)))[1:50],
                    rownames(TCGA.prScores.ess.pred.mat))
anova.list <- lapply(topTFs,function(TF){
  temp <- TCGA.prScores.ess.pred.mat[TF,]
  temp.df <- cbind.data.frame(score=temp,label=as.factor(orderedLabels))
  anova.res <- aov(score~label,data=temp.df)
  p <- summary(anova.res)[[1]][['Pr(>F)']][1]
})
anova.list <- unlist(anova.list)
names(anova.list) <- topTFs
anova.list[which(anova.list<1E-5)]
# intersect(names(anova.list[which(anova.list<1E-3)]),names(sigTFs))
# sigTFs[intersect(names(anova.list[which(anova.list<1E-3)]),names(sigTFs))]
p.adjust(anova.list,"bonferroni")

library(gplots)
TF <- "NEUROG2"
plotTFSubtypeScores <- function(TF,yupper,ylower) {
  temp <- TCGA.prScores.ess.pred.mat[TF,]
  lbl <- as.factor(orderedLabels)
  lbl <- factor(lbl,levels(lbl)[c(1,3,4,2)])
  temp.df <- cbind.data.frame(score=temp,label=lbl)
  #levels(temp.df$label) <- subtypes
  p <- ggplot(temp.df, aes(label,score,fill = label))
  p + geom_boxplot(aes(label,score,fill = label)) +
    geom_point(aes(x = label), shape = 21, position = 
                 position_jitterdodge(jitter.width = 0.5, jitter.height=0.4, 
                                      dodge.width=0.9)) +
    scale_fill_manual(values=c('brown3','cornflowerblue','mediumseagreen','goldenrod1')) +
    theme(panel.background = element_blank()) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA,size = 2)) + ylim(ylower,yupper) + ggtitle(TF)
}

# for (TF in names(anova.list[which(anova.list<1E-10)])) {
#   plotTFSubtypeScores(TF)
# }
#xxx <- lapply(names(anova.list[which(anova.list<1E-10)]),plotTFSubtypeScores)
plotTFSubtypeScores("HOXA3",1,-1)
plotTFSubtypeScores("NFE2",2,-2)
plotTFSubtypeScores("FOSL1",2,-2)
plotTFSubtypeScores("MYBL2",0.25,-1.8)
plotTFSubtypeScores("FOXF1",1,-1)

TFoI <- "YY1"
plotTFRes(TFoI)
cor(ccle.essentiality.test[TFoI,],
    ccle.essentiality.test.pred[TFoI,])

# plotTFSubtypeScores("RFX2")
# plotTFSubtypeScores("ELK4")
# # ETV6 seems interesting
# plotTFSubtypeScores("ETV6",50,-0)
# # IRF9 seems interesting
# plotTFSubtypeScores("IRF9",500,-500)
# # BACH1 seems to be a Neural coregulatory signature TF
# plotTFSubtypeScores("BACH1",5,-5)
# # NFIC is (?) a tumor suppressor
# plotTFSubtypeScores("NFIC",5000,-5000)
# # GABPA is a tumor suppressor
# plotTFSubtypeScores("GABPA",500,-500)
# # NFE2L2 is a tumor suppressor
# plotTFSubtypeScores("NFE2L2",5000,-5000)
# # CREB1 seems to be a Proneural coregulatory signature TF
# plotTFSubtypeScores("CREB1",5000,-5000)
# plotTFSubtypeScores("FOXJ1",500,-500)
# plotTFSubtypeScores("HOXA2")
# plotTFSubtypeScores("KLF1")
# plotTFSubtypeScores("POU2F2")
# plotTFSubtypeScores("THRA")
# plotTFSubtypeScores("PAX8")
# plotTFSubtypeScores("MAF")
# plotTFSubtypeScores("PBX1")
# plotTFSubtypeScores("JUN")
# plotTFSubtypeScores("ELK4")
# plotTFSubtypeScores("FOXF1")
# plotTFSubtypeScores("HOXC10")
# plotTFSubtypeScores("GFI1")
# plotTFSubtypeScores("HOXD3")
# plotTFSubtypeScores("HLX")
# plotTFSubtypeScores("PAX7")
# plotTFSubtypeScores("BATF")
# plotTFSubtypeScores("TCF4")
# plotTFSubtypeScores("MYBL2")
# plotTFSubtypeScores("ZNF410")
# plotTFSubtypeScores("IRF2")
# plotTFSubtypeScores("E2F8")
# plotTFSubtypeScores("LMX1B")
# plotTFSubtypeScores("CRX")
# plotTFSubtypeScores("ONECUT1")
# plotTFSubtypeScores("FOSL2")
# plotTFSubtypeScores("ERG")
# plotTFSubtypeScores("TCF12")
# plotTFSubtypeScores("ETV2")
# plotTFSubtypeScores("IRX4")
# plotTFSubtypeScores("GATA4")
# plotTFSubtypeScores("ZKSCAN3")
# plotTFSubtypeScores("FOXM1")
# plotTFSubtypeScores("ETV6")
# plotTFSubtypeScores("CREM")
# plotTFSubtypeScores("MYOD1")
# plotTFSubtypeScores("NFKB2")
# plotTFSubtypeScores("GBX2")
# plotTFSubtypeScores("MGA")
# plotTFSubtypeScores("FOXF1")

TCGA.clinical.data <- read.delim("data/GBM.clin.merged.txt",header=F,row.names=1,
                                 check.names=F)
TCGA.clinical.samples <- toupper(as.character(TCGA.clinical.data["patient.bcr_patient_barcode",]))
names(TCGA.clinical.samples) <- as.character(TCGA.clinical.data["patient.bcr_patient_barcode",])
TCGA.samplesWithExprs <- sapply(names(labels[1:544]),function(x){
  ele <- strsplit(x,".",fixed=T)[[1]][1:3]
  paste(as.list(ele),collapse="-")
})
availTCGASamples <- intersect(TCGA.samplesWithExprs,TCGA.clinical.samples)
availTCGAExprs <- expr[,1:544][,names(TCGA.samplesWithExprs)[match(availTCGASamples,
                                                                   TCGA.samplesWithExprs)]]
availTCGALabels <- labels[1:544][colnames(availTCGAExprs)]
availTCGASurvival <- TCGA.clinical.data['patient.days_to_death',
                                        match(availTCGASamples,TCGA.clinical.samples)]
cmpl <- complete.cases(t(availTCGASurvival))
availTCGAExprs <- availTCGAExprs[,cmpl]
annot.df <- cbind(label=(availTCGALabels[cmpl]),
                  survival=as.numeric(t(availTCGASurvival[cmpl])))

# JUN and MYBL2 appear to be interesting!
TFoI <- "MYBL2"
subtype <- "Proneural"
qntl <- 0.25
scoresVec <- TCGA.prScores.ess.pred.mat[TFoI,intersect(colnames(TCGA.prScores.ess.pred.mat),
                                                       rownames(annot.df)[annot.df[,"label"]==subtype])]
scores.qntl <- quantile(scoresVec,seq(0,1,0.25))
samplesToPlot <- names(scoresVec[scoresVec < scores.qntl[5]])
plot(as.numeric(availTCGAExprs[TFoI,samplesToPlot]),
     log2(as.numeric(annot.df[samplesToPlot,"survival"])),pch=20)
cor(as.numeric(availTCGAExprs[TFoI,samplesToPlot]),
    log2(as.numeric(annot.df[samplesToPlot,"survival"])),method="spearman")

# Prepare data for Kaplan-Meier curves in Prism
survival.splits <- lapply(subtypes,function(s){
  samples <- rownames(annot.df)[annot.df[,"label"]==s]
  exprsVec <- data.matrix(availTCGAExprs)[TFoI,samples]
  names(exprsVec) <- samples
  exprs.qntl <- quantile(exprsVec,c(0.5,0.5))
  TFoI.high.samples <- names(exprsVec)[exprsVec >= exprs.qntl[2]]
  TFoI.low.samples <- names(exprsVec)[exprsVec < exprs.qntl[1]]
  survival.data.highExprs <- as.numeric(annot.df[TFoI.high.samples,"survival"])
  survival.data.lowExprs <- as.numeric(annot.df[TFoI.low.samples,"survival"])
  fn.high <- paste0("data/perturbation/",s,"_",TFoI,"highExprs.survival.txt")
  write.table(survival.data.highExprs,fn.high,sep='\t',row.names=F,col.names=F,quote=F)
  fn.low <- paste0("data/perturbation/",s,"_",TFoI,"lowExprs.survival.txt")
  write.table(survival.data.lowExprs,fn.low,sep='\t',row.names=F,col.names=F,quote=F)
})

s <- "Mesenchymal"
plot(as.numeric(availTCGAExprs["NFE2",rownames(annot.df)[annot.df[,"label"]==s]]),
     log2(as.numeric(annot.df[annot.df[,"label"]==s,"survival"])),pch=20)
cor(as.numeric(availTCGAExprs["NFE2",rownames(annot.df)[annot.df[,"label"]==s]]),
    log2(as.numeric(annot.df[annot.df[,"label"]==s,"survival"])),method="spearman")

cor(as.numeric(availTCGAExprs["ERG",rownames(annot.df)[annot.df[,"label"]=="Classical"]]),
    log2(as.numeric(annot.df[annot.df[,"label"]=="Classical","survival"])),method="spearman")

plot(as.numeric(availTCGAExprs["TFAP2A",rownames(annot.df)[annot.df[,"label"]=="Proneural"]]),
     log2(as.numeric(annot.df[annot.df[,"label"]=="Proneural","survival"])))

plot(as.numeric(availTCGAExprs["TFAP2A",rownames(annot.df)[annot.df[,"label"]=="Proneural"]]),
     log2(as.numeric(annot.df[annot.df[,"label"]=="Proneural","survival"])))

# MGA is interesting-ish.
TF <- "MGA"
subtype <- "Classical"
cor(as.numeric(availTCGAExprs[TF,rownames(annot.df)[annot.df[,"label"]==subtype]]),
    log2(as.numeric(annot.df[annot.df[,"label"]==subtype,"survival"])),method="spearman")
plot(as.numeric(availTCGAExprs[TF,rownames(annot.df)[annot.df[,"label"]==subtype]]),
     log2(as.numeric(annot.df[annot.df[,"label"]==subtype,"survival"])),pch=16,
     xlim=c(min(availTCGAExprs[TF,]),max(availTCGAExprs[TF,])),
     ylim=c(log2(min(as.numeric(annot.df[,"survival"]))),
            log2(max(as.numeric(annot.df[,"survival"])))))

saveTFScoresSurvival <- function(TF,subtype) {
  x <- as.numeric(availTCGAExprs[TF,rownames(annot.df)[annot.df[,"label"]==subtype]])
  y <- log2(as.numeric(annot.df[annot.df[,"label"]==subtype,"survival"]))
  fn <- paste0("analyses/0_0.1/",TF,"_",subtype,"_score_vs_survival.txt")
  write.table(cbind(x,y),fn,row.names = F,col.names = F,quote=F,sep='\t')
}
lapply(subtypes,saveTFScoresSurvival,TF="MYBL2")

# FOSL2 is interesting-ish.
TF <- "FOSL2"
subtype <- "Neural"
cor(as.numeric(availTCGAExprs[TF,rownames(annot.df)[annot.df[,"label"]==subtype]]),
    log2(as.numeric(annot.df[annot.df[,"label"]==subtype,"survival"])),method="spearman")
plot(as.numeric(availTCGAExprs[TF,rownames(annot.df)[annot.df[,"label"]==subtype]]),
     log2(as.numeric(annot.df[annot.df[,"label"]==subtype,"survival"])),pch=16,
     xlim=c(min(availTCGAExprs[TF,]),max(availTCGAExprs[TF,])),
     ylim=c(log2(min(as.numeric(annot.df[,"survival"]))),
            log2(max(as.numeric(annot.df[,"survival"])))))

# TCF4 is interesting.
TF <- "TCF4"
subtype <- "Mesenchymal"
cor(as.numeric(availTCGAExprs[TF,rownames(annot.df)[annot.df[,"label"]==subtype]]),
    log2(as.numeric(annot.df[annot.df[,"label"]==subtype,"survival"])),method="spearman")
plot(as.numeric(availTCGAExprs[TF,rownames(annot.df)[annot.df[,"label"]==subtype]]),
     log2(as.numeric(annot.df[annot.df[,"label"]==subtype,"survival"])),pch=16,
     xlim=c(min(availTCGAExprs[TF,]),max(availTCGAExprs[TF,])),
     ylim=c(log2(min(as.numeric(annot.df[,"survival"]))),
            log2(max(as.numeric(annot.df[,"survival"])))))

# GFI1 is interesting.
TF <- "GFI1"
subtype <- "Mesenchymal"
cor(as.numeric(availTCGAExprs[TF,rownames(annot.df)[annot.df[,"label"]==subtype]]),
    log2(as.numeric(annot.df[annot.df[,"label"]==subtype,"survival"])),method="spearman")
plot(as.numeric(availTCGAExprs[TF,rownames(annot.df)[annot.df[,"label"]==subtype]]),
     log2(as.numeric(annot.df[annot.df[,"label"]==subtype,"survival"])),pch=16)

# FOXF1 is interesting.
TF <- "FOXF1"
subtype <- "Mesenchymal"
cor.test(as.numeric(availTCGAExprs[TF,rownames(annot.df)[annot.df[,"label"]==subtype]]),
    log2(as.numeric(annot.df[annot.df[,"label"]==subtype,"survival"])),method="spearman")
plot(as.numeric(availTCGAExprs[TF,rownames(annot.df)[annot.df[,"label"]==subtype]]),
     log2(as.numeric(annot.df[annot.df[,"label"]==subtype,"survival"])),pch=16)


# Get results from RANSAC masking
# Plot figure 5d
getRANSACTFScoresSurvival <- function(TF,subtype) {
  fn <- paste0("analyses/0_0.1/",TF,"_",subtype,"_corWithSurvival_RANSAC.txt")
  temp <- data.matrix(read.delim(fn,header=F,check.names=F))
  data.frame(x=temp[,1],y=temp[,2])
}
RANSAC_TFSurvival <- lapply(subtypes,getRANSACTFScoresSurvival,TF='MYBL2')
RANSAC_regression <- lapply(RANSAC_TFSurvival,function(df){
  lm(y~x,data=df)
})

names(RANSAC_TFSurvival) <- names(RANSAC_regression) <- subtypes

availTCGAExprs.bySubtype <- lapply(subtypes,function(s){
  availTCGAExprs[,rownames(annot.df)[annot.df[,"label"]==s]]
})
names(availTCGAExprs.bySubtype) <- subtypes

TF <- "MYBL2"
subtype <- 'Mesenchymal'
subtypePalettes <- c('Reds','Blues','Greens','Oranges')
subtypeColors <- c("brown3","cornflowerblue","mediumseagreen","goldenrod1" )
names(subtypeColors) <- subtypes
crp <- colorRampPalette(c("white",subtypePalettes[subtype]))
names(subtypePalettes) <- subtypes
xlower <- min(as.numeric(availTCGAExprs[TF,]))
xupper <- max(as.numeric(availTCGAExprs[TF,]))
ylower <- min(log2(as.numeric(annot.df[,"survival"])))
yupper <- max(log2(as.numeric(annot.df[,"survival"])))
data.df <- data.frame(x=as.numeric(availTCGAExprs.bySubtype[[subtype]][TF,]),
                      y=log2(as.numeric(annot.df[annot.df[,"label"]==subtype,"survival"])))
ggplot(data.df,aes(x=x,y=y)) + stat_density_2d(aes(fill = ..density..),
                                                               geom="raster",
                                                               contour=F) +
  scale_fill_gradient(low = "white",high=subtypeColors[subtype]) + 
  scale_x_continuous(limits=c(xlower-1,xupper+1),breaks=seq(floor(xlower),ceiling(xupper),2)) +
  scale_y_continuous(limits=c(ylower-1,yupper+1),breaks=seq(floor(ylower),ceiling(yupper),2)) +
  theme(legend.position = 'none',panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.background = element_blank(),
        panel.border = element_rect(fill = NA,size = 2)) + geom_point(color='black') + 
  geom_abline(intercept = RANSAC_regression[[subtype]]$coefficients[1],
              slope = RANSAC_regression[[subtype]]$coefficients[2],
              size = 1.25,linetype = 'dashed')



subtype <- "Mesenchymal"
library(spatstat)
subtypeColors <- c("brown3","cornflowerblue","mediumseagreen","goldenrod1")
names(subtypeColors) <- subtypes
#pdf(file = paste0("analyses/0_0.1/",TFoI,"_",tgt,"_",subtype,".pdf"),width = 10,height = 4)
par(mar = c(5, 5, 1, 1))
pppo=ppp(x=as.numeric(availTCGAExprs.bySubtype[[subtype]][TF,]),
         y=log2(as.numeric(annot.df[annot.df[,"label"]==subtype,"survival"])),
         window = owin(c(min(availTCGAExprs[TF,]),max(availTCGAExprs[TF,])),
                       c(log2(min(as.numeric(annot.df[,"survival"]))),
                         log2(max(as.numeric(annot.df[,"survival"]))))),
         xrange = c(min(availTCGAExprs[TF,]),max(availTCGAExprs[TF,])),
         yrange = c(log2(min(as.numeric(annot.df[,"survival"]))),
                    log2(max(as.numeric(annot.df[,"survival"])))))
den=density(pppo,kernel="gaussian",edge=T,diggle=T,adjust=0.6)
plot(den,main='Survival vs TF expression',col=colorRampPalette(c("white",subtypeColors[subtype]))(10),
     xlim=c(min(availTCGAExprs[TF,]),
            max(availTCGAExprs[TF,])),
     ylim=c(log2(min(as.numeric(annot.df[,"survival"]))),
            log2(max(as.numeric(annot.df[,"survival"])))),xlab=TF,ylab='log2 survival (days)')
points(availTCGAExprs.bySubtype[[subtype]][TF,],
       log2(as.numeric(annot.df[annot.df[,"label"]==subtype,"survival"])),
       xlim = (c(min(availTCGAExprs[TF,]),max(availTCGAExprs[TF,]))),
       ylim=c(log2(min(as.numeric(annot.df[,"survival"]))),
              log2(max(as.numeric(annot.df[,"survival"])))),pch=20)

subtype <- "Mesenchymal"

# HLX is interesting.
TF <- "HLX"
subtype <- "Mesenchymal"
cor(as.numeric(availTCGAExprs[TF,rownames(annot.df)[annot.df[,"label"]==subtype]]),
    log2(as.numeric(annot.df[annot.df[,"label"]==subtype,"survival"])),method="spearman")
plot(as.numeric(availTCGAExprs[TF,rownames(annot.df)[annot.df[,"label"]==subtype]]),
     log2(as.numeric(annot.df[annot.df[,"label"]==subtype,"survival"])),pch=16,
     xlim=c(min(availTCGAExprs[TF,]),max(availTCGAExprs[TF,])),
     ylim=c(log2(min(as.numeric(annot.df[,"survival"]))),
            log2(max(as.numeric(annot.df[,"survival"])))))




# Cross-check with clinical data (survival etc.)

TCGA.clinical.data <- read.delim("data/GBM.clin.merged.txt",header=F,row.names=1,
                                 check.names=F)
TCGA.clinical.samples <- toupper(as.character(TCGA.clinical.data["patient.bcr_patient_barcode",]))
names(TCGA.clinical.samples) <- as.character(TCGA.clinical.data["patient.bcr_patient_barcode",])
TCGA.samplesWithExprs <- sapply(names(labels[1:544]),function(x){
  ele <- strsplit(x,".",fixed=T)[[1]][1:3]
  paste(as.list(ele),collapse="-")
})
availTCGASamples <- intersect(TCGA.samplesWithExprs,TCGA.clinical.samples)
availTCGAExprs <- expr[,1:544][,names(TCGA.samplesWithExprs)[match(availTCGASamples,
                                                                   TCGA.samplesWithExprs)]]
availTCGALabels <- labels[1:544][colnames(availTCGAExprs)]
availTCGASurvival <- TCGA.clinical.data['patient.days_to_death',
                                        match(availTCGASamples,TCGA.clinical.samples)]
cmpl <- complete.cases(t(availTCGASurvival))
availTCGAExprs <- availTCGAExprs[,cmpl]
annot.df <- cbind(label=(availTCGALabels[cmpl]),
                  survival=as.numeric(t(availTCGASurvival[cmpl])))
plot(as.numeric(availTCGAExprs["TFAP2A",rownames(annot.df)[annot.df[,"label"]=="Proneural"]]),
     log2(as.numeric(annot.df[annot.df[,"label"]=="Proneural","survival"])),pch=20)
cor(as.numeric(availTCGAExprs["TFAP2A",rownames(annot.df)[annot.df[,"label"]=="Mesenchymal"]]),
     log2(as.numeric(annot.df[annot.df[,"label"]=="Mesenchymal","survival"])),method="spearman")

cor(as.numeric(availTCGAExprs["ERG",rownames(annot.df)[annot.df[,"label"]=="Classical"]]),
     log2(as.numeric(annot.df[annot.df[,"label"]=="Classical","survival"])),method="spearman")

plot(as.numeric(availTCGAExprs["TFAP2A",rownames(annot.df)[annot.df[,"label"]=="Proneural"]]),
    log2(as.numeric(annot.df[annot.df[,"label"]=="Proneural","survival"])))

plot(as.numeric(availTCGAExprs["TFAP2A",rownames(annot.df)[annot.df[,"label"]=="Proneural"]]),
     log2(as.numeric(annot.df[annot.df[,"label"]=="Proneural","survival"])))

