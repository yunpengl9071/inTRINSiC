rm(list=ls())

setwd('~/Data/multilayerNetwork/')
source("src/coreFunctions.R")

subtypes <- c("Classical","Neural","Proneural","Mesenchymal")
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

# lambda <- "0_0.1"
# capacities.list <- lapply(subtypes,function(x){
#   fn <- paste0("data/regrOutput/stability/",x,"_TFGeneCapacities_ri_nomi__lambda",lambda,"_rs_896454.txt")
#   temp <- data.matrix(read.delim(fn,header=T,row.names=1,check.names=F))
#   constid <- which(rownames(temp)=="const")
#   temp[-constid,]
# })
# names(capacities.list) <- subtypes

unionNet <- read.table("data/piq/TFs_for_genes.txt",sep="\t",header=F)
allTFs <- unique(as.character(unionNet[,2]))
allGenes <- unique(as.character(unionNet[,1]))
unionNet.byGene <- split(unionNet,f=unionNet[,1])
regulators <- lapply(unionNet.byGene,function(x){
  as.character(x[,2])
})

allTFInteractions.tuples <- read.delim("data/knownTFInteractionTuples.BIOGRID.txt",
                                       header=F,check.names = F)
allTFInteractions.tuples <- as.character(allTFInteractions.tuples[,1])


l2list <- c(0,0.001,0.01,0.1,1,10)
lambdas <- sapply(l2list,function(x){paste0(as.character(x),"_0")})

# Nonlinear models
capacities.list.byLambda <- lapply(lambdas,function(lambda){
  cap.list <- lapply(subtypes,function(x){
    fn <- paste0("data/regrOutput/paramSweep/nonlinear/",x,
                 "_TFGeneCapacities_ri_nomi__lambda",
                 lambda,".txt")
    temp <- data.matrix(read.delim(fn,header=T,row.names=1,check.names=F))
    constid <- which(rownames(temp)=="const")
    temp[-constid,]
  })
  names(cap.list) <- subtypes
  cat(".")
  cap.list
})
names(capacities.list.byLambda) <- lambdas

# log2F.threshold <- 5
# log2F.threshold.list <- lapply(capacities.list.byLambda,function(capmat){
#   sapply(capmat,function(x){
#     sum(abs(log2(x))>log2F.threshold)/sum(x!=1)
#   })
# })
# 
source("src/computeTFCor.R")

# geneTFCor.list.byLambda <- lapply(l2list,function(x){
#   lambda <- paste0(as.character(x),"_0")
#   capmats <- capacities.list.byLambda[[lambda]]
#   cat(lambda)
#   cat("\n")
#   computeTFCor(capmats,lambda,linear=F,log2FThres=log2F.threshold)
# })
# names(geneTFCor.list.byLambda) <- lambdas

# Should try RANSAC in scikit-learn!!
# In RANSAC there is no need to put a prior threshold on F values
# test.RANSAC <- data.matrix(read.delim("data/Classical_cormat_0.1_0.txt",header=F))
# lambda <- "0.1_0"
# capacities.list <- lapply(subtypes,function(x){
#   fn <- paste0("data/regrOutput/",x,"_TFGeneCapacities_ri_nomi__lambda",lambda,".txt")
#   temp <- data.matrix(read.delim(fn,header=T,row.names=1,check.names=F))
#   constid <- which(rownames(temp)=="const")
#   temp[-constid,]
# })
# names(capacities.list) <- subtypes
# 
# rownames(test.RANSAC) <- colnames(test.RANSAC) <- rownames(capacities.list[[1]])
# 
# examplePairs <- list()
# exampleTuples <- c()
# for (i in 1:length(allTFs)) {
#   for (j in 1:length(allTFs)) {
#     corVec <- test.RANSAC[i,j]
#     if (!is.na(corVec) & abs(corVec) > 0.6) {
#       examplePairs[[paste0(rownames(test.RANSAC)[i],"_",rownames(test.RANSAC)[j])]] <- corVec
#       exampleTuples[[paste0(rownames(test.RANSAC)[i],"_",rownames(test.RANSAC)[j])]] <- paste0(rownames(test.RANSAC)[i],"_",rownames(test.RANSAC)[j])
#     }
#   }
# }
# examplePairs <- do.call(rbind.data.frame,examplePairs)
# colnames(examplePairs) <- 'Classical'
# rownames(examplePairs) <- exampleTuples

# library(doParallel)
# library(foreach)
# cl <- parallel::makeCluster(6,outfile="")
# doParallel::registerDoParallel(cl)
# 
# geneTFCor.list.byLambda <- foreach(i = c(1:length(lambdas))) %dopar% {
#   lambda <- lambdas[i]
#   #capmats <- capacities.list.byLambda[[lambda]]
#   cat(lambda)
#   cat("\n")
#   for (subtype in subtypes) {
#     #Sys.setenv(PATH="/Users/yunpengl/anaconda2/bin/python")
#     cmd <- paste0("/Users/yunpengl/anaconda2/bin/python src/computeTFCor_RANSAC.py ", subtype, " ", lambda)
#     system(cmd)
#   }
#   #computeTFCor(capmats,lambda,linear=F,log2FThres = log2F.threshold)
# }
# stopCluster(cl)

# log2F.threshold <- 5
# geneTFCor.list.byLambda <- lapply(lambdas,function(lambda){
#   capmats <- lapply(c(1:length(subtypes)),function(i){
#     fn <- paste0("analyses/",subtypes[i],"_TF_TF_cormat_DFNet_log2FThres_",
#                  as.character(log2F.threshold),"lambda_",as.character(lambda),".txt")
#     geneTF.coreg.mat <- data.matrix(read.delim(fn,header=T,row.names=1,check.names=F))
#   })
#   names(capmats) <- subtypes
#   cat(".")
#   capmats
# })
# names(geneTFCor.list.byLambda) <- lambdas

rs <- 7792275
#log2F.threshold <- 5
geneTFCor.list.byLambda <- lapply(lambdas,function(lambda){
  capmats <- lapply(c(1:length(subtypes)),function(i){
    fn <- paste0("data/",subtypes[i],"_cormat_",as.character(lambda),"_rs_",as.character(rs),".txt")
    geneTF.coreg.mat <- data.matrix(read.delim(fn,header=F,check.names=F))
    geneTF.coreg.mat[is.nan(geneTF.coreg.mat)] <- 0
    rownames(geneTF.coreg.mat) <- colnames(geneTF.coreg.mat) <- rownames(capacities.list.byLambda[[1]][[1]])
    geneTF.coreg.mat
  })
  names(capmats) <- subtypes
  cat(".")
  capmats
})
names(geneTFCor.list.byLambda) <- lambdas

library(igraph)
allTFPairs <- c()
count <- 1
for (TF1 in allTFs) {
  for (TF2 in allTFs) {
    allTFPairs[count] <- paste0(TF1,"_",TF2)
    count <- count + 1
  }
}
TFCor.tuples.lists.byLambda <- lapply(geneTFCor.list.byLambda,function(geneRegs.TFs.coreg){
  cat("\n")
  TFCor.tuples.list <- lapply(geneRegs.TFs.coreg,function(gtc){
    g <- graph.adjacency(gtc,mode=c("directed"),weighted = T)
    el <- get.edgelist(g)
    weights <- get.edge.attribute(graph = g,name = "weight")
    el1 <- apply(el[,1:2],1,function(x){
      paste0(x[1],"_",x[2])
    })
    el2 <- apply(el[,1:2],1,function(x){
      paste0(x[2],"_",x[1])
    })
    names(weights) <- el1
    cat(".")
    allweights <- rep(0,length(allTFPairs))
    names(allweights) <- allTFPairs
    allweights[el1] <- weights
    allweights[el2] <- weights
    allweights
  })
})
names(TFCor.tuples.lists.byLambda) <- lambdas

temp <- lapply(c(1:6),function(i){
  l <- lambdas[[i]]
  for (x in 1:length(subtypes)) {
    subtype <- subtypes[x]
    tlist <- TFCor.tuples.lists.byLambda[[i]][[x]]
    fn <- paste0("data/TFCor/",subtype,"_TFCorTuples_",l,"_rs_",as.character(rs),".txt")
    write.table(data.frame(tlist),fn,sep="\t",col.names=F,quote=F)
  }
})

precision <- c()
recall <- c()
auc <- c()
fisher.p <- c()
library(zoo)
for (lambda in lambdas) {
  TFCor.tuples.list <- TFCor.tuples.lists.byLambda[[lambda]]
  percent.seq <- seq(0.01,1,0.01)
  commonTFCor.tuples <- Reduce("intersect",lapply(TFCor.tuples.list,names))
  commonTFCor.tuples.cor <- lapply(TFCor.tuples.list,function(x){
    x[commonTFCor.tuples]
  })
  commonTFCor.tuples.cor <- data.matrix(do.call(cbind.data.frame,commonTFCor.tuples.cor))
  commonTFCor.tuples.avgCor <- rowMeans(commonTFCor.tuples.cor)
  commonTFCor.tuples.avgCor <- rev(sort(commonTFCor.tuples.avgCor))
  prec <- c()
  rec <- c()
  TFCor.tuples.hits.list <- c()
  for (p in percent.seq) {
    hits <- intersect(names(commonTFCor.tuples.avgCor[1:ceiling(p*length(commonTFCor.tuples.avgCor))]),
                      allTFInteractions.tuples)
    prec <- c(prec,length(hits)/ceiling(p*length(commonTFCor.tuples.avgCor)))
    rec <- c(rec,length(hits)/length(allTFInteractions.tuples))
    TFCor.tuples.hits.list <- c(TFCor.tuples.hits.list,
                                length(hits)/ceiling(p*length(commonTFCor.tuples.avgCor)))
  }
  plot(rec,prec,type="l",lwd=3)
  # How to tell if the following QQ-like plot shows significant
  # deviation from random subsetting of the list of known interactors?
  # plot(percent.seq,TFCor.tuples.hits.list,type="l",lwd=3,
  #      xlab="percent top correlated TF pairs",
  #      ylab="proportion of TF pairs with known interactions")
  validpct <- which(complete.cases(TFCor.tuples.hits.list))
  TFCor.tuples.hits.list <- TFCor.tuples.hits.list[validpct]
  p = percent.seq[validpct][which(TFCor.tuples.hits.list==max(TFCor.tuples.hits.list))[1]]
  #p = percent.seq[which(TFCor.tuples.hits.list==max(TFCor.tuples.hits.list))[1]]
  #p <- 1
  hits <- intersect(names(commonTFCor.tuples.avgCor[1:ceiling(p*length(commonTFCor.tuples.avgCor))]),
                    allTFInteractions.tuples)
  x = matrix(c(length(hits), # predicted pairs that are also in 'truth' set
               ceiling(p*length(commonTFCor.tuples.avgCor)), # all predicted pairs
               length(allTFInteractions.tuples), # all 'true' pairs
               length(allTFs)^2),  # pairs to predict from
             2,2)
  # x = matrix(c(length(hits), # predicted pairs that are also in 'truth' set
  #              ceiling(p*length(commonTFCor.tuples.avgCor)), # all predicted pairs
  #              length(intersect(names(commonTFCor.tuples.avgCor),allTFInteractions.tuples)), # all 'true' pairs
  #              length(commonTFCor.tuples.avgCor)),  # pairs to predict from
  #            2,2)
  print(fisher.test(x))
  precision[lambda] <- length(hits)/ceiling(p*length(commonTFCor.tuples.avgCor))
  recall[lambda] <- length(hits)/length(allTFInteractions.tuples)
  #recall[lambda] <- length(hits)/length(commonTFCor.tuples.avgCor)
  AUC <- sum(diff(rec[complete.cases(prec)])*rollmean(prec[complete.cases(prec)],2))
  print(AUC)
  auc[lambda] <- AUC
  fisher.p[lambda] <- fisher.test((x))$p.value
}
# For linear models only
# auc["10_0"] <- 0
# fisher.p["10_0"] <- 1
# lambda = 0.1 seems the best?
plot(auc*(-log(fisher.p)),type="l",lwd=3)
plot(auc,type="l",lwd=3)
stats.df <- cbind(fisher.p,auc,-log10(fisher.p)*auc)
colnames(stats.df) <- c("fisher.p","auc","-log10(fisher.p)*auc")
write.table(stats.df,paste0("analyses/interactionEnrichmentSummary_RANSAC","_rs_",
                            as.character(rs),".txt"),
            sep="\t",quote=F)
#lambda <- lambdas[5]
#plot(precision,recall,xlim=c(0,0.05),ylim=c(0,0.02),pch=16)
#lines(precision[order(precision)],recall[order(precision)])














##########################################################################
# Linear models
capacities.list.byLambda <- lapply(lambdas,function(lambda){
  cap.list <- lapply(subtypes,function(x){
    fn <- paste0("data/regrOutput/paramSweep/linear/",x,
                 "_TFGeneCapacities_ri_nomi__lambda",
                 lambda,".txt")
    temp <- data.matrix(read.delim(fn,header=T,row.names=1,check.names=F))
    constid <- which(rownames(temp)=="const")
    temp[-constid,]
  })
  names(cap.list) <- subtypes
  cat(".")
  cap.list
})
names(capacities.list.byLambda) <- lambdas

source("src/computeTFCor.R")

# library(doParallel)
# library(foreach)
# cl <- parallel::makeCluster(6,outfile="")
# doParallel::registerDoParallel(cl)
# 
# geneTFCor.list.byLambda <- foreach(i = c(1:length(lambdas))) %dopar% {
#   lambda <- lambdas[i]
#   FcutoffVec <- log2F.threshold.list[[i]]
#   capmats <- capacities.list.byLambda[[lambda]]
#   cat(lambda)
#   cat("\n")
#   computeTFCor(capmats,lambda,linear=T,FcutoffVec=log2F.threshold)
# }
# stopCluster(cl)

# Update: calculate TF interaction pair enrichment with the 
# optimal lambdas for expression prediction (see prism files)
# in each type of model, and match the top percent of F values
# taken using the appropriate F cutoffs
# For linear, lambda2 = 0.001 looks good
# For nonlinear, lambda2 = 0.01 looks good
# Match F cutoff percentage to that used in nonlinear lambda2=0.01
# lambda <- lambdas[2]
# FcutoffVec <- log2F.threshold.list$`0.01_0`
# capmats <- capacities.list.byLambda[[lambda]]
# geneTFCor.list.byLambda <- list(computeTFCor(capmats,lambda,linear=T,FcutoffVec=FcutoffVec))
# l <- lambdas[[2]]
# for (x in 1:length(subtypes)) {
#   subtype <- subtypes[x]
#   tlist <- TFCor.tuples.lists.byLambda[[1]][[x]]
#   fn <- paste0("data/TFCor/",subtype,"_linear_TFCorTuples_",l,".txt")
#   write.table(data.frame(tlist),fn,sep="\t",col.names=F,quote=F)
# }
# geneTFCor.list.byLambda <- lapply(l2list,function(x){
#   lambda <- paste0(as.character(x),"_0")
#   capmats <- capacities.list.byLambda[[lambda]]
#   cat(lambda)
#   cat("\n")
#   computeTFCor(capmats,lambda,linear=T)
# })
# names(geneTFCor.list.byLambda) <- lambdas

#log2F.threshold <- 0
#log2F.threshold <- 5
geneTFCor.list.byLambda <- lapply(lambdas,function(lambda){
  capmats <- lapply(c(1:length(subtypes)),function(i){
    fn <- paste0("data/",subtypes[i],"_cormat_",as.character(lambda),"_lin.txt")
    geneTF.coreg.mat <- data.matrix(read.delim(fn,header=F,check.names=F))
    geneTF.coreg.mat[is.nan(geneTF.coreg.mat)] <- 0
    rownames(geneTF.coreg.mat) <- colnames(geneTF.coreg.mat) <- rownames(capacities.list.byLambda[[1]][[1]])
    geneTF.coreg.mat
  })
  names(capmats) <- subtypes
  cat(".")
  capmats
})
names(geneTFCor.list.byLambda) <- lambdas

library(igraph)

TFCor.tuples.lists.byLambda <- lapply(geneTFCor.list.byLambda,function(geneRegs.TFs.coreg){
  cat("\n")
  TFCor.tuples.list <- lapply(geneRegs.TFs.coreg,function(gtc){
    g <- graph.adjacency(gtc,mode=c("directed"),weighted = T)
    el <- get.edgelist(g)
    weights <- get.edge.attribute(graph = g,name = "weight")
    el1 <- apply(el[,1:2],1,function(x){
      paste0(x[1],"_",x[2])
    })
    el2 <- apply(el[,1:2],1,function(x){
      paste0(x[2],"_",x[1])
    })
    names(weights) <- el1
    cat(".")
    allweights <- rep(0,length(allTFPairs))
    names(allweights) <- allTFPairs
    allweights[el1] <- weights
    allweights[el2] <- weights
    allweights
  })
})
names(TFCor.tuples.lists.byLambda) <- lambdas

temp <- lapply(c(1:6),function(i){
  l <- lambdas[[i]]
  for (x in 1:length(subtypes)) {
    subtype <- subtypes[x]
    tlist <- TFCor.tuples.lists.byLambda[[i]][[x]]
    fn <- paste0("data/TFCor/",subtype,"_linear_TFCorTuples_",l,".txt")
    write.table(data.frame(tlist),fn,sep="\t",col.names=F,quote=F)
  }
})

precision <- c()
recall <- c()
auc <- c()
fisher.p <- c()
library(zoo)
for (lambda in lambdas) {
  TFCor.tuples.list <- TFCor.tuples.lists.byLambda[[lambda]]
  percent.seq <- seq(0.01,1,0.01)
  commonTFCor.tuples <- Reduce("intersect",lapply(TFCor.tuples.list,names))
  commonTFCor.tuples.cor <- lapply(TFCor.tuples.list,function(x){
    x[commonTFCor.tuples]
  })
  commonTFCor.tuples.cor <- data.matrix(do.call(cbind.data.frame,commonTFCor.tuples.cor))
  commonTFCor.tuples.avgCor <- rowMeans(commonTFCor.tuples.cor)
  commonTFCor.tuples.avgCor <- rev(sort(commonTFCor.tuples.avgCor))
  prec <- c()
  rec <- c()
  TFCor.tuples.hits.list <- c()
  for (p in percent.seq) {
    hits <- intersect(names(commonTFCor.tuples.avgCor[1:round(p*length(commonTFCor.tuples.avgCor))]),
                      allTFInteractions.tuples)
    prec <- c(prec,length(hits)/round(p*length(commonTFCor.tuples.avgCor)))
    rec <- c(rec,length(hits)/length(allTFInteractions.tuples))
    TFCor.tuples.hits.list <- c(TFCor.tuples.hits.list,
                                length(hits)/round(p*length(commonTFCor.tuples.avgCor)))
  }
  plot(rec,prec,type="l",lwd=3)
  # How to tell if the following QQ-like plot shows significant
  # deviation from random subsetting of the list of known interactors?
  # plot(percent.seq,TFCor.tuples.hits.list,type="l",lwd=3,
  #      xlab="percent top correlated TF pairs",
  #      ylab="proportion of TF pairs with known interactions")
  validpct <- which(complete.cases(TFCor.tuples.hits.list))
  TFCor.tuples.hits.list <- TFCor.tuples.hits.list[validpct]
  p = percent.seq[validpct][which(TFCor.tuples.hits.list==max(TFCor.tuples.hits.list))[1]]
  #p = percent.seq[which(TFCor.tuples.hits.list==max(TFCor.tuples.hits.list))[1]]
  #p <- 1
  hits <- intersect(names(commonTFCor.tuples.avgCor[1:round(p*length(commonTFCor.tuples.avgCor))]),
                    allTFInteractions.tuples)
  x = matrix(c(length(hits), # predicted pairs that are also in 'truth' set
               round(p*length(commonTFCor.tuples.avgCor)), # all predicted pairs
               length(allTFInteractions.tuples), # all 'true' pairs
               length(allTFs)^2),  # pairs to predict from
             2,2)
  # x = matrix(c(length(hits), # predicted pairs that are also in 'truth' set
  #              round(p*length(commonTFCor.tuples.avgCor)), # all predicted pairs
  #              length(intersect(names(commonTFCor.tuples.avgCor),allTFInteractions.tuples)), # all 'true' pairs
  #              length(commonTFCor.tuples.avgCor)),  # pairs to predict from
  #            2,2)
  print(fisher.test(x))
  precision[lambda] <- length(hits)/round(p*length(commonTFCor.tuples.avgCor))
  recall[lambda] <- length(hits)/length(allTFInteractions.tuples)
  AUC <- sum(diff(rec[complete.cases(prec)])*rollmean(prec[complete.cases(prec)],2))
  print(AUC)
  auc[lambda] <- AUC
  fisher.p[lambda] <- fisher.test((x))$p.value
}
# For linear models only
# auc["10_0"] <- 0
# fisher.p["10_0"] <- 1
# lambda = 0.1 seems the best?
plot(auc*(-log(fisher.p)),type="l",lwd=3)
stats.df <- cbind(fisher.p,auc,-log10(fisher.p)*auc)
colnames(stats.df) <- c("fisher.p","auc","-log10(fisher.p)*auc")
write.table(stats.df,"analyses/interactionEnrichmentSummary_linear_RANSAC.txt",
            sep="\t",quote=F)




# TFCor.tuples.lists.byLambda <- lapply(geneTFCor.list.byLambda,function(geneRegs.TFs.coreg){
#   TFCor.tuples.list <- lapply(geneRegs.TFs.coreg,function(x){
#     TFCor.tuples <- rep(0,sum(x!=0))
#     nonZeroIdx <- which(x!=0)
#     nr <- nrow(x)
#     tfs <- rownames(x)
#     count <- 0
#     cat("\n")
#     for (idx in nonZeroIdx) {
#       count <- count + 1
#       if (count %% 1000 == 0) {
#         cat(".")
#       }
#       temp <- idx %% nr
#       if (temp == 0) {
#         rowInd <- nr
#         colInd <- idx %/% nr
#       } else {
#         rowInd <- temp
#         colInd <- idx %/% nr + 1
#       }
#       if (rowInd != colInd) {
#         tfi <- tfs[rowInd]
#         tfj <- tfs[colInd]
#         tup <- paste0(tfi,"_",tfj)
#         TFCor.tuples[tup] <- x[rowInd,colInd]
#       }
#     }
#     cat("\n")
#     TFCor.tuples
#   })
#   cat("\n")
#   lapply(TFCor.tuples.list,function(x){
#     x[x!=0]
#   })
# })
# 
# 
