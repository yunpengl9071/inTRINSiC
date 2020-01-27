#rm(list=ls())

setwd('~/Data/multilayerNetwork/')
source("src/coreFunctions.R")

subtypes <- c("Classical","Neural","Proneural","Mesenchymal")
lambda <- "0_0.1"
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

lambda <- "0.1_0"
capacities.list <- lapply(subtypes,function(x){
  fn <- paste0("data/regrOutput/paramSweep/nonlinear/",x,"_TFGeneCapacities_ri_nomi__lambda",lambda,".txt")
  temp <- data.matrix(read.delim(fn,header=T,row.names=1,check.names=F))
  constid <- which(rownames(temp)=="const")
  temp[-constid,]
})
names(capacities.list) <- subtypes


capacities.log2.list <- lapply(capacities.list,log2)


library(gplots)
myheatmapfun <- colorRampPalette((c(rgb(15,61,111,maxColorValue=255),
                                    rgb(125,182,212,maxColorValue=255),
                                    rgb(255,255,255,maxColorValue=255),
                                    rgb(234,145,118,maxColorValue=255),
                                    rgb(102,3,32,maxColorValue=255))))
#myheatmapfun <- colorRampPalette((c("white","darkblue")))

heatmap.data <- capacities.log2.list[[1]]
cap.sums <- Reduce("+",capacities.log2.list)
cap.means <- cap.sums/4
heatmap.data <- heatmap.data[rowMeans(abs(cap.means))>0.15,colMeans(abs(cap.means))>0.25]
# library(biclust)
# xxx <- biclust(2^heatmap.data,method=BCSpectral())
#cRow <- hclust(as.dist(1-cor(heatmap.data)))
#cCol <- hclust(as.dist(1-cor(t(heatmap.data))))
cRow <- kmeans(heatmap.data,centers=5)$cluster
cCol <- kmeans(t(heatmap.data),centers=10)$cluster
cRow <- names(sort(cRow))
cCol <- names(sort(cCol))
dev.off()
hmRes <- heatmap.2(heatmap.data,
                   #Colv=NA,
                   #Rowv=NA,
                   scale="none",
                   density.info="none", trace="none",
                   col=myheatmapfun(75),
                   reorderfun = function(d,w) { d },
                   #cexCol = 0.8,
                   #srtCol = 45,
                   #cexRow = 0.8,
                   #rowsep = c(1:4),
                   #colsep = c(1:4),
                   #sepwidth = c(0.02,0.02),
                   #offsetCol = -0.5,
                   #RowSideColors = rscVec,
                   dendrogram="none",
                   key.title=NA,
                   key.xlab=NA,
                   lhei=c(1,8),
                   breaks=seq(-2,2,length.out=76),
                   labRow = "",
                   labCol = "",
                   margins=c(5,5)
)
library(pheatmap)
heatmap.data <- capacities.log2.list[[4]]
cap.sums <- Reduce("+",capacities.log2.list)
cap.means <- cap.sums/4
heatmap.data <- heatmap.data[rowMeans(abs(cap.means))>0.15,colMeans(abs(cap.means))>0.25]
xxx <- pheatmap(heatmap.data,breaks=seq(-10,10,length.out=76),
         clustering_method = "ward.D",clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = myheatmapfun(75))

# rowOrder <- xxx$tree_row$labels
# colOrder <- xxx$tree_col$labels
# heatmap.data <- capacities.log2.list[[2]]
# heatmap.data <- heatmap.data[rowOrder,colOrder]
# xxx <- pheatmap(heatmap.data,breaks=seq(-10,10,length.out=76),cluster_rows = F,
#                 cluster_cols = F,
#                 color = myheatmapfun(75))

# A F value is significant only when it induces gene expression changes
# of 2 fold or greater
abslog2FCutoff <- 1
# capacities.log2.el.list <- lapply(subtypes,function(x){
#   fn <- paste0(x,'_TFGeneCapacities_ri_nomi_DHSFANTOM_lambda_',lambda,'_edgelist.txt')
#   capellist <- read.delim(fn,header=F,check.names=F)
#   capellist[,3] <- log2(as.numeric(as.character(capellist[,3])))
#   capellist
# })
# names(capacities.log2.el.list) <- subtypes
# capacities.log2.el.sigF.list <- lapply(capacities.log2.el.list,function(x){
#   sigFRows <- which(abs(as.numeric(as.character(x[,3]))) > abslog2FCutoff)
#   x[sigFRows,]
# })
# temp <- lapply(capacities.log2.list,function(x){
#   absx <- abs(x)
#   sigEL <- list()
#   for (gene in colnames(x)) {
#     for (TF in rownames(x)) {
#       if (absx[TF,gene] > abslog2FCutoff) {
#         edge <- c(TF,gene,as.character(absx[TF,gene]))
#         sigEL <- c(sigEL,edge)
#       }
#     }
#   }
# })

unionNet <- read.table("data/piq/TFs_for_genes.txt",sep="\t",header=F)
allTFs <- unique(as.character(unionNet[,2]))
allGenes <- unique(as.character(unionNet[,1]))
unionNet.byGene <- split(unionNet,f=unionNet[,1])
regulators <- lapply(unionNet.byGene,function(x){
  as.character(x[,2])
})
unionNet.byTF <- split(unionNet,f=unionNet[,2])


# 0. Utility functions

library(vioplot)
compareExprs <- function(x) {
  datalist<- list("Classical"=data.matrix(expr)[x,names(commonLabels)[commonLabels=="Classical"]],
                  "Neural"=data.matrix(expr)[x,names(commonLabels)[commonLabels=="Neural"]],
                  "Proneural"=data.matrix(expr)[x,names(commonLabels)[commonLabels=="Proneural"]],
                  "Mesenchymal"=data.matrix(expr)[x,names(commonLabels)[commonLabels=="Mesenchymal"]])
  ylim.upper <- sapply(datalist,function(x){
    range(x)[2]
  })
  ylim.upper <- max(ylim.upper)
  vioplot(datalist[[1]],datalist[[2]],
          datalist[[3]],datalist[[4]],
          names=names(datalist),col="grey",ylim=c(2,ylim.upper))
}
plotTFTgt <- function(TF,tgt,subtype) {
  TF.exprs <- data.matrix(expr)[TF,names(commonLabels)[commonLabels==subtype]]
  tgt.exprs <- data.matrix(expr)[tgt,names(commonLabels)[commonLabels==subtype]]
  plot(TF.exprs,tgt.exprs,xlab=TF,ylab=tgt,main=subtype,pch=16)
}
exprs.anova <- function(x) {
  dataVec<- c(data.matrix(expr)[x,names(commonLabels)[commonLabels=="Classical"]],
              data.matrix(expr)[x,names(commonLabels)[commonLabels=="Neural"]],
              data.matrix(expr)[x,names(commonLabels)[commonLabels=="Proneural"]],
              data.matrix(expr)[x,names(commonLabels)[commonLabels=="Mesenchymal"]])
  labelVec <- c(rep("Classical",sum(commonLabels=="Classical")),
                rep("Neural",sum(commonLabels=="Neural")),
                rep("Proneural",sum(commonLabels=="Proneural")),
                rep("Mesenchymal",sum(commonLabels=="Mesenchymal")))
  exprs.df <- cbind.data.frame(labelVec,dataVec)
  summary(aov(dataVec ~ labelVec, data = exprs.df))[[1]][1,"Pr(>F)"]
}
model.pred <- function(tfvec, mivec, fvecmi, fvectf, const) {
  if (length(fvecmi)==0 | sum(fvecmi)==length(fvecmi)) {
    temp <- log2(const) + sum(log2((1 + fvectf * (tfvec^2)) / (1 + tfvec^2)))
    temp
  } else {
    temp <- log2(const) + sum(0.9 * log2((1 + fvectf * (tfvec)^2) / (1 + tfvec^2))) + sum(0.1 * log2((1 + fvecmi * (mivec)^2) / (1 + mivec^2)))
    temp
  }
}

colcol <- c(rep("red",sum(commonLabels=="Classical")),
            rep("blue",sum(commonLabels=="Neural")),
            rep("green",sum(commonLabels=="Proneural")),
            rep("orange",sum(commonLabels=="Mesenchymal")))

# 1. Compare how many edges unique to each network were selected by
# the regression algorithm - look at F values above a cutoff

log2Fcutoff <- 1
sigFValues.list <- lapply(capacities.list,function(cap){
  abs.log2.cap <- abs(log2(cap))
  allTFs <- rownames(cap)
  allGenes <- colnames(cap)
  sigEdges <- lapply(1:nrow(cap),function(i){
    sigFGenes <- allGenes[which(abs.log2.cap[i,] > log2Fcutoff)]
    sapply(sigFGenes,function(sfg){
      paste0(allTFs[i],"_",sfg)
    })
  })
  cat(".")
  unlist(sigEdges)
})
DHS.only.net <- read.table("data/piq/DHS_0.01_only_edges.txt",sep="\t",header=F)
FANTOM5.only.net <- read.table("data/piq/FANTOM5_0.1_only_edges.txt",sep="\t",header=F)
inters.net <- read.table("data/piq/DHS_0.01_FANTOM5_0.1_inters_edges.txt",sep="\t",header=F)
DHS.only.edgeset <- unique(sapply(1:nrow(DHS.only.net),function(x){
  paste0(DHS.only.net[x,1],"_",DHS.only.net[x,2])
}))
FANTOM5.only.edgeset <- unique(sapply(1:nrow(FANTOM5.only.net),function(x){
  paste0(FANTOM5.only.net[x,1],"_",FANTOM5.only.net[x,2])
}))
inters.edgeset <- unique(sapply(1:nrow(inters.net),function(x){
  paste0(inters.net[x,1],"_",inters.net[x,2])
}))
# Report some coverage stats
inters.list <- lapply(sigFValues.list,function(x){
  c(length(intersect(DHS.only.edgeset,x))/length(DHS.only.edgeset),
    length(intersect(FANTOM5.only.edgeset,x))/length(FANTOM5.only.edgeset),
    length(intersect(inters.edgeset,x))/length(inters.edgeset))
})
inters.table <- do.call(cbind.data.frame,inters.list)
rownames(inters.table) <- c("DHS only","FANTOM5 only","intersection")
write.table(inters.table,"analyses/0_0.1/DHS_FANTOM5_coveragePercentage.txt",
            sep="\t",quote=F)

# 2. Cross-validation.

set.seed(972831654)

CV.TVSamples <- lapply(c(1:length(sampleByLabels)),function(x){
  nsamples.subtype <- length(sampleByLabels[[x]])
  nsamp <- round(nsamples.subtype/CV.fold)
  sampleNames.shuffled <- sample((sampleByLabels[[x]]),size=nsamples.subtype,
                                 replace = F)
  temp <- list()
  for (i in c(1:CV.fold)) {
    start <- (i-1)*nsamp + 1
    end <- min(i*nsamp,nsamples.subtype)
    trainingSamples <- sampleNames.shuffled[-c(start:end)]
    validationSamples <- sampleNames.shuffled[c(start:end)]
    temp[[i]] <- list(trn=trainingSamples,vln=validationSamples)
  }
  temp
})

capmat.CV.list <- lapply(1:CV.fold,function(x){
  cat(".")
  capmat.list <- lapply(1:length(subtypes),function(y){
    fn <- paste0("data/regrOutput/CV_Ning/Result1509/0.1/",subtypes[y],"_TFGeneCapacities_ri_nomi__lambda0.1_0_CV_",
                 as.character(x),".txt")
    data.matrix(read.csv(fn,sep="\t",check.names=F))
  })
})

CV.error.list <- lapply(c(1:length(subtypes)),function(x){
  exprs.est <- lapply(1:CV.fold,function(i){
    validationSamples <- CV.TVSamples[[x]][[i]]$vln
    tfmat <- expr[allTFs,names(validationSamples)]
    capmat <- capmat.CV.list[[i]][[x]]
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
    colnames(exprEstimates) <- names(validationSamples)
    rownames(exprEstimates) <- colnames(capmat)
    cat(".")
    exprEstimates
  })
  cat("\n")
  exprs.est
})

tr.error.list <- lapply(c(1:length(subtypes)),function(x){
  exprs.est <- lapply(1:CV.fold,function(i){
    trainingSamples <- CV.TVSamples[[x]][[i]]$trn
    tfmat <- expr[allTFs,names(trainingSamples)]
    capmat <- capmat.CV.list[[i]][[x]]
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
    colnames(exprEstimates) <- names(trainingSamples)
    rownames(exprEstimates) <- colnames(capmat)
    cat(".")
    exprEstimates
  })
  cat("\n")
  exprs.est
})

# Use symmetric mean absolute percentage error (sMAPE)
CV.error.mat <- matrix(0,length(subtypes),CV.fold)
for (i in 1:length(subtypes)) {
  for (j in 1:CV.fold) {
    exprs.obs <- expr[rownames(CV.error.list[[i]][[j]]),
                      colnames(CV.error.list[[i]][[j]])]
    temp <- (abs(CV.error.list[[i]][[j]]) + abs(data.matrix(exprs.obs)))/2
    CV.error.mat[i,j] <- median(rowMeans(abs(CV.error.list[[i]][[j]]-data.matrix(exprs.obs))/temp))
  }
}
tr.error.mat <- matrix(0,length(subtypes),CV.fold)
for (i in 1:length(subtypes)) {
  for (j in 1:CV.fold) {
    exprs.obs <- expr[rownames(tr.error.list[[i]][[j]]),
                      colnames(tr.error.list[[i]][[j]])]
    temp <- (abs(tr.error.list[[i]][[j]]) + abs(data.matrix(exprs.obs)))/2
    tr.error.mat[i,j] <- median(rowMeans(abs(tr.error.list[[i]][[j]]-data.matrix(exprs.obs))/temp))
  }
}
write.table(t(CV.error.mat),"analyses/0_0.1/lambda_0_0.1_CV_error_mat.txt",row.names=F,
            col.names=F,sep="\t",quote=F)
write.table(t(tr.error.mat),"analyses/0_0.1/lambda_0_0.1_tr_error_mat.txt",row.names=F,
            col.names=F,sep="\t",quote=F)


# Then look at inter-subtype prediction errors
# Try using genes whose expression differ among the subtypes (ANOVA)?
# Or use expression signatures?
# diff.genes <- names(tg.anova.p.bonf)[tg.anova.p.bonf < 1E-4]
# subtypeCentroids <- read.delim("subtype_centroids.txt",header=T,row.names=1)
# diff.genes <- intersect(rownames(subtypeCentroids),colnames(capacities.list[[1]]))

capPr2mRNA.list <- lapply(subtypes,function(x){
  #fn <- paste0("data/regrOutput/stability/",x,"_TFGeneCapacities_ri_nomi__lambda",lambda,"_rs_896454.txt")
  fn <- paste0("data/regrOutput/paramSweep/nonlinear/",x,"_TFGeneCapacities_ri_nomi__lambda",lambda,".txt")
  temp <- data.matrix(read.delim(fn,header=T,row.names=1,check.names=F))
  #temp[,diff.genes]
})
names(capPr2mRNA.list) <- subtypes
predExprs.list <- list()
for (i in c(1:length(subtypes))) {
  predExprs.list[[i]] <- list()
  for (j in c(1:length(subtypes))) {
    regSubtype <- subtypes[i]
    type <- subtypes[j]
    idx <- names(commonLabels)[commonLabels==type]
    tfmat <- expr[allTFs,idx]
    capmat <- capPr2mRNA.list[[regSubtype]]
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
    predExprs.list[[i]][[j]] <- exprEstimates
  }
  cat("\n")
}
for (i in c(1:length(subtypes))) {
  for (j in c(1:length(subtypes))) {
    write.table(predExprs.list[[i]][[j]],
                paste0("analyses/0_0.1/crossSubtypeExprsPred_",subtypes[i],"_predicts_",subtypes[j],".txt"),
                sep="\t",quote=F)
  }
}

subtypeCentroids <- read.delim("data/subtype_centroids.txt",header=T,row.names=1)


sampleOrderedBySubtype <- names(c(sampleByLabels[[1]],
                                  sampleByLabels[[2]],
                                  sampleByLabels[[3]],
                                  sampleByLabels[[4]]))
tg.anova.p <- lapply(colnames(capacities.list[[1]]),function(x){
  subtypeVec <- c(rep("Classical",length(sampleByLabels$Classical)),
                  rep("Neural",length(sampleByLabels$Neural)),
                  rep("Proneural",length(sampleByLabels$Proneural)),
                  rep("Mesenchymal",length(sampleByLabels$Mesenchymal)))
  exprsVec <- as.numeric(expr[x,sampleOrderedBySubtype])
  exprs.df <- cbind.data.frame(subtypeVec,exprsVec)
  summary(aov(exprsVec ~ subtypeVec, data = exprs.df))[[1]][1,"Pr(>F)"]
})
names(tg.anova.p) <- colnames(capacities.list[[1]])
tg.anova.p <- unlist(tg.anova.p)
tg.anova.p.bonf <- p.adjust(tg.anova.p,method = "bonferroni")
subtypeMeanExprs <- lapply(subtypes,function(x){
  subtypeExprs <- commonExpr[,names(commonLabels)[commonLabels==x]]
  rowMeans(subtypeExprs)
})
subtypeMeanExprs <- do.call(cbind.data.frame,subtypeMeanExprs)
subtypeMeanExprs <- data.matrix(subtypeMeanExprs)
colnames(subtypeMeanExprs) <- subtypes
leaveOneSD <- lapply(1:length(subtypes),function(x){
  LOOExprs <- subtypeMeanExprs[,-x]
  apply(LOOExprs,1,sd)
})
leaveOneSD <- do.call(cbind.data.frame,leaveOneSD)
leaveOneSD <- data.matrix(leaveOneSD)
leaveOneSD.subtypeSig <- leaveOneSD[rownames(subtypeCentroids),]
exprsSigSubtype <- apply(leaveOneSD.subtypeSig,1,function(x){
  id <- which(x==min(x))[1]
  subtypes[id]
})
exprsSigSubtype.list <- lapply(subtypes,function(x){
  intersect(names(exprsSigSubtype)[exprsSigSubtype==x],colnames(capacities.list[[1]]))
})
tg.meanCV <- apply(subtypeMeanExprs,1,function(x){
  sd(x)/abs(mean(x))
})[colnames(capacities.list[[1]])]
tg.meanCV <- tg.meanCV[tg.anova.p.bonf<0.0001]
diff.genes <- names(rev(sort(tg.meanCV)))[1:500]



predError.mat <- matrix(0,length(subtypes),length(subtypes))
for (i in c(1:length(subtypes))) {
  for (j in c(1:length(subtypes))) {
    expr.est <- predExprs.list[[i]][[j]][diff.genes,]
    expr.obs <- data.matrix(expr[rownames(expr.est),colnames(expr.est)])
    temp <- (abs(expr.est) + abs(expr.obs))/2
    est.error <- mean(abs(expr.est-expr.obs)/temp)
    predError.mat[i,j] <- est.error
  }
}
# Try different ways of scoring prediction accuracy
# predAccuracy.mat <- matrix(0,length(subtypes),length(subtypes))
# for (i in c(1:length(subtypes))) {
#   for (j in c(1:length(subtypes))) {
#     expr.est <- predExprs.list[[i]][[j]][diff.genes,]
#     expr.obs <- data.matrix(expr[rownames(expr.est),colnames(expr.est)])
#     temp <- (abs(expr.est) + abs(expr.obs))/2
#     est.error.mat <- (abs(expr.est-expr.obs)/temp)
#     predAccuracy.mat[i,j] <- (sum(rowMeans(est.error.mat)<0.25))/nrow(est.error.mat)
#   }
# }
write.table((predError.mat),"analyses/0_0.1/lambda_0_0.1_crossSubtype_error_mat.txt",row.names=F,
            col.names=F,sep="\t",quote=F)

library(gplots)
# myheatmapfun <- colorRampPalette((c(rgb(231,224,219,maxColorValue=255),
#                                     rgb(234,145,118,maxColorValue=255),
#                                     rgb(102,3,32,maxColorValue=255))))
myheatmapfun <- colorRampPalette((c("white","darkblue")))

heatmap.data <- 1-predError.mat
dev.off()
hmRes <- heatmap.2(t(scale(t(heatmap.data[complete.cases(heatmap.data),]))),
                   Colv=NA,
                   Rowv=NA,
                   scale="none",
                   density.info="none", trace="none",
                   col=myheatmapfun(75),
                   reorderfun = function(d,w) { d },
                   cexCol = 0.8,
                   #srtCol = 45,
                   cexRow = 0.8,
                   rowsep = c(1:4),
                   colsep = c(1:4),
                   sepwidth = c(0.02,0.02),
                   #offsetCol = -0.5,
                   #RowSideColors = rscVec,
                   dendrogram="none",
                   key.title=NA,
                   key.xlab=NA,
                   lhei=c(1,8),
                   breaks=seq(-1.5,1.5,length.out=76),
                   labRow = subtypes,
                   labCol = subtypes,
                   margins=c(10,10)
)



# 3. Whether / how does F values correlate with correlation coefficients?
# Try Spearman.
binaryRegMat <- matrix(0,length(allTFs),length(regulators))
rownames(binaryRegMat) <- allTFs
colnames(binaryRegMat) <- names(regulators)
for (gene in names(regulators)) {
  binaryRegMat[regulators[[gene]],gene] <- 1
}

capmat.corr <- lapply(subtypes,function(x){
  subtypeExprs <- commonExpr[,names(commonLabels)[commonLabels==x]]
  cor(t(subtypeExprs[allTFs,]),t(subtypeExprs[allGenes,]),method = "spearman")
})
capmat.vec <- lapply(capacities.list,function(x){
  temp <- x[allTFs,names(regulators)][binaryRegMat!=0]
  as.numeric(temp)
})
capmat.corr.vec <- lapply(capmat.corr,function(x){
  temp <- x[allTFs,names(regulators)][binaryRegMat!=0]
  as.numeric(temp)
})
for (i in c(1:length(capmat.vec))) {
  fn <- paste0("analyses/0_0.1/FvsCorr_DFNet_",subtypes[i],".pdf")
  pdf(fn,width = 6,height = 5)
  plot(log2(capmat.vec[[i]]),capmat.corr.vec[[i]],pch=".",
       xlab="log2 F values",
       ylab="Spearman's correlation coefficient",
       main=subtypes[i])
  dev.off()
}

# 4. How does distribution of F vary with in-degrees of targets?

inDeg.list <- lapply(capacities.list,function(x){
  colSums(x!=1)
})
outDeg.list <- lapply(capacities.list,function(x){
  rowSums(x!=1)
})
inDeg.bins <- seq(0,342,1)
FvalByDin.list <- lapply(1:length(capacities.list),function(x){
  Fval.bins <- lapply(c(1:(length(inDeg.bins)-1)),function(i){
    tgtsInBin <- names(which(inDeg.list[[x]] >= inDeg.bins[i] & inDeg.list[[x]] < inDeg.bins[i+1]))
    Fvec <- capacities.list[[x]][,tgtsInBin]
    histdata <- hist(log2(as.numeric(Fvec[Fvec!=1])),breaks=seq(-70,25,1))
    cat(".")
    histdata$counts
  })
  binnedFvalMat <- do.call(rbind.data.frame,Fval.bins)
  binnedFvalMat
})

library(gplots)
myheatmapfun <- colorRampPalette((c(rgb(15,61,111,maxColorValue=255),
                                    rgb(125,182,212,maxColorValue=255),
                                    rgb(231,224,219,maxColorValue=255),
                                    rgb(234,145,118,maxColorValue=255),
                                    rgb(102,3,32,maxColorValue=255))))
myheatmapfun <- colorRampPalette((c(rgb(231,224,219,maxColorValue=255),
                                    rgb(234,145,118,maxColorValue=255),
                                    rgb(102,3,32,maxColorValue=255))))
heatmap.data <- data.matrix(FvalByDin.list[[4]])[,46:95]
dev.off()
library(pheatmap)
pheatmap(heatmap.data,cluster_rows = F,cluster_cols = F,
         labels_row = "",labels_col = "", breaks=seq(0,75,length.out=76),
         color = myheatmapfun(75))
hmRes <- heatmap.2(heatmap.data[complete.cases(heatmap.data),],
                   Colv=NA,
                   Rowv=NA,
                   scale="none",
                   density.info="none", trace="none",
                   col=myheatmapfun(75),
                   reorderfun = function(d,w) { d },
                   #cexCol = 0.3,
                   #srtCol = 45,
                   #cexRow = 1,
                   #offsetCol = -0.5,
                   #RowSideColors = rscVec,
                   dendrogram="none",
                   key.title=NA,
                   key.xlab=NA,
                   lhei=c(1,8),
                   breaks=seq(0,75,length.out=76),
                   labRow = as.character(inDeg.bins[2:length(inDeg.bins)]),
                   labCol = as.character(seq(-25,25,1)[2:51])
                   #margins=c(10,10)
)


# 4.1 How does F vary with the range of expression?
# Look in a target- and regulator-centric way respectively
# Should take degree into consideration too

exprsRanges.list <- lapply(subtypes,function(x){
  subtypeExprs <- commonExpr[,names(commonLabels)[commonLabels==x]]
  exprsRanges <- apply(subtypeExprs,1,function(y){
    temp <- range(y)
    abs(temp[1]-temp[2])
  })
  exprsRanges
})

# Calculate median of the absolute F values of regulators for each gene in each subtype
avgFByTarget.list <- lapply(subtypes,function(x){
  avgFByTarget <- sapply(colnames(capacities.log2.list[[x]]),function(y){
    if (length(regulators[[y]])>0) {
      median(abs((capacities.log2.list[[x]])[regulators[[y]],y]))
    } else {
      0
    }
  })
  avgFByTarget
})

# Calculate median of the absolute F values of regulators for each TF in each subtype
avgFByTF.list <- lapply(subtypes,function(x){
  avgFByTF <- sapply(rownames(capacities.log2.list[[x]]),function(y){
    if (T) {
      median(abs((capacities.log2.list[[x]])[y,capacities.log2.list[[x]][y,]!=0]))
    } else {
      0
    }
  })
  avgFByTF
})

colorFun <- colorRampPalette(c("gray","darkred"),alpha=0.7)
colorVec <- colorFun(max(unlist(inDeg.list)))
tgtGeneInDegreeColVec <- colorVec[inDeg.list[[1]]]
colorVec.TF <- colorFun(max(unlist(outDeg.list)))
TFOutDegreeColVec <- colorVec.TF[outDeg.list[[1]]]

# TF F vs outdegree looks interesting.
for (i in c(1:length(subtypes))) {
  fn1 <- paste0("analyses/0_0.1/",subtypes[i],"_DFNet_abslog2FvsTgtInDegree.pdf")
  fn2 <- paste0("analyses/0_0.1/",subtypes[i],"_DFNet_abslog2FvsTFExprsRange.pdf")
  pdf(file=fn1,width = 6,height = 5)
  subtype <- subtypes[i]
  plot(exprsRanges.list[[i]][colnames(capacities.list[[i]])],
       avgFByTarget.list[[i]][colnames(capacities.list[[i]])],col=tgtGeneInDegreeColVec,
       pch=16,
       xlab="target gene expression range",
       ylab="median TF regulatory capacity (|log2(F)|)",
       main=subtype)
  dev.off()
  pdf(file=fn2,width = 6,height = 5)
  subtype <- subtypes[i]
  plot(exprsRanges.list[[i]][names(avgFByTF.list[[i]])],
       avgFByTF.list[[i]],col=TFOutDegreeColVec,
       pch=16,
       xlab="TF expression range",
       ylab="median TF regulatory capacity (|log2(F)|)",
       main=subtype)
  dev.off()
}
#dev.off()
# Also look at average F vs average TF expression (per Rick's suggestion 10/25/17)
TF.avgExprs.list <- lapply(subtypes,function(x){
  exprs.subtype <- commonExpr[,names(commonLabels)[commonLabels==x]]
  TF.avgExprs <- rowMeans(exprs.subtype[allTFs,])
  TF.avgExprs
})
TF.CV.list <- lapply(subtypes,function(x){
  exprs.subtype <- commonExpr[,names(commonLabels)[commonLabels==x]]
  TF.avgExprs <- rowMeans(exprs.subtype[allTFs,])
  TF.SD <- apply(exprs.subtype[allTFs,],1,sd)
  TF.SD/TF.avgExprs
})
for (i in c(1:length(subtypes))) {
  fn1 <- paste0("analyses/0_0.1/",subtypes[i],"_DFNet_abslog2FvsTFAvgExprs.pdf")
  pdf(file=fn1,width = 6,height = 5)
  subtype <- subtypes[i]
  plot(TF.avgExprs.list[[i]][names(avgFByTF.list[[i]])],
       avgFByTF.list[[i]],col=TFOutDegreeColVec,
       pch=16,
       xlab="TF mean expression",
       ylab="median TF regulatory capacity (|log2(F)|)",
       ylim=c(0,2),
       main=subtype)
  dev.off()
}
for (i in c(1:length(subtypes))) {
  fn1 <- paste0("analyses/0_0.1/",subtypes[i],"_DFNet_abslog2FvsTFexprsCV.pdf")
  pdf(file=fn1,width = 6,height = 5)
  subtype <- subtypes[i]
  plot(TF.CV.list[[i]][names(avgFByTF.list[[i]])],
       avgFByTF.list[[i]],col="black",
       pch=16,
       xlab="TF expression variability (CV)",
       ylab="median TF regulatory capacity (|log2(F)|)",
       ylim=c(0,1.5),
       main=subtype)
  dev.off()
}

extremeExamples <- lapply((1:length(subtypes)),function(i){
  left <- names(avgFByTF.list[[i]])[avgFByTF.list[[i]]>0.8]
  right <- names(TF.CV.list[[i]])[TF.CV.list[[i]]>0.2]
  list(left=left,right=right)
})

# 5. Are signature TF-gene pairs also significant in label-shuffled data?
# 01-03-17: no shuffled data available yet. Leave significance part for later
# First get significant signature TF-gene pairs using the 'out of range' criterion
niter <- 30
capmat.unshuffled <- lapply(capacities.list,function(x){
  x[allTFs,]
})
# # for (subtype in cl) {
# #   fn <- paste0(subtype,"_TFGeneCapacities_ri_nomi_randomLabels_",paste0("_lambda",lambda),
# #                "_iter_0.txt")
# #   tfm.mat <- read.delim(fn,header = T,row.names = 1,sep="\t",check.names = F)
# #   capmat.unshuffled[[subtype]] <- data.matrix(tfm.mat[allTFs,])
# # }
lambda <- "0.1_0"
capmat.shuffled <- list()
for (i in c(1:(niter))) {
  capmat <- list()
  for (subtype in subtypes) {
    fn <- paste0("data/regrOutput/labelShuffling/0.1/",subtype,"_TFGeneCapacities_ri_nomi_",paste0("_lambda",lambda),
                 "_randomShuffle_",as.character(i),".txt")
    tfm.mat <- read.delim(fn,header = T,row.names = 1,sep="\t",check.names = F)
    if (!("ZNF8" %in% rownames(tfm.mat))) {
      tempRow <- rep(1,ncol(tfm.mat))
      TFs <- intersect(rownames(tfm.mat),c(allTFs,'const'))
      temp <- rbind(tfm.mat,tempRow)
      rownames(temp) <- c(TFs,"ZNF8")
      tfm.mat <- temp
    }
    geneRegs <- data.matrix(tfm.mat[allTFs,])
    capmat[[subtype]] <- geneRegs
  }
  cat(".")
  capmat.shuffled[[i]] <- capmat
}

nsd <- 5

#source("src/extractFSignatures.R")
source("src/obtainSignatures.R")
# 
FValueSig.res <- obtainSignatures(nsd=5,capmat=capmat.unshuffled,adjPThres = 1E-10)


# For signature analyses see 'extractFSignatures.R' for now...
fn <- "analyses/0_0.1/domSigSubtype.filtered_percentCutoff0.000.multipleTargets5.withSCExprsTop3000.txt"
# write.table(sigTFs,
#             fn,
#             col.names=F,sep="\t",quote=F)
sigTFs <- (read.delim(fn,header=F,check.names=F))
sigTFNames <- as.character(sigTFs[,1])
sigTFs <- as.character(sigTFs[,2])
names(sigTFs) <- sigTFNames


# 6. Check TF partnering / clusters
# For top correlated TF-TF pairs, do they share motifs? Are they more closely connected
# in the protein network? 
# 
# Need to reduce zero-inflation & correct for 2-D outliers
geneRegs.all <- lapply(capmat.unshuffled,log2)
geneRegs.mean <- Reduce("+",geneRegs.all)
geneRegs.mean <- geneRegs.mean / 4
geneRegs.var <- lapply((geneRegs.all),function(x){(x-geneRegs.mean)^2})
geneRegs.var <- Reduce("+",geneRegs.var) / 4
geneRegs.sd <- geneRegs.var^0.5
geneRegs.TF <- lapply(geneRegs.all,function(x){
  x[allTFs,]
})
geneRegs.TF.sd <- geneRegs.sd[allTFs,]
meanTFCV <- (geneRegs.TF.sd)/abs(geneRegs.mean[allTFs,])
meanTFCV <- t(apply(meanTFCV,1,function(x){
  y <- x
  y[is.na(y)] <- 0
  y
}))

# TF-TF correlations are pre-computed!

log2F.threshold <- 5

lambda <- "0.1_0"
allTFInteractions.tuples <- read.delim("data/knownTFInteractionTuples.BIOGRID.txt",
                                       header=F,check.names = F)
allTFInteractions.tuples <- as.character(allTFInteractions.tuples[,1])

geneTFCor.list <- lapply(c(1:length(subtypes)),function(i){
  fn <- paste0("data/",subtypes[i],"_cormat_",as.character(lambda),"_rs_7792275.txt")
  geneTF.coreg.mat <- data.matrix(read.delim(fn,header=F,check.names=F))
  geneTF.coreg.mat[is.nan(geneTF.coreg.mat)] <- 0
  rownames(geneTF.coreg.mat) <- colnames(geneTF.coreg.mat) <- rownames(capacities.list[[1]])
  geneTF.coreg.mat
})
names(geneTFCor.list) <- subtypes
rs <- 7792275

geneTFCor.list <- lapply(c(1:length(subtypes)),function(i){
  cat(".")
  fn <- paste0("data/",subtypes[i],"_cormat_",as.character(lambda),"_rs_",
               as.character(rs),".txt")
  geneTF.coreg.mat <- data.matrix(read.delim(fn,header=F,check.names=F))
  geneTF.coreg.mat[is.nan(geneTF.coreg.mat)] <- 0
  rownames(geneTF.coreg.mat) <- colnames(geneTF.coreg.mat) <- rownames(capacities.list[[1]])
  geneTF.coreg.mat[lower.tri(geneTF.coreg.mat)] <- t(geneTF.coreg.mat)[lower.tri(geneTF.coreg.mat)]
  geneTF.coreg.mat
})
names(geneTFCor.list) <- subtypes


library(corrplot)

myheatmapfun <- colorRampPalette((c(rgb(15,61,111,maxColorValue=255),
                                    rgb(125,182,212,maxColorValue=255),
                                    rgb(231,224,219,maxColorValue=255),
                                    rgb(234,145,118,maxColorValue=255),
                                    rgb(102,3,32,maxColorValue=255))))

heatmap.data <- geneTFCor.list[[1]]
#corrplot(heatmap.data,type="upper")
clust.res <- hclust(dist(heatmap.data),method = "average")
#clust.res <- kmeans(heatmap.data,centers = 5)
heatmap.data <- geneTFCor.list[[4]]
pheatmap(heatmap.data[clust.res$order,clust.res$order],
         show_rownames = F,show_colnames = F,color=myheatmapfun(75),
         cluster_rows = F,cluster_cols = F,cellwidth = 1,cellheight = 1)

library(igraph)
library(igraph)
allTFPairs <- c()
count <- 1
for (TF1 in allTFs) {
  for (TF2 in allTFs) {
    allTFPairs[count] <- paste0(TF1,"_",TF2)
    count <- count + 1
  }
}

TFCor.tuples.list <- lapply(geneTFCor.list,function(gtc){
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
names(TFCor.tuples.list) <- subtypes

# TFCor.tuples.list <- lapply(geneTFCor.list,function(gtc){
#   g <- graph.adjacency(gtc,mode=c("directed"),weighted = T)
#   el <- get.edgelist(g)
#   weights <- get.edge.attribute(graph = g,name = "weight")
#   el <- apply(el[,1:2],1,function(x){
#     paste0(x[1],"_",x[2])
#   })
#   names(weights) <- el
#   cat(".")
#   weights
# })

#TFCor.tuples.list <- TFCor.tuples.lists.byLambda[[lambda]]
percent.seq <- seq(0.01,1,0.01)
commonTFCor.tuples <- Reduce("intersect",lapply(TFCor.tuples.list,names))
commonTFCor.tuples.cor <- lapply(TFCor.tuples.list,function(x){
  x[commonTFCor.tuples]
})
commonTFCor.tuples.cor <- data.matrix(do.call(cbind.data.frame,commonTFCor.tuples.cor))
commonTFCor.tuples.avgCor <- rowMeans(commonTFCor.tuples.cor)
commonTFCor.tuples.avgCor <- rev(sort(commonTFCor.tuples.avgCor))
# Output top correlated pairs (0.01) for visualization
qThreshold <- 0.05
visPairs <- commonTFCor.tuples.avgCor[1:ceiling(qThreshold*length(commonTFCor.tuples.avgCor))]
hits <- intersect(names(commonTFCor.tuples.avgCor[1:ceiling(qThreshold*length(commonTFCor.tuples.avgCor))]),
                  allTFInteractions.tuples)
visPairs <- visPairs[hits]
network.data <- lapply(names(visPairs),function(p){
  temp <- strsplit(p,"_",fixed=T)[[1]]
  TF1 <- temp[1]
  TF2 <- temp[2]
  col <- "black"
  p.rev <- paste0(TF2,"_",TF1)
  if (p %in% hits | p.rev %in% hits) {
    col <- "red"
  }
  c(TF1,TF2,col)
})
network.data <- do.call(rbind.data.frame,network.data)
write.table(network.data,"analyses/0_0.1/topCommonCorrelatedTFPairs_BIOGRID.txt",
            sep='\t',col.names=F,row.names=F,quote=F)

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
write.table(cbind(rec,prec),"analyses/0_0.1/rec_prec_0.1_maxrec1.txt",
            sep='\t',row.names=F,col.names=c("rec","prec"),quote = F)
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
precision <- length(hits)/ceiling(p*length(commonTFCor.tuples.avgCor))
recall <- length(hits)/length(allTFInteractions.tuples)
#recall[lambda] <- length(hits)/length(commonTFCor.tuples.avgCor)
AUC <- sum(diff(rec[complete.cases(prec)])*rollmean(prec[complete.cases(prec)],2))
print(AUC)
auc <- AUC
fisher.p <- fisher.test((x))$p.value
# AUC <- sum(diff(rec)*rollmean(prec,2))
# print(AUC)
# auc[lambda] <- AUC
# fisher.p[lambda] <- fisher.test((x))$p.value
# 
# 

# Use random sampling of TF-TF pairs to show enrichment
sampleSize <- ceiling(p*length(commonTFCor.tuples.avgCor))
allTFPairs <- c()
count <- 1
for (TF1 in allTFs) {
  for (TF2 in allTFs) {
    allTFPairs[count] <- paste0(TF1,"_",TF2)
    count <- count + 1
  }
}
set.seed(12345678)
niter <- 2000
random.draw.hits <- sapply(1:niter,function(i){
  TF.pairs <- sample(allTFPairs,sampleSize,replace = F)
  length(intersect(TF.pairs,allTFInteractions.tuples))/sampleSize
})
hist(random.draw.hits,xlab="Precision",breaks = 50, xlim=c(0,0.03))
sum(random.draw.hits>precision)/niter
abline(v=precision,lwd=1.5)

tfs <- rownames(geneRegs.TFs.coreg[[1]])
geneRegs.TFs.coreg <- geneTFCor.list
examplePairs <- list()
exampleTuples <- c()
for (i in 1:length(allTFs)) {
  for (j in 1:length(allTFs)) {
    corVec <- sapply(geneRegs.TFs.coreg,function(x){
      x[i,j]
    })
    idx <- which(abs(corVec)==max(abs(corVec)))[1]
    if (max(abs(corVec)) > 0.5 & max(abs(corVec[-idx])) < 0.1) {
      examplePairs[[paste0(tfs[i],"_",tfs[j])]] <- corVec
      exampleTuples[[paste0(tfs[i],"_",tfs[j])]] <- paste0(tfs[i],"_",tfs[j])
    }
  }
}
examplePairs <- do.call(rbind.data.frame,examplePairs)
colnames(examplePairs) <- subtypes
rownames(examplePairs) <- exampleTuples

TF1 <- "MAX"
TF2 <- "MYC"
subtype <- 4
ind <- ((capacities.log2.list[[subtype]][TF1,])>-2 & (capacities.log2.list[[subtype]][TF2,])>-2 &
          (capacities.log2.list[[subtype]][TF1,])<2 & (capacities.log2.list[[subtype]][TF2,])<2 &
          (capacities.log2.list[[subtype]][TF1,])!=0 & (capacities.log2.list[[subtype]][TF2,])!=0)

cor(capacities.log2.list[[subtype]][TF1,ind],
    capacities.log2.list[[subtype]][TF2,ind])
pointCol <- rep("black",(sum(ind)))
pointCol[capacities.log2.list[[subtype]][TF2,ind] < -0.5] <- "salmon"
plot(capacities.log2.list[[subtype]][TF1,ind],
     capacities.log2.list[[subtype]][TF2,ind],pch=20,xlim=c(-2,2),ylim=c(-2,2),col=pointCol)

library(spatstat)
subtypeColors <- c("brown3","cornflowerblue","mediumseagreen","goldenrod1")
names(subtypeColors) <- subtypes
#pdf(file = paste0("analyses/0_0.1/",TFoI,"_",tgt,"_",subtype,".pdf"),width = 10,height = 4)
par(mar = c(2, 2, 2, 2))
pppo=ppp(x=capacities.log2.list[[subtype]][TF1,ind],
         y=capacities.log2.list[[subtype]][TF2,ind],
         window = owin(c(-2,2),c(-2,2)))
den=density(pppo,kernel="gaussian",edge=T,diggle=T,adjust=0.4)
plot(den,main='TF vs target',col=colorRampPalette(c("white",subtypeColors[subtype]))(10),xlim=c(-2,2),
     ylim=c(-2,2),xlab=TF1,ylab=TF2)
points(capacities.log2.list[[subtype]][TF1,ind],
       capacities.log2.list[[subtype]][TF2,ind],
       xlim=c(-2,2),ylim=c(-2,2),pch=20,col="black")

# library(mvoutlier)
# geneRegs.TFs.coreg <- lapply(geneRegs.TF, function(x){
#   SDs <- apply(x,1,sd)
#   temp <- x[SDs!=0,]
#   temp <- temp[,colSums(abs(temp))!=0]
#   cormat <- matrix(0,nrow(temp),nrow(temp))
#   rownames(cormat) <- rownames(temp)
#   colnames(cormat) <- rownames(temp)
#   count <- 0
#   for (i in rownames(cormat)) {
#     count <- count + 1
#     if (count %% 10 == 0) {
#       cat(".")
#     }
#     for (j in colnames(cormat)) {
#       if (i!=j) {
#         tfi <- temp[i,]
#         tfj <- temp[j,]
#         tgts <- which(abs(tfi)>log2F.threshold & abs(tfj)>log2F.threshold)
#         if (length(tgts)>3) {
#           df <- cbind(tfi[tgts],tfj[tgts])
#           pcout.res <- pcout(df)
#           pcout.scores <- pcout.res$wfinal
#           # Kick out two most likely outliers
#           pcout.lowScores <- order(pcout.scores)[1:2]
#           # df <- df[pcout(df)$wfinal01 > 0,]
#           df <- df[-pcout.lowScores,]
#           if (nrow(df) > 3) {
#             cormat[i,j] <- cor(df[,1],df[,2],method="pearson")
#           } else {
#             cormat[i,j] <- 0
#           }
#         } else {
#           cormat[i,j] <- 0
#         }
#       }
#     }
#   }
#   cat("\n")
#   cormat
# })
# #dump("geneRegs.TFs.coreg","geneRegs.TFs.coreg._DFNet_030418.R")
# #source("geneRegs.TFs.coreg._DFNet.R")
# for (i in c(1:length(subtypes))) {
#   fn <- paste0("analyses/0_0.1/",subtypes[i],"_TF_TF_cormat_DFNet_log2F_2.txt")
#   write.table(geneRegs.TFs.coreg[[i]],fn,sep="\t",col.names=T,row.names=T,quote=F)
# }
# geneRegs.TFs.coreg <- lapply(subtypes,function(x){
#   fn <- paste0(x,"_TF_TF_cormat_DFNet.txt")
#   data.matrix(read.delim(fn,header=T,row.names=1,check.names=F))
# })

#geneRegs.TFs.coreg.kmeans <- kmeans(geneRegs.TFs.coreg[[2]],centers = 10)
#rowOrder <- order(geneRegs.TFs.coreg.kmeans$cluster)
library(gplots)
myheatmapfun <- colorRampPalette((c(rgb(15,61,111,maxColorValue=255),
                                    rgb(125,182,212,maxColorValue=255),
                                    rgb(231,224,219,maxColorValue=255),
                                    rgb(234,145,118,maxColorValue=255),
                                    rgb(102,3,32,maxColorValue=255))))
heatmap.data <- geneTFCor.list[[1]]
goodCors <- which(sapply(1:nrow(heatmap.data),function(x){
  mean(abs(heatmap.data[x,]))
})>0.1)
goodCors <- 1:nrow(heatmap.data)
heatmap.data <- heatmap.data[goodCors,goodCors]
d <- dist(heatmap.data,method = "euclidean")
fit <- hclust(as.dist(1-cor(heatmap.data)),method = "ward.D")
kmeans.res <- kmeans(heatmap.data,centers = 10)$cluster
kmeans.order <- order(kmeans.res)
#myclustfun <- function(x){hclust(x,method="complete")}
#heatmap.data <- heatmap.data[xxx,xxx]
#heatmap.data[abs(heatmap.data)<0.3] <- 0
dev.off()
hmRes <- heatmap.2(heatmap.data[kmeans.order,kmeans.order],
                   Colv=NA,
                   Rowv=NA,
                   #symm=T,revC=T,
                   scale="none",
                   density.info="none", trace="none",
                   col=myheatmapfun(75),
                   hclustfun = myclustfun,
                   reorderfun = function(d,w) { d },
                   cexCol = 0.8,
                   #srtCol = 45,
                   cexRow = 0.8,
                   #offsetCol = -0.5,
                   #RowSideColors = rscVec,
                   dendrogram="none",
                   key.title=NA,
                   key.xlab=NA,
                   #lhei=c(1,8),
                   breaks=seq(-0.5,0.5,length.out=76),
                   labRow = "",
                   labCol = ""
)
goodTFsInOrder <- rownames(heatmap.data[fit$order,fit$order])[hmRes$rowInd]
heatmap.data <- geneRegs.TFs.coreg[[1]]
heatmap.data <- heatmap.data[goodTFsInOrder,goodTFsInOrder]
dev.off()
hmRes <- heatmap.2(heatmap.data,
                   Colv=NA,
                   Rowv=NA,
                   #symm=T,revC=T,
                   scale="none",
                   density.info="none", trace="none",
                   col=myheatmapfun(75),
                   hclustfun = myclustfun,
                   #reorderfun = function(d,w) { d },
                   cexCol = 0.8,
                   #srtCol = 45,
                   cexRow = 0.8,
                   #offsetCol = -0.5,
                   #RowSideColors = rscVec,
                   dendrogram="none",
                   key.title=NA,
                   key.xlab=NA,
                   #lhei=c(1,8),
                   breaks=seq(-0.3,0.3,length.out=76),
                   labRow = "",
                   labCol = ""
)

# Try circos in R for visualization
# 010717 - not now.
# library(circlize)
# col_fun <- colorRamp2(c(-1,-0.5,0,0.5,1),
#                       c(rgb(15,61,111,maxColorValue=255),
#                         rgb(125,182,212,maxColorValue=255),
#                         rgb(231,224,219,maxColorValue=255),
#                         rgb(234,145,118,maxColorValue=255),
#                         rgb(102,3,32,maxColorValue=255)))
# chordDiagram(geneRegs.TFs.coreg[[1]],col=col_fun(geneRegs.TFs.coreg[[1]]),
#              grid.col = NA, grid.border = "black",
#              link.largest.ontop = TRUE,
#              preAllocateTracks = list(
#                list(track.height = 0.02)
#              ))

# Create data files for circos plots of TFs
# chromosomes are now TFs
sigDiffTFs <- FValueSig.res$sigTFs.anova
#sigDiffTFs <- names(sort(domSigSubtype))
# Reconsider sig. diff. TFs: take TFs that came out from filtering the 
# SD reduction matrix for signature discovery?
nonSigDiffTFs <- setdiff(rownames(geneRegs.TFs.coreg[[1]]),sigDiffTFs)
allChr <- c(allTFs)
chr.list <- lapply(1:nrow(geneRegs.TFs.coreg[[1]]),function(x){
  c("chr","-",paste0("hs",as.character(x)),allChr[x],
    0,1000,allChr[x])
})
chr.df <- do.call(rbind.data.frame,chr.list)
colnames(chr.df) <- ""
write.table(chr.df,"analyses/circos/DFNet_karyotype.TFs.txt",sep=" ",col.names=F,row.names=F,quote=F)
# Prepare links
chrNameLookUp <- as.character(chr.df[,3])
names(chrNameLookUp) <- as.character(chr.df[,7])
links.list <- lapply(geneRegs.TFs.coreg,function(x){
  linkList <- list()
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      if (abs(x[i,j])>0.8) {
        nodeA <- chrNameLookUp[rownames(x)[i]]
        nodeB <- chrNameLookUp[rownames(x)[j]]
        linkList <- c(linkList,list(c(nodeA,0,1000,nodeB,0,1000)))
      }
    }
  }
  linkList <- do.call(rbind.data.frame,linkList)
  colnames(linkList) <- ""
  linkList
})
lapply(c(1:length(subtypes)),function(x){
  fn <- paste0("analyses/0_0.1/circos/DFNet_TF_coreg_new_links_TFs_",subtypes[x],".txt")
  write.table(links.list[[x]],fn,sep=" ",row.names=F,col.names=F,quote=F)
})

# Also derive signature matrix for TF-TF correlations using the LOSD method
geneRegs.coreg.mean <- Reduce("+",geneRegs.TFs.coreg)
geneRegs.coreg.mean <- geneRegs.coreg.mean / 4
geneRegs.coreg.var <- lapply((geneRegs.TFs.coreg),function(x){(x-geneRegs.coreg.mean)^2})
geneRegs.coreg.var <- Reduce("+",geneRegs.coreg.var) / 4
geneRegs.coreg.sd <- geneRegs.coreg.var^0.5

leaveOneSD <- lapply(c(1:length(subtypes)),function(x){
  currList <- geneRegs.TFs.coreg[-x]
  LOM <- Reduce("+",currList)
  LOM <- LOM / (length(subtypes)-1)
  LOV <- lapply(currList,function(y){(y-LOM)^2})
  LOV <- Reduce("+",LOV) / (length(subtypes)-1)
  LOSD <- LOV^0.5
  LOSD
})

overallSD <- geneRegs.coreg.sd

sigTFCorSubtype <- matrix(0,nrow(geneRegs.TFs.coreg[[1]]),ncol(geneRegs.TFs.coreg[[1]]))

for (i in c(1:nrow(geneRegs.TFs.coreg[[1]]))) {
  for (j in c(1:ncol(geneRegs.TFs.coreg[[1]]))) {
    Fval <- sapply(leaveOneSD,function(x){
      x[i,j]
    })
    coregs.allSubtypes <- sapply(c(1:length(subtypes)),function(k){
      geneRegs.TFs.coreg[[k]][i,j]
    })
    if (overallSD[i,j] > 0 & max(abs(coregs.allSubtypes)) > 0.5) {
      sigTFCorSubtype[i,j] <- which(Fval==min(Fval))[1]
    }
  }
}
rownames(sigTFCorSubtype) <- rownames(geneRegs.TFs.coreg[[1]])
colnames(sigTFCorSubtype) <- colnames(geneRegs.TFs.coreg[[1]])
# table(sigTFCorSubtype[names(domSigSubtype),])
# sigTFCorSubtype.sigTFs <- sigTFCorSubtype[names(domSigSubtype),]
# sigTFCorSubtype.sigTFs.druggable <- sigTFCorSubtype.sigTFs[,druggableTFs]
# write.table(sigTFCorSubtype,"TFTFPairSigMat.txt",sep="\t",quote=F)
# write.table(sigTFCorSubtype.sigTFs.druggable,"TFTFPairSigMat.sigTFs.druggable.txt",sep="\t",quote=F)
# sigTFCorSubtype.sigTFs.ordered <- t(apply(sigTFCorSubtype.sigTFs,1,function(x){
#   prefix <- rep("0",sum(x==0))
#   suffix <- sort(table(x[x!=0]))
#   suffix <- unlist(lapply(1:length(suffix),function(i){
#     rep(names(suffix)[i],suffix[i])
#   }))
#   as.numeric(c(prefix,suffix))
# }))

# Do top correlated TF pairs interact in the STRING network?
# Inspect direct connections

# protNetwork.full <- read.delim("../../fullNetworkPrEL.txt",header=F,check.names=F)
library(igraph)
# protNetwork.full <- read.delim("data/BIOGRID-ORGANISM-Homo_sapiens-3.4.158.tab2.txt",header=T,check.names=F)
# protNetwork.full <- cbind(as.character(protNetwork.full$`Official Symbol Interactor A`),
#                           as.character(protNetwork.full$`Official Symbol Interactor B`))
# protNetwork.full <- unique(protNetwork.full)
# protNetwork.graph <- graph.edgelist(protNetwork.full,directed=T)
# TFsInProtNet <- intersect(as.character(names(V(protNetwork.graph))),allTFs)
# protNetwork.graph.TF <- induced.subgraph(protNetwork.graph,vids = TFsInProtNet)
# # col1 <- as.character(protNetwork.full[,1])
# # colId1 <- lapply(allTFs,function(x){
# #   grep(x,col1,fixed=T)
# # })
# # colId1 <- unlist(colId1)
# # colId1 <- colId1[complete.cases(colId1)]
# # col2 <- as.character(protNetwork.full[,2])
# # colId2 <- lapply(allTFs,function(x){
# #   grep(x,col2,fixed=T)
# # })
# # colId2 <- unlist(colId2)
# # colId2 <- colId2[complete.cases(colId2)]
# allTFInteractions <- get.edgelist(protNetwork.graph.TF)
# allTFInteractions.tuples <- apply(allTFInteractions,1,function(x){
#   paste0(x[1],"_",x[2])
# })
# write.table(allTFInteractions.tuples,"data/knownTFInteractionTuples.BIOGRID.txt",sep="\t",row.names=F,col.names=F,quote=F)
# TFCor.tuples.list <- lapply(geneRegs.TFs.coreg,function(x){
#   TFCor.tuples <- rep(0,sum(x!=0))
#   nonZeroIdx <- which(x!=0)
#   nr <- nrow(x)
#   tfs <- rownames(x)
#   for (idx in nonZeroIdx) {
#     temp <- idx %% nr
#     if (temp == 0) {
#       rowInd <- nr
#       colInd <- idx %/% nr
#     } else {
#       rowInd <- temp
#       colInd <- idx %/% nr + 1
#     }
#     if (rowInd != colInd) {
#       tfi <- tfs[rowInd]
#       tfj <- tfs[colInd]
#       tup <- paste0(tfi,"_",tfj)
#       TFCor.tuples[tup] <- x[rowInd,colInd]
#     }
#   }
#   cat(".")
#   TFCor.tuples
# })
# TFCor.tuples.list <- lapply(TFCor.tuples.list,function(x){
#   x[x!=0]
# })
# For now only look at positively correlated TF pairs
# TFCor.tuples.list <- lapply(TFCor.tuples.list,function(x){
#   temp <- rev(sort(x))
#   temp[temp > 0]
# })
percent.seq <- seq(0,0.1,0.001)
# commonTFCor.tuples <- Reduce("intersect",lapply(TFCor.tuples.list,names))
commonTFCor.tuples <- Reduce("intersect",lapply(TFCor.tuples.list,names))
commonTFCor.tuples.cor <- lapply(TFCor.tuples.list,function(x){
  x[commonTFCor.tuples]
})
commonTFCor.tuples.cor <- data.matrix(do.call(cbind.data.frame,commonTFCor.tuples.cor))
commonTFCor.tuples.avgCor <- rowMeans(commonTFCor.tuples.cor)
commonTFCor.tuples.avgCor <- rev(sort(commonTFCor.tuples.avgCor))
# commonTFCor.tuples.list <- lapply(percent.seq,function(x){
#   subList <- lapply(TFCor.tuples.list,function(y){
#     names(y[1:round(x*length(y))])
#   })
#   commonTuples <- Reduce("intersect",subList)
# })
# TFCor.tuples.hits.list <- lapply(TFCor.tuples.list,function(x){
#   inters.list <- sapply(percent.seq,function(p){
#     topTFCor.tuples <- x[1:round(p*length(x))]
#     length(intersect(names(topTFCor.tuples),allTFInteractions.tuples))/round(p*length(x))
#   })
# })
# plot(percent.seq,TFCor.tuples.hits.list[[1]],type="l",lwd=3)
TFCor.tuples.hits.list <- sapply(percent.seq,function(p){
  hits <- intersect(names(commonTFCor.tuples.avgCor[1:round(p*length(commonTFCor.tuples.avgCor))]),
                    allTFInteractions.tuples)
  length(hits)/round(p*length(commonTFCor.tuples.avgCor))
})

# How to tell if the following QQ-like plot shows significant
# deviation from random subsetting of the list of known interactors?
plot(percent.seq,TFCor.tuples.hits.list,type="l",lwd=3,
     xlab="percent top correlated TF pairs",
     ylab="proportion of TF pairs with known interactions")

p = 0.02
hits <- intersect(names(commonTFCor.tuples.avgCor[1:round(p*length(commonTFCor.tuples.avgCor))]),
                  allTFInteractions.tuples)
x = matrix(c(length(hits),
             round(p*length(commonTFCor.tuples.avgCor)),
             length(intersect(names(commonTFCor.tuples.avgCor),allTFInteractions.tuples)),
             length(commonTFCor.tuples.avgCor)),2,2)
fisher.test(x) # p value is 0.006047

# The above analysis did not take into consideration indirect interactions
# Highly correlated TF pairs might also act in complexes where direct interactions were
# not observed or because they both interact with the same third TF
# To obtain this information, first build a TF-TF intxn network from BIOGRID,
# then look at shortest paths of length 2 or less, i.e. for each TF pair in the list
# of top correlations at each cutoff, see if their shortest path length in the BIOGRID
# network is at most 2.
library(igraph)
shortestPL <- distances(protNetwork.graph.TF)
availTFs <- rownames(shortestPL)
commonTFCor.tuples.dist <- sapply(names(commonTFCor.tuples.avgCor),function(x){
  temp <- strsplit(x,"_",fixed=T)[[1]]
  TF1 <- temp[1]
  TF2 <- temp[2]
  if (TF1 %in% availTFs & TF2 %in% availTFs) {
    shortestPL[TF1,TF2] <= 2
  } else {
    F
  }
})

percent.seq <- seq(0.001,1,0.01)

TFCor.tuples.hits.list <- sapply(percent.seq,function(p){
  topCor <- commonTFCor.tuples.dist[1:round(p*length(commonTFCor.tuples.avgCor))]
  hits <- names(topCor)[topCor==T]
  length(hits)/round(p*length(commonTFCor.tuples.avgCor))
})

# How to tell if the following QQ-like plot shows significant
# deviation from random subsetting of the list of known interactors?
plot(percent.seq,TFCor.tuples.hits.list,type="l",lwd=3,
     xlab="percent top correlated TF pairs",
     ylab="proportion of TF pairs with known interactions")

p = 0.005

topCor <- commonTFCor.tuples.dist[1:round(p*length(commonTFCor.tuples.avgCor))]
hits <- names(topCor)[topCor==T]
x = matrix(c(length(hits),
             round(p*length(commonTFCor.tuples.avgCor)),
             sum(commonTFCor.tuples.dist),
             length(commonTFCor.tuples.avgCor)),2,2)
fisher.test(x) # p value is 0.006047


# 7. Enumerate motifs
library(igraph)
for (subtype in subtypes) {
  tg <- graph.data.frame(capacities.log2.el.sigF.list[[subtype]],directed=T)
  res <- maximal.cliques(as.undirected(tg),min=3,max=3)
  tgv <- V(tg)$name
  res <- lapply(res,function(x){tgv[x]})
  res <- do.call(cbind.data.frame, res)
  rownames(res) <- c()
  colnames(res) <- c()
  write.table(t(res),paste0(subtype,"_3cliques.sigFlog2_1.txt"),
              sep="\t",col.names=F,row.names=F,quote=F)
}

# 8. Correlation with scRNA-seq data

scExprs <- read.delim("data/GSE84465_GBM_All_data.csv",
                      sep=" ",header = T, row.names=1, check.names=F,quote = "\"")
scExprs <- data.matrix(scExprs)
scExprs <- scExprs[rowSums(scExprs)!=0,]
expressedGenes <- which(apply(scExprs,1,function(x){
  sum(x>0) > 0.1*ncol(scExprs)
}))
scExprs <- (scExprs)[expressedGenes,]
scExprs <- log2(scExprs+1)

geneUniv <- rownames(scExprs)

mapped_gene_TFs <- read.delim("data/mapped_sc_gene_TFs.txt",header=F,check.names=F)

TFList <- unique(c(as.character(mapped_gene_TFs[,2])))
tgtGenes <- unique(as.character(mapped_gene_TFs[,1]))
withAllRegs <- lapply(tgtGenes,function(x){
  if (length(intersect(regulators[[x]],geneUniv))>=0.85*length(regulators[[x]]) &
      length(regulators[[x]]) >= 3) {
    T
  } else {
    F
  }
})
tgtGenesWithAllRegs <- tgtGenes[unlist(withAllRegs)]
length(tgtGenesWithAllRegs)

capacities.sc.list <- lapply(subtypes,function(x){
  fn <- paste0("data/regrOutput/scExprs/",x,"_TFGeneCapacities_ri_nomi__lambda0.1_0_sc.txt")
  capacities.sc <- read.delim(fn,sep="\t",check.names=F)
  data.matrix(capacities.sc)
})
names(capacities.sc.list) <- subtypes

lambda <- "0.1_0"
capacities.list <- lapply(subtypes,function(x){
  fn <- paste0("data/regrOutput/paramSweep/nonlinear/",x,"_TFGeneCapacities_ri_nomi__lambda",lambda,".txt")
  temp <- data.matrix(read.delim(fn,header=T,row.names=1,check.names=F))
  #constid <- which(rownames(temp)=="const")
  #temp[-constid,]
})
names(capacities.list) <- subtypes

corWithBulk.list <- lapply(subtypes,function(x){
  capmat.bulk <- log2(capacities.list[[x]])
  capmat.sc <- log2(capacities.sc.list[[x]])
  TFs <- intersect(rownames(capmat.sc),rownames(capmat.bulk))
  tgts <- intersect(colnames(capmat.sc),colnames(capmat.bulk))
  sapply(c(1:length(TFs)),function(y){
    cor(capmat.bulk[TFs,tgts][y,],capmat.sc[TFs,tgts][y,],method="pearson")
  })
})
names(corWithBulk.list) <- subtypes

# Simulate background median distribution of correlations
# using random permutations

niter <- 1000

corWithBulk.bg.list <- lapply(subtypes,function(x){
  capmat.bulk <- log2(capacities.list[[x]])
  capmat.sc <- log2(capacities.sc.list[[x]])
  cat('\n')
  perm.medians <- sapply(c(1:niter),function(i){
    cat(".")
    TFs <- intersect(rownames(capmat.sc),rownames(capmat.bulk))
    tgts <- intersect(colnames(capmat.sc),colnames(capmat.bulk))
    ntgts <- length(tgts)
    allcor <- sapply(c(1:length(TFs)),function(y){
      capvec.sc.perturb <- capmat.sc[TFs,tgts][y,]
      capvec.sc.perturb <- sample(capvec.sc.perturb,ntgts,replace = F)
      cor(capmat.bulk[TFs,tgts][y,],capvec.sc.perturb,method="pearson")
    })
    median(allcor)
  }) 
})
names(corWithBulk.list) <- subtypes

library(ggplot2)
library(ggridges)

colorVec <- c("brown3","cornflowerblue","mediumseagreen","goldenrod1")
data.df <- data.frame("Correlation"=c(corWithBulk.list[[1]],
                                                 corWithBulk.list[[2]],
                                                 corWithBulk.list[[3]],
                                                 corWithBulk.list[[4]]),
                                 "Subtype"=factor(c(rep("Classical",length(corWithBulk.list[[1]])),
                                                    rep("Neural",length(corWithBulk.list[[2]])),
                                                    rep("Proneural",length(corWithBulk.list[[3]])),
                                                    rep("Mesenchymal",length(corWithBulk.list[[4]]))),
                                                  levels=(c("Mesenchymal","Proneural","Neural","Classical"))))
colvec <- colorVec
#colvec[as.numeric(sigTFs[TFoI])] <- subtypeColors[as.numeric(sigTFs[TFoI])]
ggplot(data.df, aes(x = Correlation, y = Subtype,fill=Subtype)) + geom_density_ridges2() +
  scale_fill_manual(values=rev(colvec)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

colorVec <- c("brown3","cornflowerblue","mediumseagreen","goldenrod1")
data.df <- data.frame("Correlation"=c(corWithBulk.bg.list[[1]],
                                      corWithBulk.bg.list[[2]],
                                      corWithBulk.bg.list[[3]],
                                      corWithBulk.bg.list[[4]]),
                      "Subtype"=factor(c(rep("Classical",length(corWithBulk.bg.list[[1]])),
                                         rep("Neural",length(corWithBulk.bg.list[[2]])),
                                         rep("Proneural",length(corWithBulk.bg.list[[3]])),
                                         rep("Mesenchymal",length(corWithBulk.bg.list[[4]]))),
                                       levels=(c("Mesenchymal","Proneural","Neural","Classical"))))
colvec <- colorVec
#colvec[as.numeric(sigTFs[TFoI])] <- subtypeColors[as.numeric(sigTFs[TFoI])]
ggplot(data.df, aes(x = Correlation, y = Subtype,fill=Subtype)) + geom_density_ridges2() +
  scale_fill_manual(values=rev(colvec)) + xlim(c(-0.6,0.6)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

for (i in c(1:length(subtypes))) {
  fn <- paste0("analyses/0_0.1/",subtypes[i],"_scRegrResCorWithBulk.pdf")
  pdf(fn,width=6,height=5)
  hist(corWithBulk.list[[i]][abs(corWithBulk.list[[i]])<0.9999999],20,col=colorVec[i],main=subtypes[i],
       xlab="Pearson correlation",freq=F)
  dev.off()
}
for (i in c(1:length(subtypes))) {
  print(wilcox.test(corWithBulk.list[[i]][abs(corWithBulk.list[[i]])<0.9999999]))
}
# hist(corWithBulk.list[[1]][abs(corWithBulk.list[[1]])!=1],20,col="grey",main="Classical",
#      xlab="Avana gene effect")
# hist(corWithBulk.list[[2]][abs(corWithBulk.list[[2]])!=1],20,col="grey",main="Neural",
#      xlab="Avana gene effect")
# hist(corWithBulk.list[[3]][abs(corWithBulk.list[[3]])!=1],20,col="grey",main="Proneural",
#      xlab="Avana gene effect")
# hist(corWithBulk.list[[4]][abs(corWithBulk.list[[4]])!=1],20,col="grey",main="Mesenchymal",
#      xlab="Avana gene effect")

vioplot(corWithBulk.list[[1]][round(abs(corWithBulk.list[[1]]),5)!=1],
        corWithBulk.list[[2]][round(abs(corWithBulk.list[[2]]),5)!=1],
        corWithBulk.list[[3]][round(abs(corWithBulk.list[[3]]),5)!=1],
        corWithBulk.list[[4]][round(abs(corWithBulk.list[[4]]),5)!=1],
        col="grey")

for (subtype in subtypes) {
  fn <- paste0("analyses/0_0.1/",subtype,"_scTFFVecCorWithBulk.txt")
  write.table(corWithBulk.list[[subtype]][abs(round(corWithBulk.list[[subtype]],5))<1],
              fn,sep="\t",row.names=F,col.names=F,quote=F)
}

# Latest update: to explain whether differential partnering
# may explain the altered behavior of TFs, look at the correlation
# of F values specifically at the differentially regulated genes
# in that signature subtype and compare with that in other subtypes
# for a given TF. 

# Separate signature target genes into F-low and F-high (relative to
# overall average) ones.

for (sigTF in names(sigTFs)) {
  sigSubtype.sigTF <- as.numeric(sigTFs[sigTF])
  TFoI <- sigTF
  print(TFoI)
  tgtGenes <- which(sigSubtype[TFoI,]==sigSubtype.sigTF)
  heatmap.data <- lapply(capacities.log2.list,function(x){
    x[TFoI,names(tgtGenes)]
  })
  heatmap.data <- data.matrix(do.call(rbind.data.frame,heatmap.data))
  rownames(heatmap.data) <- subtypes
  colnames(heatmap.data) <- names(tgtGenes)
  sig.F.low <- apply(heatmap.data,2,function(x){
    x[sigSubtype.sigTF] < mean(x)
  })
  sig.F.high <- !sig.F.low
  heatmap.F.low <- heatmap.data[,sig.F.low,drop=F]
  heatmap.F.high <- heatmap.data[,sig.F.high,drop=F]
  write.table(heatmap.data,paste0("analyses/0_0.1/",TFoI,"_sigTgts.txt"),sep='\t',
              quote=F)
  if (ncol(heatmap.F.low) >= 5) {
    write.table(heatmap.F.low,paste0("analyses/0_0.1/",TFoI,"_sigTgts.F.low.txt"),sep='\t',
                quote=F)
    write.table(match(colnames(heatmap.F.low),colnames(capacities.log2.list[[1]]))-1,
                paste0("analyses/0_0.1/",TFoI,"_sigTgts.F.low.idx.txt"),sep="\t",
                quote=F,row.names=F,col.names=F)
    system(paste0("/Users/yunpengl/anaconda2/bin/python src/computeTFCor_RANSAC_wrapper.py 0.1_0 ",
                  sigTF," ",as.character(match(sigTF,rownames(capacities.list[[1]])))," ",
                  paste0("analyses/0_0.1/",TFoI,"_sigTgts.F.low.idx.txt "),"low"))
  }
  if (ncol(heatmap.F.high) >= 5) {
    write.table(heatmap.F.high,paste0("analyses/0_0.1/",TFoI,"_sigTgts.F.high.txt"),sep='\t',
                quote=F)
    write.table(match(colnames(heatmap.F.high),colnames(capacities.log2.list[[1]]))-1,
                paste0("analyses/0_0.1/",TFoI,"_sigTgts.F.high.idx.txt"),sep="\t",
                quote=F,row.names=F,col.names=F)
    system(paste0("/Users/yunpengl/anaconda2/bin/python src/computeTFCor_RANSAC_wrapper.py 0.1_0 ",
                  sigTF," ",as.character(match(sigTF,rownames(capacities.list[[1]])))," ",
                  paste0("analyses/0_0.1/",TFoI,"_sigTgts.F.high.idx.txt "),"high"))
  }
}

corOverSigTgt.low.list <- list()
corOverSigTgt.high.list <- list()

for (sigTF in names(sigTFs)) {
  sigSubtype.sigTF <- as.numeric(sigTFs[sigTF])
  TFoI <- sigTF
  print(TFoI)
  tgtGenes <- which(sigSubtype[TFoI,]==sigSubtype.sigTF)
  heatmap.data <- lapply(capacities.log2.list,function(x){
    x[TFoI,names(tgtGenes)]
  })
  heatmap.data <- data.matrix(do.call(rbind.data.frame,heatmap.data))
  rownames(heatmap.data) <- subtypes
  colnames(heatmap.data) <- names(tgtGenes)
  sig.F.low <- apply(heatmap.data,2,function(x){
    x[sigSubtype.sigTF] < mean(x)
  })
  sig.F.high <- !sig.F.low
  heatmap.F.low <- heatmap.data[,sig.F.low,drop=F]
  heatmap.F.high <- heatmap.data[,sig.F.high,drop=F]
  if (ncol(heatmap.F.low) >= 5) {
    corOverSigTgt.low <- read.delim(paste0('analyses/0_0.1/',
                                           sigTF,'_sigTgtCormat_low_0.1_0.txt'),
                                    sep="\t",header=F)
    corOverSigTgt.low <- data.matrix(corOverSigTgt.low)
    corOverSigTgt.low.list[[sigTF]] <- corOverSigTgt.low
  }
  if (ncol(heatmap.F.high) >= 5) {
    corOverSigTgt.high <- read.delim(paste0('analyses/0_0.1/',
                                           sigTF,'_sigTgtCormat_high_0.1_0.txt'),
                                    sep="\t",header=F)
    corOverSigTgt.high <- data.matrix(corOverSigTgt.high)
    corOverSigTgt.high.list[[sigTF]] <- corOverSigTgt.high
  }
}

corOverSigTgt.low.list.cleaned <- lapply(names(corOverSigTgt.low.list),function(x){
  cormat <- corOverSigTgt.low.list[[x]]
  subtype <- as.numeric(sigTFs[x])
  colnames(cormat) <- rownames(capacities.list[[1]])
  cormat[is.nan(cormat)] <- 0
  cormat <- cormat[,colSums(cormat)!=0]
  colMinusOneSDReduction <- apply(cormat,2,function(y){
    sd(y)-sd(y[-subtype])
  })
  cormat <- cormat[,rev(order(colMinusOneSDReduction))[1:max(c(100,ncol(cormat)))]]
  cormat
})
names(corOverSigTgt.low.list.cleaned) <- names(corOverSigTgt.low.list)

corOverSigTgt.high.list.cleaned <- lapply(names(corOverSigTgt.high.list),function(x){
  cormat <- corOverSigTgt.high.list[[x]]
  subtype <- as.numeric(sigTFs[x])
  colnames(cormat) <- rownames(capacities.list[[1]])
  cormat[is.nan(cormat)] <- 0
  cormat <- cormat[,colSums(cormat)!=0]
  colMinusOneSDReduction <- apply(cormat,2,function(y){
    sd(y)-sd(y[-subtype])
  })
  cormat <- cormat[,rev(order(colMinusOneSDReduction))[1:max(c(150,ncol(cormat)))]]
  cormat
})
names(corOverSigTgt.high.list.cleaned) <- names(corOverSigTgt.high.list)

TFoI <- "FOXO3"
pheatmap(corOverSigTgt.high.list.cleaned[[TFoI]],
         cluster_rows = F)
boxplot(corOverSigTgt.low.list.cleaned[[TFoI]][1,],
        corOverSigTgt.low.list.cleaned[[TFoI]][2,],
        corOverSigTgt.low.list.cleaned[[TFoI]][3,],
        corOverSigTgt.low.list.cleaned[[TFoI]][4,],col="gray")

boxplot(corOverSigTgt.high.list.cleaned[[TFoI]][1,],
        corOverSigTgt.high.list.cleaned[[TFoI]][2,],
        corOverSigTgt.high.list.cleaned[[TFoI]][3,],
        corOverSigTgt.high.list.cleaned[[TFoI]][4,],col="gray")

library(ggplot2)
library(ggridges)

subtypeColors <- c("brown3","cornflowerblue","mediumseagreen","goldenrod1" )
for (TFoI in names(corOverSigTgt.low.list.cleaned)) {
  outfn <- paste0("analyses/0_0.1/sigTFCorAtDiffRegGenes/",TFoI,"corOverSigTgt_low.pdf")
  data.df <- data.frame("Correlation"=c(corOverSigTgt.low.list.cleaned[[TFoI]][1,],
                                        corOverSigTgt.low.list.cleaned[[TFoI]][2,],
                                        corOverSigTgt.low.list.cleaned[[TFoI]][3,],
                                        corOverSigTgt.low.list.cleaned[[TFoI]][4,]),
                        "Subtype"=factor(c(rep("Classical",ncol(corOverSigTgt.low.list.cleaned[[TFoI]])),
                                    rep("Neural",ncol(corOverSigTgt.low.list.cleaned[[TFoI]])),
                                    rep("Proneural",ncol(corOverSigTgt.low.list.cleaned[[TFoI]])),
                                    rep("Mesenchymal",ncol(corOverSigTgt.low.list.cleaned[[TFoI]]))),
                                    levels=c("Mesenchymal","Proneural","Neural","Classical")))
  colvec <- rep("gray",4)
  colvec[as.numeric(sigTFs[TFoI])] <- subtypeColors[as.numeric(sigTFs[TFoI])]
  ggplot(data.df, aes(x = Correlation, y = Subtype,fill=Subtype,height=..scaled..)) + geom_density_ridges2(scale=1.25,stat='density') +
    scale_fill_manual(values=rev(colvec)) +
    theme(panel.background = element_blank()) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + ggtitle(TFoI) + xlim(c(-1,1)) + scale_y_discrete(expand = c(0.01, 0))
  ggsave(file = outfn, width = 10, height = 4.5)
}

subtypeColors <- c("brown3","cornflowerblue","mediumseagreen","goldenrod1" )
for (TFoI in names(corOverSigTgt.high.list.cleaned)) {
  outfn <- paste0("analyses/0_0.1/sigTFCorAtDiffRegGenes/",TFoI,"corOverSigTgt_high.pdf")
  data.df <- data.frame("Correlation"=c(corOverSigTgt.high.list.cleaned[[TFoI]][1,],
                                        corOverSigTgt.high.list.cleaned[[TFoI]][2,],
                                        corOverSigTgt.high.list.cleaned[[TFoI]][3,],
                                        corOverSigTgt.high.list.cleaned[[TFoI]][4,]),
                        "Subtype"=factor(c(rep("Classical",ncol(corOverSigTgt.high.list.cleaned[[TFoI]])),
                                           rep("Neural",ncol(corOverSigTgt.high.list.cleaned[[TFoI]])),
                                           rep("Proneural",ncol(corOverSigTgt.high.list.cleaned[[TFoI]])),
                                           rep("Mesenchymal",ncol(corOverSigTgt.high.list.cleaned[[TFoI]]))),
                                         levels=c("Mesenchymal","Proneural","Neural","Classical")))
  colvec <- rep("gray",4)
  colvec[as.numeric(sigTFs[TFoI])] <- subtypeColors[as.numeric(sigTFs[TFoI])]
  ggplot(data.df, aes(x = Correlation, y = Subtype,fill=Subtype,height=..scaled..)) + geom_density_ridges2(scale=1.25,stat='density') +
    scale_fill_manual(values=rev(colvec)) +
    theme(panel.background = element_blank()) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + ggtitle(TFoI) + xlim(c(-1,1)) + scale_y_discrete(expand = c(0.01, 0))
  ggsave(file = outfn, width = 10, height = 4.5)
}

# Now extract signatures again focusing on the differentially regulated
# target genes.

sigTFsigTgtCorSig.low.list <- lapply(names(corOverSigTgt.low.list.cleaned),function(sigTF){
  cormat <- corOverSigTgt.low.list.cleaned[[sigTF]]
  leaveOneSDReduction <- apply(cormat,2,function(x){
    SD.all <- sd(x)
    sapply(c(1:length(subtypes)),function(i) {
      SD.all - sd(x[-i])
    })
  })
  leaveOneSDReduction <- leaveOneSDReduction[,which(apply(cormat,2,function(x){
    max(abs(x)) > 0.5
  }))]
  sigSubtypeVec <- apply(leaveOneSDReduction,2,function(x){
    which(x==max(x))[1]
  })
  sigSubtypeCounts <- table(sigSubtypeVec)
  #names(sigSubtypeCounts) <- subtypes
  sigSubtype <- names(sigSubtypeCounts)[which(sigSubtypeCounts==max(sigSubtypeCounts))[1]]
  sigSubtypePercentage <- sigSubtypeCounts[sigSubtype]/sum(sigSubtypeCounts)
  list(subtypes[as.numeric(sigSubtype)],sigSubtypePercentage)
})
sigTFsigTgtCorSig.low.df <- data.frame(sigSubtype=unlist(sapply(sigTFsigTgtCorSig.low.list,
                                                                function(x){
                                                                  if (length(x[[1]])>0) {
                                                                    x[[1]]
                                                                  } else {
                                                                    "NA"
                                                                  }
                                                                  })),
                                       sigSubtypePercentage=unlist(sapply(sigTFsigTgtCorSig.low.list,
                                                                          function(x){
                                                                            if (length(x[[2]])>0) {
                                                                              x[[2]]
                                                                            } else {
                                                                              "NA"
                                                                            }
                                                                            })))

sigTFsigTgtCorSig.high.list <- lapply(names(corOverSigTgt.high.list.cleaned),function(sigTF){
  cormat <- corOverSigTgt.high.list.cleaned[[sigTF]]
  leaveOneSDReduction <- apply(cormat,2,function(x){
    SD.all <- sd(x)
    sapply(c(1:length(subtypes)),function(i) {
      SD.all - sd(x[-i])
    })
  })
  leaveOneSDReduction <- leaveOneSDReduction[,which(apply(cormat,2,function(x){
    max(abs(x)) > 0.5
  }))]
  sigSubtypeVec <- apply(leaveOneSDReduction,2,function(x){
    which(x==max(x))[1]
  })
  sigSubtypeCounts <- table(sigSubtypeVec)
  #names(sigSubtypeCounts) <- subtypes
  sigSubtype <- names(sigSubtypeCounts)[which(sigSubtypeCounts==max(sigSubtypeCounts))[1]]
  sigSubtypePercentage <- sigSubtypeCounts[sigSubtype]/sum(sigSubtypeCounts)
  list(subtypes[as.numeric(sigSubtype)],sigSubtypePercentage)
})
sigTFsigTgtCorSig.high.df <- data.frame(sigSubtype=unlist(sapply(sigTFsigTgtCorSig.high.list,
                                                                function(x){
                                                                  if (length(x[[1]])>0) {
                                                                    x[[1]]
                                                                  } else {
                                                                    "NA"
                                                                  }
                                                                })),
                                       sigSubtypePercentage=unlist(sapply(sigTFsigTgtCorSig.high.list,
                                                                          function(x){
                                                                            if (length(x[[2]])>0) {
                                                                              x[[2]]
                                                                            } else {
                                                                              "NA"
                                                                            }
                                                                          })))

sigTFCorSubtype <- matrix(0,nrow(geneRegs.TFs.coreg[[1]]),ncol(geneRegs.TFs.coreg[[1]]))

for (i in c(1:nrow(geneRegs.TFs.coreg[[1]]))) {
  for (j in c(1:ncol(geneRegs.TFs.coreg[[1]]))) {
    Fval <- sapply(leaveOneSD,function(x){
      x[i,j]
    })
    coregs.allSubtypes <- sapply(c(1:length(subtypes)),function(k){
      geneRegs.TFs.coreg[[k]][i,j]
    })
    if (overallSD[i,j] > 0 & max(abs(coregs.allSubtypes)) > 0.5) {
      sigTFCorSubtype[i,j] <- which(Fval==min(Fval))[1]
    }
  }
}
rownames(sigTFCorSubtype) <- rownames(geneRegs.TFs.coreg[[1]])
colnames(sigTFCorSubtype) <- colnames(geneRegs.TFs.coreg[[1]])
