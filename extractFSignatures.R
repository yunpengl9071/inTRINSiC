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

lambda <- "0.1_0"
capacities.list <- lapply(subtypes,function(x){
  fn <- paste0("data/regrOutput/paramSweep/nonlinear/",x,"_TFGeneCapacities_ri_nomi__lambda",lambda,".txt")
  temp <- data.matrix(read.delim(fn,header=T,row.names=1,check.names=F))
  constid <- which(rownames(temp)=="const")
  temp[-constid,]
})
names(capacities.list) <- subtypes

# # Use Neural results from rerun with different random seed
# temp <- data.matrix(read.delim("Neural_TFGeneCapacities_ri_nomi__lambda0_0_newRNSD.txt"))
# constid <- which(rownames(temp)=="const")
# temp <- temp[-constid,]
# 
# # 03/10/18: check Neural results against a re-run
# capacities.neural.new <- read.delim("Neural_TFGeneCapacities_ri_nomi_DHSFANTOM_new__lambda0_0.txt",
#                                     header=T,row.names=1,check.names=F)
# constid <- which(rownames(capacities.neural.new)=="const")
# capacities.neural.new <- data.matrix(capacities.neural.new[-constid,])
# # The two matrices are exactly the same.
# # There might be some systematic bias in the neural data (or the other subtypes)
# # Or - could it be a bad random seed?
# diff <- capacities.neural.new - capacities.list$Neural
# 
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
  allexprs <- data.matrix(expr)[tgt,]
  allexprs.TF <- data.matrix(expr)[TF,]
  llim <- min(allexprs)
  ulim <- max(allexprs)
  llimy <- min(allexprs.TF)
  ulimy <- max(allexprs.TF)
  plot(TF.exprs,tgt.exprs,xlab=TF,ylab=tgt,main=subtype,pch=16,ylim=c(llim-0.1,ulim+0.1),
       xlim=c(llimy-0.1,ulimy+0.1))
}
plotTFTgt <- function(TF,tgt,subtype,xl=NA,yl=NA) {
  TF.exprs <- data.matrix(expr)[TF,names(commonLabels)[commonLabels==subtype]]
  tgt.exprs <- data.matrix(expr)[tgt,names(commonLabels)[commonLabels==subtype]]
  if (is.na(xl)) {
    allexprs <- data.matrix(expr)[tgt,]
    llim <- min(allexprs)
    ulim <- max(allexprs)
  } else {
    llim <- xl[1]
    ulim <- xl[2]
  }
  if (is.na(yl)) {
    allexprs.TF <- data.matrix(expr)[TF,]
    llimy <- min(allexprs.TF)
    ulimy <- max(allexprs.TF)
  } else {
    llimy <- yl[1]
    ulimy <- yl[2]
  }
  plot(TF.exprs,tgt.exprs,xlab=TF,ylab=tgt,main=subtype,pch=16,ylim=c(llim-0.1,ulim+0.1),
       xlim=c(llimy-0.1,ulimy+0.1))
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

lambda <- "0.1_0"
capacities.list <- lapply(subtypes,function(x){
  fn <- paste0("data/regrOutput/paramSweep/nonlinear/",x,"_TFGeneCapacities_ri_nomi__lambda",lambda,".txt")
  temp <- data.matrix(read.delim(fn,header=T,row.names=1,check.names=F))
  constid <- which(rownames(temp)=="const")
  temp[-constid,]
})
names(capacities.list) <- subtypes


capacities.log2.list <- lapply(capacities.list,log2)

niter <- 30
capmat.unshuffled <- lapply(capacities.list,function(x){
  x[allTFs,allGenes]
})

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

# Signature mining method 1: a signature TF-gene pair should only show up as
# subtype-specific in one subtype

# Signature mining method 2: a signature TF-gene pair should display an F value in one subtype
# that is significantly different from that in the other three subtypes - this makes more sense.
# Should start with background-corrected (through label shuffling) subtype-specific F matrices,
# then use the SD reduction method.
# Also filter the final results by log2 fold change in the F values.

subtypeSpecificMats.bgAdj <- FValueSig.res$subtypeSpecificFMats
names(subtypeSpecificMats.bgAdj) <- subtypes
temp <- lapply(subtypes,function(x){
  fn <- paste0("data/",x,"_",lambda,"_subtypeSpecificMats.bgAdj.txt")
  write.table(subtypeSpecificMats.bgAdj[[x]],fn,sep="\t",quote=F)
})


subtypeSpecificMats.bgAdj <- lapply(subtypes,function(x){
  fn <- paste0("data/",x,"_",lambda,"_subtypeSpecificMats.bgAdj.txt")
  data.matrix(read.delim(fn,sep="\t",header=T,row.names=1,check.names=F))
})
names(subtypeSpecificMats.bgAdj) <- subtypes
geneRegs.all <- lapply(subtypeSpecificMats.bgAdj,log2)
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
geneRegs.TF.mean <- geneRegs.mean[allTFs,]
leaveOneSD <- lapply(c(1:length(subtypes)),function(x){
  currList <- geneRegs.TF[-x]
  LOM <- Reduce("+",currList)
  LOM <- LOM / (length(subtypes)-1)
  LOV <- lapply(currList,function(y){(y-LOM)^2})
  LOV <- Reduce("+",LOV) / (length(subtypes)-1)
  LOSD <- LOV^0.5
  LOSD
})
overallSD <- geneRegs.TF.sd
# Use the amount of reduction in variance as measurement for strength of signature
# Also set a cutoff for log2 change in F values
# Seems like the oddity in the number of Neural signature TF-gene pairs decreases
# as the threshold becomes more stringent. Points to potential bias in data or poor
# regression performance.
log2FDiffThreshold <- 5
sigSubtype <- matrix(0,nrow(geneRegs.TF[[1]]),ncol(geneRegs.TF[[1]]))
SDReduction.mat <- matrix(0,nrow(geneRegs.TF[[1]]),ncol(geneRegs.TF[[1]]))
for (i in c(1:nrow(geneRegs.TF[[1]]))) {
  for (j in c(1:ncol(geneRegs.TF[[1]]))) {
    Fval <- sapply(leaveOneSD,function(x){
      x[i,j]
    })
    if (overallSD[i,j] > 0) {
      Fvec <- sapply(geneRegs.TF,function(x){
        x[i,j]
      })
      ss <- which(Fval==min(Fval))[1]
      log2FDiff <- abs(Fvec[ss]-mean(Fvec[-ss]))
      if (log2FDiff > log2FDiffThreshold) {
        SDReduction.mat[i,j] <- overallSD[i,j]-Fval[ss]
        sigSubtype[i,j] <- ss
      }
    }
  }
}
rownames(sigSubtype) <- rownames(geneRegs.TF[[1]])
colnames(sigSubtype) <- colnames(geneRegs.TF[[1]])
write.table(sigSubtype,paste0("data/",lambda,"_sigSubtype.mat.txt"),
            sep="\t",quote=F)
sigSubtype <- data.matrix(read.delim(paste0("data/",lambda,"_sigSubtype.mat.txt"),
                         header=T,row.names = 1,check.names=F))

# Filter signature subtype matrix by percentage of rewired targets as well
# as out-degree distribution of TFs
outDegrees <- sapply(unionNet.byTF,nrow)
outDegrees <- outDegrees[allTFs]
sigOutDegrees <- (outDegrees >= 50)
percentSigCutoff <- 0.000
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
GSE57872.all <- read.delim("data/GSE57872_GBM_data_matrix.txt",header=T,row.names=1)
GSE57872.all <- data.matrix(GSE57872.all)
GSE57872.topExprsGenes <- rev(order(rowMeans(GSE57872.all)))[1:3000]

# write.table(domSigSubtype.filtered[domSigSubtype.filtered!="0"][domSigSubtype.filtered.multipleTgts],
#             "domSigSubtype.filtered_percentCutoff0.005.multipleTarget.txt",sep="\t",quote=F)
sigTFs <- domSigSubtype.filtered[intersect(rownames(GSE57872.all)[GSE57872.topExprsGenes],
                                 names(domSigSubtype.filtered[domSigSubtype.filtered!="0"][domSigSubtype.filtered.multipleTgts]))]
# '0.005' is outdated in the fn below
fn <- "analyses/0_0.1/domSigSubtype.filtered_percentCutoff0.000.multipleTargets5.withSCExprsTop3000.txt"
write.table(sigTFs,
            fn,
            col.names=F,sep="\t",quote=F)
sigTFs <- (read.delim(fn,header=F,check.names=F))
sigTFNames <- as.character(sigTFs[,1])
sigTFs <- as.character(sigTFs[,2])
names(sigTFs) <- sigTFNames

# Plot a heatmap of signature participation for TFs that show single cell expression
sigSubtype.sigTFs <- sigSubtype[names(sigTFs),]

sigSubtype.participation <- apply(sigSubtype.sigTFs,1,function(x){
  temp <- x[x!='0']
  tally <- table(x)
  for (i in c('1','2','3','4')) {
    if (!(i %in% names(tally))) {
      tally[i] <- 0
    }
  }
  tally <- tally[c('1','2','3','4')]
  tally/sum(tally)
})
rownames(sigSubtype.participation) <- subtypes
colnames(sigSubtype.participation) <- names(sigTFs)

myheatmapfun <- colorRampPalette((c(rgb(255,255,255,maxColorValue=255),
                                    rgb(234,145,118,maxColorValue=255),
                                    rgb(102,3,32,maxColorValue=255))))
library(pheatmap)
xxx <- pheatmap(sigSubtype.participation,breaks=seq(0.2,0.8,length.out=76),
                clustering_method = "ward.D",cluster_rows = F,
                clustering_distance_cols = "correlation",
                color = viridis(75),dendrogram=F,fontsize=12)
subtype <- 'Classical'
for (subtype in subtypes) {
  sigSubtype.participation.subtype <- sigSubtype.participation[,
                                                               names(sigTFs)
                                                               [sigTFs==as.character(match(subtype,
                                                                                           subtypes))]]
  pheatmap(sigSubtype.participation.subtype,breaks=seq(0.2,0.6,length.out=76),
           clustering_method = "ward.D",cluster_rows = F,
           clustering_distance_cols = "correlation",
           color = viridis(75),dendrogram=F,fontsize=12,
           treeheight_row = 0, treeheight_col = 0,
           cellwidth = 20,cellheight = 20)
}

# Example plots of signatures

TFoI <- "MXI1"
# TFoI <- "MEIS1"
tgtGenes <- which(sigSubtype[TFoI,]==as.numeric(domSigSubtype.filtered[TFoI]))
heatmap.data <- lapply(capacities.log2.list,function(x){
  x[TFoI,names(tgtGenes)]
})
heatmap.data <- data.matrix(do.call(rbind.data.frame,heatmap.data))
rownames(heatmap.data) <- subtypes
colnames(heatmap.data) <- names(tgtGenes)
# write.table(heatmap.data,paste0("analyses/0_0.1/",TFoI,"_sigTgts.txt"),sep='\t',
#             quote=F)
library(pheatmap)
xxx <- pheatmap(heatmap.data,breaks=seq(-10,10,length.out=76),
                clustering_method = "ward.D",cluster_rows = F,
                clustering_distance_cols = "euclidean",
                color = myheatmapfun(75),dendrogram=F)

heatmap.data <- data.matrix(do.call(rbind.data.frame,heatmap.data))
rownames(heatmap.data) <- subtypes
colnames(heatmap.data) <- names(tgtGenes)
#dev.off()
hmRes <- heatmap.2(heatmap.data,
                   #Colv=NA,
                   Rowv=NA,
                   scale="none",
                   density.info="none", trace="none",
                   col=myheatmapfun(75),
                   reorderfun = function(d,w) { d },
                   cexCol = 0.5,
                   srtCol = 45,
                   cexRow = 0.7,
                   offsetCol = -0.5,
                   #RowSideColors = rscVec,
                   dendrogram="none",
                   key.title=NA,
                   key.xlab=NA,
                   breaks=seq(-3,3,length.out=76),
                   margins=c(10,10),
                   lhei=c(1,8))

tgt <- "ITGA5"
plotTFTgt(TFoI,tgt,subtype = "Classical")
plotTFTgt(TFoI,tgt,subtype = "Neural")
plotTFTgt(TFoI,tgt,subtype = "Proneural")
plotTFTgt(TFoI,tgt,subtype = "Mesenchymal")

# Use density coutours to make trends more visible.

library(ggplot2)
subtypeColors <- c("brown3","cornflowerblue","mediumseagreen","goldenrod1")
names(subtypeColors) <- subtypes
subtype <- subtypes[4]
TFvsTgt.exprs <- cbind.data.frame(data.matrix(expr)[TFoI,names(labels)[labels==subtype]],
                                  data.matrix(expr)[tgt,names(labels)[labels==subtype]])
fn <- paste0("analyses/0_0.1/",TFoI,"_",tgt,"_",subtype,"_exprs.txt")
write.table(TFvsTgt.exprs,fn,sep="\t",col.names=F,row.names=F,quote=F)

# Get RANSAC inliers and regression
getRANSACTFGeneExprs <- function(TFoI,tgt,subtype) {
  fn <- paste0("analyses/0_0.1/",TFoI,"_",tgt,"_",subtype,"_exprs_RANSAC.txt")
  temp <- data.matrix(read.delim(fn,header=F,check.names=F))
  data.frame(x=temp[,1],y=temp[,2])
}
# fn <- paste0("analyses/0_0.1/",TFoI,"_",tgt,"_exprs_RANSAC.txt")
# TFvsTgt.exprs.inliers <- getRANSACTFGeneExprs(TFoI,tgt,subtypes[1])
# TFvsTgt.regression <- lm(y~x,data=TFvsTgt.exprs.inliers)
RANSAC_TFTgt <- lapply(subtypes,getRANSACTFGeneExprs,TF=TFoI,tgt=tgt)
RANSAC_regression <- lapply(RANSAC_TFTgt,function(df){
  lm(y~x,data=df)
})
names(RANSAC_TFTgt) <- subtypes
names(RANSAC_regression) <- subtypes
#colnames(TFvsTgt.exprs.inliers) <- c(TFoI,tgt)
# ggplot(TFvsTgt.exprs, aes(x=TFoI,y=tgt)) + geom_bin2d(bins=20) + theme_bw() +
#   scale_fill_gradientn(limits=c(0,10), breaks=seq(0, 40, by=1), colours=c("white","pink"))
# ggplot(TFvsTgt.exprs, aes(x=TFoI,y=tgt)) + geom_density_2d(bins=10) + theme_bw() 
library(spatstat)
#pdf(file = paste0("analyses/0_0.1/",TFoI,"_",tgt,"_",subtype,".pdf"),width = 10,height = 4)
# par(mar = c(1, 1, 1, 1))
# pppo=ppp(x=TFvsTgt.exprs[,1],y=TFvsTgt.exprs[,2],window = owin(c(6,12),c(5,11)))
# den=density(pppo,kernel="gaussian",edge=T,diggle=T,adjust=0.38)
# plot(den,main='TF vs target',col=colorRampPalette(c("white",subtypeColors[subtype]))(10),xlim=c(6,12),
#      ylim=c(5,11),xlab=TFoI,ylab=tgt)
# points(TFvsTgt.exprs,xlim=c(6,12),ylim=c(5,11),pch=20,col="black")
#dev.off()

subtype <- subtypes[4]
TFvsTgt.exprs <- cbind.data.frame(data.matrix(expr)[TFoI,names(labels)[labels==subtype]],
                                  data.matrix(expr)[tgt,names(labels)[labels==subtype]])
colnames(TFvsTgt.exprs) <- c(TFoI,tgt)
subtypePalettes <- c('Reds','Blues','Greens','Oranges')
crp <- colorRampPalette(c("white",subtypeColors[subtype]))
names(subtypePalettes) <- subtypes
xlower <- min(as.numeric(expr[TFoI,]))
xupper <- max(as.numeric(expr[TFoI,]))
ylower <- min((as.numeric(expr[tgt,])))
yupper <- max((as.numeric(expr[tgt,])))
# data.df <- data.frame(x=as.numeric(availTCGAExprs.bySubtype[[subtype]][TF,]),
#                       y=log2(as.numeric(annot.df[annot.df[,"label"]==subtype,"survival"])))
ggplot(TFvsTgt.exprs,aes(x=TFvsTgt.exprs[,TFoI],y=TFvsTgt.exprs[,tgt])) + stat_density_2d(aes(fill = ..density..),
                                               geom="raster",
                                               contour=F) +
  scale_fill_gradient(low = "white",high=subtypeColors[subtype]) + 
  scale_x_continuous(limits=c(xlower-0.25,xupper+0.25),breaks=seq(floor(xlower),ceiling(xupper),2)) +
  scale_y_continuous(limits=c(ylower-0.25,yupper+0.25),breaks=seq(floor(ylower),ceiling(yupper),2)) +
  theme(legend.position = 'none',panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.background = element_blank(),
        panel.border = element_rect(fill = NA,size = 2)) + geom_point(color='black') + 
  geom_abline(intercept = RANSAC_regression[[subtype]]$coefficients[1],
              slope = RANSAC_regression[[subtype]]$coefficients[2],
              linetype = 'dashed',
              size = 1.25)

for (sigTF in names(sigTFs)) {
  tgtGenes <- which(sigSubtype[sigTF,]==as.numeric(domSigSubtype.filtered[sigTF]))
  fn <- paste0("analyses/0_0.1/",sigTF,"_sigTgtGenes.txt")
  write.table(names(tgtGenes),fn,sep="\t",row.names=F,col.names=F,quote=F)
}

# For the signature target genes of these signature TFs, are they physically close to each other?
geneAnnotByChr <- read.delim("~/MyGoogleDrive/HemannLab/cloudData/RNA-seq/151016_sorted_CBAB_reseq/geneAnnotTable.txt",
                             header=T,row.names=1)
tgtGenePos <- lapply(names(sigTFs),function(sigTF){
  tgtGenes <- which(sigSubtype[sigTF,]==as.numeric(domSigSubtype.filtered[sigTF]))
  loc.list <- lapply(names(tgtGenes),function(tg){
    idx <- match(tg,geneAnnotByChr$geneSymbol)
    as.character(geneAnnotByChr[idx,c(1,2)])
  })
  loc.df <- do.call(rbind.data.frame,loc.list)
  colnames(loc.df) <- c("chr","start")
  loc.df[complete.cases(loc.df),]
  loc.df[order(loc.df$chr),]
})
for (sigTF in names(sigTFs)) {
  tgtGenes <- which(sigSubtype[sigTF,]==as.numeric(domSigSubtype.filtered[sigTF]))
  fn <- paste0(sigTF,"_sigTgtGenes.txt")
  write.table(names(tgtGenes),fn,sep="\t",row.names=F,col.names=F,quote=F)
}

sigSubtype.tgtGeneUniv <- colnames(sigSubtype)[colSums(data.matrix(sigSubtype))!=0]

idmapper <- read.delim("data/hg19.IDMapper.txt",header=F)
sigSubtype.tgtGeneUniv <- as.character(idmapper[match(sigSubtype.tgtGeneUniv,
                                                      as.character(idmapper$V2)),"V1"])
write.table(sigSubtype.tgtGeneUniv,"sigSubtype.tgtGeneUniv.txt",
            sep="\t",col.names=F,row.names=F,quote=F)

sigTFs <- names(domSigSubtype.filtered[intersect(rownames(GSE57872.all),
                                                 names(domSigSubtype.filtered[domSigSubtype.filtered!="0"][domSigSubtype.filtered.multipleTgts]))])


# Make signature matrices for Prism.
heatmap.data <- sigSubtype.filtered[names(sigTFs),]
domSigSubtype <- apply(heatmap.data,1,function(x){
  temp <- table(x)
  temp <- temp[names(temp) != "0"]
  names(temp)[which(temp==max(temp))][1]
})
sigTFColors <- rep("red",length(domSigSubtype))
sigTFColors[domSigSubtype=="2"] <- "blue"
sigTFColors[domSigSubtype=="3"] <- "green"
sigTFColors[domSigSubtype=="4"] <- "orange"
names(sigTFColors) <- names(domSigSubtype)
signatureMat <- apply(heatmap.data,1,function(x){
  subtypeCounts <- table(x)
  subtypeCounts <- subtypeCounts[names(subtypeCounts)!="0"]
  subtypeCounts <- sort(subtypeCounts)
  temp <- lapply(names(subtypeCounts),function(i){
    c(x[x==i])
  })
  temp <- c(x[x=="0"],unlist(temp))
  temp
})
signatureMat <- t(signatureMat)
colnames(signatureMat) <- colnames(capmat[[1]])
signatureMat <- signatureMat[,colSums(signatureMat)>0]
signatureMat <- signatureMat[order(signatureMat[,ncol(signatureMat)]),]
# rowOrder <- c("POU3F2","TCF12","ATF4","EGR1","ZEB1",
#               "SOX9","TCF4")
signatureMat <- signatureMat
signatureMat.data <- lapply(1:nrow(signatureMat),function(x){
  temp <- table(signatureMat[x,])
  temp <- temp[names(temp)!="0"]
  for (i in as.character(1:length(subtypes))) {
    if (!(i %in% names(temp))) {
      temp[i] <- 0
    }
  }
  temp[as.character(1:length(subtypes))]
})
names(signatureMat.data) <- rownames(signatureMat)
signatureMat.data <- do.call(rbind.data.frame,signatureMat.data)
colnames(signatureMat.data) <- subtypes
signatureMat.data <- t(apply(data.matrix(signatureMat.data),1,function(x){
  x/sum(x)
}))
rownames(signatureMat.data) <- rownames(signatureMat)
signatureMat.data <- signatureMat.data[order(apply(signatureMat.data,1,max)),]
# Need to segregate into four subtypes for visual clarity.
sigMatList <- list()
for (i in 1:length(subtypes)) {
  signatureMat.data.subtype <- signatureMat.data[names(sigTFs)[sigTFs==as.character(i)],]
  signatureMat.data.subtype <- signatureMat.data.subtype[order(apply(signatureMat.data.subtype,1,max)),]
  write.table(signatureMat.data.subtype,paste0("analyses/0_0.1/F_signatureMat_",subtypes[i],"_Prism.txt"),sep="\t",
              quote=F)
  sigMatList[[i]] <- signatureMat.data.subtype
}
sigMatList <- do.call(rbind.data.frame,sigMatList)
write.table(sigMatList,"sigMats_Prism_orderedBySubtype.txt",quote=F,sep="\t")
write.table(signatureMat.data,"analyses/0_0.1/F_signatureMat_Prism.txt",sep="\t",
            quote=F)

dev.off()
hmRes <- heatmap.2(signatureMat,
                   Colv=NA,
                   Rowv=NA,
                   scale="none",
                   density.info="none", trace="none",
                   col=c("white","red","navyblue","darkgreen","orange"),
                   reorderfun = function(d,w) { d },
                   cexCol = 0.3,
                   srtCol = 45,
                   cexRow = 1,
                   offsetCol = -0.5,
                   rowsep = c(1:nrow(signatureMat)),
                   sepcolor = "white",
                   #RowSideColors = rscVec,
                   dendrogram="none",
                   key.title=NA,
                   key.xlab=NA,
                   lhei=c(1,8),
                   #margins=c(10,10),
                   labCol="")

# For correlation signatures, sometimes the correlation structure is more complex
# than a single cluster of correlated F values.
# Try calculating correlation in combination with clustering for TF pairs that have
# lots of non-zero correlations, or simply set a rectangular gate.
# Also kick out F value pairs where one of them is zero. Those make less sense.
lambda <- "0.1_0"
capacities.list <- lapply(subtypes,function(x){
  fn <- paste0("data/regrOutput/paramSweep/nonlinear/",x,"_TFGeneCapacities_ri_nomi__lambda",lambda,".txt")
  temp <- data.matrix(read.delim(fn,header=T,row.names=1,check.names=F))
  #constid <- which(rownames(temp)=="const")
  #temp[-constid,]
})
names(capacities.list) <- subtypes
capacities.log2.list <- lapply(capacities.list,log2)

TF1 <- "POU1F1"
TF2 <- "SMAD4"
subtype <- 2
ind <- ((capacities.log2.list[[subtype]][TF1,])>-5 & (capacities.log2.list[[subtype]][TF2,])>-5 &
          (capacities.log2.list[[subtype]][TF1,])<5 & (capacities.log2.list[[subtype]][TF2,])<5 &
          (capacities.log2.list[[subtype]][TF1,])!=0 & (capacities.log2.list[[subtype]][TF2,])!=0)

cor(capacities.log2.list[[subtype]][TF1,ind],
     capacities.log2.list[[subtype]][TF2,ind])
plot(capacities.log2.list[[subtype]][TF1,ind],
    capacities.log2.list[[subtype]][TF2,ind],pch=20,xlim=c(-5,5),ylim=c(-5,5))

for (subtype in 1:4) {
  TF1 <- "ETV7"
  TF2 <- "ARNT"
  
  ind <- ((capacities.log2.list[[subtype]][TF1,])!=0 & (capacities.log2.list[[subtype]][TF2,])!=0)
  
  cor(capacities.log2.list[[subtype]][TF1,ind],
      capacities.log2.list[[subtype]][TF2,ind])
  plot(capacities.log2.list[[subtype]][TF1,ind],
       capacities.log2.list[[subtype]][TF2,ind],pch=20,xlim=c(-2,2),ylim=c(-2,2))
  
}
TF1 <- "TFAP2C"
TF2 <- "SIX1"

#log2F.threshold <- 0
# library(mvoutlier)
# rectCutoff <- 5
# geneRegs.TFs.coreg <- lapply(capacities.log2.list, function(x){
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
#       if (i>j) {
#         tfi <- temp[i,]
#         tfj <- temp[j,]
#         tgts <- which(tfi != 0 & tfj != 0)
#         # tgts <- which(abs(tfi)>log2F.threshold & abs(tfj)>log2F.threshold)
#         if (length(tgts) > 5) {
#           df <- data.matrix(cbind(tfi[tgts],tfj[tgts]))
#           if (length(tgts) <= 10) {
#             cormat[i,j] <- cor(df[,1],df[,2],method="pearson")
#           } else {
#             # Or, use kmeans clustering
#             # Empirically set ncluster to 3 and choose largest correlation
#             #kmeans.res <- kmeans(df,centers = 3,iter.max = 30)
#             #hclust.res <- hclust(dist((df)),method = "ward.D2")
#             #hclust.res <- cutree(hclust.res,k=5)
#             #print(table(kmeans.res$cluster))
#             # kmeans.res.list <- sapply((c(1:3)),function(cid){
#             #   df.cluster <- df[(kmeans.res$cluster)==cid,,drop=F]
#             #   #print(nrow(df.cluster))
#             #   if (nrow(df.cluster) > 10) {
#             #     cor(df.cluster[,1],df.cluster[,2],method="spearman")
#             #   } else {
#             #     0
#             #   }
#             # })
#             # kmeans.res.len <- sapply(c(1:3),function(cid){
#             #   sum((kmeans.res$cluster)==cid)
#             # })
#             # cormat[i,j] <- kmeans.res.list[which(kmeans.res.len==max(kmeans.res.len))[1]]
#             gateIdx <- which(abs(df[,1]) < rectCutoff & abs(df[,2]) < rectCutoff)
#             if (length(gateIdx) > 10) {
#               cormat[i,j] <- cor(df[gateIdx,1],df[gateIdx,2],method="pearson")
#             } else {
#               cormat[i,j] <- 0
#             }
#           }
#           
#           # pcout.res <- pcout(df)
#           # pcout.scores <- pcout.res$wfinal
#           # # Kick out two most likely outliers
#           # pcout.lowScores <- order(pcout.scores)[1:2]
#           # # df <- df[pcout(df)$wfinal01 > 0,]
#           # df <- df[-pcout.lowScores,]
#           # if (nrow(df) > 3) {
#           #   cormat[i,j] <- cor(df[,1],df[,2],method="pearson")
#           # } else {
#           #   cormat[i,j] <- 0
#           # }
#         } else {
#           cormat[i,j] <- 0
#         }
#       }
#     }
#   }
#   cat("\n")
#   cormat
# })
# 
# for (i in c(1:length(subtypes))) {
#   fn <- paste0("analyses/0_0.1/",subtypes[i],"_TF_TF_cormat_lambda_0_0.1_rectGate.txt")
#   write.table(geneRegs.TFs.coreg[[i]],fn,sep="\t",col.names=T,row.names=T,quote=F)
# }

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
      examplePairs[[paste0(allTFs[i],"_",allTFs[j])]] <- corVec
      exampleTuples[[paste0(allTFs[i],"_",allTFs[j])]] <- paste0(allTFs[i],"_",allTFs[j])
    }
  }
}
examplePairs <- do.call(rbind.data.frame,examplePairs)
colnames(examplePairs) <- subtypes
rownames(examplePairs) <- exampleTuples

TF1 <- "MYB"
TF2 <- "MSX2"
subtype <- 4
ind <- ((capacities.log2.list[[subtype]][TF1,])>-2 & (capacities.log2.list[[subtype]][TF2,])>-2 &
          (capacities.log2.list[[subtype]][TF1,])<2 & (capacities.log2.list[[subtype]][TF2,])<2 &
          (capacities.log2.list[[subtype]][TF1,])!=0 & (capacities.log2.list[[subtype]][TF2,])!=0)

cor(capacities.log2.list[[subtype]][TF1,ind],
    capacities.log2.list[[subtype]][TF2,ind])
plot(capacities.log2.list[[subtype]][TF1,ind],
     capacities.log2.list[[subtype]][TF2,ind],pch=20,xlim=c(-5,5),ylim=c(-5,5))

library(spatstat)
subtypeColors <- c("brown3","cornflowerblue","mediumseagreen","goldenrod1")
names(subtypeColors) <- subtypes
#pdf(file = paste0("analyses/0_0.1/",TFoI,"_",tgt,"_",subtype,".pdf"),width = 10,height = 4)
par(mar = c(1, 1, 1, 1))
pppo=ppp(x=capacities.log2.list[[subtype]][TF1,ind],
         y=capacities.log2.list[[subtype]][TF2,ind],
         window = owin(c(-2,2),c(-2,2)))
den=density(pppo,kernel="gaussian",edge=T,diggle=T,adjust=0.4)
plot(den,main='TF vs target',col=colorRampPalette(c("white",subtypeColors[subtype]))(10),xlim=c(-2,2),
     ylim=c(-2,2),xlab=TF1,ylab=TF2)
points(capacities.log2.list[[subtype]][TF1,ind],
       capacities.log2.list[[subtype]][TF2,ind],
       xlim=c(-2,2),ylim=c(-2,2),pch=20,col="black")

TF1 <- "MYB"
TF2 <- "MSX2"
subtype <- subtypes[3]
for (subtype in subtypes) {
  TFCoreg <- cbind.data.frame(capacities.log2.list[[subtype]][TF1,],
                              capacities.log2.list[[subtype]][TF2,])
  ind <- ((capacities.log2.list[[subtype]][TF1,])!=0 & 
            (capacities.log2.list[[subtype]][TF2,])!=0)
  TFCoreg <- TFCoreg[ind,]
  fn <- paste0("analyses/0_0.1/",TF1,"_",TF2,"_",subtype,"_coreg.txt")
  write.table(TFCoreg,fn,sep="\t",
              col.names=F,row.names=F,quote=F)
}

# Get RANSAC inliers to plot trend line.
getRANSACTFCoreg <- function(TF1,TF2,subtype) {
  fn <- paste0("analyses/0_0.1/",TF1,"_",TF2,"_",subtype,"_coreg_RANSAC.txt")
  temp <- data.matrix(read.delim(fn,header=F,check.names=F))
  data.frame(x=temp[,1],y=temp[,2])
}
RANSAC_TFCoreg <- lapply(subtypes,getRANSACTFCoreg,TF1=TF1,TF2=TF2)
RANSAC_regression <- lapply(RANSAC_TFCoreg,function(df){
  lm(y~x,data=df)
})
names(RANSAC_TFCoreg) <- subtypes
names(RANSAC_regression) <- subtypes

library(spatstat)

subtype <- subtypes[4]

for (subtype in subtypes) {
  TFCoreg <- cbind.data.frame(capacities.log2.list[[subtype]][TF1,],
                              capacities.log2.list[[subtype]][TF2,])
  ind <- ((capacities.log2.list[[subtype]][TF1,])!=0 & 
            (capacities.log2.list[[subtype]][TF2,])!=0)
  TFCoreg <- TFCoreg[ind,]
  colnames(TFCoreg) <- c(TF1,TF2)
  subtypePalettes <- c('Reds','Blues','Greens','Oranges')
  crp <- colorRampPalette(c("white",subtypeColors[subtype]))
  names(subtypePalettes) <- subtypes
  xlower <- -4
  xupper <- 4
  ylower <- -4
  yupper <- 4
  #xlower <- min(as.numeric(capacities.log2.list[[subtype]][TF1,]))
  #xupper <- max(as.numeric(expr[TFoI,]))
  #ylower <- min((as.numeric(expr[tgt,])))
  #yupper <- max((as.numeric(expr[tgt,])))
  # data.df <- data.frame(x=as.numeric(availTCGAExprs.bySubtype[[subtype]][TF,]),
  #                       y=log2(as.numeric(annot.df[annot.df[,"label"]==subtype,"survival"])))
  ggplot(TFCoreg,aes(x=TFCoreg[,TF1],y=TFCoreg[,TF2])) + stat_density_2d(aes(fill = ..density..),
                                                                         geom="raster",
                                                                         contour=F) +
    scale_fill_gradient(low = "white",high=subtypeColors[subtype]) + 
    scale_x_continuous(limits=c(xlower-1,xupper+1),breaks=seq(floor(xlower),ceiling(xupper),2)) +
    scale_y_continuous(limits=c(ylower-1,yupper+1),breaks=seq(floor(ylower),ceiling(yupper),2)) +
    theme(legend.position = 'none',panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),panel.background = element_blank(),
          panel.border = element_rect(fill = NA,size = 2)) + geom_point(color='black',size=1) + 
    geom_abline(intercept = RANSAC_regression[[subtype]]$coefficients[1],
                slope = RANSAC_regression[[subtype]]$coefficients[2],
                linetype = 'dashed',
                size = 1.25)
  
}



# geneRegs.TFs.coreg <- lapply(subtypes,function(x){
#   fn <- paste0("analyses/",x,"_TF_TF_cormat_lambda_0_0_rectGate_100918.txt")
#   data.matrix(read.delim(fn,header=T,row.names=1,check.names=F))
# })
# 
# TFCor.tuples.list <- lapply(geneRegs.TFs.coreg,function(x){
#   TFCor.tuples <- rep(0,sum(x!=0))
#   nonZeroIdx <- which(x!=0)
#   nr <- nrow(x)
#   tfs <- rownames(x)
#   count <- 1
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
#       names(TFCor.tuples)[count] <- tup
#       TFCor.tuples[tup] <- x[rowInd,colInd]
#     }
#     count <- count + 1
#     if (count %% 1000 == 0) {
#       cat(".")
#     }
#   }
#   cat("\n")
#   TFCor.tuples
# })
# TFCor.tuples.list <- lapply(TFCor.tuples.list,function(x){
#   x[x!=0]
# })
# 
# library(igraph)
# protNetwork.full <- read.delim("data/BIOGRID-ORGANISM-Homo_sapiens-3.4.158.tab2.txt",header=T,check.names=F)
# protNetwork.full <- cbind(as.character(protNetwork.full$`Official Symbol Interactor A`),
#                           as.character(protNetwork.full$`Official Symbol Interactor B`))
# protNetwork.full <- unique(protNetwork.full)
# protNetwork.graph <- graph.edgelist(protNetwork.full,directed=T)
# TFsInProtNet <- intersect(as.character(names(V(protNetwork.graph))),allTFs)
# protNetwork.graph.TF <- induced.subgraph(protNetwork.graph,vids = TFsInProtNet)
# allTFInteractions <- get.edgelist(protNetwork.graph.TF)
# allTFInteractions.tuples <- apply(allTFInteractions,1,function(x){
#   paste0(x[1],"_",x[2])
# })
# percent.seq <- seq(0,0.5,0.001)
# commonTFCor.tuples <- Reduce("intersect",lapply(TFCor.tuples.list,names))
# commonTFCor.tuples.cor <- lapply(TFCor.tuples.list,function(x){
#   x[commonTFCor.tuples]
# })
# commonTFCor.tuples.cor <- data.matrix(do.call(cbind.data.frame,commonTFCor.tuples.cor))
# colnames(commonTFCor.tuples.cor) <- subtypes
# write.table(commonTFCor.tuples.cor,"analyses/commonTFCor.tuples.cor.100918.txt",
#             sep="\t",quote=F)
# commonTFCor.tuples.cor <- data.matrix(read.delim("analyses/commonTFCor.tuples.cor.100918.txt",
#                                      header=T,row.names=1,check.names=F))
# commonTFCor.tuples.avgCor <- rowMeans(commonTFCor.tuples.cor)
# commonTFCor.tuples.avgCor <- rev(sort(commonTFCor.tuples.avgCor))
# TFCor.tuples.hits.list <- sapply(percent.seq,function(p){
#   hits <- intersect(names(commonTFCor.tuples.avgCor[1:round(p*length(commonTFCor.tuples.avgCor))]),
#                     allTFInteractions.tuples)
#   length(hits)/round(p*length(commonTFCor.tuples.avgCor))
# })
# plot(percent.seq,TFCor.tuples.hits.list,type="l",lwd=3,
#      xlab="percent top correlated TF pairs",
#      ylab="proportion of TF pairs with known interactions")
# p = 0.008
# hits <- intersect(names(commonTFCor.tuples.avgCor[1:round(p*length(commonTFCor.tuples.avgCor))]),
#                   allTFInteractions.tuples)
# x = matrix(c(length(hits),
#              round(p*length(commonTFCor.tuples.avgCor)),
#              length(intersect(names(commonTFCor.tuples.avgCor),allTFInteractions.tuples)),
#              length(commonTFCor.tuples.avgCor)),2,2)
# fisher.test(x) # p value is 3.618E-5

# Now look for correlation signatures and compare with F value signatures

log2F.threshold <- 2

lambda <- "0.1_0"
allTFInteractions.tuples <- read.delim("data/knownTFInteractionTuples.BIOGRID.txt",
                                       header=F,check.names = F)
allTFInteractions.tuples <- as.character(allTFInteractions.tuples[,1])
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

geneRegs.TFs.coreg <- geneTFCor.list

# TFoI.topIntxn.bySubtype <- lapply(subtypes,function(s){
#   corVec <- geneRegs.TFs.coreg[[s]][TFoI,]
#   topCor <- corVec[corVec>0.7]
# })
# TFoI.allInteractors <- Reduce("union",lapply(TFoI.topIntxn.bySubtype,names))
# TFoI.topIntxn.bySubtype <- lapply(subtypes,function(s){
#   geneRegs.TFs.coreg[[s]][TFoI,TFoI.allInteractors]+
#     geneRegs.TFs.coreg[[s]][TFoI.allInteractors,TFoI]
# })
# TFoI.corMat <- do.call(rbind.data.frame,TFoI.topIntxn.bySubtype)
# rownames(TFoI.corMat) <- subtypes
# colnames(TFoI.corMat) <- TFoI.allInteractors
# write.table(TFoI.corMat,paste0("analyses/0_0.1/",TFoI,"topCor_allSubtypes.txt"),
#             sep="\t",quote=F)
# myheatmapfun <- colorRampPalette((c(rgb(15,61,111,maxColorValue=255),
#                                     rgb(125,182,212,maxColorValue=255),
#                                     rgb(255,255,255,maxColorValue=255),
#                                     rgb(234,145,118,maxColorValue=255),
#                                     rgb(102,3,32,maxColorValue=255))))
# heatmap.data <- data.matrix(TFoI.corMat)
# dev.off()
# library(pheatmap)
# pheatmap(heatmap.data,cluster_rows = F,cluster_cols = T,
#          breaks=seq(-1,1,length.out=76),
#          color = myheatmapfun(75))


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

# Update: extract correlation signature the same way as
# for the regulatory signature.

cormat.unshuffled <- geneTFCor.list
nshuffle <- 30
lambda <- '0.1_0'
cormat.shuffled <- lapply(1:nshuffle,function(i){
  cat("\n")
  temp <- lapply(subtypes,function(subtype){
    cat(".")
    fn <- paste0('data/TFCor/shuffleLabels/',subtype,'_cormat_',lambda,'_randomShuffle_',
                 as.character(i),'.txt')
    geneTF.coreg.mat <- data.matrix(read.delim(fn,header=F,check.names=F))
    if (nrow(geneTF.coreg.mat)<length(allTFs)) {
      geneTF.coreg.mat <- rbind(geneTF.coreg.mat,rep(0,length(allTFs)-1))
      geneTF.coreg.mat <- cbind(geneTF.coreg.mat,rep(0,length(allTFs)))
      allTFNames <- rownames(capacities.list[[1]])
      missingIdx <- which(allTFNames == "ZNF8") 
      rownames(geneTF.coreg.mat) <- colnames(geneTF.coreg.mat) <- c(allTFNames[-missingIdx],"ZNF8")
      geneTF.coreg.mat <- geneTF.coreg.mat[allTFNames,allTFNames]
    }
    geneTF.coreg.mat[is.nan(geneTF.coreg.mat)] <- 0
    rownames(geneTF.coreg.mat) <- colnames(geneTF.coreg.mat) <- rownames(capacities.list[[1]])
    geneTF.coreg.mat[lower.tri(geneTF.coreg.mat)] <- t(geneTF.coreg.mat)[lower.tri(geneTF.coreg.mat)]
    geneTF.coreg.mat
    #data.matrix(read.delim(fn,header=F,check.names=F))
  })
  names(temp) <- subtypes
  temp
})

nsd <- 5
niter <- 30
source('src/obtainCorSignatures.R')
res <- obtainCorSignatures(nsd,cormat.unshuffled,shuffled=F,sigTFs=NULL,adjPThres=1E-10)
subtypeSpecificCorMats <- res$subtypeSpecificCorMats

names(subtypeSpecificCorMats) <- subtypes
geneCoregs.all <- subtypeSpecificCorMats
geneCoregs.mean <- Reduce("+",geneCoregs.all)
geneCoregs.mean <- geneCoregs.mean / 4
geneCoregs.var <- lapply((geneCoregs.all),function(x){(x-geneCoregs.mean)^2})
geneCoregs.var <- Reduce("+",geneCoregs.var) / 4
geneCoregs.sd <- geneCoregs.var^0.5
geneCoregs.TF <- lapply(geneCoregs.all,function(x){
  x[allTFs,]
})
geneCoregs.TF.sd <- geneCoregs.sd[allTFs,]
meanTFCV <- (geneCoregs.TF.sd)/abs(geneCoregs.mean[allTFs,])
meanTFCV <- t(apply(meanTFCV,1,function(x){
  y <- x
  y[is.na(y)] <- 0
  y
}))
geneCoregs.TF.mean <- geneCoregs.mean[allTFs,]
leaveOneSD <- lapply(c(1:length(subtypes)),function(x){
  currList <- geneCoregs.TF[-x]
  LOM <- Reduce("+",currList)
  LOM <- LOM / (length(subtypes)-1)
  LOV <- lapply(currList,function(y){(y-LOM)^2})
  LOV <- Reduce("+",LOV) / (length(subtypes)-1)
  LOSD <- LOV^0.5
  LOSD
})
overallSD <- geneCoregs.TF.sd
# Use the amount of reduction in variance as measurement for strength of signature
# Also set a cutoff for log2 change in F values
# Seems like the oddity in the number of Neural signature TF-gene pairs decreases
# as the threshold becomes more stringent. Points to potential bias in data or poor
# regression performance.
#log2FDiffThreshold <- 5
corDiffThreshold <- 0.2
sigCorSubtype <- matrix(0,nrow(geneCoregs.TF[[1]]),ncol(geneCoregs.TF[[1]]))
SDReduction.mat <- matrix(0,nrow(geneCoregs.TF[[1]]),ncol(geneCoregs.TF[[1]]))
for (i in c(1:nrow(geneCoregs.TF[[1]]))) {
  for (j in c(1:ncol(geneCoregs.TF[[1]]))) {
    Fval <- sapply(leaveOneSD,function(x){
      x[i,j]
    })
    if (overallSD[i,j] > 0) {
      Fvec <- sapply(geneCoregs.TF,function(x){
        x[i,j]
      })
      ss <- which(Fval==min(Fval))[1]
      corDiff <- abs(Fvec[ss]-mean(Fvec[-ss]))
      if (corDiff > corDiffThreshold) {
        SDReduction.mat[i,j] <- overallSD[i,j]-Fval[ss]
        sigCorSubtype[i,j] <- ss
      }
    }
  }
}
rownames(sigCorSubtype) <- rownames(geneCoregs.TF[[1]])
colnames(sigCorSubtype) <- colnames(geneCoregs.TF[[1]])

outDegrees <- sapply(unionNet.byTF,nrow)
outDegrees <- outDegrees[allTFs]
sigOutDegrees <- (outDegrees >= 50)
percentSigCutoff <- 0.000
sigCorSubtype.filtered <- sigCorSubtype[sigOutDegrees,]
domSigCorSubtype.filtered <- sapply(rownames(sigCorSubtype.filtered),function(x){
  temp <- table(sigCorSubtype.filtered[x,])
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
# domSigCorSubtype.filtered.multipleTgts <- sapply(names(domSigCorSubtype.filtered[domSigCorSubtype.filtered!="0"]),function(x){
#   tgts <- tgtGenes <- which(sigCorSubtype[x,]==as.numeric(domSigCorSubtype.filtered[x]))
#   length(tgts) > 5
# })
GSE57872.all <- read.delim("data/GSE57872_GBM_data_matrix.txt",header=T,row.names=1)
GSE57872.all <- data.matrix(GSE57872.all)
GSE57872.topExprsGenes <- rev(order(rowMeans(GSE57872.all)))[1:3000]

# write.table(domSigSubtype.filtered[domSigSubtype.filtered!="0"][domSigSubtype.filtered.multipleTgts],
#             "domSigSubtype.filtered_percentCutoff0.005.multipleTarget.txt",sep="\t",quote=F)
sigCorTFs <- domSigCorSubtype.filtered[intersect(rownames(GSE57872.all)[GSE57872.topExprsGenes],
                                           names(domSigCorSubtype.filtered[domSigCorSubtype.filtered!="0"]))]
sigTFs[intersect(names(sigCorTFs),names(sigTFs))]
sigCorTFs[intersect(names(sigCorTFs),names(sigTFs))]
sum(sigTFs[intersect(names(sigCorTFs),names(sigTFs))]==
      sigCorTFs[intersect(names(sigCorTFs),names(sigTFs))])
consistentTFs <- names(which(sigTFs[intersect(names(sigCorTFs),names(sigTFs))]==
      sigCorTFs[intersect(names(sigCorTFs),names(sigTFs))]))
sigTFs[consistentTFs]
sigCorTFs[consistentTFs]

# overallSD <- geneRegs.coreg.sd
# # Generate correlation signatures for Prism
# coregCutoffs <- c(0.1,0.15,0.2,0.25,0.3,0.35,0.4)
# res <- c()
# coregCutoff <- c(0.2)
# for (coregCutoff in coregCutoff) {
#   sigTFCorSubtype <- matrix(0,nrow(geneRegs.TFs.coreg[[1]]),ncol(geneRegs.TFs.coreg[[1]]))
#   for (i in c(1:nrow(geneRegs.TFs.coreg[[1]]))) {
#     for (j in c(1:ncol(geneRegs.TFs.coreg[[1]]))) {
#       Fval <- sapply(leaveOneSD,function(x){
#         x[i,j]
#       })
#       coregs.allSubtypes <- sapply(c(1:length(subtypes)),function(k){
#         geneRegs.TFs.coreg[[k]][i,j]
#       })
#       if (overallSD[i,j] > 0 & max(abs(coregs.allSubtypes)) > coregCutoff) {
#         sigTFCorSubtype[i,j] <- which(Fval==min(Fval))[1]
#       }
#     }
#   }
#   write.table(sigTFCorSubtype,"analyses/0_0.1/sigTFCorSubtype.mat.txt",sep='\t',quote=F)
#   rownames(sigTFCorSubtype) <- rownames(geneRegs.TFs.coreg[[1]])
#   colnames(sigTFCorSubtype) <- colnames(geneRegs.TFs.coreg[[1]])
#   table(sigTFCorSubtype[names(domSigSubtype),])
#   sigTFCorSubtype.sigTFs <- sigTFCorSubtype[names(domSigSubtype),]
#   sigTFCorSubtype.sigTFs.ordered <- t(apply(sigTFCorSubtype.sigTFs,1,function(x){
#     prefix <- rep(0,sum(x==0))
#     suffix <- sort(table(x[x!=0]))
#     suffix <- unlist(lapply(1:length(suffix),function(i){
#       rep(names(suffix)[i],suffix[i])
#     }))
#     as.numeric(c(prefix,suffix))
#   }))
#   sigTFCorSubtype.sigTFs.ordered <- sigTFCorSubtype.sigTFs.ordered[rownames(signatureMat),]
#   sigTFCorSubtype.sigTFs.ordered <- sigTFCorSubtype.sigTFs.ordered[,colSums(sigTFCorSubtype.sigTFs.ordered)!=0]
#   cat(".")
#   res <- c(res,sum(sigTFCorSubtype.sigTFs.ordered[,ncol(sigTFCorSubtype.sigTFs.ordered)]==
#         signatureMat[,ncol(signatureMat)]))
# }

sigCorTFSubtype <- sigCorSubtype[names(sigCorTFs),]
cor.sigSubtype.participation <- apply(sigCorTFSubtype,1,function(x){
  temp <- x[x!='0']
  tally <- table(x)
  for (i in c('1','2','3','4')) {
    if (!(i %in% names(tally))) {
      tally[i] <- 0
    }
  }
  tally <- tally[c('1','2','3','4')]
  tally/sum(tally)
})
rownames(cor.sigSubtype.participation) <- subtypes
colnames(cor.sigSubtype.participation) <- names(sigCorTFs)

myheatmapfun <- colorRampPalette((c(rgb(255,255,255,maxColorValue=255),
                                    rgb(234,145,118,maxColorValue=255),
                                    rgb(102,3,32,maxColorValue=255))))
library(pheatmap)
library(viridis)
xxx <- pheatmap(cor.sigSubtype.participation,breaks=seq(0.2,0.8,length.out=76),
                clustering_method = "ward.D",cluster_rows = F,
                clustering_distance_cols = "correlation",
                color = magma(75),dendrogram=F,fontsize=12)
subtype <- 'Classical'
for (subtype in subtypes) {
  cor.sigSubtype.participation.subtype <- cor.sigSubtype.participation[,
                                                               names(sigCorTFs)
                                                               [sigCorTFs==as.character(match(subtype,
                                                                                           subtypes))]]
  pheatmap(cor.sigSubtype.participation.subtype,breaks=seq(0.2,0.8,length.out=76),
           clustering_method = "ward.D",cluster_rows = F,
           clustering_distance_cols = "correlation",
           color = magma(75),dendrogram=F,fontsize=12,
           treeheight_row = 0, treeheight_col = 0,
           cellwidth = 20,cellheight = 20)
}

# Plot examples of rewired interactions
TFoI <- "MYC"
TFoI.sigCorSubtype <- 4
TFoI.sigIntxn.bySubtype <- lapply(subtypes,function(s){
  sigIntxn <- colnames(sigTFCorSubtype)[sigTFCorSubtype[TFoI,]==TFoI.sigCorSubtype]
  corVec <- geneRegs.TFs.coreg[[s]][TFoI,sigIntxn]
})
TFoI.corMat <- do.call(rbind.data.frame,TFoI.sigIntxn.bySubtype)
rownames(TFoI.corMat) <- subtypes
colnames(TFoI.corMat) <- colnames(sigTFCorSubtype)[sigTFCorSubtype[TFoI,]==
                                                     TFoI.sigCorSubtype]
TFoI.corMat <- data.matrix(TFoI.corMat)
TFoI.corMat <- TFoI.corMat[,apply(TFoI.corMat,2,function(x){max(abs(x))>0.8})]
write.table(TFoI.corMat,paste0("analyses/0_0.1/",TFoI,"sigCor_allSubtypes.txt"),
            sep="\t",quote=F)
myheatmapfun <- colorRampPalette((c(rgb(15,61,111,maxColorValue=255),
                                    rgb(125,182,212,maxColorValue=255),
                                    rgb(255,255,255,maxColorValue=255),
                                    rgb(234,145,118,maxColorValue=255),
                                    rgb(102,3,32,maxColorValue=255))))
heatmap.data <- data.matrix(TFoI.corMat)
dev.off()
library(pheatmap)
pheatmap(heatmap.data,cluster_rows = F,cluster_cols = T,
         breaks=seq(-1,1,length.out=76),
         color = myheatmapfun(75))

# Find examples to show

signatureMat.data.cor <- lapply(1:nrow(sigTFCorSubtype.sigTFs.ordered),function(x){
  temp <- table(sigTFCorSubtype.sigTFs.ordered[x,])
  temp <- temp[names(temp)!="0"]
  for (i in as.character(1:length(subtypes))) {
    if (!(i %in% names(temp))) {
      temp[i] <- 0
    }
  }
  temp[as.character(1:length(subtypes))]
})
names(signatureMat.data.cor) <- rownames(sigTFCorSubtype.sigTFs.ordered)
signatureMat.data.cor <- do.call(rbind.data.frame,signatureMat.data.cor)
colnames(signatureMat.data.cor) <- subtypes
signatureMat.data.cor <- t(apply(data.matrix(signatureMat.data.cor),1,function(x){
  x/sum(x)
}))
rownames(signatureMat.data.cor) <- rownames(sigTFCorSubtype.sigTFs.ordered)
sigTFCor.domSubtypes <- apply(sigTFCorSubtype.sigTFs.ordered,1,function(x){
  temp <- table(x)
  temp <- temp[names(temp)!="0"]
  names(temp)[temp==max(temp)][1]
})
signatureMat.data.cor <- signatureMat.data.cor[rownames(signatureMat.data),]
write.table(signatureMat.data.cor,"analyses/0_0.1/cor_signatureMat_Prism.txt",sep="\t",
            quote=F)
sigCorMatList <- list()
for (i in 1:length(subtypes)) {
  signatureMat.data.cor.subtype <- signatureMat.data.cor[names(sigTFCor.domSubtypes)[sigTFCor.domSubtypes==as.character(i)],]
  signatureMat.data.cor.subtype <- signatureMat.data.cor.subtype[order(apply(signatureMat.data.cor.subtype,1,max)),]
  write.table(signatureMat.data.cor.subtype,paste0("analyses/0_0.1/cor_signatureMat_",subtypes[i],"_Prism.txt"),sep="\t",
              quote=F)
  sigCorMatList[[i]] <- signatureMat.data.cor.subtype
}
sigCorMatList <- do.call(rbind.data.frame,sigCorMatList)
# sigCorMatList <- sigCorMatList[rownames(sigMatList),]
write.table(sigCorMatList,"analyses/0_0.1/sigCorMats_Prism_orderedBySubtype.txt",quote=F,sep="\t")


domSig.F <- apply(signatureMat.data,1,function(x){
  which(x==max(x))[1]
})
domSig.cor <- apply(signatureMat.data.cor,1,function(x){
  which(x==max(x))[1]
})
consistentSig <- names(which(domSig.cor==domSig.F))
write.table(signatureMat.data[consistentSig,],"analyses/0_0.1/F_signatureMatConstt_Prism.txt",sep="\t",
            quote=F)
write.table(signatureMat.data.cor[consistentSig,],"analyses/0_0.1/cor_signatureMatConstt_Prism.txt",sep="\t",
            quote=F)

write.table(signatureMat,"analyses/0_0.1/F_signatureMat.txt",sep="\t",
            quote=F)
write.table(sigTFCorSubtype.sigTFs.ordered,"analyses/0_0.1/cor_signatureMat.txt",sep="\t",
            quote=F)

# domSigSubtype <- apply(heatmap.data,1,function(x){
#   temp <- table(x)
#   temp <- temp[names(temp) != "0"]
#   names(temp)[which(temp==max(temp))][1]
# })
# sigTFColors <- rep("red",length(domSigSubtype))
# sigTFColors[domSigSubtype=="2"] <- "blue"
# sigTFColors[domSigSubtype=="3"] <- "green"
# sigTFColors[domSigSubtype=="4"] <- "orange"
# names(sigTFColors) <- names(domSigSubtype)
# 

dev.off()
hmRes <- heatmap.2(sigTFCorSubtype.sigTFs.ordered,
                   Colv=NA,
                   Rowv=NA,
                   scale="none",
                   density.info="none", trace="none",
                   col=c("white","red","navyblue","darkgreen","orange"),
                   reorderfun = function(d,w) { d },
                   cexCol = 0.3,
                   srtCol = 45,
                   cexRow = 1,
                   offsetCol = -0.5,
                   rowsep = c(1:nrow(sigTFCorSubtype.sigTFs.ordered)),
                   sepcolor = "white",
                   #RowSideColors = rscVec,
                   dendrogram="none",
                   key.title=NA,
                   key.xlab=NA,
                   lhei=c(1,8),
                   #margins=c(10,10),
                   labCol="")

# Prepare circos data
geneRegs.TFs.coreg <- lapply(geneRegs.TFs.coreg,function(x){
  constid <- which(rownames(x)=="const")
  if (length(constid)) {
    x[-constid,]
  } else {
    x
  }
})
# geneRegs.TFs.coreg.topTFs <- lapply(geneRegs.TFs.coreg,function(x){
#   rowMeanAbs <- apply(x,1,function(y){mean(abs(y))})
#   x[rowMeanAbs > 0.05,rowMeanAbs > 0.05]
# })
# topCorTFs <- Reduce("intersect",lapply(geneRegs.TFs.coreg.topTFs,rownames))
# 
#sigDiffTFs <- rownames(signatureMat)
sigDiffTFs <- names(sigCorTFs)[order(sigCorTFs)]
# nonSigDiffTFs <- intersect(setdiff(rownames(geneRegs.TFs.coreg[[1]]),sigDiffTFs),topCorTFs)
interactorTFs <- lapply(geneRegs.TFs.coreg,function(x){
  submat <- x[sigDiffTFs,]
  colnames(submat)[which(apply(submat,2,function(y){mean(abs(y))>0.5}))]
})
interactorTFs <- Reduce("union",interactorTFs)
allChr <- c(sigDiffTFs,setdiff(interactorTFs,sigDiffTFs))
geneRegs.TFs.coreg.topTFs <- lapply(geneRegs.TFs.coreg,function(x){
  x[sigDiffTFs,allChr]
})
chr.list <- lapply(1:length(allChr),function(x){
  c("chr","-",paste0("hs",as.character(x)),allChr[x],
    0,10000,allChr[x])
})
chr.df <- do.call(rbind.data.frame,chr.list)
colnames(chr.df) <- ""
write.table(chr.df,"analyses/0_0.1/circos/karyotype.TFs.txt",sep=" ",col.names=F,row.names=F,quote=F)
# Prepare links
chrNameLookUp <- as.character(chr.df[,3])
names(chrNameLookUp) <- as.character(chr.df[,7])
chr.color <- rep("black",nrow(chr.df))
names(chr.color) <- names(chrNameLookUp)
chr.color[names(sigCorTFs)[which(sigCorTFs=="1")]] <- "red"
chr.color[names(sigCorTFs)[which(sigCorTFs=="2")]] <- "blue"
chr.color[names(sigCorTFs)[which(sigCorTFs=="3")]] <- "green"
chr.color[names(sigCorTFs)[which(sigCorTFs=="4")]] <- "orange"
chr.list <- lapply(1:length(allChr),function(x){
  c("chr","-",paste0("hs",as.character(x)),allChr[x],
    0,10000,chr.color[x])
})
chr.col.df <- do.call(rbind.data.frame,chr.list)
colnames(chr.col.df) <- ""
# The original 'analyses/circos/..' file was overwritten by 0.1 data!!!
write.table(chr.col.df,"analyses/0_0.1/circos/karyotype.TFs.colored.txt",sep=" ",col.names=F,row.names=F,quote=F)

links.list <- lapply(geneRegs.TFs.coreg.topTFs,function(x){
  linkList <- list()
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      if (abs(x[i,j])>0.9) {
        nodeA <- chrNameLookUp[rownames(x)[i]]
        nodeB <- chrNameLookUp[colnames(x)[j]]
        linkList <- c(linkList,list(c(nodeA,0,10000,nodeB,0,10000)))
      }
    }
  }
  linkList <- do.call(rbind.data.frame,linkList)
  colnames(linkList) <- ""
  linkList
})
#(rownames(x)[i] %in% sigDiffTFs || rownames(x)[j] %in% interactorTFs)
lapply(c(1:length(subtypes)),function(x){
  fn <- paste0("analyses/0_0.1/circos/TF_coreg_links_sigDiffCorTFs_",subtypes[x],".txt")
  write.table(links.list[[x]],fn,sep=" ",row.names=F,col.names=F,quote=F)
})

# Need to zoom in onto specific TFs to show 'local rewiring'.
links.list <- lapply(geneRegs.TFs.coreg.topTFs,function(x){
  linkList <- list()
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      if (abs(x[i,j])>0.9 & rownames(x)[i]=="MXI1") {
        nodeA <- chrNameLookUp[rownames(x)[i]]
        nodeB <- chrNameLookUp[colnames(x)[j]]
        linkList <- c(linkList,list(c(nodeA,0,10000,nodeB,0,10000)))
      }
    }
  }
  linkList <- do.call(rbind.data.frame,linkList)
  colnames(linkList) <- ""
  linkList
})
#(rownames(x)[i] %in% sigDiffTFs || rownames(x)[j] %in% interactorTFs)
lapply(c(1:length(subtypes)),function(x){
  fn <- paste0("analyses/0_0.1/circos/MXI1_TF_coreg_links_sigDiffTFs_",subtypes[x],".txt")
  write.table(links.list[[x]],fn,sep=" ",row.names=F,col.names=F,quote=F)
})

# chr.color <- rep("black",nrow(chr.df))
# names(chr.color) <- names(chrNameLookUp)
# chr.color[names(sigTFs)[which(sigTFs=="1")]] <- "red"
# chr.color[names(sigTFs)[which(sigTFs=="2")]] <- "blue"
# chr.color[names(sigTFs)[which(sigTFs=="3")]] <- "green"
# chr.color[names(sigTFs)[which(sigTFs=="4")]] <- "orange"
# 
# chr.color.df <- data.frame(rep("chr",length(chrNameLookUp)),
#                            rep("-",length(chrNameLookUp)),
#                            (chrNameLookUp),
#                            chr.color)
# write.table(chr.color.df,"analyses/circos/chr.color.txt",sep=" ",
#             row.names=F,col.names=F,quote=F)

# Now see how many differentially regulated target genes are also expression signatures -
# this shows whether differential expression landscape can be explained by differential
# TF regulatory activity.

subtypeCentroids <- read.delim("data/subtype_centroids.txt",header=T,row.names=1)

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
colnames(leaveOneSD) <- subtypes
#leaveOneSD.subtypeSig <- leaveOneSD[rownames(subtypeCentroids),]
leaveOneSD.subtypeSig <- leaveOneSD[colnames(capacities.log2.list[[1]]),]
exprsSigSubtype <- apply(leaveOneSD.subtypeSig,1,function(x){
  id <- which(x==min(x))[1]
  subtypes[id]
})
leaveOneSDReduction.subtypeSig <- apply(leaveOneSD.subtypeSig,1,function(x){
  id <- which(x==min(x))[1]
  mean(x[-id]) - x[id]
})
topExprsSigGenes <- names(leaveOneSDReduction.subtypeSig)[rev(order(leaveOneSDReduction.subtypeSig))[1:2000]]
exprsSigSubtype <- exprsSigSubtype[topExprsSigGenes]
exprsSigSubtype.list <- lapply(subtypes,function(x){
  intersect(names(exprsSigSubtype)[exprsSigSubtype==x],colnames(capacities.log2.list[[1]]))
})
names(exprsSigSubtype.list) <- subtypes
exprsSigSubtype.domSubtype <- lapply(exprsSigSubtype.list,function(x){
  domSubtype <- sapply(x,function(tgt){
    counts <- table(sigSubtype.filtered[,tgt])
    if (length(counts)>1) {
      counts <- counts[names(counts)!="0"]
      names(counts)[counts==max(counts)][1] 
    } else {
      if (names(counts)!="0") {
        names(counts)
      } else {
        "0"
      }
    }
  })
})
cstExprsRegTgts <- lapply(1:length(exprsSigSubtype.domSubtype),function(x){
  names(exprsSigSubtype.domSubtype[[x]])[which(exprsSigSubtype.domSubtype[[x]]==as.character(x))]
})
names(cstExprsRegTgts) <- subtypes
cstExprsRegTgts.TFs <- lapply(c(1:length(cstExprsRegTgts)),function(x){
  domSigTFs <- lapply(cstExprsRegTgts[[x]],function(tgt){
    rownames(sigSubtype.filtered)[sigSubtype.filtered[,tgt]==as.character(x)]
  })
  unique(unlist(domSigTFs))
})
names(cstExprsRegTgts.TFs) <- subtypes

# Update on rewired target identification: if a gene is part of the expression
# signature, as long as its F profile shows any signature TF-gene pair that
# is consistent with the subtype label, it is considered rewired.

exprsSigSubtype.rewiredTgts <- lapply(1:length(exprsSigSubtype.list),function(x){
  subtypeLabel <- as.character(x)
  isExprsRegSig <- sapply(exprsSigSubtype.list[[x]],function(tgt){
    subtypeLabel %in% as.character(sigSubtype.filtered[,tgt])
  })
  names(exprsSigSubtype)[(which(isExprsRegSig))]
})
exprsSigSubtype.rewiredTgtTFs <- lapply(c(1:length(cstExprsRegTgts)),function(x){
  domSigTFs <- lapply(exprsSigSubtype.rewiredTgts[[x]],function(tgt){
    if (tgt %in% colnames(sigSubtype.filtered)) {
      rownames(sigSubtype.filtered)[sigSubtype.filtered[,tgt]==as.character(x)]
    } else {
      return()
    }
  })
  unique(unlist(domSigTFs))
})
names(exprsSigSubtype.rewiredTgtTFs) <- subtypes
# Few or none of the TFs of rewired expression signature genes are expression signatures
# themselves, meaning that the differential regulation of these targets are not explained
# by differential expression of TFs.
exprsSigSubtype.rewiredTgtTFs.inExprsSig <- lapply(names(exprsSigSubtype.rewiredTgtTFs),function(x){
  intersect(exprsSigSubtype.rewiredTgtTFs[[x]],exprsSigSubtype.list[[x]])
})
# These are the percentages for the above statement.
sapply(exprsSigSubtype.rewiredTgtTFs.inExprsSig,length)/sapply(exprsSigSubtype.rewiredTgtTFs,length)

exprsSigSubtype.rewiredTgtTFs.byTgt <- lapply(c(1:length(cstExprsRegTgts)),function(x){
  rewiredTgtTFs <- lapply(exprsSigSubtype.rewiredTgts[[x]],function(tgt){
    rownames(sigSubtype.filtered)[sigSubtype.filtered[,tgt]==as.character(x)]
  })
  names(rewiredTgtTFs) <- exprsSigSubtype.rewiredTgts[[x]]
  exprsSigTFs <- lapply(names(rewiredTgtTFs),function(n){
    dst <- rewiredTgtTFs[[n]]
    dst.DE <- length(intersect(dst,exprsSigSubtype.list[[x]]))/length(dst)
  })
})
names(exprsSigSubtype.rewiredTgtTFs.byTgt) <- subtypes
boxplot(unlist(exprsSigSubtype.rewiredTgtTFs.byTgt[[1]]),
        unlist(exprsSigSubtype.rewiredTgtTFs.byTgt[[2]]),
        unlist(exprsSigSubtype.rewiredTgtTFs.byTgt[[3]]),
        unlist(exprsSigSubtype.rewiredTgtTFs.byTgt[[4]]))

rewiredTgts.domSigSubtype <- apply(sigSubtype.filtered,2,function(x){
  temp <- table(x)
  if(length(temp)>1) {
    temp <- temp[names(temp)!="0"]
    names(temp)[temp==max(temp)][1]
  } else {
    "0"
  }
})

rewiredTgt.enrichment.list <- lapply(1:length(cstExprsRegTgts),function(x){
  m <- sum(rewiredTgts.domSigSubtype==as.character(x)) # regulatory signature genes
  y <- length(exprsSigSubtype.rewiredTgts[[x]]) # both reg-sig and exprs-sig
  n <- ncol(sigSubtype.filtered) - m # non-regulatory signature genes
  k <- length(exprsSigSubtype.domSubtype[[x]]) # exprs-sig
  hyperG.test <- phyper(y-1,m,n,k,lower.tail = F)
  list(m,y,n,k,hyperG.test)
})

# What about the signature genes that did not undergo txn rewiring?
# (These are genes whose column in the subtype signature matrix did not contain
# any subtype label)
# Did they show diff. expression because of that of their regulators?
# How should this be demonstrated?
# One possibility: take the regulatory vectors for a given nonrewired gene across
# all four subtypes and calculate correlation among them.
# The guess would be that the correlation will be higher in nonrewired genes than in
# the rewired ones.
# In addition, the average correlation between target genes and its regulators
# across subtypes should be higher in nonrewired ones than in rewired ones if the
# above hypothesis is true.

nonRewiredExprsSig.list <- lapply(c(1:length(subtypes)),function(x){
  rewiredExprsSigs <- exprsSigSubtype.rewiredTgts[[x]]
  allExprsSigs <- (exprsSigSubtype.list[[x]])
  setdiff(allExprsSigs,rewiredExprsSigs)
})

nonRewiredExprsSig.FCor.list <- lapply(c(1:length(subtypes)),function(x){
  nonRewiredExprsSigs <- nonRewiredExprsSig.list[[x]]
  averageFCor <- sapply(nonRewiredExprsSigs,function(nres){
    regs <- regulators[[nres]]
    FVecs <- lapply(subtypes,function(subtype){
      capacities.log2.list[[subtype]][regs,nres]
    })
    FVecs.mat <- do.call(cbind.data.frame,FVecs)
    FVecs.cormat <- cor(data.matrix(FVecs.mat))
    avgCor <- mean(FVecs.cormat[upper.tri(FVecs.cormat)])
    avgCor
  })
  averageFCor
})

rewiredExprsSig.FCor.list <- lapply(c(1:length(subtypes)),function(x){
  rewiredExprsSigs <- exprsSigSubtype.rewiredTgts[[x]]
  averageFCor <- sapply(rewiredExprsSigs,function(res){
    regs <- regulators[[res]]
    FVecs <- lapply(subtypes,function(subtype){
      capacities.log2.list[[subtype]][regs,res]
    })
    FVecs.mat <- do.call(cbind.data.frame,FVecs)
    FVecs.cormat <- cor(data.matrix(FVecs.mat))
    avgCor <- mean(FVecs.cormat[upper.tri(FVecs.cormat)])
    avgCor
  })
  averageFCor
})

nonRewiredExprsSig.exprsCorWithRegs.list <- lapply(c(1:length(subtypes)),function(x){
  nonRewiredExprsSigs <- nonRewiredExprsSig.list[[x]]
  averageExprsCorWithRegs <- sapply(nonRewiredExprsSigs,function(nres){
    regs <- regulators[[nres]]
    nres.subtypeMeanExprs <- subtypeMeanExprs[nres,]
    regs.subtypeMeanExprs <- subtypeMeanExprs[regs,,drop=F]
    allCor <- sapply(c(1:nrow(regs.subtypeMeanExprs)),function(idx){
      cor(nres.subtypeMeanExprs,regs.subtypeMeanExprs[idx,])
    })
    avgCor <- mean(abs(allCor[complete.cases(allCor)]))
    avgCor
  })
  averageExprsCorWithRegs
})


nonRewiredExprsSigs.all <- unlist(nonRewiredExprsSig.list)
nonRewiredExprsSig.exprsCorWithRegs.all <- sapply(nonRewiredExprsSigs.all,function(nres){
  regs <- regulators[[nres]]
  nres.subtypeMeanExprs <- subtypeMeanExprs[nres,]
  regs.subtypeMeanExprs <- subtypeMeanExprs[regs,,drop=F]
  allCor <- sapply(c(1:nrow(regs.subtypeMeanExprs)),function(idx){
    cor(nres.subtypeMeanExprs,regs.subtypeMeanExprs[idx,])
  })
  avgCor <- mean(abs(allCor[complete.cases(allCor)]))
  avgCor
})

rewiredExprsSigs.all <- unlist(exprsSigSubtype.rewiredTgts)
rewiredExprsSig.exprsCorWithRegs.all <- sapply(rewiredExprsSigs.all,function(nres){
  regs <- regulators[[nres]]
  nres.subtypeMeanExprs <- subtypeMeanExprs[nres,]
  regs.subtypeMeanExprs <- subtypeMeanExprs[regs,,drop=F]
  allCor <- sapply(c(1:nrow(regs.subtypeMeanExprs)),function(idx){
    cor(nres.subtypeMeanExprs,regs.subtypeMeanExprs[idx,])
  })
  avgCor <- mean(abs(allCor[complete.cases(allCor)]))
  avgCor
})

data.df <- data.frame("Correlation"=c(rewiredExprsSig.exprsCorWithRegs.all,
                                      nonRewiredExprsSig.exprsCorWithRegs.all),
                      "Subgroup"=factor(c(rep("With differential regulation",
                                              length(rewiredExprsSig.exprsCorWithRegs.all)),
                                          rep("Without differential regulation",
                                              length(nonRewiredExprsSig.exprsCorWithRegs.all))),
                                        levels=c("With differential regulation",
                                                 "Without differential regulation")))
colvec <- c("gray","gray")
colvec[1] <- "purple"
# ggplot(data.df, aes(x = Correlation, y = Subgroup,fill=Subgroup)) + geom_density_ridges2() +
#   scale_fill_manual(values=(colvec)) +
#   theme(panel.background = element_blank()) +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) + ggtitle(subtype) + xlim(c(0,1))
ggplot(data.df, aes(x=Correlation, fill=Subgroup)) + geom_density(alpha=.5) +
  scale_fill_manual(values=(colvec)) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + ggtitle(subtype) + xlim(c(0.2,0.8))

rewiredExprsSig.exprsCorWithRegs.list <- lapply(c(1:length(subtypes)),function(x){
  rewiredExprsSigs <- exprsSigSubtype.rewiredTgts[[x]]
  averageExprsCorWithRegs <- sapply(rewiredExprsSigs,function(res){
    regs <- regulators[[res]]
    res.subtypeMeanExprs <- subtypeMeanExprs[res,]
    regs.subtypeMeanExprs <- subtypeMeanExprs[regs,,drop=F]
    allCor <- sapply(c(1:nrow(regs.subtypeMeanExprs)),function(idx){
      cor(res.subtypeMeanExprs,regs.subtypeMeanExprs[idx,])
    })
    avgCor <- mean(abs(allCor[complete.cases(allCor)]))
    avgCor
  })
  averageExprsCorWithRegs
})

for (i in c(1:length(subtypes))) {
  fn <- paste0("analyses/0_0.1/",subtypes[i],"nonRewiredExprsSig.FCor.txt")
  write.table(nonRewiredExprsSig.FCor.list[[i]][complete.cases(nonRewiredExprsSig.FCor.list[[i]])],
              fn,sep="\t",row.names=F,col.names=F,quote=F)
  fn <- paste0("analyses/0_0.1/",subtypes[i],"rewiredExprsSig.FCor.txt")
  write.table(rewiredExprsSig.FCor.list[[i]][complete.cases(rewiredExprsSig.FCor.list[[i]])],
              fn,sep="\t",row.names=F,col.names=F,quote=F)
  fn <- paste0("analyses/0_0.1/",subtypes[i],"nonRewiredExprsSig.exprsCorWithRegs.txt")
  write.table(nonRewiredExprsSig.exprsCorWithRegs.list[[i]][complete.cases(nonRewiredExprsSig.exprsCorWithRegs.list[[i]])],
              fn,sep="\t",row.names=F,col.names=F,quote=F)
  fn <- paste0("analyses/0_0.1/",subtypes[i],"rewiredExprsSig.exprsCorWithRegs.txt")
  write.table(rewiredExprsSig.exprsCorWithRegs.list[[i]][complete.cases(rewiredExprsSig.exprsCorWithRegs.list[[i]])],
              fn,sep="\t",row.names=F,col.names=F,quote=F)
}

library(ggplot2)
library(ggridges)

subtypeColors <- c("brown3","cornflowerblue","mediumseagreen","goldenrod1" )

for (i in 1:length(subtypes)) {
  subtype <- subtypes[i]
  outfn <- paste0("analyses/0_0.1/",subtype,"corWithRegExprs.pdf")
  data.df <- data.frame("Correlation"=c(rewiredExprsSig.exprsCorWithRegs.list[[i]],
                                        nonRewiredExprsSig.exprsCorWithRegs.list[[i]]),
                        "Subgroup"=factor(c(rep("With differential regulation",
                                                            length(rewiredExprsSig.exprsCorWithRegs.list[[i]])),
                                           rep("Without differential regulation",
                                               length(nonRewiredExprsSig.exprsCorWithRegs.list[[i]]))),
                                         levels=c("With differential regulation",
                                                  "Without differential regulation")))
  colvec <- c("gray","gray")
  colvec[1] <- subtypeColors[i]
  ggplot(data.df, aes(x = Correlation, y = Subgroup,fill=Subgroup)) + geom_density_ridges2() +
    scale_fill_manual(values=(colvec)) +
    theme(panel.background = element_blank()) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + ggtitle(subtype) + xlim(c(0,1))
  ggsave(file = outfn, width = 10, height = 4.5)
}

x1 <- ecdf(rewiredExprsSig.exprsCorWithRegs.list[[1]])
x2 <- ecdf(nonRewiredExprsSig.exprsCorWithRegs.list[[1]])
# first plot
plot(x, y1)
# second plot
par(new = TRUE)
plot(x, y2, axes = FALSE, xlab = "", ylab = "")

for (i in 1:length(subtypes)) {
  subtype <- subtypes[i]
  outfn <- paste0("analyses/0_0.1/",subtype,"corWithRegExprsHistOverlay.pdf")
  data.df <- data.frame("Correlation"=c(rewiredExprsSig.exprsCorWithRegs.list[[i]],
                                        nonRewiredExprsSig.exprsCorWithRegs.list[[i]]),
                        "Subgroup"=factor(c(rep("With differential regulation",
                                                length(rewiredExprsSig.exprsCorWithRegs.list[[i]])),
                                            rep("Without differential regulation",
                                                length(nonRewiredExprsSig.exprsCorWithRegs.list[[i]]))),
                                          levels=c("With differential regulation",
                                                   "Without differential regulation")))
  colvec <- c("gray","gray")
  colvec[1] <- subtypeColors[i]
  ggplot(data.df, aes(x=Correlation, fill=Subgroup)) + geom_density(alpha=.3) +
    scale_fill_manual(values=(colvec)) +
    theme(panel.background = element_blank()) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + ggtitle(subtype) + xlim(c(0.2,0.8))
  ggsave(file = outfn, width = 10, height = 4.5)
}
# The non-rewired expression signature genes do not have TFs that show signature
# behavior with these target genes.
nonRewiredExprsSigTF.percent.list <- lapply(c(1:length(subtypes)),function(x){
  regTFs.list <- unique(unlist(lapply(nonRewiredExprsSig.list[[x]],function(tgt){
    rownames(sigSubtype.filtered)[sigSubtype.filtered[,tgt]==as.character(x)]
  })))
  length(intersect(regTFs.list,names(exprsSigSubtype.domSubtype[[x]])))/length(regTFs.list)
})
# Are any of the expression signatures TFs themselves?
# (If so do they regulate the genes within the same signature set?)
# exprsSigTFs.list <- lapply(c(1:length(cl)),function(x){
#   intersect(cstExprsRegTgts[[x]],rownames(geneRegs.TF[[1]]))
# })
exprsSigTFs.list <- lapply(c(1:length(subtypes)),function(x){
  intersect(exprsSigSubtype.list[[x]],rownames(capacities.log2.list[[1]]))
})
# Find out whether TFs within expression signature regulate
# target genes showing both signature expression AND signature regulation
exprsSigTFs.reg.list <- lapply(c(1:length(subtypes)),function(x){
  intersect(exprsSigTFs.list[[x]],cstExprsRegTgts.TFs[[x]])
})
# Updated version
exprsSigTFs.reg.list <- lapply(c(1:length(subtypes)),function(x){
  intersect(exprsSigTFs.list[[x]],exprsSigSubtype.rewiredTgtTFs[[x]])
})


# What about looking at signatures in a target-centric way?




# Look at correlation with scRNA-seq data
# First need to check how many target genes in the scRNA-seq data have
# all of their regulators covered in the dataset
scExprs <- read.delim("data/GSE84465_GBM_All_data.csv",
                      sep=" ",header = T, row.names=1, check.names=F,quote = "\"")
scExprs <- data.matrix(scExprs)
scExprs <- scExprs[rowSums(scExprs)!=0,]
expressedGenes <- which(apply(scExprs,1,function(x){
  sum(x>0) > 0.1*ncol(scExprs)
}))
geneUniv <- rownames(scExprs)[expressedGenes]
gene_TFs <- read.table("data/piq/TFs_for_genes.txt",sep="\t",header=F)

mappableEdges <- lapply(1:nrow(gene_TFs),function(x){
  if (!(gene_TFs[x,1] %in% geneUniv) | !(gene_TFs[x,2] %in% geneUniv)) {
    F
  } else {
    T
  }
})

mapped_gene_TFs <- gene_TFs[unlist(mappableEdges),]

TFList <- unique(c(as.character(mapped_gene_TFs[,2])))
tgtGenes <- unique(as.character(mapped_gene_TFs[,1]))
withAllRegs <- lapply(tgtGenes,function(x){
  if (length(intersect(regulators[[x]],geneUniv))==length(regulators[[x]])) {
    T
  } else {
    F
  }
})
tgtGenesWithAllRegs <- tgtGenes[unlist(withAllRegs)]
write.table(tgtGenesWithAllRegs,"data/sc_targetGeneList.txt",sep="\t",
            row.names = F,col.names = F,quote = F)
allTFs.mappable <- lapply(tgtGenesWithAllRegs,function(x){
  regulators[[x]]
})
allTFs.mappable <- unique(unlist(allTFs.mappable))
TFList <- allTFs.mappable



