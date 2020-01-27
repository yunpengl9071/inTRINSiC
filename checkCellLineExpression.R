rm(list=ls())

setwd('~/Data/multilayerNetwork/')
centroids <- read.table("data/subtype_centroids.txt",sep="\t",header=T)
sig.genes <- as.character(centroids$Gene.Symbol)
# Try using centroids for actual expression profiles instead...
expr <- read.csv("data/hgMatrix.txt",sep="\t")
labels <- read.table("data/class_labels.txt",sep="\t",header=F)
labels <- as.character(labels[,1])
colvec <- rep("Classical",ncol(expr))
colvec[labels=="Neural"] <- "Neural"
colvec[labels=="Proneural"] <- "Proneural"
colvec[labels=="Mesenchymal"] <- "Mesenchymal"
names(labels) <- colnames(expr)

# tumor_data <- read.table("data/hgmatrix.txt",sep="\t",
#                          header=T,row.names=1)
# tumor_data <- as.matrix(tumor_data)
# tumor_data <- sweep(tumor_data,1, apply(tumor_data,1,mean,na.rm=T))
# tumor_data <- t(apply(tumor_data,1,function(x){x/sd(x)}))
# 
# tumor_data_sg <- data.matrix(tumor_data_sg[sig.genes,])
# 
# class_labels <- apply(t(scale(t(data.matrix(expr[sig.genes,])))), 2, function(x){
#   best_label <- ""
#   mindist <- Inf
#   for (label in c("Proneural","Neural","Classical","Mesenchymal")) {
#     #cendist <- 1-cor(x,centroids[,label],method = "pearson")
#     cendist <- sum((x-centroids[,label])^2)
#     if (cendist < mindist) {
#       best_label <- label
#       mindist <- cendist
#     }
#   }
#   return(best_label)
# })

# ccle.exprs <- read.delim("~/Data/RegevLab/CCLE_exprs_cleaned_avgProbes.txt",
#                          header=T,row.names=1,check.names=F)
# ccle.cns.exprs <- data.matrix(ccle.exprs[,grep("CENTRAL_NERVOUS",colnames(ccle.exprs),fixed=T)])
# write.table(ccle.cns.exprs,"CCLE.CNS.exprs.txt",sep="\t",quote=F)
ccle.cns.exprs <- data.matrix(read.delim("data/CCLE.CNS.exprs.txt",header=T,row.names=1,
                                         check.names=F))
avail.sigGenes <- intersect(rownames(ccle.cns.exprs),rownames(expr))

# Visualize correlation distance to GBM patient subtype centroids

library(Rtsne)

rownames(centroids) <- centroids$Gene.Symbol
centroids <- data.matrix(centroids[,2:5])
#combined.exprs <- cbind(centroids[avail.sigGenes,],ccle.cns.exprs[avail.sigGenes,])
combined.exprs <- cbind(t((t(expr[avail.sigGenes,]))),t((t(ccle.cns.exprs[avail.sigGenes,]))))
# Perform quantile normalization
# install.packages(pkgs="~/Downloads/preprocessCore_1.34.0.tar.gz",repos=NULL,type="source")
# library(preprocessCore)
#combined.exprs.qn <- normalize.quantiles(data.matrix(combined.exprs))
# Try limma's batch effect removal tool
library(limma)
combined.exprs.rbe <- removeBatchEffect(data.matrix(combined.exprs)[intersect(rownames(combined.exprs),
                                                                  sig.genes),],
                                        batch=c(rep("A",ncol(expr)),rep("B",ncol(ccle.cns.exprs))),
                                        covariates=NULL)
# intersGenes <- intersect(rownames(combined.exprs),
#                          sig.genes)
# Looks much better with RBE!!! (Is this actually good?)
dist.mat <- dist(t(combined.exprs.rbe))
#mads <- rev(order(apply(combined.exprs.rbe,1,mad)))
# tsne.clustMat <- Rtsne(t(combined.exprs.rbe),is_distance = F,
#                        perplexity=15,initial_dims = 500,pca=F)
tsne.clustMat.dist <- Rtsne(dist.mat,is_distance = T,
                       perplexity=15,pca=F)
ccle.exprs.rbe <- combined.exprs.rbe[,colnames(ccle.cns.exprs)]
ccle.exprs.rbe.scaled <- t(scale(t(ccle.exprs.rbe)))
ccle.availSigGenes <- intersect(rownames(ccle.exprs.rbe.scaled),rownames(centroids))
ccle.class.labels <- apply(ccle.exprs.rbe.scaled[ccle.availSigGenes,], 2, function(x){
  best_label <- ""
  mindist <- Inf
  for (label in c("Proneural","Neural","Classical","Mesenchymal")) {
    #cendist <- 1-cor(x,centroids[,label],method = "pearson")
    cendist <- sum((x-centroids[ccle.availSigGenes,label])^2)
    if (cendist < mindist) {
      best_label <- label
      mindist <- cendist
    }
  }
  return(best_label)
})

cellLineCol <- rep("red",length(ccle.class.labels))
cellLineCol[ccle.class.labels == "Neural"] <- "blue"
cellLineCol[ccle.class.labels == "Proneural"] <- "green"
cellLineCol[ccle.class.labels == "Mesenchymal"] <- "orange"
colvec <- rep("red",ncol(expr))
colvec[labels=="Neural"] <- "blue"
colvec[labels=="Proneural"] <- "green"
colvec[labels=="Mesenchymal"] <- "orange"
colvec.all <- c(colvec,cellLineCol)

# Try diffusion map for visualization
library(diffusionMap)
combined.exprs.rbe.diffuse <- diffuse(dist(t(combined.exprs.rbe)))
plot(combined.exprs.rbe.diffuse$X[,1:2],pch=20,col=colvec.all)

subtypeColors <- c("brown3","cornflowerblue","mediumseagreen","goldenrod1" )

library(gplots)
ggplot(data = data.frame(tsne.clustMat.dist$Y), aes(x = tsne.clustMat.dist$Y[,1],
                                                    y = tsne.clustMat.dist$Y[,2],
                                                    shape = as.factor(c(rep("TCGA",ncol(expr)),
                                                                        rep("CCLE",ncol(ccle.cns.exprs)))),
                                                    color = as.factor(c(labels,ccle.class.labels)))) + 
  geom_point(size = 3) + scale_shape_manual(values=c(13, 16)) +
  scale_color_manual(values=c('brown3','goldenrod1','cornflowerblue','mediumseagreen')) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA,size = 2),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20)) + xlim(-35,35) + ylim(-40,20)
  #geom_text(size = 2, color = 'black')

ggplot(data = data.frame(tsne.clustMat.dist$Y[1:544,]), aes(x = tsne.clustMat.dist$Y[1:544,1],
                                                    y = tsne.clustMat.dist$Y[1:544,2],
                                                    shape = as.factor(c(rep("TCGA",ncol(expr)))),
                                                    color = as.factor(c(labels)))) + 
  geom_point(size = 3) + scale_shape_manual(values=c(16)) +
  scale_color_manual(values=c('brown3','goldenrod1','cornflowerblue','mediumseagreen')) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA,size = 2),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20)) + xlim(-35,35) + ylim(-40,20)
#geom_text(size = 2, color = 'black')

avana.scores <- read.csv('data/gene_effect.csv',row.names = 1,check.names = F)
avana.celllinenames <- read.csv('data/DepMap-2018q3-celllines.csv',row.names = 1,check.names=F)
ccle.C <- match(names(ccle.class.labels)[ccle.class.labels=="Classical"],avana.celllinenames$CCLE_Name)
ccle.N <- match(names(ccle.class.labels)[ccle.class.labels=="Neural"],avana.celllinenames$CCLE_Name)
ccle.P <- match(names(ccle.class.labels)[ccle.class.labels=="Proneural"],avana.celllinenames$CCLE_Name)
ccle.M <- match(names(ccle.class.labels)[ccle.class.labels=="Mesenchymal"],avana.celllinenames$CCLE_Name)
ccle.C <- rownames(avana.celllinenames)[ccle.C]
ccle.N <- rownames(avana.celllinenames)[ccle.N]
ccle.P <- rownames(avana.celllinenames)[ccle.P]
ccle.M <- rownames(avana.celllinenames)[ccle.M]
TFoI <- "STAT5A "
TFoI.C <- avana.scores[ccle.C,grep(TFoI,colnames(avana.scores))]
TFoI.C <- TFoI.C[complete.cases(TFoI.C)]
TFoI.N <- avana.scores[ccle.N,grep(TFoI,colnames(avana.scores))]
TFoI.N <- TFoI.N[complete.cases(TFoI.N)]
TFoI.P <- avana.scores[ccle.P,grep(TFoI,colnames(avana.scores))]
TFoI.P <- TFoI.P[complete.cases(TFoI.P)]
TFoI.M <- avana.scores[ccle.M,grep(TFoI,colnames(avana.scores))]
TFoI.M <- TFoI.M[complete.cases(TFoI.M)]

# TCGA.exprs.rbe <- combined.exprs.rbe[,colnames(expr)]
# TCGA.exprs.rbe.scaled <- t(scale(t(TCGA.exprs.rbe)))
# TCGA.centroids <- lapply(c("Proneural","Neural","Classical","Mesenchymal"),function(x){
#   TCGA.subtype.samples <- names(labels)[labels==x]
#   TCGA.subtype.centroid <- rowMeans(data.matrix(TCGA.exprs.rbe.scaled[ccle.availSigGenes,
#                                                                       TCGA.subtype.samples]))
# })
# TCGA.centroids <- data.matrix(do.call(cbind.data.frame,TCGA.centroids))
# colnames(TCGA.centroids) <- c("Proneural","Neural","Classical","Mesenchymal")

centroids.avail <- centroids[rownames(ccle.exprs.rbe.scaled),]
ccle.class.scores <- apply(ccle.exprs.rbe.scaled,2,function(x){
  dists <- apply(centroids.avail,2,function(y){
    #1-cor(x,y,method = "spearman")
    1/sqrt(sum((x-y)^2))
  })
  dists
})
#ccle.class.scores <- t(scale(t(ccle.class.scores)))

ccle.Classical <- names(ccle.class.labels)[ccle.class.labels=="Classical"]
ccle.Classical.scores <- apply(ccle.class.scores[,ccle.Classical],2,function(x){
  x["Classical"]
})
ccle.Classical.scores <- (ccle.Classical.scores - min(ccle.Classical.scores))/
  (max(ccle.Classical.scores)-min(ccle.Classical.scores))
ccle.Classical.scores <- rev(sort(ccle.Classical.scores))
write.table(ccle.Classical.scores,"analyses/ccle.Classical.scores.txt",
            sep="\t",row.names=gsub("_CENTRAL_NERVOUS_SYSTEM","",
                                    names(ccle.Classical.scores),fixed=T),
            quote=F)

ccle.Neural <- names(ccle.class.labels)[ccle.class.labels=="Neural"]
ccle.Neural.scores <- apply(ccle.class.scores[,ccle.Neural],2,function(x){
  x["Neural"]
})
ccle.Neural.scores <- (ccle.Neural.scores - min(ccle.Neural.scores))/
  (max(ccle.Neural.scores)-min(ccle.Neural.scores))
ccle.Neural.scores <- rev(sort(ccle.Neural.scores))
write.table(ccle.Neural.scores,"analyses/ccle.Neural.scores.txt",
            sep="\t",row.names=gsub("_CENTRAL_NERVOUS_SYSTEM","",
                                    names(ccle.Neural.scores),fixed=T),
            quote=F)

ccle.Proneural <- names(ccle.class.labels)[ccle.class.labels=="Proneural"]
ccle.Proneural.scores <- apply(ccle.class.scores[,ccle.Proneural,drop=F],2,function(x){
  x[1]
})
ccle.Proneural.scores <- (ccle.Proneural.scores - min(ccle.Proneural.scores))/
  (max(ccle.Proneural.scores)-min(ccle.Proneural.scores))
ccle.Proneural.scores <- rev(sort(ccle.Proneural.scores))
write.table(ccle.Proneural.scores,"analyses/ccle.Proneural.scores.txt",
            sep="\t",row.names=gsub("_CENTRAL_NERVOUS_SYSTEM","",
                                    names(ccle.Proneural.scores),fixed=T),
            quote=F)

ccle.Mesenchymal <- names(ccle.class.labels)[ccle.class.labels=="Mesenchymal"]
ccle.Mesenchymal.scores <- apply(ccle.class.scores[,ccle.Mesenchymal],2,function(x){
  x["Mesenchymal"]
})
ccle.Mesenchymal.scores <- (ccle.Mesenchymal.scores - min(ccle.Mesenchymal.scores))/
  (max(ccle.Mesenchymal.scores)-min(ccle.Mesenchymal.scores))
ccle.Mesenchymal.scores <- rev(sort(ccle.Mesenchymal.scores))
write.table(ccle.Mesenchymal.scores,"analyses/ccle.Mesenchymal.scores.txt",
            sep="\t",row.names=gsub("_CENTRAL_NERVOUS_SYSTEM","",
                                    names(ccle.Mesenchymal.scores),fixed=T),
            quote=F)

cellLineCol <- rep("Classical",length(ccle.class.labels))
cellLineCol[ccle.class.labels == "Neural"] <- "Neural"
cellLineCol[ccle.class.labels == "Proneural"] <- "Proneural"
cellLineCol[ccle.class.labels == "Mesenchymal"] <- "Mesenchymal"

#plot(tsne.clustMat.dist$Y,pch=16,col=c(colvec,cellLineCol))

cellLineNames <- gsub("_CENTRAL_NERVOUS_SYSTEM","",colnames(ccle.cns.exprs),fixed=T)
textLabels <- (c(rep("",ncol(expr)),cellLineNames))
library(ggplot2)
library(Rtsne)
dist.mat <- dist(((t(combined.exprs.rbe))))
tsne.clustMat.dist <- Rtsne(dist.mat,is_distance = T,
                            perplexity=10,pca=F)

ggplot(data = data.frame(tsne.clustMat.dist$Y), aes(x = tsne.clustMat.dist$Y[,1],
                                                    y = tsne.clustMat.dist$Y[,2],
                                            shape = as.factor(c(rep("TCGA",ncol(expr)),
                                                      rep("CCLE",ncol(ccle.cns.exprs)))),
                                            color = as.factor(c(labels,ccle.class.labels)),
                                            label = as.factor(textLabels))) + 
  geom_point(size = 3) + scale_shape_manual(values=c(7, 16)) +
  scale_color_manual(values=c('red','orange','blue','green')) +
  geom_text(size = 2, color = 'black')


library(gplots)
myheatmapfun <- colorRampPalette((c(rgb(15,61,111,maxColorValue=255),
                                    rgb(125,182,212,maxColorValue=255),
                                    rgb(231,224,219,maxColorValue=255),
                                    rgb(234,145,118,maxColorValue=255),
                                    rgb(102,3,32,maxColorValue=255))))

heatmap.data <- combined.exprs.rbe
csc <- c(colvec,cellLineCol)
csc[csc=="Classical"] <- "red"
csc[csc=="Neural"] <- "blue"
csc[csc=="Proneural"] <- "green"
csc[csc=="Mesenchymal"] <- "orange"
labcol <- rep("",length(csc))
labcol[(length(colvec)+1):length(labcol)] <- cellLineNames
dev.off()
colOrder <- kmeans(t(combined.exprs.rbe),centers=4)$cluster
hmRes <- heatmap.2(heatmap.data[,order(colOrder)],
                   Colv=NA,
                   #Rowv=NA,
                   #symm=T,revC=T,
                   scale="none",
                   density.info="none", trace="none",
                   col=myheatmapfun(75),
                   #hclustfun = myclustfun,
                   #reorderfun = function(d,w) { d },
                   cexCol = 0.8,
                   #srtCol = 45,
                   cexRow = 0.8,
                   #offsetCol = -0.5,
                   ColSideColors = csc,
                   dendrogram="none",
                   key.title=NA,
                   key.xlab=NA,
                   #lhei=c(1,8),
                   breaks=seq(-1.5,1.5,length.out=76),
                   labRow = "",
                   labCol = labcol[order(colOrder)]
)

dev.off()
csc <- c(colvec)
csc[csc=="Classical"] <- "red"
csc[csc=="Neural"] <- "blue"
csc[csc=="Proneural"] <- "green"
csc[csc=="Mesenchymal"] <- "orange"
hmRes <- heatmap.2(data.matrix(expr[avail.sigGenes,]),
                   #Colv=NA,
                   #Rowv=NA,
                   #symm=T,revC=T,
                   scale="row",
                   density.info="none", trace="none",
                   col=myheatmapfun(75),
                   #hclustfun = myclustfun,
                   #reorderfun = function(d,w) { d },
                   cexCol = 0.8,
                   #srtCol = 45,
                   cexRow = 0.8,
                   #offsetCol = -0.5,
                   ColSideColors = csc,
                   dendrogram="none",
                   key.title=NA,
                   key.xlab=NA,
                   #lhei=c(1,8),
                   breaks=seq(-1.5,1.5,length.out=76),
                   labRow = "",
                   labCol = ""
)

# plot(tsne.clustMat.dist$Y[match(colnames(ccle.cns.exprs),colnames(combined.exprs.rbe)),],
#      pch=16,col=cellLineCol,
#      xlab="dimension 1", ylab = "dimension 2")
# cellLineNames <- gsub("_CENTRAL_NERVOUS_SYSTEM","",colnames(ccle.cns.exprs),fixed=T)
# text(tsne.clustMat.dist$Y[match(colnames(ccle.cns.exprs),colnames(combined.exprs.rbe)),1],
#      tsne.clustMat.dist$Y[match(colnames(ccle.cns.exprs),colnames(combined.exprs.rbe)),2]+0.75,
#      cellLineNames,cex=0.5)
# plot(tsne.clustMat$Y,pch=16,
#      col=c(colvec,rep("black",ncol(ccle.cns.exprs))),
#      xlab="dimension 1", ylab = "dimension 2")
# colvec2 <- c(colvec,rep("black",ncol(ccle.cns.exprs)))
# colvec2[topnames.idx] <- "purple"
# shapevec <- rep(16,length(colvec2))
# shapevec[topnames.idx] <- 12
# plot(tsne.clustMat$Y,pch=shapevec,
#      col=colvec2,
#      xlab="dimension 1", ylab = "dimension 2")

#text(tsne.clustMat.dist$Y,colnames(combined.exprs),cex=0.3)
#text(tsne.clustMat$Y[topnames.idx,],colnames(combined.exprs)[topnames.idx],cex=1)
# MDS
# library(MASS)
# d <- as.dist(1-cor((combined.exprs.rbe)))
# #fit <- isoMDS(d, k=2)
# fit <- cmdscale(d, k = 2)
# #x <- fit$points[,1]
# #y <- fit$points[,2]
# plot(fit, xlab="Coordinate 1", ylab="Coordinate 2",
#      main="metric MDS",	pch=16, col = c(colvec,rep("black",ncol(ccle.cns.exprs))))

# Now for each cell line try classifying into the most likely subtype.
# Reward closeness to one centroid while far from the rest (on average). (Not doing this...)
# # Use random forest!
# library(randomForest)
# 
# rf <- randomForest((t(combined.exprs.rbe[intersect(sig.genes,avail.sigGenes),colnames(expr)])),
#                    as.factor(labels),importance = T,
#                    keep.forest = T)
# 
# classRes <- predict(rf,t(combined.exprs.rbe[intersect(sig.genes,avail.sigGenes),
#                                             colnames(ccle.cns.exprs)]))


# centroids <- lapply(names(table(labels)),function(x){
#   rowMeans(combined.exprs.rbe[,colnames(expr)][,names(labels)[labels==x]])
# })
# subtypeDist <- apply((combined.exprs.rbe[,colnames(ccle.cns.exprs)]),2,function(x){
#   sims <- sapply(centroids,function(y){
#     cor(x,y)
#   })
#   sims
# })
# subtypeDist <- t(subtypeDist)
# colnames(subtypeDist) <- names(table(labels))
# classRes <- apply(subtypeDist,1,function(x){
#   colnames(subtypeDist)[which(x==max(x))]
# })

# plot(tsne.clustMat.dist$Y[match(colnames(ccle.cns.exprs),colnames(combined.exprs.rbe)),],
#      pch=16,col=cellLineCol,
#      xlab="dimension 1", ylab = "dimension 2")
# cellLineNames <- gsub("_CENTRAL_NERVOUS_SYSTEM","",colnames(ccle.cns.exprs),fixed=T)
# text(tsne.clustMat.dist$Y[match(colnames(ccle.cns.exprs),colnames(combined.exprs.rbe)),1],
#      tsne.clustMat.dist$Y[match(colnames(ccle.cns.exprs),colnames(combined.exprs.rbe)),2]+0.75,
#      cellLineNames,cex=0.5)
# names(classRes) <- cellLineNames
#write.table(classRes,"CCLE_classLabels.txt",sep="\t",quote=F)

# Maybe should classify the CCLE lines using the same way TCGA samples were
# classified: nearest centroid method.

# sigGenes.inters <- intersect(sig.genes,avail.sigGenes)
# ccle.cns.zscores <- t((t(combined.exprs.rbe)))[,545:613]
# centroids.avail <- centroids[sigGenes.inters,]
# ccle.cns.sig.zscores <- ccle.cns.zscores[sigGenes.inters,]

# classes <- c("Proneural","Neural","Classical","Mesenchymal")
# ccle.class.labels <- apply(ccle.cns.sig.zscores,2,function(x){
#   dists <- apply(centroids.avail,2,function(y){
#     sqrt(sum((x-y)^2))
#     #1-cor(x,y)
#   })
#   classes[which(dists==min(dists))]
# })
# 
# ccle.class.labels <- as.character(classRes)
# names(ccle.class.labels) <- names(classRes)
