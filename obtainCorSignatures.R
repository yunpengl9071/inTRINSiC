# If shuffled==F, provide unshuffled matrix to capmat
# If shuffled==T, provide subtype-specific matrices from previous shuffling to capmat

obtainCorSignatures <- function(nsd,cormat,shuffled=F,sigTFs=NULL,adjPThres=0.0001) {
  if (!shuffled) {
    cat("Generating subtype-specific correlation value matrix...\n")
    nonzeros <- lapply(cormat,function(x){which(x!=0)})
    nonzeros <- Reduce("union",nonzeros)
    subtypeSpecific.cor.mat.list <- lapply(c(1:length(subtypes)),function(x){
      cat(".\n")
      subtypeSpecific.cor.mat <- matrix(1,nrow(cormat[[1]]),ncol(cormat[[1]]))
      allIdx <- nonzeros
      # Also filter F values
      for (idx in allIdx) {
        rand.cor.list <- sapply(c(1:niter),function(y){
          cormat.shuffled[[y]][[x]][idx]
        })
        rf.sd <- sd((rand.cor.list))
        rf.mean <- mean((rand.cor.list))
        rf.upper <- rf.mean + nsd * rf.sd
        rf.lower <- rf.mean - nsd * rf.sd
        corval <- (cormat.unshuffled[[x]][idx])
        #log2FVal <- log2(Fval)
        if ((corval) > rf.upper | (corval) < rf.lower) {
          subtypeSpecific.cor.mat[idx] <- corval
        }
      }
      subtypeSpecific.cor.mat
    }) 
    subtypeSpecific.cor.mat.list <- lapply(subtypeSpecific.cor.mat.list,function(x){
      temp <- x
      rownames(temp) <- rownames(cormat[[1]])
      colnames(temp) <- colnames(cormat[[1]])
      temp
    })
    geneRegs.all <- subtypeSpecific.cor.mat.list
  } else {
    geneRegs.all <- cormat
  }
  cat("\nDetermining significantly different TFs using ANOVA...\n")
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
  TFs.anova.p <- lapply(allTFs,function(x){
    subtypeVec <- c(rep("Classical",ncol(meanTFCV)),
                    rep("Neural",ncol(meanTFCV)),
                    rep("Proneural",ncol(meanTFCV)),
                    rep("Mesenchymal",ncol(meanTFCV)))
    corVec <- c(geneRegs.all[[1]][x,],
              geneRegs.all[[2]][x,],
              geneRegs.all[[3]][x,],
              geneRegs.all[[4]][x,])
    TF.df <- cbind.data.frame(subtypeVec,corVec)
    summary(aov(corVec ~ subtypeVec, data = TF.df))[[1]][1,"Pr(>F)"]
  })
  names(TFs.anova.p) <- allTFs
  TFs.anova.p <- unlist(TFs.anova.p)
  sigDiffTFs <- names(which(p.adjust(TFs.anova.p,method = "BH")<adjPThres))
  sigSubtype <- matrix(0,nrow(geneRegs.TF[[1]]),ncol(geneRegs.TF[[1]]))
  if (T) {
    cat("\nGenerating TF-TF-subtype signature tuples using the leave-one-out SD reduction method...\n")
    sigDiffCorMat <- geneRegs.TF[[1]][,]
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
    SDReduction.mat <- matrix(0,nrow(geneRegs.TF[[1]]),ncol(geneRegs.TF[[1]]))
    for (i in c(1:nrow(geneRegs.TF[[1]]))) {
      for (j in c(1:ncol(geneRegs.TF[[1]]))) {
        Fval <- sapply(leaveOneSD,function(x){
          x[i,j]
        })
        if (overallSD[i,j] > 0) {
          sigSubtype[i,j] <- which(Fval==min(Fval))[1]
          SDReduction.mat[i,j] <- overallSD[i,j]-Fval[sigSubtype[i,j]]
        }
      }
    }
    rownames(sigSubtype) <- rownames(geneRegs.TF[[1]])
    colnames(sigSubtype) <- colnames(geneRegs.TF[[1]])
  }
  if (!shuffled) {
    res <- list(subtypeSpecificCorMats=subtypeSpecific.cor.mat.list,
                sigTFs.anova=sigDiffTFs,
                SDReductionMat=SDReduction.mat,
                sigSubtypeMat=sigSubtype)
  } else {
    res <- list(subtypeSpecificCorMats=cormat,
                sigTFs.anova=sigDiffTFs,
                SDReductionMat=SDReduction.mat,
                sigSubtypeMat=sigSubtype)
  }
  res
}