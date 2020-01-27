crossValidation <- function(rnseed,
                            fnprefix="data/regrOutput/CV/",
                            fnsuffix="_TFGeneCapacities_ri_nomi__lambda0_0_CV_",
                            CV.fold=5,
                            method="RMSE") {
  # set.seed(rnseed)
  # CV.TVSamples <- lapply(c(1:length(sampleByLabels)),function(x){
  #   nsamples.subtype <- length(sampleByLabels[[x]])
  #   nsamp <- round(nsamples.subtype/CV.fold)
  #   sampleNames.shuffled <- sample((sampleByLabels[[x]]),size=nsamples.subtype,
  #                                  replace = F)
  #   temp <- list()
  #   for (i in c(1:CV.fold)) {
  #     start <- (i-1)*nsamp + 1
  #     end <- min(i*nsamp,nsamples.subtype)
  #     trainingSamples <- sampleNames.shuffled[-c(start:end)]
  #     validationSamples <- sampleNames.shuffled[c(start:end)]
  #     temp[[i]] <- list(trn=trainingSamples,vln=validationSamples)
  #   }
  #   temp
  # })
  capmat.CV.list <- lapply(1:CV.fold,function(x){
    cat(".")
    capmat.list <- lapply(1:length(subtypes),function(y){
      #cat(y)
      fn <- paste0(fnprefix,subtypes[y],fnsuffix,
                   as.character(x),".txt")
      temp <- data.matrix(read.csv(fn,sep="\t",check.names=F))
      if (!("ZNF8" %in% rownames(temp))) {
        tempRow <- rep(1,ncol(temp))
        TFs <- rownames(temp)
        temp <- rbind(temp,tempRow)
        rownames(temp) <- c(TFs,"ZNF8")
        temp
      } else {
        temp
      }
    })
  })
  cat("\n")
  CV.error.list <- lapply(c(1:length(subtypes)),function(x){
    exprs.est <- lapply(1:CV.fold,function(i){
      fn <- paste0(fnprefix,subtypes[x],fnsuffix,
                   as.character(i),".vlnExprs.txt")
      #file.exists(fn)
      if (F){
        data.matrix(read.delim(fn,header=T,row.names=1,check.names=F))
      } else {
        validationSamples <- CV.TVSamples[[x]][[i]]$vln
        tfmat <- expr[allTFs,names(validationSamples)]
        capmat <- capmat.CV.list[[i]][[x]]
        exprEstimates <- lapply(colnames(capmat),function(z){
          regs <- regulators[[z]]
          tfsubmat <- tfmat[regs,,drop=F]
          capvec <- capmat[regs,z]
          const <- capmat["const",z]
          # cat(x)
          # cat('\n')
          #cat(z)
          exprsEst <- as.numeric(apply(tfsubmat,2,model.pred,
                                       mivec=c(),fvecmi=c(),fvectf=capvec,const=const))
        })
        exprEstimates <- do.call(rbind.data.frame,exprEstimates)
        exprEstimates <- data.matrix(exprEstimates)
        colnames(exprEstimates) <- names(validationSamples)
        rownames(exprEstimates) <- colnames(capmat)
        cat(".")
        write.table(exprEstimates,fn,sep="\t",quote=F)
        exprEstimates
      }
    })
    cat("\n")
    exprs.est 
  })
  tr.error.list <- lapply(c(1:length(subtypes)),function(x){
    exprs.est <- lapply(1:CV.fold,function(i){
      fn <- paste0(fnprefix,subtypes[x],fnsuffix,
                   as.character(i),".trnExprs.txt")
      #file.exists(fn)
      if (F){
        data.matrix(read.delim(fn,header=T,row.names=1,check.names=F))
      } else {
        trainingSamples <- CV.TVSamples[[x]][[i]]$trn
        tfmat <- expr[allTFs,names(trainingSamples)]
        capmat <- capmat.CV.list[[i]][[x]]
        exprEstimates <- lapply(colnames(capmat),function(z){
          regs <- regulators[[z]]
          tfsubmat <- tfmat[regs,,drop=F]
          capvec <- capmat[regs,z]
          const <- capmat["const",z]
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
        write.table(exprEstimates,fn,sep="\t",quote=F)
        exprEstimates
      }
    })
    cat("\n")
    exprs.est
  })
  # for (i in c(1:length(subtypes))) {
  #   fn <- paste0(fnprefix,subtypes[y],fnsuffix,
  #                as.character(x),".trnExprs.txt")
  #   write.table(CV.error.list,fn,sep="\t",quote=F)
  # }
  # Use symmetric mean absolute percentage error (sMAPE)
  # also try RMSE
  if (method == "sMAPE") {
    CV.error.mat <- matrix(0,length(subtypes),CV.fold)
    for (i in 1:length(subtypes)) {
      for (j in 1:CV.fold) {
        exprs.obs <- expr[rownames(CV.error.list[[i]][[j]]),
                          colnames(CV.error.list[[i]][[j]])]
        temp <- (abs(CV.error.list[[i]][[j]]) + abs(data.matrix(exprs.obs)))/2
        CV.error.mat[i,j] <- mean(rowMeans(abs(CV.error.list[[i]][[j]]-data.matrix(exprs.obs))/temp))
      }
    }
    tr.error.mat <- matrix(0,length(subtypes),CV.fold)
    for (i in 1:length(subtypes)) {
      for (j in 1:CV.fold) {
        exprs.obs <- expr[rownames(tr.error.list[[i]][[j]]),
                          colnames(tr.error.list[[i]][[j]])]
        temp <- (abs(tr.error.list[[i]][[j]]) + abs(data.matrix(exprs.obs)))/2
        tr.error.mat[i,j] <- mean(rowMeans(abs(tr.error.list[[i]][[j]]-data.matrix(exprs.obs))/temp))
      }
    }
  } else if (method == "RMSE") {
    CV.error.mat <- matrix(0,length(subtypes),CV.fold)
    for (i in 1:length(subtypes)) {
      for (j in 1:CV.fold) {
        exprs.obs <- expr[rownames(CV.error.list[[i]][[j]]),
                          colnames(CV.error.list[[i]][[j]])]
        CV.error.mat[i,j] <- mean((rowMeans((CV.error.list[[i]][[j]]-data.matrix(exprs.obs))^2))^0.5)
      }
    }
    tr.error.mat <- matrix(0,length(subtypes),CV.fold)
    for (i in 1:length(subtypes)) {
      for (j in 1:CV.fold) {
        exprs.obs <- expr[rownames(tr.error.list[[i]][[j]]),
                          colnames(tr.error.list[[i]][[j]])]
        tr.error.mat[i,j] <- mean((rowMeans((tr.error.list[[i]][[j]]-data.matrix(exprs.obs))^2))^0.5)
      }
    }
  } else if (method == "log2ratio") {
    CV.error.mat <- matrix(0,length(subtypes),CV.fold)
    for (i in 1:length(subtypes)) {
      for (j in 1:CV.fold) {
        exprs.obs <- expr[rownames(CV.error.list[[i]][[j]]),
                          colnames(CV.error.list[[i]][[j]])]
        temp <- abs(log2(CV.error.list[[i]][[j]]/(data.matrix(exprs.obs))))
        temp <- temp[complete.cases(temp),]
        CV.error.mat[i,j] <- median(rowMeans(temp))
      }
    }
    tr.error.mat <- matrix(0,length(subtypes),CV.fold)
    for (i in 1:length(subtypes)) {
      for (j in 1:CV.fold) {
        exprs.obs <- expr[rownames(tr.error.list[[i]][[j]]),
                          colnames(tr.error.list[[i]][[j]])]
        temp <- abs(log2(tr.error.list[[i]][[j]]/(data.matrix(exprs.obs))))
        temp <- temp[complete.cases(temp),]
        tr.error.mat[i,j] <- median(rowMeans(temp))
        #tr.error.mat[i,j] <- mean((rowMeans(log2((tr.error.list[[i]][[j]]+1E-5)/(data.matrix(exprs.obs)+1E-5)))))
      }
    }
  }
  list(CV.error.mat,tr.error.mat)
  # write.table(t(CV.error.mat),"CV/CV_error_mat.txt",row.names=F,
  #             col.names=F,sep="\t",quote=F)
  # write.table(t(tr.error.mat),"CV/tr_error_mat.txt",row.names=F,
  #             col.names=F,sep="\t",quote=F)
}



# test <- lapply(5:10,function(x){
#   lapply(1:1,function(i){
#     lapply(1:5, function(x){
#       cat(x)
#     })
#     cat("\n")
#     cat(x)
#     cat("\n")
#   })
# })
