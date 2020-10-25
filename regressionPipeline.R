rm(list=ls())

# Set the library paths on the server where you run this script
.libPaths('/ahg/regevdata/projects/NHLeukemia/rlib/')

setwd("/ahg/regevdata/projects/txnRegModeling/regression/")

args = commandArgs(trailingOnly = T)
# You can run each subtype/cell type using a command line argument
subtypes <- args[1]
#subtypes <- c("Classical","Proneural","Mesenchymal")
cat("Preparing input...\n\n")

# Load utility functions
source("coreFunctions.R")

# Read in expression matrix in the format of CSV/TSV with both row names and col names
expr <- read.delim("data/hgMatrix_CGGA.txt",sep=" ",header=T,row.names=1,check.names=F)
expr <- data.matrix(expr)

# expr <- read.csv("hgMatrix_CGGA.txt",sep="\t")
# expr <- data.matrix(expr)

# We need a list of gene-to-TF mapping from the backbone, binary network
# constructed previously (from ATAC-seq or DNase-seq data etc.)

gene_TFs <- read.table("TFs_for_genes.txt",sep="\t",header=F)
TFList <- unique(c(as.character(gene_TFs[,2])))
TFList <- intersect(TFList,rownames(expr))
tgtGenes <- unique(as.character(gene_TFs[,1]))
tgtGenes <- intersect(tgtGenes, rownames(expr))
options(stringsAsFactors = F)
unionNet.TFs <- sapply(1:nrow(gene_TFs),function(i){
  (as.character(gene_TFs[i,1]) %in% tgtGenes) &
    (as.character(gene_TFs[i,2]) %in% TFList)
})
gene_TFs <- gene_TFs[unionNet.TFs,]
tgtGenes <- unique(as.character(gene_TFs[,1]))
write.table(tgtGenes,"targetGeneList.txt",sep="\t",
            row.names = F,col.names = F,quote = F)

# Lambda controls regularization strength in the nonlinear regression model
# Example: L2 strength is 0.1, so input is 0_0.1
# Currently the regression pipeline does not support L1 regression since
# it is harder to compute the gradient function.
lambdas <- args[2]

# Read in class labels: first column is sample name matching expression matrix
# and second column is subtype/cell type names
labels <- read.table("class_labels_CGGA.txt",sep="\t",header=F)
labels <- as.character(labels[,1])
names(labels) <- colnames(expr)
commonSamples <- colnames(expr)
commonLabels.all <- labels[commonSamples]
commonLabels <- commonLabels.all
commonExpr <- expr[,commonSamples]
expr.curr <- expr
sampleByLabels <- lapply(names(table(commonLabels.all)),function(x){
  commonLabels.all[commonLabels.all==x]
})
names(sampleByLabels) <- names(table(commonLabels.all))

# Generate training and validation subsets for expression data
# or randomly shuffle labels to estimate significance for
# subtype-specific regulation

set.seed(972831654)

currentIter <- c(1)

# The third argument provides the mode of regression: 
# vanilla, CV or shuffle sample labels
regrMod <- args[3]
CV.fold <- 1
rmd <- 0

if (length(grep("CV",regrMod,fixed=T))) {
  rmd <- 1
  CV.fold <- as.numeric(args[4])
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
  save("CV.TVSamples",file="CV.TVSamples.R")
  currentIter <- c(as.numeric(args[5]))
} else if (length(grep("ShuffleLabels",args[3],fixed=T))) {
  rmd <- 2
  nshuffle <- as.numeric(args[4])
  labelList <- lapply(c(1:nshuffle),function(i){
    rn <- sample(commonLabels,length(commonLabels),replace = F)
    temp <- commonLabels
    temp[names(commonLabels)] <- as.character(rn)
    temp
  })
  #labelList <- c(list(commonLabels),labelList)
  currentIter <- c(as.numeric(args[5]))
}

for (j in currentIter) {
  cat(paste0("\n",args[3]," #",as.character(j),"\n"))
  if (rmd == 1) {
    commonSamples <- lapply(CV.TVSamples,function(x){
      x[[j]]$trn
    })
    commonSamples <- unlist(commonSamples)
    commonLabels <- commonLabels.all[names(commonSamples)]
    commonExpr <- expr[,names(commonSamples)]
  } else if (rmd == 2) {
    commonSamples <- labelList[[j]]
    commonLabels <- commonSamples
    commonExpr <- expr[,names(commonSamples)]
  } else {
    commonExpr <- expr[,commonSamples]
  }
  #commonLabels <- commonLabels.all[names(commonSamples)]
  #commonExpr <- expr[,names(commonSamples)]
  #commonmiExpr <- miExpr[,names(commonSamples)]
  expr.curr <- commonExpr
  iteration  <- 0
  maxit <- 1
  while (iteration < maxit) {
    #cat(paste0("\niteration: ",iteration,"\n"))
    cat("\nAssigning files for regression...\n\n")
    if (F) {
      for (type in subtypes) {
        if (rmd == 1 & CV.fold > 1) {
          system(paste0("mkdir ",type))
          system(paste0("mkdir ",type,"/CV_",as.character(j)))
          system(paste0("mkdir ",type,"/CV_",as.character(j),"/genes"))
        } else if (rmd == 2) {
          system(paste0("mkdir ",type))
          system(paste0("mkdir ",type,"/randShfl_",as.character(j)))
          system(paste0("mkdir ",type,"/randShfl_",as.character(j),"/genes"))
        } else {
          system(paste0("mkdir ",type))
          system(paste0("mkdir ",type,"/genes"))
        }
        subtype <- names(commonLabels)[commonLabels==type]
        subtypeExpr <- commonExpr[,subtype]
        # subtypemiExpr <- commonmiExpr[miRNAList,subtype]
        subtypeTFExpr <- expr[TFList,subtype]
        if (rmd == 1 & CV.fold > 1) {
          write.table(subtypeExpr,paste0(type,"/CV_",as.character(j),"/gene_expr.txt"),
                      sep="\t",quote=F)
          write.table(subtypeTFExpr,paste0(type,"/CV_",as.character(j),"/TF_expr.txt"),
                      sep="\t",quote=F)
        } else if (rmd == 2) {
          write.table(subtypeExpr,paste0(type,"/randShfl_",as.character(j),"/gene_expr.txt"),
                      sep="\t",quote=F)
          write.table(subtypeTFExpr,paste0(type,"/randShfl_",as.character(j),"/TF_expr.txt"),
                      sep="\t",quote=F)
        } else {
          write.table(subtypeExpr,paste0(type,"/gene_expr.txt"),
                      sep="\t",quote=F)
          write.table(subtypeTFExpr,paste0(type,"/TF_expr.txt"),
                      sep="\t",quote=F)
        }
        typeLabel <- subtype
        count <- 0
        if (rmd == 1 & CV.fold > 1) {
          for (gene in tgtGenes) {
            count <- count + 1
            fp <- file.path(mainDir=paste(getwd(),"/",type,"/CV_",as.character(j),"/","genes",sep=""),gene)
            if (!dir.exists(fp)) {
              dir.create(fp)
              TFs <- as.character(gene_TFs[gene_TFs[,1]==gene,2])
              TFs <- match(TFs,TFList)
              write.table(TFs,paste(getwd(),"/",type,"/CV_",as.character(j),"/","genes/",gene,"/","TFs.txt",sep=""),
                          sep="\t",row.names=F,col.names=F,quote=F) 
            }
            if (count %% 100 == 0) {
              cat(".")
            }
          } 
        } else if (rmd == 2) {
          for (gene in tgtGenes) {
            count <- count + 1
            fp <- file.path(mainDir=paste(getwd(),"/",type,"/randShfl_",as.character(j),"/genes",sep=""),gene)
            if (!dir.exists(fp)) {
              dir.create(fp)
              TFs <- as.character(gene_TFs[gene_TFs[,1]==gene,2])
              TFs <- match(TFs,TFList)
              write.table(TFs,paste(getwd(),"/",type,"/randShfl_",as.character(j),"/","genes/",gene,"/","TFs.txt",sep=""),
                          sep="\t",row.names=F,col.names=F,quote=F)
            }
            if (count %% 100 == 0) {
              cat(".")
            }
          } 
        } else {
          for (gene in tgtGenes) {
            count <- count + 1
            fp <- file.path(mainDir=paste(getwd(),"/",type,"/","genes",sep=""),gene)
            if (!dir.exists(fp)) {
              dir.create(fp)
              TFs <- as.character(gene_TFs[gene_TFs[,1]==gene,2])
              TFs <- match(TFs,TFList)
              write.table(TFs,paste(getwd(),"/",type,"/","genes/",gene,"/","TFs.txt",sep=""),
                          sep="\t",row.names=F,col.names=F,quote=F) 
            }
            if (count %% 100 == 0) {
              cat(".")
            }
          } 
        }
        count <- 0
        if (rmd == 1 & CV.fold > 1) {
          for (gene in tgtGenes) {
            count <- count + 1
            exprVec <- commonExpr[gene,subtype]
            exprfn <- paste(getwd(),"/",type,"/CV_",as.character(j),"/","genes/",gene,"/","exprs.txt",sep="")
            if (!file.exists(exprfn)) {
              write.table(exprVec,exprfn,sep="\t",row.names=F,col.names=F,quote=F)
            }
            if (count %% 100 == 0) {
              cat(".")
            }
          }
        } else if (rmd == 2) {
          for (gene in tgtGenes) {
            count <- count + 1
            exprVec <- commonExpr[gene,subtype]
            exprfn <- paste(getwd(),"/",type,"/randShfl_",as.character(j),"/","genes/",gene,"/","exprs.txt",sep="")
            if (!file.exists(exprfn)) {
              write.table(exprVec,exprfn,sep="\t",row.names=F,col.names=F,quote=F)
            }
            if (count %% 100 == 0) {
              cat(".")
            }
          }
        } else {
          for (gene in tgtGenes) {
            count <- count + 1
            exprVec <- commonExpr[gene,subtype]
            exprfn <- paste(getwd(),"/",type,"/","genes/",gene,"/","exprs.txt",sep="")
            if (!file.exists(exprfn)) {
              write.table(exprVec,exprfn,sep="\t",row.names=F,col.names=F,quote=F)
            }
            if (count %% 100 == 0) {
              cat(".")
            }
          }
        }
      }
    }
    
    # 2) Run regression on each training set for cross-validation
    cat("\nRunning regression...\n\n")
    if (T) {
      for (l in lambdas) {
        lambda2 <- strsplit(l,"_")[[1]][2]
        temp <- lapply(subtypes,function(s){
          # Here we need to use a Matlab compiler to make matlab executables
          # so that the number of instances submitted is not limited by licensing
          # The script run_batchOpt.sh will handle the compiled matlab script batchOpt, which
          # will need to be reconfigured for the specific datasets used.
          # Refer to this website for compiling Matlab runtimes:
          # https://www.bu.edu/tech/support/research/software-and-programming/common-languages/matlab/standalone/
          cmd <- paste0('./run_batchOpt.sh /broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2017a/ ',
                        s,' ',as.character(j),' ',as.character(rmd), ' ', lambda2)
          # cmd <- paste0('matlab -nodesktop -nosplash -nodisplay -r "batchOpt_l2(\'',
          #               s,'\',\'',as.character(j),'\',',as.character(rmd),',',
          #               lambda2,');quit"')
          system(cmd)
        })
      }
    }
    
    # 3) Collect results into network files
    if (T) {
      cat("\nCollecting regression results...\n")
      for (type in subtypes) {
      	for (l in lambdas) {
	        lambda <- strsplit(l,"_",fixed=T)[[1]]
	        lambda <- paste0(lambda[2],"_",lambda[1])
      	  if (rmd > 0 & CV.fold > 1) {
      	    # The script integrateNetworks_pipeline.py takes the individual outputs from the 
      	    # regression pipeline and assembles an adjacency list with weights (F values)
      	    # for the subtype-/cell type-specific networks.
      	    cmd <- paste0("python integrateNetworks_pipeline.py ",type," ",lambda," ",
      	                  as.character(rmd)," ",as.character(currentIter))
      	  } else {
      	    cmd <- paste0("python integrateNetworks_pipeline.py ",type," ",lambda," ",
      	                  as.character(rmd))
      	  }
          system(cmd)
	      }
      }
    }
    
    # 4) Compute capacities from network files generated above and compile into adjacency matrices
    if (T) {
      cat("\nComputing capacities...\n\n")
      allTFs <- read.delim("TFs_for_genes.txt",header = F, check.names = F)
      allGenes <- unique(as.character(allTFs[,1]))
      allTFs <- unique(c(as.character(allTFs[,2])))
      for (subtype in subtypes) {
        if (rmd == 1 & CV.fold > 1) {
          fpath <- paste0(subtype,"/CV_",as.character(j),"/")
        } else if (rmd == 2) {
          fpath <- paste0(subtype,"/randShfl_",as.character(j),"/")
        } else {
          fpath <- paste0(subtype,"/")
        }
        for (l in lambdas){
	  lambda <- strsplit(l,"_",fixed=T)[[1]]
	  lambda <- paste0(lambda[2],"_",lambda[1])
          tfm <- read.delim(paste0(fpath,subtype,"_TF_gene_ri_nomi_",lambda,".txt"),sep="\t",
                            header=F,check.names = F)
          tfm.tgt <- unique(c(as.character(tfm[,2])))
	  tfm.tgt <- gsub("_g","",tfm.tgt,fixed=T)
          tfm.tgt <- unique(c(tfm.tgt,tgtGenes))
          tfm.reg <- unique(c(as.character(tfm[,1])))
	  tfm.reg <- gsub("_t","",tfm.reg,fixed=T)
          tfm.reg <- unique(c(tfm.reg,TFList))
          tfm.mat <- matrix(1,length(tfm.reg),length(tfm.tgt))
          rownames(tfm.mat) <- tfm.reg
          colnames(tfm.mat) <- tfm.tgt
          for (i in c(1:nrow(tfm))) {
            reg <- gsub("_t","",tfm[i,1],fixed=T)
            tgt <- gsub("_g","",tfm[i,2],fixed=T)
            tfm.mat[reg,tgt] <- as.numeric(as.character(tfm[i,3]))
          }
          if (rmd == 1 & CV.fold > 1) {
            write.table(tfm.mat,paste0(subtype,"_TFGeneCapacities_ri_nomi_",paste0("_lambda",lambda),
                                       "_CV_",as.character(j),".txt"),
                        sep="\t",col.names=T,row.names=T,quote=F)
          } else if (rmd == 2) {
            write.table(tfm.mat,paste0(subtype,"_TFGeneCapacities_ri_nomi_",paste0("_lambda",lambda),
                                       "_randomShuffle_",as.character(j),".txt"),
                        sep="\t",col.names=T,row.names=T,quote=F)
          } else {
            write.table(tfm.mat,paste0(subtype,"_TFGeneCapacities_ri_nomi_",paste0("_lambda",lambda),
                                       ".txt"),
                        sep="\t",col.names=T,row.names=T,quote=F)
          }
          cat(".")
        }
        cat("\n")
      }
    }
    # The iteration counter is for future purposes only
    iteration <- iteration + 1
  }
}












