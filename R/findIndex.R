#' findIndex
#'
#' Find a reasonable set of genes (one mode and at least 25% nonzero values) to use for simulation.
#' @inheritParams scDD
#' 
#'  @importClassesFrom Biobase ExpressionSet
#'  
#'  @importMethodsFrom Biobase exprs
#'  
#'  @importMethodsFrom Biobase featureNames
#'  
#'  @importMethodsFrom Biobase sampleNames
#'
#' @importFrom mclust Mclust mclustBIC
#' 
#' @importFrom BiocParallel bplapply
#'
#' @return Vector of indices for a reasonable set of genes that can be used for simulation.


findIndex <- function(SCdat){
  
# separate ExpressionSet by condition  
Dataset1 <- SCdat[,SCdat$condition==1]
Dataset2 <- SCdat[,SCdat$condition==2]
  
### Find zero-percent
zeropercent <- matrix(data=0, nrow=dim(Dataset1)[1], ncol=2)
rownames(zeropercent) <- featureNames(Dataset1)
zeropercent[,1] <- apply(exprs(Dataset1), 1, function(x) sum(x==0)/length(x))
zeropercent[,2] <- apply(exprs(Dataset2), 1, function(x) sum(x==0)/length(x))

### log-transformation before Mclust
logdata1 <- t(apply(exprs(Dataset1), 1, log))
logdata2 <- t(apply(exprs(Dataset2), 1, log))
logdata1[logdata1==-Inf] <- 0
logdata2[logdata2==-Inf] <- 0
rownames(logdata1) <- rownames(logdata2) <- featureNames(Dataset1)

### Select the genes with less than 75% expression values are 0 and with 1 mode
nonzeroindex <- which((zeropercent[,1]<0.75)&(zeropercent[,2]<0.75))
NumofMode <- matrix(data=0,nrow=nrow(exprs(Dataset1[nonzeroindex,])),ncol=2)
rownames(NumofMode) <- featureNames(Dataset1)[nonzeroindex]

a <- logdata1[nonzeroindex, ]
b <- logdata2[nonzeroindex, ]

mcmat <- function(ROW, MAT){
  y <- MAT[ROW,]
  mclust::Mclust(y[y!=0])$G
}
NumofMode[,1] <- unlist(BiocParallel::bplapply(1:nrow(a), mcmat, MAT=a))
NumofMode[,2] <- unlist(BiocParallel::bplapply(1:nrow(b), mcmat, MAT=b))
  
onemodeindex <- which(NumofMode[,2]==1 & NumofMode[,1]==1)
index <- nonzeroindex[onemodeindex]
return(index)
}
