#' testKS
#'
#' Function to perform KS test
#'
#' @inheritParams singleCellSimu 
#' 
#' @param dat Matrix of single-cell RNA-seq data with genes in rows and samples 
#' in columns.
#'
#' @param condition Vector containing the indicator of which condition each 
#' sample 
#'  (in the columns of \code{dat}) belongs to.
#' 
#' @param DEIndex Vector containing the row numbers of the DE genes
#'
#' @param inclZero Logical indicating whether to include zero in the test of 
#' different distributions
#'
#' @export 
#' 
#' @references Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, 
#' Kendziorski C. A statistical approach for identifying differential 
#' distributions
#' in single-cell RNA-seq experiments. Genome Biology. 2016 Oct 25;17(1):222. 
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-
#' 1077-y}
#'
#' @return List object containing the significant gene indices, their adjusted 
#' p-values, and (if DE genes are supplied)
#' the power and fdr.
#'   
#' @examples 
#' 
#' # load toy simulated example ExpressionSet to find KS genes
#' 
#' data(scDatExSim)
#' 
#' 
#' # load SingleCellExperiment package to facilitate subset operations
#' 
#' library(SingleCellExperiment)
#' 
#' 
#' # check that this object is a member of the ExpressionSet class
#' # and that it contains 200 samples and 30 genes
#' 
#' class(scDatExSim)
#' show(scDatExSim)
#' 
#' # perform KS test and obtain adjusted p-values
#' RES_KS <- testKS(normcounts(scDatExSim), scDatExSim$condition, inclZero=FALSE,
#'                  numDE=20, DEIndex=1:20)


testKS <- function(dat, condition, inclZero=TRUE, numDE=NULL, DEIndex){
  onegene <- function(x){
    x1 <- x[condition==unique(condition)[1]]
    x2 <- x[condition==unique(condition)[2]]
    
    if (!inclZero){
      x1 <- (x1[x1>0])
      x2 <- (x2[x2>0])
    }
    suppressWarnings(ks.test(x1,x2, exact=FALSE))$p.value
  }
  
  ks.pval.unadj <- apply(dat, 1, function(x) onegene(x) )
  ks.pval <- p.adjust(ks.pval.unadj, method="BH")
  sig_genes_ks <- which(ks.pval < 0.05)
  names(sig_genes_ks) <- rownames(dat)[sig_genes_ks]
  names(ks.pval) <- rownames(dat)
  
  if(!is.null(numDE)){
    power <- sum(sig_genes_ks %in% DEIndex) / numDE
    fdr <- length(sig_genes_ks[!(sig_genes_ks %in% DEIndex)]) / 
      length(sig_genes_ks)
    return(list(genes=sig_genes_ks, p=ks.pval, p.unadj=ks.pval.unadj, 
                power=power, fdr=fdr))
  }
  return(list(genes=sig_genes_ks, p=ks.pval, p.unadj=ks.pval.unadj))
}
