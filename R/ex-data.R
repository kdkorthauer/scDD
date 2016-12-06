#' Toy example data list 
#'
#' Toy example data list (one item for each of two conditions) for 100 genes to illustrate how to use the function
#'  \code{\link{preprocess}}.  
#'
#' @name scDatExList
#'
#' @docType data
#'
#' @usage data(scDatExList)
#'
#' @format A list of two matrices (one for each of two conditions) labeled "C1" and "C2".  Each matrix contains data for 100 genes
#'   and a variable number of samples (78 in C1 and 64 in C2).
#'   
#' @references Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, Kendziorski C. A statistical approach for identifying differential distributions
#' in single-cell RNA-seq experiments. Genome Biology. 2016 Oct 25;17(1):222. \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1077-y}
#'
#' @keywords datasets
#'
#'
NULL


#' Toy example data
#'
#' Toy example data in \code{\link[Biobase]{ExpressionSet}} format for 500 genes to illustrate how to generate simulated data from example data
#'  using \code{\link{simulateSet}}.  
#'
#' @name scDatEx
#'
#' @docType data
#'
#' @usage data(scDatEx)
#'
#' @format An object of class \code{\link[Biobase]{ExpressionSet}} containing data for 500 genes for 142 samples
#'  (78 from condition 1 and 64 from condition 2).  Condition labels (1 or 2) are stored in the phenoData slot. 
#'   
#' @references Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, Kendziorski C. A statistical approach for identifying differential distributions
#' in single-cell RNA-seq experiments. Genome Biology. 2016 Oct 25;17(1):222. \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1077-y}   
#'
#' @keywords datasets
#'
#'
NULL


#' Toy example of simulated data
#'
#' Toy example data in \code{\link[Biobase]{ExpressionSet}} format for 500 genes to illustrate how to generate simulated data from example data
#'  using \code{\link{simulateSet}}.  Contains 5 genes from each category (DE, DP, DM, DB, EE, and EP).
#'
#' @name scDatExSim
#'
#' @docType data
#'
#' @usage data(scDatExSim)
#'
#' @format An object of class \code{\link[Biobase]{ExpressionSet}} containing data for 30 genes for 200 samples
#'  (100 from condition 1 and 100 from condition 2).  Condition labels (1 or 2) are stored in the phenoData slot.  
#'  Row names of the assayData slot contain the two letter category label that the gene was simulated from (e.g.  'EE', 'DB', ...)
#'  along with the row number (1-30).
#'   
#' @references Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, Kendziorski C. A statistical approach for identifying differential distributions
#' in single-cell RNA-seq experiments. Genome Biology. 2016 Oct 25;17(1):222. \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1077-y}
#'
#' @keywords datasets
#'
#'
NULL
