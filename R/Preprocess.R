#' preprocess
#'
#' Function to preprocess a list of data matrices to (1) combine 
#' into one matrix and (2)
#' only keep genes with a certain number of nonzero entries.
#'
#' 
#' @param DataList A list object where each item contains a matrix of 
#'  data with rows 
#'  designating genes and samples designating columns.  
#'  The name of the list objects represents
#'  their condition.
#'
#' @param ConditionNames Character vector of length 1 or 2 
#' which contains the name(s) of the
#'  conditions (the names of the items in \code{DataList}) to be processed.
#' 
#' @param zero.thresh A numeric value between 0 and 1 that represents
#'  the maximum proportion
#'  of zeroes per gene allowable in the processed dataset
#'  
#' @param scran_norm Logical indicating whether or not to normalize the data
#'  using scran Normalization from \code{scran}
#'  
#' @param median_norm Logical indicating whether or not to normalize the data
#'  using Median Normalization from \code{EBSeq}
#' 
#' @export
#' 
#' @importFrom EBSeq MedianNorm
#' 
#' @importFrom EBSeq GetNormalizedMat
#' 
#' @importFrom scran computeSumFactors
#'
#' @references Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, 
#' Kendziorski C. A statistical approach for identifying differential 
#' distributions
#' in single-cell RNA-seq experiments. Genome Biology. 2016 Oct 25;17(1):222.
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-
#' 1077-y}
#'
#' @return pe_mat Processed matrix with genes in rows and samples in columns. 
#'  Column names 
#'  indicate condition.
#'  
#' @examples 
#'  
#'  # load toy example data list
#'  
#'  data(scDatExList)
#'  
#'  
#'  # check that the data is formated as a list of 2 matrices 
#'  # (one for each of 2 conditions), 
#'  # that each matrix has 100 rows (one for each gene), 
#'  # and that the number of columns in 
#'  # each matrix corresponds to the number of samples in 
#'  # each condition (78 and 64, respectively)
#'  
#'  str(scDatExList)
#'  
#'  
#'  # get the names of the conditions to pass to the preprocess function
#'  
#'  condition.names <- names(scDatExList)
#'  
#'  
#'  # apply the preprocess function to reformat the data into one 
#'  # data matrix with 100 rows and 78+64=142 columns
#'  # set the zero.thresh argument to 1 so that genes are filtered out if they
#'  # are all zero set the median_norm argument to FALSE to return raw data
#'  
#'  scDatExMat <- preprocess(scDatExList, ConditionNames=condition.names, 
#'                           zero.thresh=1)
#'  
#'  
#'  # apply the preprocess function again, but this time threshold on the 
#'  # proportion of zeroes and apply scran normalization
#'  # set the zero.thresh argument to 0.75 so that genes with more than 75% 
#'  # zeroes are filtered out 
#'  # set the scran_norm argument to TRUE to return scran normalized counts
#'  
#'  normDat.scran <- preprocess(scDatExList, 
#'                                     ConditionNames=condition.names, 
#'                                     zero.thresh=0.90, scran_norm=TRUE)
#'  
#'  # set the median_norm argument to TRUE to return Median normalized counts
#'  
#'  normDat.median <- preprocess(scDatExList, 
#'                                     ConditionNames=condition.names, 
#'                                     zero.thresh=0.90, median_norm=TRUE)

preprocess <- function(DataList, ConditionNames, 
                       zero.thresh=0.9,
                       scran_norm=FALSE,
                       median_norm=FALSE){
  if(!is.list(DataList)){
    stop("Input data must be a list of data matrices, 
         where each item is a different condition")
  }
  
  if(median_norm & scran_norm){
    stop(paste0("Specified conflicting options for normalization. Only one ",
                " of scran_norm and median_norm can be set to TRUE"))
  }
  
  if(length(ConditionNames)==1){
    these <- which(names(DataList) %in% ConditionNames)
    pe_mat <- DataList[[these]]
    condition <- rep(1, ncol(pe_mat))
    
    # threshold on zero percentage
    pz1 <- apply(pe_mat, 1, function(x) sum(x==0)/length(x))
    pe_mat <- pe_mat[pz1 <= zero.thresh ,]
    
    # median normalize
    if (median_norm){ 
      message("Performing Median Normalization...")
      libsize <- MedianNorm(pe_mat)
      pe_mat <- GetNormalizedMat(pe_mat, libsize)
    }else if(scran_norm){
      message("Performing scran Normalization...")
      libsize <-  computeSumFactors(pe_mat)
      pe_mat <- GetNormalizedMat(pe_mat, libsize)
    }
    
  }else if(length(ConditionNames)==2){
    these <- which(names(DataList) %in% ConditionNames)
    pe_mat <- list(DataList[[these[1]]], DataList[[these[2]]])
    condition <- c(rep(1, ncol(pe_mat[[1]])), rep(2, ncol(pe_mat[[2]])))
    
    # threshold on zero percentage
    pz1 <- apply(pe_mat[[1]], 1, function(x) sum(x==0)/length(x))
    pz2 <- apply(pe_mat[[2]], 1, function(x) sum(x==0)/length(x))
    pe_mat <- cbind(pe_mat[[1]], pe_mat[[2]])
    pe_mat <- pe_mat[pz1 <= zero.thresh & pz2 <= zero.thresh,]
    
    # normalize
    if (median_norm){ 
      message("Performing Median Normalization")
      libsize <- MedianNorm(pe_mat)
      pe_mat <- GetNormalizedMat(pe_mat, libsize)
    }else if(scran_norm){
      message("Performing scran Normalization")
      libsize <-  computeSumFactors(pe_mat)
      pe_mat <- GetNormalizedMat(pe_mat, libsize)
    }
    
  }else{ stop("Only 1 or 2 conditions supported")}
  return(pe_mat)
  }
