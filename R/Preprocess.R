#' preprocess
#'
#' Function to preprocess SingleCellExperiment object (1) to
#' only keep genes with a certain number of nonzero entries,
#' and (2) optionally apply a normalization procedure.
#'
#' @param SCdat An object of class \code{SingleCellExperiment} that contains 
#'  single-cell expression and metadata. The \code{assays} 
#'   slot contains a named list of matrices, where the normalized counts are 
#'   housed in the one named \code{normcounts}, and unnormalized counts are 
#'   stored in the one names \code{counts}. If either \code{scran_norm} or 
#'   \code{median_norm} is set to TRUE, the \code{normcounts} slot will be 
#'   created from the \code{counts} slot. The counts and normalized counts 
#'   matrices should have one
#'    row for each gene and one sample for each column.  
#'   The \code{colData} slot should contain a data.frame with one row per 
#'   sample and columns that contain metadata for each sample.  This data.frame
#'   should contain a variable that represents biological condition, which is 
#'   in the form of numeric values (either 1 or 2) that indicates which 
#'   condition each sample belongs to (in the same order as the columns of 
#'   \code{normcounts}).  Optional additional metadata about each cell can also
#'   be contained in this data.frame, and additional information about the 
#'   experiment can be contained in the \code{metadata} slot as a list.
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
#' @param condition A character object that contains the name of the column in 
#' \code{colData} that represents 
#'  the biological group or condition of interest (e.g. treatment versus 
#'  control).  Note that this variable should only contain two 
#'  possible values since \code{scDD} can currently only handle two-group 
#'  comparisons.  The default option assumes that there
#'  is a column named "condition" that contains this variable. 
#' 
#' @export
#' 
#' @importFrom EBSeq MedianNorm
#' 
#' @importFrom EBSeq GetNormalizedMat
#' 
#' @importFrom scran computeSumFactors
#' 
#' @import SummarizedExperiment
#'
#' @import SingleCellExperiment
#' 
#'
#' @references Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, 
#' Kendziorski C. A statistical approach for identifying differential 
#' distributions
#' in single-cell RNA-seq experiments. Genome Biology. 2016 Oct 25;17(1):222.
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-
#' 1077-y}
#'
#' @return An object of class \code{SingleCellExperiment} with genes removed if
#' they have more than \code{zero.thresh} zeroes, and the \code{normcounts} 
#' assay added if either \code{scran_norm} or \code{median_norm} is set to TRUE
#' and only \code{counts} is provided. If \code{normcounts} already exists and
#' either \code{scran_norm} or \code{median_norm} is set to TRUE, then the new
#' normalized counts are placed in the \code{normcounts} assay slot, and the
#' original values are moved to a new slot called \code{normcounts-orig}.
#'  
#' @examples 
#'  
#'  # load toy example SingleCellExperiment object
#'  
#'  data(scDatEx)
#'  
#'  # apply the preprocess function to filter out genes if they have more than
#'  # 75% zero
#'  
#'  scDatEx <- preprocess(scDatEx, zero.thresh=0.75)
#'  
#'  # apply the preprocess function again, but this time threshold on the 
#'  # proportion of zeroes and apply scran normalization
#'  # set the zero.thresh argument to 0.75 so that genes with more than 75% 
#'  # zeroes are filtered out 
#'  # set the scran_norm argument to TRUE to return scran normalized counts
#'  
#'  scDatEx.scran <- preprocess(scDatEx, zero.thresh=0.75, scran_norm=TRUE)
#'  
#'  # set the median_norm argument to TRUE to return Median normalized counts
#'  
#'  scDatEx.median <- preprocess(scDatEx, zero.thresh=0.75, median_norm=TRUE)

preprocess <- function(SCdat,
                       condition="condition",
                       zero.thresh=0.9,
                       scran_norm=FALSE,
                       median_norm=FALSE){
  
  # check whether SCdat is a member of the SingleCellExperiment class
  if(!("SingleCellExperiment" %in% class(SCdat))){
    stop("Please provide a valid 'SingleCellExperiment' object.")
  }
  
  # check that condition inputs are valid
  if (length(unique(colData(SCdat)[[condition]])) != 2 | 
      length(colData(SCdat)[[condition]]) != dim(SCdat)[2]){
    stop("Error: Please specify valid condition labels.")
  }
  
  if(median_norm & scran_norm){
    stop(paste0("Specified conflicting options for normalization. Only one ",
                " of scran_norm and median_norm can be set to TRUE"))
  }
  
  if (sum(median_norm, scran_norm) == 1 &
        !("counts" %in% assayNames(SCdat)) ){
    stop(paste0("If median or scran norm is specified, the input object ",
                "must contain raw counts in 'counts'."))
  }
  
  if (sum(median_norm, scran_norm) == 1 &
       "normcounts" %in% assayNames(SCdat) ){
    warning(paste0("median or scran norm is specified and the 'normcounts' ",
                "assay already exists; replacing 'normcounts' in output ",
                "with the specified normalization method. Original contents ",
                "of 'normcounts' are now in 'normcounts-orig'."))
    assay(SCdat, "normcounts-orig") <- normcounts(SCdat)
  }
    
  if (sum(median_norm, scran_norm) == 0 &
     !("normcounts" %in% assayNames(SCdat))){
    stop(paste0("If median or scran norm is not specified, the input object ",
                " must contain normalized counts in 'normcounts'"))
  }
  
  
  if (sum(median_norm, scran_norm) == 0){
    pe_mat <- normcounts(SCdat)
  }else{
    pe_mat <- assay(SCdat, "counts")
  }
  
  # reference category/condition - the first listed one
  ref <- unique(colData(SCdat)[[condition]])[1]
  cond <- colData(SCdat)[[condition]]
 
  # threshold on zero percentage (within each condition)
  pz1 <- apply(pe_mat[,cond==ref], 1, function(x) sum(x==0)/length(x))
  pz2 <- apply(pe_mat[,cond!=ref], 1, function(x) sum(x==0)/length(x))
  SCdat <- SCdat[pz1 <= zero.thresh & pz2 <= zero.thresh,]
  
  # normalize
  if (median_norm){ 
    message("Performing Median Normalization")
    libsize <- MedianNorm(assay(SCdat, "counts"))
    normcounts(SCdat) <- GetNormalizedMat(assay(SCdat, "counts"),
                                          libsize)
  }else if(scran_norm){
    message("Performing scran Normalization")
    SCdat <-  computeSumFactors(SCdat, assay.type="counts")
    normcounts(SCdat) <- GetNormalizedMat(assay(SCdat, "counts"),
                                          sizeFactors(SCdat))
  }
  
  return(SCdat)
}
