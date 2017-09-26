#' results
#' 
#' extract results objects after running scDD analysis
#' 
#' @details Convenient helper function to extract the results (gene 
#' classifications, pvalues, and clustering information). Results 
#' data.frames/matrices are stored in the 
#' \code{metadata} slot and can also be accessed without the help of this 
#' convenience function by calling \code{metadata(SCdat)}.
#'   
#' @param SCdat An object of class \code{SingleCellExperiment} that contains 
#' normalized single-cell expression and metadata, and the output of the
#' \code{scDD} function.
#' 
#' @param type A character variable specifying which output is desired, 
#'  with possible values "Genes", "Zhat.c1", "Zhat.c2", and 
#'  "Zhat.overall".  The default value is "Genes", which contains a 
#'  a data frame with nine columns: 
#' gene name (matches rownames of SCdat), permutation p-value for testing of 
#' independence of 
#'  condition membership with clustering, Benjamini-Hochberg adjusted version 
#'  of the previous column, p-value for test of difference in dropout rate
#'   (only for non-DD genes), 
#'  Benjamini-Hochberg adjusted version of the previous column, name of the 
#'  DD (DE, DP, DM, DB) pattern or DZ (otherwise NS = not significant), the 
#'  number of clusters identified overall, the number of clusters identified in 
#'  condition 1 alone, and the number of clusters identified in condition 
#'  2 alone.
#'  
#'  If \code{type} is "Zhat.c1", then a \code{matrix} is returned
#'  that contains the fitted cluster memberships (partition estimates Z)
#'  for each sample (cluster number given by 1,2,3,...) in columns and
#'  gene in rows only for condition 1.  The same information is returned 
#'  only for condition 2, and for the overall clustering, when \code{type}
#'  is set to "Zhat.c2" or "Zhat.overall", respectively.  
#'   Zeroes, which are not involved in the clustering, are
#'   labeled as zero.    
#'   
#' @return A \code{data.frame} which contains either the gene classification 
#'  and p-value results, or cluster membership information, as detailed in the
#'  description of the \code{type} input parameter.
#'  
#' @export
#'  
#' @examples 
#'  
#' # load toy simulated example SingleCellExperiment object to find DD genes  
#' data(scDatExSim)
#' 
#' # set arguments to pass to scDD function
#' 
#' prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
#' 
#' # call the scDD function to perform permutations and classify DD genes
#' 
#' scDatExSim <- scDD(scDatExSim, prior_param=prior_param, testZeroes=FALSE)
#' 
#' # extract main results object
#' 
#' RES <- results(scDatExSim)

results <- function(SCdat, type=c("Genes", "Zhat.c1" , "Zhat.c2", 
                                  "Zhat.combined")){
  type <- match.arg(type)
  return(metadata(SCdat)[[type]])
}