
#' luOutlier
#'
#' Count the number of clusters with at least \code{min.size} samples
#' 
#' @details Function to obtain a count of the number of clusters that is robust to outliers.  Requires at least \code{min.size} samples to be considered
#'  in the robust count.
#' 
#' @param x Numeric vector of cluster membership (1st item (named \code{class}) in list returned by \code{\link{mclustRestricted}})
#' 
#' @param min.size Numeric value for the minimum number of samples a cluster must have to be considered in the robust count.  Default is 3.
#'   
#' @references Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, Kendziorski C. A statistical approach for identifying differential distributions
#' in single-cell RNA-seq experiments. Genome Biology. 2016 Oct 25;17(1):222. \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1077-y}
#'   
#' @return The robust count of the number of unique clusters excluding those with less than \code{min.size} samples.

luOutlier <- function(x, min.size=3){
  return(sum(table(x)>=min.size))
}


#' lu
#'
#' Shortcut for \code{length(unique())}
#' 
#' @details Function to obtain a count of the number of clusters 
#' 
#' @param x Numeric vector of cluster membership (1st item (named \code{class}) in list returned by \code{\link{mclustRestricted}})
#'
#' @references Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, Kendziorski C. A statistical approach for identifying differential distributions
#' in single-cell RNA-seq experiments. Genome Biology. 2016 Oct 25;17(1):222. \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1077-y}
#'
#' @return The count of the number of unique clusters.

lu <- function(x){
  return(length(unique(x)))
}


#' findOutliers
#'
#' Find the clusters that are considered outliers
#' 
#' @details Function to obtain a count of the number of clusters that is robust to outliers.  Requires at least \code{min.size} samples to be considered
#'  in the robust count.
#' 
#' @param clustering Numeric vector of cluster membership (1st item (named \code{class}) in list returned by \code{\link{mclustRestricted}})
#' 
#' @param min.size Numeric value for the minimum number of samples a cluster must have to be considered in the robust count.  Default is 3.
#'   
#' @references Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, Kendziorski C. A statistical approach for identifying differential distributions
#' in single-cell RNA-seq experiments. Genome Biology. 2016 Oct 25;17(1):222. \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1077-y}
#'
#' @return The robust count of the number of unique clusters excluding those with less than \code{min.size} samples.

findOutliers <- function(clustering, min.size=3){
  good <- names(table(clustering))[table(clustering)>=min.size]
  return(as.numeric(good))
}
