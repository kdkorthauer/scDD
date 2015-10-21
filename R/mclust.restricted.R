#' mclustRestricted
#'
#' Function to determine how many normal mixture components are present.
#' 
#' @details Robust to detecting multiple components that are close together by enforcing that the distance between two clusters 
#'  of appreciable size (at least 4 samples) be at least 2 standard deviations and not have variances that differ by more than 
#'  a ratio of 3 for clusters of small size (10% or less of the sample). 
#' 
#' @param y Numeric vector of values to fit to a normal mixture model with Mclust.
#' 
#' @param restrict Logical indicating whether or not to enforce the restriction on cluster separation by at least 2 standard deviations
#'   and ratio of largest to smallest variance less than 3 when small clusters are present (containing 10% or less of the sample).
#'   If False, then Mclust results as is are returned.
#'   
#'  
#' @importFrom mclust Mclust
#' 
#' @return List object with (1) vector of cluster membership, (2) cluster means, (3) cluster variances, (4) number of model parameters,
#'  (5) sample size, (6) BIC of selected model, and (6) loglikelihood of selected model. 



mclustRestricted <- function(y, restrict=TRUE){
  mc <- Mclust(y, warn=FALSE, modelNames=c("V"), G=1:5)	
  cl <- mc$classification
  comps <- mc$G
  
  if(comps > length(unique(cl))){
    mc <- Mclust(y, warn=FALSE, modelNames=c("V"), G=1:(comps-1))	
    cl <- mc$classification
    comps <- mc$G
  }
  
  if (restrict) { 
    if (comps > 1){
      compmeans <- as.numeric(by(y, cl, mean))
      meandiff <- diff(compmeans)/max(sqrt(mc$parameters$variance$sigmasq))
      if( min(diff(compmeans)) < 1 ) {
        meandiff <- diff(compmeans)/sd(y)
      }
      
      max.cl <- which.max(mc$parameters$variance$sigmasq)
      min.cl <- which.min(mc$parameters$variance$sigmasq)
      vardiff <- mc$parameters$variance$sigmasq[max.cl] / mc$parameters$variance$sigmasq[min.cl]
      nmin <- table(cl)[min.cl]
      
      tries <- 0	
      
      if (length(vardiff)==0){
        vardiff <- 1
      }
      if (length(nmin)==0){
        nmin <- Inf
      }
      if (length(meandiff)==0){
        meandiff <- Inf
      }
      
      # enforce that clusters have to be separated by at least 2 sds 
      while( (min(meandiff) < 2 | (vardiff > 3 & nmin < 0.10*length(y))) & min(table(cl))>2 & tries <=4 & sum(is.na(mc$class))==0){
        comps_old <- comps
        comps <- comps-1 
        mc <- Mclust(y, warn=FALSE, modelNames=c("V"), G=1:comps)
        cl <- mc$classification
        comps <- mc$G
        
        tries <- tries + 1
        if (comps > 1 & comps-comps_old==1){
          compmeans <- as.numeric(by(y, cl, mean))
          meandiff <- diff(compmeans)/max(sqrt(mc$parameters$variance$sigmasq))
         
          max.cl <- which.max(mc$parameters$variance$sigmasq)
          min.cl <- which.min(mc$parameters$variance$sigmasq)
          vardiff <- mc$parameters$variance$sigmasq[max.cl] / mc$parameters$variance$sigmasq[min.cl]
          nmin <- table(cl)[min.cl]
        }else{
          meandiff <- 10
          bicdiff <- -Inf
          nmin <- 100
        }
      }
    }
  }
  
  # enforce clusters to have increasing means
  cl2 <- cl
  clust.order <- 1
  if (comps > 1){
    clust.order <- order(mc$parameters$mean)
    nms.clust <- names(mc$parameters$mean[clust.order])
    for (i in 1:comps){
      cl2[cl==i] <- as.numeric(nms.clust[i])
    }
  }
  return(list(class=cl2, mean=mc$parameters$mean[clust.order], var=mc$parameters$variance$sigmasq[clust.order], 
              df=mc$df, n=mc$n, bic=mc$bic, loglik=mc$loglik, model=mc$modelName))
}
