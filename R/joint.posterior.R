#' jointPosterior
#'
#' Function to obtain the normalized joint posterior of the data and partition. 
#' 
#' @details Calculates the normalized joint posterior of the data and partition
#'  under the Product Partition Model formulation
#'   of the Dirichlet Process Mixture model.
#' 
#' @param y Numeric data vector for one gene (log-transformed non-zeroes)
#' 
#' @param mcobj Object returned by \code{\link{mclustRestricted}}
#' 
#' @param alpha Value for the Dirichlet concentration parameter
#' 
#' @param m0 Prior mean value for generating distribution of cluster means
#' 
#' @param s0 Prior precision value for generating distribution of cluster means
#' 
#' @param a0 Prior shape parameter value for the generating distribution of 
#' cluster precision
#' 
#' @param b0 Prior scale parameter value for the generating distribution of 
#' cluster precision
#' 
#' @references Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, 
#' Kendziorski C. A statistical approach for identifying differential 
#' distributions
#' in single-cell RNA-seq experiments. Genome Biology. 2016 Oct 25;17(1):222. 
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-
#' 1077-y}
#'
#' @return log joint posterior value


jointPosterior <- function(y, mcobj, alpha, m0, s0, a0, b0){
  
  r <- length(unique(mcobj$class))
  nk <- NULL
  sumyk <- NULL
  sumyk2 <- NULL
  for (k in 1:r){ 
    nk <- c(nk, sum(mcobj$class==k)) 
    sumyk <- c(sumyk, sum(y[mcobj$class==k]))
    sumyk2 <- c(sumyk2, sum(y[mcobj$class==k]^2))
  }
  
  sk <- s0 + nk
  mk <- (s0*m0 + sumyk) / sk
  ak <- a0 + nk
  bk <- b0 + sumyk2 + s0*m0^2 - sk*mk^2
  J <- sum(nk)
  
  
  # including prior for partition	
  logpost <- 0
  for (k in 1:r){ 
    logpost <- logpost + lgamma(nk[k]) + lgamma(ak[k]/2) - 
      (ak[k]/2)*log(bk[k]/2) - 0.5*log(sk[k]) 
  }
  logpost <- logpost + r*log(alpha) + lgamma(alpha) - J*lgamma(a0/2) + 
    (J*a0/2)*log(b0/2) + J*0.5*log(s0) - (J/2)*log(2*pi)- lgamma(alpha + J)
  
  return(logpost)
}


