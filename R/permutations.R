#' permMclustCov
#'
#' Function to  obtain bayes factor for permutations of one gene's residuals
#' 
#' @details Obtains bayes factor numerator for data vector \code{y} 
#' representing one gene 
#' 
#' @param y Numeric data vector for one gene
#' 
#' @inheritParams jointPosterior 
#' 
#' @inheritParams scDD
#' 
#' @param nperms Number of permutations of residuals to evaulate
#' 
#' @param condition Vector of condition indicators for each sample
#' 
#' @param C Matrix of confounder variables, where there is one row for each
#'  sample and one column for each covariate.
#' 
#' @param remove.zeroes Logical indicating whether zeroes need to be removed 
#' from \code{y}
#' 
#' @param log.transf Logical indicating whether the data is in the raw scale 
#' (if so, will be log-transformed)
#' 
#' @param restrict Logical indicating whether to perform restricted Mclust 
#' clustering where close-together clusters are joined.
#' 
#' @param ref one of two possible values in condition; represents the 
#' referent category.
#' 
#' @importFrom BiocParallel bplapply
#' 
#' @references Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, 
#' Kendziorski C. A statistical approach for identifying differential 
#' distributions
#' in single-cell RNA-seq experiments. Genome Biology. 2016 Oct 25;17(1):222. 
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-
#' 1077-y}
#'
#' @return Bayes factor numerator for the current permutation

permMclustCov <- function(y, nperms, C, condition, remove.zeroes=TRUE, 
                          log.transf=TRUE, restrict=FALSE, 
                          alpha, m0, s0, a0, b0, ref, min.size){
  orig.y <- y
  
  if(remove.zeroes & log.transf){
    C <- C[y>0]
    cond <- condition[y>0]
    y <- log(y[y>0])
  }else if(log.transf){
    C <- C[y>0]
    cond <- condition
    y <- log(y)
  }else if(remove.zeroes){
    C <- C[y>0]
    cond <- condition[y>0]
    y <- y[y>0]
  }
  
  # fit reduced model
  X <- model.matrix(~ C)
  model <- lm(y ~ C)
  params <- summary(model)$coef[,1]
  
  
  getPerm <- function(model, X, params, orig.y, condition, restrict, 
                      remove.zeroes, log.transf){
    new.resid <- sample(residuals(model), replace=FALSE)
    new.y0 <- as.vector(exp(X %*% params + new.resid))
    new.y <- orig.y
    new.y[new.y > 0] <- new.y0
    
    y1 <- new.y[condition==ref]
    y2 <- new.y[condition!=ref]
    if(remove.zeroes & log.transf){
      y1 <- log(y1[y1>0])
      y2 <- log(y2[y2>0])
    }else if(log.transf){
      y1 <- log(y1)
      y2 <- log(y2)
    }else if(remove.zeroes){
      y1 <- y1[y1>0]
      y2 <- y2[y2>0]
    }
    
    oa <- mclustRestricted(c(y1,y2), restrict=restrict, min.size=min.size)
    c1 <- mclustRestricted(y1, restrict=restrict, min.size=min.size)
    c2 <- mclustRestricted(y2, restrict=restrict, min.size=min.size)
    
    bf.p <- jointPosterior(y1, c1, alpha, m0, s0, a0, b0) + 
      jointPosterior(y2, c2, alpha, m0, s0, a0, b0) -
      jointPosterior(c(y1,y2), oa, alpha, m0, s0, a0, b0)
    return(bf.p)
  }

  bf.p <- unlist(bplapply(1:nperms, function(x) 
    getPerm(model, X, params, orig.y, condition, 
            restrict, remove.zeroes, log.transf)))
  
  return(bf.p)
}


#' permMclust
#'
#' Function to  obtain bayes factor numerator for permutations of one gene
#' 
#' @details Obtains bayes factor numerator for data vector
#'  \code{y} representing one gene 
#' 
#' @inheritParams jointPosterior
#' 
#' @inheritParams scDD
#' 
#' @param y Numeric data vector for one gene
#' 
#' @param nperms Number of permutations of residuals to evaulate
#' 
#' @param condition Vector of condition indicators for each sample
#' 
#' @param remove.zeroes Logical indicating whether 
#' zeroes need to be removed from \code{y}
#' 
#' @param log.transf Logical indicating whether the data is 
#' in the raw scale (if so, will be log-transformed)
#' 
#' @param ref one of two possible values in condition; 
#' represents the referent category.
#' 
#' @param restrict Logical indicating whether to perform restricted Mclust
#'  clustering where close-together clusters are joined.
#' 
#' @importFrom BiocParallel bplapply
#' 
#' @references Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R,
#'  Kendziorski C. A statistical approach for identifying differential 
#'  distributions
#' in single-cell RNA-seq experiments. Genome Biology. 2016 Oct 25;17(1):222. 
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-
#' 1077-y}
#'
#' @return Bayes factor numerator for the current permutation

permMclust <- function(y, nperms, condition, remove.zeroes=TRUE, 
                       log.transf=TRUE, restrict=FALSE, 
                       alpha, m0, s0, a0, b0, ref, min.size){
  
  y.orig <- y
  
  getPerm <- function(y, condition, restrict, remove.zeroes, log.transf){
    if(remove.zeroes & log.transf){
      cond <- condition[y>0]
      y <- log(y[y>0])
    }else if(log.transf){
      cond <- condition
      y <- log(y)
    }else if(remove.zeroes){
      cond <- condition[y>0]
      y <- y[y>0]
    }
    new <- sample(1:length(y), replace=FALSE)
    new.y <- y[new]
    y1 <- new.y[cond==ref]
    y2 <- new.y[cond!=ref]
    
    c1 <- mclustRestricted(y1, restrict=restrict, min.size=min.size)
    c2 <- mclustRestricted(y2, restrict=restrict, min.size=min.size)
    
    bf.p <- jointPosterior(y1, c1, alpha, m0, s0, a0, b0) + 
      jointPosterior(y2, c2, alpha, m0, s0, a0, b0) 
    return(bf.p)
  }
  
  bf.p <- unlist(bplapply(1:nperms, function(x) 
    getPerm(y.orig, condition, restrict, remove.zeroes, log.transf)))
  
  return(bf.p)
}

#' permMclustGene
#'
#' Function to obtain bayes factor for all permutations of one gene 
#' (not parallelized; to be used when parallelizing over Genes)
#' 
#' @details Obtains bayes factor for data vector \code{y} representing one gene 
#' 
#' @param y Numeric data vector for one gene
#'
#' @param ref one of two possible values in condition;
#'  represents the referent category.
#' 
#' @inheritParams jointPosterior 
#' 
#' @inheritParams scDD
#' 
#' @inheritParams permMclustCov
#' 
#' @importFrom BiocParallel bplapply
#' 
#' @references Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R,
#'  Kendziorski C. A statistical approach for identifying differential
#'   distributions
#' in single-cell RNA-seq experiments. Genome Biology. 2016 Oct 25;17(1):222. 
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-
#' 1077-y}
#'
#' @return Bayes factor numerator for the current permutation
permMclustGene <- function(y, adjust.perms, nperms, condition, 
                           remove.zeroes=TRUE, log.transf=TRUE, restrict=TRUE, 
                           alpha, m0, s0, a0, b0, C, ref, min.size){
  orig.y <- y
  
  if(remove.zeroes & log.transf){
    C <- C[y>0]
    cond <- condition[y>0]
    y <- log(y[y>0])
  }else if(log.transf){
    C <- C[y>0]
    cond <- condition
    y <- log(y)
  }else if(remove.zeroes){
    C <- C[y>0]
    cond <- condition[y>0]
    y <- y[y>0]
  }
  
  if (adjust.perms){
    # fit reduced model
    X <- model.matrix(~ C)
    model <- lm(y ~ C)
    params <- summary(model)$coef[,1]
  }
  
  # for loop to go through all nperms permutations (not 
  # parallelized since this function is to be
  # used when there are multiple genes running at a time)
  bf.p <- rep(NA, nperms)
  for (b in 1:nperms){
    if (adjust.perms){
      new.resid <- sample(residuals(model), replace=FALSE)
      new.y0 <- as.vector(exp(X %*% params + new.resid))
      new.y <- orig.y
      new.y[new.y > 0] <- new.y0
      y1 <- new.y[condition==ref]
      y2 <- new.y[condition!=ref]
      
      if(remove.zeroes & log.transf){
        y1 <- log(y1[y1>0])
        y2 <- log(y2[y2>0])
      }else if(log.transf){
        y1 <- log(y1)
        y2 <- log(y2)
      }else if(remove.zeroes){
        y1 <- y1[y1>0]
        y2 <- y2[y2>0]
      }
      
      oa <- mclustRestricted(c(y1,y2), restrict=restrict, min.size=min.size)
      c1 <- mclustRestricted(y1, restrict=restrict, min.size=min.size)
      c2 <- mclustRestricted(y2, restrict=restrict, min.size=min.size)
      
      bf.p[b] <- jointPosterior(y1, c1, alpha, m0, s0, a0, b0) + 
                 jointPosterior(y2, c2, alpha, m0, s0, a0, b0) -
                 jointPosterior(c(y1,y2), oa, alpha, m0, s0, a0, b0)
    }else{
      # don't adjust permutations for covariates; don't 
      # need to calculate denominator
      # of BF since it will be the same in observed and permuted 
      # sets with no residual adjustment
      new <- sample(1:length(y), replace=FALSE)
      new.y <- y[new]
      y1 <- new.y[cond==ref]
      y2 <- new.y[cond!=ref]
      
      c1 <- mclustRestricted(y1, restrict=restrict, min.size=min.size)
      c2 <- mclustRestricted(y2, restrict=restrict, min.size=min.size)
      
      bf.p[b] <- jointPosterior(y1, c1, alpha, m0, s0, a0, b0) + 
                 jointPosterior(y2, c2, alpha, m0, s0, a0, b0) 
    }
  }
  return(bf.p)
}


#' permZero
#'
#' Function to generate random permutations of nonzero values.
#' 
#' @details Generates random permutations for all genes, where the zeroes are 
#' kept fixed (i.e. only permute the nonzero condition labels).
#' 
#' @param m Number of permuted sets to generate.
#' 
#' @param size Number of samples present in the dataset
#' 
#' @param zmat Matrix of indicators of whether the original data value is 
#' zero or not.  Should contain the 
#'   same number of rows and columns as original data matrix.
#' 
#' @references Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, 
#' Kendziorski C. A statistical approach for identifying differential 
#' distributions
#' in single-cell RNA-seq experiments. Genome Biology. 2016 Oct 25;17(1):222. 
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-
#' 1077-y}
#'
#' @return a list of length 'm' (nperms) where each item is a 'ngenes' by 
#' 'size' matrix

permZero <- function(m, size, zmat) { # Obtain m unique combinations of 1:size
  
  # Function to obtain a new permutation.
  newperm <- function(zeroes) {
    # Generate a permutation
    p <- sample(1:(size-length(zeroes)))
    
    # put back into zero context
    newp <- (1:size)
    if (length(zeroes)>0){
      newp[-zeroes] <- newp[-zeroes][p]
      newp[zeroes] <- zeroes
    } else {
      newp <- newp[p]
    }
    newp                   # Return this (new) permutation
  }
  
  # Obtain m unique permutations.
  permlist <- vector("list", m)
  for (perms in 1:m){
    permlist[[perms]] <- t(apply(zmat, 1, function(x) 
      newperm(if (sum(x==1)>0){ which(x==1) }else{NULL})))
  } 
  return(permlist)
}
# Returns a list of length 'm' (nperms) where each 
# item is a 'ngenes' by 'size' matrix
