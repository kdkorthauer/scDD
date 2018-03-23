#' getPosteriorParams
#' 
#' Given the observations for a single gene and its clustering information, 
#' return the calculated posterior parameters
#' 
#' @inheritParams jointPosterior
#' 
#' @return A list of posterior parameter values under the DP mixture model 
#' framework, given the data and prior parameter values.
#' 
#' @references Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, 
#' Kendziorski C. A statistical approach for identifying differential 
#' distributions
#' in single-cell RNA-seq experiments. Genome Biology. 2016 Oct 25;17(1):222. 
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-
#' 1077-y}
#' 

getPosteriorParams <- function(y, mcobj, alpha, m0, s0, a0, b0){
  
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
  
  return(list(mk=mk, sk=sk, ak=ak, bk=bk))
}

#' classifyDD
#'
#' Classify significantly DD genes into the four categories (DE, DP, DM or DB)
#'  based on posterior distributions of cluster mean parameters
#'
#' @inheritParams jointPosterior
#' 
#' @inheritParams scDD
#' 
#' @param pe_mat Matrix with genes in rows and samples in columns.  
#' Column names indicate condition.
#'
#' @param condition Vector of condition indicators (with two possible values).
#' 
#' @param sig_genes Vector of the indices of significantly DD genes 
#' (indicating the row number of \code{pe_mat})
#'  
#' @param oa List item with one item for each gene where the first element 
#' contains the cluster membership for 
#'  each nonzero sample in the overall (pooled) fit.
#'
#' @param c1 List item with one item for each gene where the first element 
#' contains the cluster membership for 
#'  each nonzero sample in condition 1 only fit
#'  
#' @param c2 List item with one item for each gene where the first element 
#' contains the cluster membership for 
#'  each nonzero sample in condition 2 only fit
#'  
#' @param log.nonzero Logical indicating whether to perform log 
#' transformation of nonzero values.
#'
#' @param ref one of two possible values in condition; 
#' represents the referent category.
#' 
#' @references Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, 
#' Kendziorski C. A statistical approach for identifying differential 
#' distributions
#' in single-cell RNA-seq experiments. Genome Biology. 2016 Oct 25;17(1):222. 
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-
#' 1077-y}
#'
#' @return cat Character vector of the same length as \code{sig_genes} that 
#' indicates which category of 
#'  DD each significant gene belongs to (DE, DP, DM, DB, or NC (no call))
#' 

classifyDD <- function(pe_mat, condition, sig_genes, oa, c1, c2, 
                       alpha, m0, s0, a0, b0, log.nonzero=TRUE, 
                      adjust.perms=FALSE, ref, min.size=3){
    ms <- min.size
    mc.oa <- oa
    mc.c1 <- c1
    mc.c2 <- c2
    
    oa <- lapply(oa, function(x) x[[1]])
    c1 <- lapply(c1, function(x) x[[1]])
    c2 <- lapply(c2, function(x) x[[1]])
    
    c.oa <- sapply(sig_genes, function(x) luOutlier(oa[[x]], min.size=ms))
    c.c1 <- sapply(sig_genes, function(x) luOutlier(c1[[x]], min.size=ms))
    c.c2 <- sapply(sig_genes, function(x) luOutlier(c2[[x]], min.size=ms))
    
    cdr <- apply(pe_mat, 2, function(x) sum(x>0)/length(x))
    
    # for each gene
    cat <- rep(NA, length(sig_genes))
    
    s <- 1
    for (g in sig_genes){
      
      y <- pe_mat[g,]
      
      if (log.nonzero){
        cond <- condition[y>0]
        cdr0 <- cdr[y>0]
        y <- log(y[y>0])
      }else{
        cond <- condition
        cdr0 <- cdr
      }
      
      if(length(unique(c1[[g]])) != max(c1[[g]])){
        missing.label <- which(! ((1:max(c1[[g]]) %in% unique(c1[[g]])) ))
        new.labels <- c1[[g]]
        for (l in (missing.label+1):max(c1[[g]])){
          new.labels[new.labels==l] <- l-1
        }		
        c1[[g]] <- new.labels
        mc.c1[[g]]$class <- c1[[g]]
        mc.c1[[g]]$mean <- mc.c1[[g]]$mean[-missing.label]
        mc.c1[[g]]$var <- mc.c1[[g]]$var[-missing.label]
      }
      if(length(unique(c2[[g]])) != max(c2[[g]])){
        missing.label <- which(! ((1:max(c2[[g]]) %in% unique(c2[[g]])) ))
        new.labels <- c2[[g]]
        for (l in (missing.label+1):max(c2[[g]])){
          new.labels[new.labels==l] <- l-1
        }		
        c2[[g]] <- new.labels
        mc.c2[[g]]$class <- c2[[g]]
        mc.c2[[g]]$mean <- mc.c2[[g]]$mean[-missing.label]
        mc.c2[[g]]$var <- mc.c2[[g]]$var[-missing.label]
      }
      if(length(unique(oa[[g]])) != max(oa[[g]])){
        missing.label <- which(! ((1:max(oa[[g]]) %in% unique(oa[[g]])) ))
        new.labels <- oa[[g]]
      for (l in (missing.label+1):max(oa[[g]])){
          new.labels[new.labels==l] <- l-1
      }		
        oa[[g]] <- new.labels
        mc.oa[[g]]$class <- oa[[g]]
        mc.oa[[g]]$mean <- mc.oa[[g]]$mean[-missing.label]
        mc.oa[[g]]$var <- mc.oa[[g]]$var[-missing.label]
      }
      
      
      
      params.c1 <- getPosteriorParams(y[cond==ref], mc.c1[[g]], 
                                      alpha, m0, s0, a0, b0)
      params.c2 <- getPosteriorParams(y[cond!=ref], mc.c2[[g]], 
                                      alpha, m0, s0, a0, b0)
      
      c.c1.no <- findOutliers(c1[[g]], min.size=ms)
      c.c2.no <- findOutliers(c2[[g]], min.size=ms)
      
      # posterior analysis of cluster distances
      num.comparisons <- length(c.c1.no)*length(c.c2.no)
      comparisons <- rep(NA, num.comparisons)
      
      t <- 1
      for (a in c.c1.no){
        for (b in c.c2.no){
          
            samp.diff <- quantile(rt(10000, df=round(params.c1$ak[a]))*
                                    sqrt(params.c1$bk[a]/
                                        (params.c1$ak[a]*params.c1$sk[a])) +
                                    params.c1$mk[a] - 
                                    rt(10000, df=round(params.c2$ak[b]))*
                                    sqrt(params.c2$bk[b]/(params.c2$ak[b]*
                                                        params.c2$sk[b])) 
                                  - params.c2$mk[b],
                                  c(0, 1))
            comparisons[t] <- 1*(  (0 < min(samp.diff) & 0 < max(samp.diff)) | 
                                  (-0 > min(samp.diff) & -0 > max(samp.diff)) )
          
          t <- t+1
        }
      }
      
      if (c.c1[s]==c.c2[s]){ # same number of clusters in each condition
        if (c.c1[s]==1){ # one cluster in each condition
            if (sum(comparisons)>0){
              cat[s] <- "DE"
            }else{
              cat[s] <- "NC"
            }
        }else{
          if (c.oa[s]==c.c1[s]){ # at least two clusters overall, same number 
                                 # within each condition as overall
            if (sum(comparisons)<=c.c1[s]){
              cat[s] <- "DP"
            }else{
              cat[s] <- "NC"
            }
          }else if(c.oa[s]>c.c1[s]){ # equal number of clusters within each 
                              #conditon (at least 2) and more than that overall
            if (sum(comparisons)> c.c1[s]*(c.c1[s]-1)){
              # Require that the multimodal DE genes show 
              # evidence of an overall shift in mean
              # in addition to having at most one pair of component overlaps
              if(adjust.perms){
                comparison <- summary(lm(y ~ cdr0 + factor(cond)))$coef[3,4]
              }else{
                comparison <- t.test(y~factor(cond))$p.value	
              }
              if (comparison < 0.01){
                cat[s] <- "DE"
              }else{
                cat[s] <- "NC"
              }
            }else{
              cat[s] <- "NC"
            }
          }else{
            cat[s] <- "NC"
          }
        }
      }else{ # not equal cluster number in each condition
        if (sum(comparisons)==length(comparisons)){
          cat[s] <- "DB"
        }else{
          cat[s] <- "DM"
        }
      }
      s <- s + 1
    }
    names(cat) <- names(sig_genes)
    return(cat)
  } 

#' testZeroes
#' 
#' Test for a difference in the proportion of zeroes between conditions 
#' for a specified set of genes
#' 
#' @details Test for a difference in the proportion of zeroes between 
#' conditions that is not explained by the 
#'   detection rate.  Utilizes Bayesian logistic regression.
#' 
#' @param dat Matrix of single cell expression data with genes in rows and 
#' samples in columns.  
#' 
#' @param cond Vector of condition labels
#' 
#' @param these vector of row numbers (gene numbers) to test for a difference
#'  in the proportion of zeroes.
#' 
#' @return Vector of FDR adjusted p-values
#' 
#'
#' @importFrom arm bayesglm  

testZeroes <- function(dat, cond, these=1:nrow(dat)){
  detection <- colSums(dat>0)/nrow(dat)
  
  onegene <- function(y, detection){
    if (sum(y==0) > 0){
      M0 <- suppressWarnings(arm::bayesglm(y>0 ~ detection, 
                                           family=binomial(link="logit"),
                                           Warning=FALSE))
      M1 <- suppressWarnings(arm::bayesglm(y>0 ~ detection + factor(cond), 
                                           family=binomial(link="logit"),
                                           Warning=FALSE))
      return(anova(M1, M0, test="Chisq")[2,5])
    }else{
      return(NA)
    }
  }
  
  pval <- unlist(bplapply(seq_along(these), 
                   function(j) onegene(y=dat[these[j],], detection=detection)))
  return(pval)
}

#' feDP
#'
#' Function to identify additional DP genes, since clustering process can be 
#' consistent within each condition 
#'   and still have differential proportion within each mode.  
#'   The Bayes factor score also tends to be small when
#'   the correct number of clusters is not correctly detected; 
#'   in that case differential proportion will manifest
#'   as a mean shift.
#'   
#' @details The Fisher's Exact test is used to test for independence of 
#' condition membership and clustering when 
#'   the clustering is the same across conditions as it is overall 
#'   (and is multimodal). When clustering within 
#'   condition is not multimodal or is different across conditions 
#'   (most often the case), an FDR-adjusted t-test
#'   is performed to detect overall mean shifts.
#'
#' @inheritParams classifyDD
#' 
#' @inheritParams scDD
#'
#' @return cat Character vector of the same length as \code{sig_genes} 
#' that indicates which nonsignificant genes by
#'  the permutation test belong to the DP category
#' 
feDP <- function(pe_mat, condition, sig_genes, oa, c1, c2, log.nonzero=TRUE, 
                 testZeroes=FALSE, adjust.perms=FALSE, min.size=3){
  if(testZeroes){
    pval.thresh <- 0.025
  }else{
    pval.thresh <- 0.05
  }
  
  oa <- lapply(oa, function(x) x[[1]])
  c1 <- lapply(c1, function(x) x[[1]])
  c2 <- lapply(c2, function(x) x[[1]])
  
  ns_genes <- (1:nrow(pe_mat))[-sig_genes]
  c.oa <- sapply(ns_genes, function(x) luOutlier(oa[[x]], min.size))
  c.c1 <- sapply(ns_genes, function(x) luOutlier(c1[[x]], min.size))
  c.c2 <- sapply(ns_genes, function(x) luOutlier(c2[[x]], min.size))
  
  cdr <- apply(pe_mat, 2, function(x) sum(x>0)/length(x))
  cat <- rep(NA, length(ns_genes))
  pval.ns <- rep(NA, length(ns_genes))
  s <- 1
  
  ref <- unique(condition)[1]
  
  for (g in ns_genes){
    y <- pe_mat[g,]
    cond <- condition[y>0]
    cdr0 <- cdr[y>0]
    y <- log(y[y>0])
    
    # add runif(-0.1,0.1) jitter if all y vals are identical
    if (length(unique(y[cond==ref]))==1 | length(unique(y[cond!=ref]))==1){
      y <- y + runif(length(y), -0.1, 0.1)
    }
    
    # detect shifts in mean (to catch DP genes with an incorrect # components)
    if(adjust.perms){
      pval.ns[s] <- summary(lm(y ~ cdr0 + factor(cond)))$coef[3,4]
    }else{
      pval.ns[s] <- t.test(y~factor(cond))$p.value	
    }
  
    # check whether clustering process is consistent within each
    # condition as overall
    if (c.c1[s]==c.c2[s] & c.c1[s]==c.oa[s] & c.oa[s]>1){
      # non-outlying cluster names
      c.oa.no <- findOutliers(oa[[g]], min.size)
      c.c1.no <- findOutliers(c1[[g]], min.size)
      c.c2.no <- findOutliers(c2[[g]], min.size)
      
      test <- 1*(fisher.test(table(oa[[g]][oa[[g]] %in% c.oa.no], 
                                   cond[oa[[g]] %in% c.oa.no]))$p.value < 0.05)
      if (test==1){
        cat[s] <- "DP"
      }else{
        cat[s] <- "NS"
      }
    }else{
      cat[s] <- "NS"
    }
  s <- s+1
  }
  
  pval.ns <- p.adjust(pval.ns, method="BH")
  cat[pval.ns < pval.thresh & cat != "DP"] <- "NC"
  cat[pval.ns < pval.thresh & c.c1 == c.c2 & c.c1 == c.oa & c.c1 > 1] <- "DP"			

  names(pval.ns) <- cat
  return(pval.ns)
}


