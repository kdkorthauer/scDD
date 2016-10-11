#' scDD
#' 
#' Find genes with differential distributions (DD) across two conditions
#' 
#' @details Find genes with differential distributions (DD) across two conditions.  Models each log-transformed gene as a Dirichlet 
#'   Process Mixture of normals and uses a permutation test to determine whether condition membership is independent of sample clustering.
#'   The FDR adjusted (Benjamini-Hochberg) permutation p-value is returned along with the classification of each significant gene 
#'   (with p-value less than 0.05 (or 0.025 if also testing for a difference in the proportion of zeroes)) into one of four categories 
#'   (DE, DP, DM, DB).  For genes that do not show significant influence, of condition on clustering, an optional test of whether the 
#'   proportion of zeroes (dropout rate) is different across conditions is performed (DZ).
#'   
#' @param SCdat An object of class \code{ExpressionSet} that contains normalized single-cell expression and metadata, where the \code{assayData} 
#'   slot contains one row for each gene and one sample for each column.  The \code{PhenoData} slot should contain a vector of numeric values
#'   (either 1 or 2) that indicates which 
#'   condition each sample belongs to (in the same order as the columns of \code{assayData}).  Optional additional metadata about the 
#'   experiment can be contained in the \code{experimentData} slot.
#' 
#' @param prior_param A list of prior parameter values to be used when modeling each gene as a mixture of DP normals.  Default 
#'    values are given that specify a vague prior distribution on the cluster-specific means and variances.
#'    
#' @param permutations The number of permutations to be used in calculating empirical p-values.  If set to zero (default),
#'   the full Bayes Factor permutation test will not be performed.  Instead, a fast procedure to identify the genes with significantly different
#'   expression distributions will be performed using the nonparametric Kolmogorov-Smirnov test, which tests the null hypothesis that 
#'   the samples are generated from the same continuous distribution.  This test will yield
#'   slightly lower power than the full permutation testing framework (this effect is more pronounced at smaller sample 
#'   sizes, and is more pronounced in the DB category), but is orders of magnitude faster.  This option
#'   is recommended when compute resources are limited.  The remaining steps of the scDD framework will remain unchanged
#'   (namely, categorizing the significant DD genes into patterns that represent the major distributional changes, 
#'   as well as the ability to visualize the results with violin plots using the \code{sideViolin} function).
#' 
#' @param testZeroes Logical indicating whether or not to test for a difference in the proportion of zeroes
#' 
#' @param adjust.perms Logical indicating whether or not to adjust the permutation tests for the sample
#'   detection rate (proportion of nonzero values).  If true, the residuals of a linear model adjusted for 
#'   detection rate are permuted, and new fitted values are obtained using these residuals.
#'   
#' @param n.cores integer number of cores to use when computing the Bayes Factor estimates for the full permutation testing 
#'  framework.  Defaults to the number of cores returned by \code{parallel::detectCores()}.  To use fewer cores, specify a number
#'  less than the number of cores on your machine.
#'  
#' @param parallelBy For the permutation test (if invoked), the manner in which to parallelize.  The default option
#'  is \code{"Genes"} which will spawn processes that divide up the genes across \code{n.cores} cores, and then loop through the permutations. 
#'  The alternate option is \code{"Permutations"} which
#'  loop through each gene and spawn processes that divide up the permutations across \code{n.cores} cores.  
#'  The default option is recommended when analyzing more genes than the number of permutations.
#' 
#' @param condition A character object that contains the name of the column in \code{phenoData} that represents 
#'  the biological group or condition of interest (e.g. treatment versus control).  Note that this variable should only contain two 
#'  possible values since \code{scDD} can currently only handle two-group comparisons.  The default option assumes that there
#'  is a column named "condition" that contains this variable. 
#' 
#' @return List with four items: the first is a data frame with nine columns: gene name (matches rownames of SCdat), permutation p-value for testing of independence of 
#'  condition membership with clustering, Benjamini-Hochberg adjusted version of the previous column, p-value for test of difference in dropout rate (only for non-DD genes), 
#'  Benjamini-Hochberg adjusted version of the previous column, name of the 
#'  DD (DE, DP, DM, DB) pattern or DZ (otherwise NS = not significant), the number of clusters identified overall, the number of clusters identified in 
#'  condition 1 alone, and the number of clusters identified in condition 2 alone. The remaining three elements are data frames (first for condition 1 and 2 combined, 
#'  then condition 1 alone, then condition 2 alone) that contains the cluster memberships for each sample (cluster 1,2,3,...) in columns and
#'  genes in rows.  Zeroes, which are not involved in the clustering, are labeled as zero.  
#'  
#' @export
#'
#' @importFrom BiocParallel bplapply  
#' 
#' @importFrom BiocParallel register
#' 
#' @importFrom BiocParallel MulticoreParam
#' 
#' @importFrom parallel detectCores
#' 
#' @import Biobase 
#'  
#' @examples 
#'  
#' # load toy simulated example ExpressionSet to find DD genes
#' 
#' data(scDatExSim)
#' 
#' 
#' # check that this object is a member of the ExpressionSet class
#' # and that it contains 200 samples and 30 genes
#' 
#' class(scDatExSim)
#' show(scDatExSim)
#' 
#' 
#' # set arguments to pass to scDD function
#' # we will perform 100 permutations on each of the 30 genes
#' 
#' prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
#' nperms <- 100
#' 
#' 
#' # call the scDD function to perform permutations, classify DD genes, and return results
#' # we won't perform the test for a difference in the proportion of zeroes since none 
#' # exists in this simulated toy example data
#' # this step will take significantly longer with more genes and/or more permutations
#' 
#' RES <- scDD(scDatExSim, prior_param=prior_param, permutations=nperms, testZeroes=FALSE)

scDD <- function(SCdat, prior_param=list(alpha=0.10, mu0=0, s0=0.01, a0=0.01, b0=0.01), permutations=0,
                  testZeroes=TRUE, adjust.perms=FALSE, n.cores=parallel::detectCores(), parallelBy=c("Genes", "Permutations"),
                  condition="condition"){
  
  # check whether SCdat is a member of the ExpressionSet class
  if(!("ExpressionSet" %in% class(SCdat))){
    stop("Please provide a valid 'ExpressionSet' object.")
  }
  
  parallelBy <- match.arg(parallelBy)
  
  # unpack prior param objects
  alpha = prior_param$alpha
  m0 = prior_param$mu0
  s0 = prior_param$s0
  a0 = prior_param$a0
  b0 = prior_param$b0
  
  # check that condition inputs are valid
  if (length(unique(phenoData(SCdat)[[condition]])) != 2 | length(phenoData(SCdat)[[condition]]) != ncol(exprs(SCdat))){
    stop("Error: Please specify valid condition labels.")
  }
  
  # reference category/condition - the first listed one
  ref <- unique(phenoData(SCdat)[[condition]])[1]
  
  # check for genes that are all (or almost all) zeroes 
  tofit <- which((rowSums(exprs(SCdat)[,phenoData(SCdat)[[condition]]==ref]>0) > 1) &
                 (rowSums(exprs(SCdat)[,phenoData(SCdat)[[condition]]!=ref]>0) > 1))
  
  if (length(tofit) < nrow(exprs(SCdat))){
    if(testZeroes){
      message("Notice: There exist genes that are all (or almost all) zero. For genes with 0 or 1 nonzero measurements per condition, only testing for DZ")    
    }else{
      message("Notice: There exist genes that are all (or almost all) zero. Skipping genes with 0 or 1 nonzero measurements per condition")    
    }
  }

  # cluster each gene in SCdat
  message("Clustering observed expression data for each gene")
  message(paste0("Setting up parallel back-end using ", n.cores, " cores" ))
  BiocParallel::register(BPPARAM = BiocParallel::MulticoreParam(workers=n.cores))
  
  oa <- c1 <- c2 <- vector("list", nrow(exprs(SCdat)[tofit,]))
  bf <- den <- comps.all <- comps.c1 <- comps.c2 <- rep(NA, nrow(exprs(SCdat)[tofit,]))
  
  if (permutations == 0){

    # function to fit one gene 
    genefit <- function(y){
      cond0 <- phenoData(SCdat)[[condition]][y>0]
      y <- log(y[y>0])
      
      oa <- mclustRestricted(y, restrict=TRUE)
      c1 <- mclustRestricted(y[cond0==ref], restrict=TRUE)
      c2 <- mclustRestricted(y[cond0!=ref], restrict=TRUE)
    
      return(list(
        oa=oa,
        c1=c1,
        c2=c2
      ))
    }
    
    out <- bplapply(1:nrow(exprs(SCdat)[tofit,]), function(x) genefit(exprs(SCdat)[tofit[x],]))
    oa <- lapply(out, function(x) x[["oa"]])
    c1 <- lapply(out, function(x) x[["c1"]])
    c2 <- lapply(out, function(x) x[["c2"]])
    rm(out); gc()
    
    comps.all <- unlist(lapply(oa, function(x) luOutlier(x$class)))
    comps.c1  <- unlist(lapply(c1, function(x) luOutlier(x$class)))
    comps.c2  <- unlist(lapply(c2, function(x) luOutlier(x$class)))
    
    message("Notice! Number of permutations is set to zero; using Kolmogorov-Smirnov to test for differences in distributions instead of the Bayes Factor permutation test")
    
    res_ks <- testKS(exprs(SCdat)[tofit,], phenoData(SCdat)[[condition]], inclZero=FALSE)
    
    if (testZeroes){
      sig <- which(res_ks$p < 0.025)
    }else{
      sig <- which(res_ks$p < 0.05)
    }
    
    pvals <- res_ks$p.unadj
    
  }else{ 

    # function to fit one gene 
    genefit <- function(y){
      cond0 <- phenoData(SCdat)[[condition]][y>0]
      y <- log(y[y>0])
      
      oa <- mclustRestricted(y, restrict=TRUE)
      c1 <- mclustRestricted(y[cond0==ref], restrict=TRUE)
      c2 <- mclustRestricted(y[cond0!=ref], restrict=TRUE)
      
      bf <- jointPosterior(y[cond0==ref], c1, alpha, m0, s0, a0, b0) + 
        jointPosterior(y[cond0!=ref], c2, alpha, m0, s0, a0, b0) 
      den <- jointPosterior(y, oa, alpha, m0, s0, a0, b0)
      return(list(
        oa=oa,
        c1=c1,
        c2=c2,
        bf=bf,
        den=den
      ))
    }
    
    out <- bplapply(1:nrow(exprs(SCdat)[tofit,]), function(x) genefit(exprs(SCdat)[tofit[x],]))
    oa <- lapply(out, function(x) x[["oa"]])
    c1 <- lapply(out, function(x) x[["c1"]])
    c2 <- lapply(out, function(x) x[["c2"]])
    bf <- unlist(lapply(out, function(x) x[["bf"]]))
    den<- unlist(lapply(out, function(x) x[["den"]]))
    rm(out); gc()
    
    comps.all <- unlist(lapply(oa, function(x) luOutlier(x$class)))
    comps.c1  <- unlist(lapply(c1, function(x) luOutlier(x$class)))
    comps.c2  <- unlist(lapply(c2, function(x) luOutlier(x$class)))
  

      # obtain Bayes Factor score numerators for each permutation
      message("Performing permutations to evaluate independence of clustering and condition for each gene")
      message(paste0("Parallelizing by ", parallelBy))
      bf.perm <- vector("list", nrow(exprs(SCdat)[tofit,]))
      names(bf.perm) <- rownames(exprs(SCdat)[tofit,])
      
      if(parallelBy=="Permutations"){
        if(adjust.perms){
          C <- apply(exprs(SCdat)[tofit,], 2, function(x) sum(x>0)/length(x))
          
          t1 <- proc.time()
          for (g in 1:nrow(exprs(SCdat)[tofit,])){
            bf.perm[[g]] <- permMclustCov(exprs(SCdat)[tofit[g],], permutations, C, phenoData(SCdat)[[condition]], remove.zeroes=TRUE, log.transf=TRUE, restrict=TRUE, 
                                          alpha, m0, s0, a0, b0, ref)
            
            if (g%%1000 == 0){
              t2 <- proc.time()
              message(paste0(g, " genes completed at ", date(), ", took ", round((t2-t1)[3]/60, 2), " minutes")) 
              t1 <- t2
            }
          }
          
        }else{
          t1 <- proc.time()
          for (g in 1:nrow(exprs(SCdat)[tofit,])){
            bf.perm[[g]] <- permMclust(exprs(SCdat[tofit[g],]), permutations, phenoData(SCdat)[[condition]], remove.zeroes=TRUE, log.transf=TRUE, restrict=TRUE, 
                                       alpha, m0, s0, a0, b0, ref)
            
            if (g%%1000 == 0){
              t2 <- proc.time()
              message(paste0(g, " genes completed at ", date(), ", took ", round((t2-t1)[3]/60, 2), " minutes")) 
              t1 <- t2
            }
          }
      }
      }else if(parallelBy=="Genes"){
        C <- apply(exprs(SCdat)[tofit,], 2, function(x) sum(x>0)/length(x))
        bf.perm <- bplapply(1:nrow(exprs(SCdat)[tofit,]), function(x) 
              permMclustGene(exprs(SCdat)[tofit[x],], adjust.perms, permutations, phenoData(SCdat)[[condition]], remove.zeroes=TRUE, log.transf=TRUE, restrict=TRUE, 
                             alpha, m0, s0, a0, b0, C, ref))
      }else{stop("Please specify either 'Permutations' or 'Genes' to parallelize by using the parallelizeBy argument")}
      
      if (adjust.perms){
        pvals <- sapply(1:nrow(exprs(SCdat)[tofit,]), function(x) sum( bf.perm[[x]] > bf[x] - den[x] ) )/(permutations)
      }else{
        pvals <- sapply(1:nrow(exprs(SCdat)[tofit,]), function(x) sum( bf.perm[[x]] > bf[x]) ) / (permutations)
      }
      
      if (testZeroes){
        sig <- which(p.adjust(pvals, method="BH") < 0.025)
      }else{
        sig <- which(p.adjust(pvals, method="BH") < 0.05)
      }
  }
  
  message("Classifying significant genes into patterns")
  dd.cats <- classifyDD(exprs(SCdat)[tofit,], phenoData(SCdat)[[condition]], sig, oa, c1, c2, alpha=alpha, m0=m0, s0=s0, a0=a0, b0=b0, log.nonzero=TRUE, ref=ref)
  
  cats <- rep("NS", nrow(exprs(SCdat)[tofit,]))
  cats[sig] <- dd.cats
  
  extraDP <- feDP(exprs(SCdat)[tofit,], phenoData(SCdat)[[condition]], sig, oa, c1, c2, log.nonzero=TRUE,
                  testZeroes=testZeroes, adjust.perms=adjust.perms)
  cats[-sig] <- names(extraDP)
  
  # classify additional genes with evidence of DD in the form of a mean shift found by 'extraDP'
  if(testZeroes){
    NCs <- which(p.adjust(pvals, method="BH") > 0.025 & cats == "NC")
  }else{
    NCs <- which(p.adjust(pvals, method="BH") > 0.05 & cats == "NC")
  }
  NC.cats <- classifyDD(exprs(SCdat)[tofit,], phenoData(SCdat)[[condition]], NCs, oa, c1, c2, alpha=alpha, m0=m0, s0=s0, a0=a0, b0=b0, log.nonzero=TRUE, ref=ref)
  cats[NCs] <- NC.cats
  
  cats.all <- pvals.all <- rep(NA, nrow(exprs(SCdat)))
  cats.all[tofit] <- cats
  pvals.all[tofit] <- pvals
   
  # zero test
  ns <- which(!(cats.all %in% c("DE", "DP", "DM", "DB")))
  pvals.z <- rep(NA, nrow(exprs(SCdat)))
  if (testZeroes){
    ztest <- testZeroes(exprs(SCdat), phenoData(SCdat)[[condition]], ns)
    pvals.z[ns] <- ztest
    cats.all[p.adjust(pvals.z, method="BH") < 0.025] <- "DZ"
    cats.all[p.adjust(pvals.z, method="BH") >= 0.025] <- "NS"
  }
  
  # build MAP objects
  MAP1 <- matrix(1, nrow=nrow(exprs(SCdat)), ncol=sum(phenoData(SCdat)[[condition]]==ref))
  MAP2 <- matrix(1, nrow=nrow(exprs(SCdat)), ncol=sum(phenoData(SCdat)[[condition]]!=ref))
  MAP <- matrix(1, nrow=nrow(exprs(SCdat)), ncol=ncol(exprs(SCdat)))
  rownames(MAP1) <- rownames(MAP2) <- rownames(MAP) <- featureNames(SCdat)
  colnames(MAP1) <- sampleNames(SCdat[,phenoData(SCdat)[[condition]]==ref])
  colnames(MAP2) <- sampleNames(SCdat[,phenoData(SCdat)[[condition]]!=ref])
  colnames(MAP) <- sampleNames(SCdat)
  MAP1[exprs(SCdat)[, phenoData(SCdat)[[condition]]==ref]==0] <- 0
  MAP2[exprs(SCdat)[, phenoData(SCdat)[[condition]]!=ref]==0] <- 0
  MAP[exprs(SCdat)==0] <- 0
  
  for (g in 1:nrow(exprs(SCdat)[tofit,])){
    MAP1[tofit[g],][exprs(SCdat[tofit[g], phenoData(SCdat)[[condition]]==ref])!=0] <- c1[[g]]$class + 1 
    MAP2[tofit[g],][exprs(SCdat[tofit[g], phenoData(SCdat)[[condition]]!=ref])!=0] <- c2[[g]]$class + 1 
    MAP[tofit[g],][exprs(SCdat[tofit[g], ])!=0] <- oa[[g]]$class + 1 
  }
  
  comps.all.ALL <- comps.c1.ALL <- comps.c2.ALL <- rep(NA, nrow(exprs(SCdat)))
  comps.all.ALL[tofit] <- comps.all
  comps.c1.ALL[tofit] <- comps.c1
  comps.c2.ALL[tofit] <- comps.c2
  
  # return...
  return(list(Genes=data.frame(gene=rownames(SCdat), nonzero.pvalue=pvals.all, nonzero.pvalue.adj=p.adjust(pvals.all, method="BH"), 
                    zero.pvalue=pvals.z, zero.pvalue.adj=p.adjust(pvals.z, method="BH"), DDcategory=cats.all, 
                    Clusters.combined=comps.all.ALL, Clusters.c1=comps.c1.ALL, Clusters.c2=comps.c2.ALL), 
          Zhat.combined=MAP, Zhat.c1=MAP1, Zhat.c2=MAP2))
}


