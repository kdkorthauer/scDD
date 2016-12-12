#' findFC
#'
#' Find the appropriate Fold Change vectors for simulation that will be use in
#'  classic differential expression case.
#' 
#' @details This code is a modified version of Sam Younkin's simulate FC 
#' function.  Major things that were changed are 
#'   (1) standard deviations are calculated only on the nonzeroes, (2) the
#'    sampling of FCs is uniform on the log scale 
#'   instead of the raw scale, and (3) the binning is done by quantiles 
#'   instead of evenly spaced along the average expression
#'   values.
#'   
#' @inheritParams singleCellSimu
#'
#' @inheritParams scDD
#' 
#' @param sd.range Numeric vector of length two which describes the interval
#'  (lower, upper) of standard deviations
#'  of fold changes to randomly select. 
#'   
#' @param N Integer value for the number of bins to divide range of fold 
#' changes for calculating standard deviations
#' 
#' @param overExpressionProb Numeric value between 0 and 1 which describes 
#' the ratio of over to under expression 
#'  values to sample.
#'  
#' @param plot.FC Logical indicating whether or not to plot the observed 
#' and simulated log2 fold changes.
#'
#' @importFrom fields stats.bin
#' 
#' @references Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, 
#' Kendziorski C. A statistical approach for identifying differential 
#' distributions
#' in single-cell RNA-seq experiments. Genome Biology. 2016 Oct 25;17(1):222. 
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-
#' 1077-y}
#'
#' @return FC.vec Return Fold Change Vectors




findFC <- function(SCdat, index, sd.range = c(1,3), N = 4, 
                   overExpressionProb = 0.5, plot.FC=FALSE, 
                   condition="condition"){

  # reference category/condition - the first listed one
  ref <- unique(colData(SCdat)[[condition]])[1]
  
# separate ExpressionSet by condition  
Dataset1 <- SCdat[,colData(SCdat)[[condition]]==ref]
Dataset2 <- SCdat[,colData(SCdat)[[condition]]!=ref]

### Find Mean Difference
x1bar <- apply(normExprs(Dataset1)[index,], 1, function(x) 
  if(sum(x!=0)>0){ mean(x[x!=0]) }else{0})
x2bar <- apply(normExprs(Dataset2)[index,], 1, function(x) 
  if(sum(x!=0)>0){ mean(x[x!=0]) }else{0})

### Find FC by MA_Plot
createFindFC <- function(x1bar, x2bar, sd.range=sd.range, N=N, 
                         overExpressionProb=overExpressionProb){
    x <- 1/2*log2(x1bar*x2bar)
    y <- log2(x2bar/x1bar)
    bindata <- fields::stats.bin(x, y, breaks=quantile(x, seq(0,1,by=1/N)))
    f <- function(x1, x2){
        if( x1 == 0 | x2 == 0 ){
            FC <- NA
        }else{
            binnum <- max(sum(1/2*log2(x1*x2) >= bindata$breaks+0.0005),1)
            sd <- bindata$stats["Std.Dev.",binnum]
            if( is.na(sd) & binnum == 1 ){
                message("Bin 1 has sd = NA.  Using sd from bin 2.")
                sd <- bindata$stats["Std.Dev.",binnum+1]
            }else if( is.na(sd) & binnum == N){
                message(paste0("Bin ", N, " has sd = NA.  Using sd from bin ",
                               N-1, "."))
                sd <- bindata$stats["Std.Dev.",N-1]
            }else if( is.na(sd) & binnum == (N+1)){
              message(paste0("Bin ", N+1, " has sd = NA.  Using sd from bin ",
                             N, "."))
              sd <- bindata$stats["Std.Dev.",N]
            }else if( is.na(sd)){
              message(paste0("Bin ", N, " has sd = NA.  Using sd from bin ",
                             N-1, "."))
              sd <- bindata$stats["Std.Dev.",N-1]
            }
            lowerFC <- sd.range[1]*sd
            upperFC <- sd.range[2]*sd
            FC.range <- 2^(c(lowerFC,upperFC))
            if( runif(1) < overExpressionProb ){
                FC <- exp(runif(1, min = log(FC.range[1]), 
                                max = log(FC.range[2])))
            }else{
                FC <- exp(runif(1, max = log(1/FC.range[1]), 
                                min = log(1/FC.range[2])))
            }
        }
        return(FC)
    }
    return(f)
}
findFC <- createFindFC(x1bar = x1bar, x2bar = x2bar, sd.range=sd.range, 
                       N=N, overExpressionProb=overExpressionProb)
FC.vec <- mapply(findFC, x1 = x1bar, x2 = x2bar)

if(plot.FC){
  par(mfrow=c(1,1))
  plot(1/2*log2(x1bar*x2bar), log2(FC.vec), pch=20, cex=0.5, col="red", 
       xlab="Mean Log2 Expression", ylab="Log2 FC", main="Observed (grey)
       and Sampled DE (red) Fold Changes")
  points(1/2*log2(x1bar*x2bar), log2(x1bar/x2bar), pch=20, cex=0.5, col="grey")
  points(1/2*log2(x1bar*x2bar), log2(FC.vec), pch=20, cex=0.5, col="red")
}
FC <- rep(NA, nrow(normExprs(Dataset1)))
FC[index] <- FC.vec
return(FC)
}
