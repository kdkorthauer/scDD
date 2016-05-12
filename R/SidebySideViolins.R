#' sideViolin
#'
#' Plots two histograms side by side with smoothed density overlay
#'
#' 
#' @param y Numeric vector of data to plot.
#'
#' @param cond Vector of condition labels corresponding to elements of \code{x}.
#' 
#' @param MAP List of MAP partition estimates with conditions as list items and samples as elements 
#'  (integer indicating which cluster each observation belongs to; zeroes belong to cluster 1)
#'
#' @param logT Logical that indicates whether to take the log(x+1) transformation.
#' 
#' @param title.gene Character vector that contains the gene name that you are plotting.
#' 
#' @param conditionLabels Character vector containing the names of the two conditions.
#' 
#' @param axes.titles Logical indicating whether or not to include axes labels on plots.  
#' 
#' @export
#' 
#' @import ggplot2
#'
#' @examples 
#' 
#' # load toy simulated example ExpressionSet to find DD genes
#' 
#' data(scDatExSim)
#' 
#' 
#' # load Biobase package to facilitate subset operations on ExpressionSet class objects
#' 
#' library(Biobase)
#' 
#' 
#' # plot side by side violin plots for Gene 1 (DE)
#' 
#' sideViolin(exprs(scDatExSim)[1,], scDatExSim$condition, title.gene=featureNames(scDatExSim)[1])
#' 
#' 
#' # plot side by side violin plots for Gene 6 (DP)
#' 
#' sideViolin(exprs(scDatExSim)[6,], scDatExSim$condition, title.gene=featureNames(scDatExSim)[6])
#' 
#' 
#' # plot side by side violin plots for Gene 11 (DM)
#' 
#' sideViolin(exprs(scDatExSim)[11,], scDatExSim$condition, title.gene=featureNames(scDatExSim)[11])
#' 
#' 
#' # plot side by side violin plots for Gene 16 (DB)
#' 
#' sideViolin(exprs(scDatExSim)[16,], scDatExSim$condition, title.gene=featureNames(scDatExSim)[16])
#' 
#' # plot side by side violin plots for Gene 21 (EP)
#' 
#' sideViolin(exprs(scDatExSim)[21,], scDatExSim$condition, title.gene=featureNames(scDatExSim)[21])
#' 
#' 
#' # plot side by side violin plots for Gene 26 (EE)
#' 
#' sideViolin(exprs(scDatExSim)[26,], scDatExSim$condition, title.gene=featureNames(scDatExSim)[26])
#' 
#' 


sideViolin <- function(y, cond, MAP=NULL, logT=TRUE, title.gene="", conditionLabels=NULL, axes.titles=TRUE){
  
  shps1 <- shps2 <- c("a", "c", "d", "e", "f")
  if(!is.null(MAP)){
    shps1 <- shps1[MAP[[1]] + 1]
    shps2 <- shps2[MAP[[2]] + 1]
  }else{
    shps1 <- rep(shps1[1], sum(cond==1))
    shps2 <- rep(shps2[1], sum(cond==2))
  }
  
  if (logT){
    shps1[y[cond==1]==0] <- "b"
    shps2[y[cond==2]==0] <- "b"
    y <- log(y+1)
  }
  
  if (length(conditionLabels)==2){
    cond[cond==1] <- conditionLabels[1]
    cond[cond==2] <- conditionLabels[2]
  }
  
  if(axes.titles){
    xlabel <- ggplot2::element_text()
    ylabel <- ggplot2::element_text()
  }else{
    xlabel <- ylabel <- ggplot2::element_blank()
  }
  
  shps=c(shps1, shps2)
  daty <- data.frame(y, cond, shps)
  
  g <- ggplot(daty, aes(factor(cond), y), aes(shape=shps))
  g + geom_jitter(alpha=0.5, color="black", position = position_jitter(width = 0.15),
                      aes(shape=shps), show_guide=FALSE) +
    geom_violin(data=daty[daty$y>0, ], alpha=0.5, aes(fill=factor(cond)), show_guide=FALSE, scale="count") + 
    ggtitle(paste0(title.gene)) +
    theme(plot.title = element_text(size=20, face="bold", vjust=2)) + 
    labs(x="Condition", y="log(EC + 1)") +
    theme(axis.text.x=element_text(size=14, vjust=0.5), 
          axis.text.y=element_text(size=14, vjust=0.5), 
          axis.title.x = xlabel,
          axis.title.y = ylabel)
}


