#' simulateSet
#'
#' Simulation of a complete dataset, where the number of each type of differential distributions and equivalent distributions is specified.
#'
#' @inheritParams singleCellSimu
#' 
#' @inheritParams scDD
#' 
#' @inheritParams findFC
#' 
#' @param nDE Number of DE genes to simulate
#' 
#' @param nDP Number of DP genes to simulate
#' 
#' @param nDM Number of DM genes to simulate
#' 
#' @param nDB Number of DB genes to simulate
#' 
#' @param nEE Number of EE genes to simulate
#' 
#' @param nEP Number of EP genes to simulate
#' 
#' @param plots Logical indicating whether or not to generate fold change and validation plots
#' 
#' @param random.seed Numeric value for a call to \code{set.seed} for reproducibility.
#' 
#' @param plot.file Character containing the file string if the plots are to be sent to a pdf instead of to the standard output.
#' 
#' @inheritParams singleCellSimu
#' 
#' @inheritParams findIndex
#' 
#' @export
#' 
#' @import Biobase 
#'
#' @return A named list of two items: the first (labeled 'Simulated_Data') is a matrix of simulated 
#'   data with \code{numSamples} columns and \code{nDE + nDP + nDM + nDB + nEE + nEP} rows 
#'   (total number of genes).  The second item (named 'FC') is a vector of the number of standard 
#'   deviations used for fold changes. For DE genes, this value is computed from the sampled fold
#'   changes obtained from \code{\link{findFC}}.  For DP, DM, DB, and EP genes, this is one of the 
#'   values in \code{modeFC}.  For EE genes, this value is \code{NA}.
#' 
#' @examples 
#' 
#' # Load toy example ExpressionSet to simulate from
#' 
#' data(scDatEx)
#' 
#' 
#' # check that this object is a member of the ExpressionSet class
#' # and that it contains 142 samples and 500 genes
#' 
#' class(scDatEx)
#' show(scDatEx)
#' 
#' 
#' # set arguments to pass to simulateSet function
#' # we will simuate 30 genes total; 5 genes of each type;
#' # and 100 samples in each of two conditions
#' 
#' nDE <- 5
#' nDP <- 5
#' nDM <- 5
#' nDB <- 5
#' nEE <- 5
#' nEP <- 5
#' numSamples <- 100
#' seed <- 816
#' 
#' 
#' # create simulated set with specified numbers of DE, DP, DM, DM, EE, and EP genes,
#' # specified number of samples, DE genes are 2 standard deviations apart, and 
#' # multimodal genes have modal distance of 4 standard deviations
#' 
#' SD <- simulateSet(scDatEx, numSamples=numSamples, nDE=nDE, nDP=nDP, nDM=nDM, nDB=nDB, 
#'                   nEE=nEE, nEP=nEP, sd.range=c(2,2), modeFC=4, plots=FALSE, 
#'                   random.seed=seed)
#'                   
#'                   
#' # convert the simulated data object returned by simulateSet into an ExpressionSet object
#' 
#' library(Biobase)   # needed to create and instance of the ExpressionSet class
#' condition <- c(rep(1, numSamples), rep(2, numSamples))
#' rownames(SD[[1]]) <- paste0(rownames(SD[[1]]), 1:nrow(SD[[1]]), sep="")
#' colnames(SD[[1]]) <- names(condition) <- paste0("Sample", 1:ncol(SD[[1]]), sep="")
#' SDExpressionSet <- ExpressionSet(assayData=SD[[1]], 
#'                      phenoData=as(data.frame(condition), "AnnotatedDataFrame"))


simulateSet <- function(SCdat, numSamples=100, nDE=250, nDP=250, nDM=250, nDB=250, 
                         nEE=5000, nEP=4000, sd.range=c(1,3), modeFC=c(2,3,4), plots=TRUE, plot.file=NULL, random.seed=284){
if(!is.null(plot.file)){
  pdf(file=paste0(plot.file))
}

print("Identifying a set of genes to simulate from...")  
index <- findIndex(SCdat)
print("Simulating DE fold changes...")  
FC <- findFC(SCdat, index, sd.range=sd.range, N=6, overExpressionProb = 0.5, plot.FC=TRUE)

constantZero <- NULL
generateZero <- "empirical"

print("Simulating individual genes...")

# pull off matrix of expression values for condition 1
Dataset1 <- exprs(SCdat[,SCdat$condition==1])

set.seed(random.seed)
### DE
SD1 <- singleCellSimu(Dataset1, Method = "DE", index, FC, modeFC, Validation = FALSE, numGenes=nDE, numDE=nDE,
                        numSamples=numSamples, generateZero=generateZero,
                        constantZero=constantZero)
Simulated_Data <- SD1[[1]]
rnms <- rep("EE", nrow(Simulated_Data))
rnms[SD1[[2]]] <- "DE"
rownames(Simulated_Data) <- rnms
pe_mat <- Simulated_Data
fcs <- SD1[[3]]

### DP
SD2 <- singleCellSimu(Dataset1, Method = "DP", index, FC, modeFC, DP = c(0.33,0.66), Validation = FALSE,
                        numGenes=nDP, numDE=nDP, numSamples=numSamples,
                        generateZero=generateZero, constantZero=constantZero)
Simulated_Data_DP <- SD2[[1]]
rnms <- rep("EE", nrow(Simulated_Data_DP))
rnms[SD2[[2]]] <- "DP"
rownames(Simulated_Data_DP) <- rnms
pe_mat <- rbind(pe_mat, Simulated_Data_DP)
fcs <- c(fcs, SD2[[3]])

### DM
SD3 <- singleCellSimu(Dataset1, Method = "DM", index, FC, modeFC, Validation = FALSE,
                        numGenes=nDM, numDE=nDM, numSamples=numSamples, 
                        generateZero=generateZero, constantZero=constantZero)
Simulated_Data_DM <- SD3[[1]]
rnms <- rep("EE", nrow(Simulated_Data_DM))
rnms[SD3[[2]]] <- "DM"
rownames(Simulated_Data_DM) <- rnms
pe_mat <- rbind(pe_mat, Simulated_Data_DM)
fcs <- c(fcs, SD3[[3]])


### DB
SD4 <- singleCellSimu(Dataset1, Method = "DB", index, FC, modeFC, DP = c(0.5,0.5), Validation = FALSE,
                        numGenes=nDB, numDE=nDB, numSamples=numSamples, 
                        generateZero=generateZero, constantZero=constantZero)
Simulated_Data_DB <- SD4[[1]]
rnms <- rep("EE", nrow(Simulated_Data_DB))
rnms[SD4[[2]]] <- "DB"
rownames(Simulated_Data_DB) <- rnms
pe_mat <- rbind(pe_mat, Simulated_Data_DB)
fcs <- c(fcs, SD4[[3]])

### EP
SD5 <- singleCellSimu(Dataset1, Method = "DP", index, FC, modeFC, DP = c(0.50,0.50), Validation = FALSE,
                        numGenes=nEP, numDE=0, numSamples=numSamples, 
                        generateZero=generateZero, constantZero=constantZero)
Simulated_Data_EP <- SD5[[1]]
rnms <- rep("EP", nrow(Simulated_Data_EP))
rnms[SD5[[2]]] <- "DP"
rownames(Simulated_Data_EP) <- rnms
pe_mat <- rbind(pe_mat, Simulated_Data_EP)
fcs <- c(fcs, SD5[[3]])

### EE
SD6<- singleCellSimu(Dataset1, Method = "DE", index, FC, modeFC, Validation = plots, 
                       numGenes=nEE, numDE=0, numSamples=numSamples, 
                       generateZero=generateZero, constantZero=constantZero)
Simulated_Data_EE <- SD6[[1]]
rnms <- rep("EE", nrow(Simulated_Data_EE))
rnms[SD6[[2]]] <- "DE"
rownames(Simulated_Data_EE) <- rnms
pe_mat <- rbind(pe_mat, Simulated_Data_EE)
fcs <- c(fcs, rep(NA, nEE))
names(fcs) <- rownames(pe_mat[(1:sum(nDE,nDP,nDM,nDB)),])

if (!is.null(plot.file)){
  dev.off()
}

# continuity correction
unifmat <- matrix(runif(nrow(pe_mat)*ncol(pe_mat)), nrow=nrow(pe_mat), ncol=ncol(pe_mat))
pe_mat2 <- pe_mat + unifmat
pe_mat2[pe_mat==0] <- 0

SD <- list(Simulated_Data=pe_mat2, FC=fcs)
print(paste0("Done! Simulated ", nDE, " DE, ", nDP, " DP, ", nDM, " DM, ", nDB, " DB, ", nEE, " EE, and ", nEP, " EP genes "))
return(SD)
}


