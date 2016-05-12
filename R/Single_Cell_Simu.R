#' calcRP
#'
#'  Calculate parameter R and P in NB distribution
#'
#' @param Emean Empirical mean
#'
#' @param Evar Empirical variance
#' 
#' @return RP Vector of two elements, first contains method of moments estimator for r and second contains 
#'   method of moments estimator for p (parameters of NB distribution)
#'   

calcRP <- function(Emean, Evar){
  RP <- rep(NA,2)
  RP[1]=max(0.01,(Emean)^2/(Evar-Emean));
  RP[2]=max(0.01,1-Emean/Evar);
  if((RP[1]==0.01)&(RP[2]==0.01))
  {
    RP[1]=1
    RP[2]=0.5
  }
  RP[2]=min(RP[2],0.99999999)
  return(RP)
}


#' calcMV
#'
#'  Calculate empirical means and variances of selected genes in a given dataset.
#'  
#' @details Calculate empirical means and variances of selected genes in a given dataset. Optionally, multiply
#'  the means and standard deviations by a fold change value, which can also vary by mean value.  If the mean is
#'  below some specified threshold \code{threshold}, use one fold change value \code{FC}.  If above the threshold, use the alternate 
#'  fold change value \code{FC.thresh}.  Estimates of mean and variance are robust to outliers.
#'
#' @param a Numeric vector of values to calculate empirical mean and variance.
#'
#' @param FC Fold change for the mean and standard deviation.  Default value is 1.
#' 
#' @param FC.thresh Alternate fold change for the mean and standard deviation when the (log nonzero) mean is above the value of \code{threshold}.  
#'  Default value is \code{FC}.
#' 
#' @param threshold Mean threshold value which dictates which fold change value to use for multiplying mean and standard deviation.
#'  Default value is Inf (so \code{FC} is always used).
#'  
#' @param include.zeroes Logical value indicating whether the zero values should be included in the calculations of the empirical means 
#'    and variances.
#' 
#' @return MV Vector of two elements, first contains the empirical mean estimate, second contains the empirical variance estimate (optionally 
#'   multiplied by a fold change).
#'   

calcMV <- function(a, FC=1, FC.thresh=NA, threshold=Inf, include.zeroes=FALSE){
  MV <- rep(NA, 2)
  if(outliers::grubbs.test(log(a+1))$p.value<0.01){
    a <- a[-which.max(a)]
  }
  
  if(!include.zeroes){
    a <- a[a!=0]
  }
  
  if (mean(log(a[a!=0])) < threshold){
    MV[1] <- mean(a)*FC
    MV[2] <- var(a)*FC^2
  }else{
    if (is.na(FC.thresh)){stop("Please specify valid value for alternate fold change FC.thresh!")}
    MV[1] <- mean(a)*FC.thresh
    MV[2] <- var(a)*FC.thresh^2
  }
  return(MV)
}

#' singleCellSimu
#'
#' Called by \code{\link{simulateSet}} to simulate a specified number of genes from 
#'   one DD category at a time.
#' 
#'
#'
#' @param Dataset1 Numeric matrix of expression values with genes in rows and samples in columns.
#'
#' @param Method Type of simulation should choose from "DE" "DP" "DM" "DB"
#'
#' @param index Reasonable set of genes for simulation
#'
#' @param FC Fold Change values for DE Simulation
#'
#' @param DP Differetial Proportion vector
#' 
#' @param modeFC Vector of values to use for fold changes between modes for DP, DM, and DB.
#'
#' @param Validation Show Validation plots or not
#' 
#' @param numGenes numeric value for the number of genes to simulate
#'
#' @param numDE numeric value for the number of genes that will differ between two conditions
#'  
#' @param generateZero Specification of how to generate the zero values.  If "\code{empirical}" (default), the
#'  observed proportion of zeroes in each gene is used for the simuated data, and the nonzeroes are simulated from a 
#'  truncated negative binomial distribution.  If "\code{simulated}", all values are simulated out of
#'  a negative binomial distribution, includling the zeroes.  If "\code{constant}", then each gene has a fixed 
#'  proportion of zeroes equal to \code{constantZero}.
#'  
#' @param constantZero Numeric value between 0 and 1 that indicates the fixed proportion of zeroes for every gene.
#'  Ignored if \code{generateZero} method is not equal to "\code{constant}". 
#' 
#' @param numSamples numeric value for the number of samples in each condition to simulate
#' 
#' @importFrom outliers grubbs.test
#'
#' @return Simulated_Data A list object where the first element contains a matrix of 
#' the simulated dataset, the second element contains the DEIndex, and the third 
#' element contains the fold change (between two conditions for DE, between two modes 
#' for DP, DM, and DB).

singleCellSimu <- function(Dataset1, Method, index, FC, modeFC, DP, Validation=FALSE, 
                             numGenes=1000, numDE=100, numSamples=100, 
                             generateZero=c("empirical", "simulated", "constant"),
                             constantZero=NULL){
  
if(length(generateZero)>1){
  generateZero <- generateZero[1]
}
### Select numGenes genes from index
samplename <- sort(sample(index,numGenes,replace=TRUE))

FC2 <- rep(1,nrow=numGenes)

f <- sample(c(rep(modeFC, ceiling(numGenes/length(modeFC)))), numGenes, replace=FALSE)
for(i in 1:numGenes){
  a <- as.numeric(Dataset1[samplename[i],])
  FC2[i] <- exp(f[i]*sqrt(log(1+var(a[a!=0])/mean(a[a!=0])^2)))
}

### Calculate means and variances for these genes
Zeropercent_Base <- as.matrix(apply(Dataset1[samplename,], 1, function(a) length(which(a == 0)) / length(a)))           
MV <- matrix(data = 0,nrow = numGenes,ncol = 4)
if(Method %in%  c("DP", "DM")){
  MV[,1:2] <- t(sapply(1:length(samplename), function(x) calcMV(Dataset1[samplename[x],], FC=1, FC.thresh=FC2[x]^(-1/2), 
                                                              threshold = 3, include.zeroes=FALSE)))
  MV[,3:4] <- t(sapply(1:length(samplename), function(x) calcMV(Dataset1[samplename[x],], FC=1, FC.thresh=FC2[x]^(-1/2), 
                                                              threshold = 3, include.zeroes=TRUE)))
}else{
  MV[,1:2] <- t(sapply(1:length(samplename), function(x) calcMV(Dataset1[samplename[x],], FC=1, include.zeroes=FALSE)))
  MV[,3:4] <- t(sapply(1:length(samplename), function(x) calcMV(Dataset1[samplename[x],], FC=1, include.zeroes=TRUE)))
}


### Calculate parameter R and P in NB distribution
if(generateZero %in% c("empirical", "constant")){
  RP <- t(apply(MV[,1:2], 1, function(x) calcRP(x[1], x[2])))
}else if (generateZero == "simulated"){
  RP <- t(apply(MV[,3:4], 1, function(x) calcRP(x[1], x[2])))
}else{
  stop("Error: Please specify a valid generateZero method!")
}


### Indentify the relationship between mean and variance
if(generateZero %in% c("empirical", "constant")){
  fit <- lm(log(MV[,2])~log(MV[,1]))
} else{
  fit <- lm(log(MV[,4])~log(MV[,3]))
}
coeff <- fit$coefficients


### Simulate data by the chosen method
Simulated_Data <- matrix(data=NA,nrow=numGenes,ncol=2*numSamples)
x <- 1:numGenes
DEIndex <- sample(x,numDE,replace=FALSE)
EEIndex <- x[!(x%in%DEIndex)]


if(Method=="DE")
{
	desim <- simuDE(Dataset1, Simulated_Data, DEIndex, samplename, Zeropercent_Base, f, FC, coeff, RP, modeFC,
	                         generateZero, constantZero)
	Simulated_Data <- desim[[1]]
  DE_FC <- desim[[2]]
} else if(Method=="DP"){
  dpsim <- simuDP(Dataset1, Simulated_Data, DEIndex, samplename, Zeropercent_Base, f, FC2, coeff, RP, modeFC, DP,
                  generateZero, constantZero)
  Simulated_Data <- dpsim[[1]]
  DE_FC <- dpsim[[2]]
}else if(Method=="DM"){
  
  dmsim <- simuDM(Dataset1, Simulated_Data, DEIndex, samplename, Zeropercent_Base, f, FC2, coeff, RP, modeFC,
                  generateZero, constantZero)
  Simulated_Data <- dmsim[[1]]
  DE_FC <- dmsim[[2]]
}else if(Method=="DB"){
  
  DBsim <- simuDB(Dataset1, Simulated_Data, DEIndex, samplename, Zeropercent_Base, f, FC2, coeff, RP, modeFC, DP,
                      generateZero, constantZero)
  Simulated_Data <- DBsim[[1]]
  DE_FC <- DBsim[[2]]
}

if(Validation==TRUE){
  validation(MV, DEIndex, Zeropercent_Base, Simulated_Data, numGenes)
}
return(list(Simulated_Data, DEIndex, DE_FC))
}
