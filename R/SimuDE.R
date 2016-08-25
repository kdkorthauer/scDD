#' simuDE
#'
#' Simulation for Classic Differentially Expressed Case.
#' 
#' @details Method called by main function \code{\link{singleCellSimu}} to simulate genes 
#'  that have different means in each condition. Not intended to be called directly by user.
#'
#' @inheritParams singleCellSimu
#'
#' @param Simulated_Data Required input empty matrix to provide structure information of 
#'  output matrix with simulated data
#'
#' @param DEIndex Index for DE genes
#'
#' @param samplename The name for genes that chosen for simulation
#'
#' @param Zeropercent_Base Zero percentage for corresponding gene expression values 
#'
#' @param coeff Relationship coefficients for Mean and Variance
#'
#' @param RP matrix for NB parameters for genes in samplename
#' 
#' @param f Fold change values (number of SDs) for each gene 
#'
#' @return Simulated_Data Simulated dataset for DE

simuDE <- function(Dataset1, Simulated_Data, DEIndex, samplename, Zeropercent_Base, f, FC, coeff, RP, modeFC,
                   generateZero, constantZero, varInflation){
  
numGenes <- nrow(Simulated_Data)
numSamples <- ncol(Simulated_Data)/2
numDE <- length(DEIndex)

### Simulate data for condition 1 first
for(i in 1:numGenes){
  if(generateZero=="empirical"){
    p <- round(numSamples*(1-Zeropercent_Base[i,1]))
  } else if (generateZero=="constant"){
    p <- round((1-constantZero)*numSamples)
  }
  
  if(generateZero %in% c("empirical", "constant")){
    temp <- matrix(data=0,nrow=p,ncol=1)
    if(is.null(varInflation)){
      for(j in 1:p){
        while(temp[j]==0){temp[j] <- rnbinom(1,RP[i,1],1-RP[i,2])}
      }
    }else{
      for(j in 1:p){
        while(temp[j]==0){temp[j] <- rnbinom(1,RP[i,3],1-RP[i,4])}
      }
    }
    ## randomly place the nonzeroes
    randp <- sample(1:numSamples, p, replace=FALSE)
    Simulated_Data[i,randp] <- temp
    if(p<numSamples)
    {Simulated_Data[i,(1:numSamples)[-randp]] <-0}
  } else{
    Simulated_Data[i,1:numSamples] <- rnbinom(numSamples,RP[i,1],1-RP[i,2])
    if(!is.null(varInflation)){
      Simulated_Data[i,1:numSamples] <- rnbinom(numSamples,RP[i,3],1-RP[i,4])
    }
  }
}

### Simulate data for condition 2 in DE
FCSimu <- FC[samplename]

f <- rep(NA, length(FCSimu))
for(i in 1:numGenes){
  a <- as.numeric(Dataset1[samplename[i],])
  if(mean(log(a[a!=0])) < 3 & FCSimu[i] < 1){
    FCSimu[i] <- 1/FCSimu[i]
  }else if(mean(log(a[a!=0])) > 7 & FCSimu[i] > 1){
    FCSimu[i] <- 1/FCSimu[i]
  }

  f[i] <- log(FCSimu[i] / sqrt(log(1+var(a[a!=0])/mean(a[a!=0])^2)))
}

# calculate empirical means and variances with fold change
MVC2 <- matrix(data=0,nrow=numGenes,ncol=2)
if(generateZero %in% c("empirical", "constant")){
  MVC2[,1:2] <- t(sapply(1:length(samplename), function(x) calcMV(Dataset1[samplename[x],], FC=FCSimu[x], include.zeroes=FALSE)))
}else{
  MVC2[,1:2] <- t(sapply(1:length(samplename), function(x) calcMV(Dataset1[samplename[x],], FC=FCSimu[x], include.zeroes=TRUE)))
}

# calculate r and p parameters of NB
RPC2 <- t(apply(MVC2[,1:2], 1, function(x) calcRP(x[1], x[2])))

# m.factor = (1 + v/m^2)^((c-1)/2)
# v.factor = ( (1 + v/m^2)^(2c) - (1 + v/m^2)^c ) /  ( (1 + v/m^2)^2 - (1+m/v^2))
if(!is.null(varInflation)){
  MVC2.Infl1 <- MVC2
  MVC2.Infl1[,1] <- MVC2[,1]*(1 + MVC2[,2]/(MVC2[,1]^2))^((varInflation[1]-1)/2) 
  MVC2.Infl1[,2] <- MVC2[,2]*( ((1 + MVC2[,2]/(MVC2[,1]^2))^(2*varInflation[1])-(1 + MVC2[,2]/(MVC2[,1]^2))^varInflation[1]) / 
                               ((1 + MVC2[,2]/(MVC2[,1]^2))^2-(1 + MVC2[,2]/(MVC2[,1]^2))) )
  
  MVC2.Infl2 <- MVC2
  MVC2.Infl2[,1] <- MVC2[,1]*(1 + MVC2[,2]/(MVC2[,1]^2))^((varInflation[2]-1)/2) 
  MVC2.Infl2[,2] <- MVC2[,2]*( ((1 + MVC2[,2]/(MVC2[,1]^2))^(2*varInflation[2])-(1 + MVC2[,2]/(MVC2[,1]^2))^varInflation[2]) / 
                                 ((1 + MVC2[,2]/(MVC2[,1]^2))^2-(1 + MVC2[,2]/(MVC2[,1]^2))) )
  
  RPC2 <- cbind(RPC2, t(apply(MVC2.Infl1[,1:2], 1, function(x) calcRP(x[1], x[2]))))
  RPC2 <- cbind(RPC2, t(apply(MVC2.Infl2[,1:2], 1, function(x) calcRP(x[1], x[2]))))
}

# now simulate genes in condition 2
for(i in 1:numGenes){
  if(generateZero=="empirical"){
    p <- round(numSamples*(1-Zeropercent_Base[i,1]))
  } else if (generateZero=="constant"){
    p <- round((1-constantZero)*numSamples)
  }
  
  if(generateZero %in% c("empirical", "constant")){
    if(i%in%DEIndex){ 
      temp <- matrix(data=0,nrow=p,ncol=1)
      if(is.null(varInflation)){
        for(j in 1:p){
          while(temp[j]==0) {temp[j] <- rnbinom(1,RPC2[i,1],1-RPC2[i,2])}
        }
      }else{
        for(j in 1:p){
          while(temp[j]==0) {temp[j] <- rnbinom(1,RPC2[i,5],1-RPC2[i,6])}
        }
      }
      ## randomly place the nonzeroes
      randp <- sample(1:numSamples, p, replace=FALSE)
      Simulated_Data[i,randp+numSamples] <- temp
      if(p<numSamples)
      {Simulated_Data[i,((numSamples+1):(2*numSamples))[-randp]] <-0}
    }else{
      temp <- matrix(data=0,nrow=p,ncol=1)
      if(is.null(varInflation)){
        for(j in 1:p){
          while(temp[j]==0){ temp[j] <- rnbinom(1,RP[i,1],1-RP[i,2])}
        }
      }else{
        for(j in 1:p){
          while(temp[j]==0){ temp[j] <- rnbinom(1,RP[i,5],1-RP[i,6])}
        }
      }
      ## randomly place the nonzeroes
      randp <- sample(1:numSamples, p, replace=FALSE)
      Simulated_Data[i,randp+numSamples] <- temp
      if(p<numSamples){
        Simulated_Data[i,((numSamples+1):(2*numSamples))[-randp]] <-0
        }
    }
  }else{
    if(i%in%DEIndex){
      if(!is.null(varInflation)){
        Simulated_Data[i,((numSamples+1):(2*numSamples))] <- rnbinom(numSamples,RPC2[i,1],1-RPC2[i,2])
      }else{
        Simulated_Data[i,((numSamples+1):(2*numSamples))] <- rnbinom(numSamples,RPC2[i,5],1-RPC2[i,6])
      }
    }else{
      if(!is.null(varInflation)){
        Simulated_Data[i,((numSamples+1):(2*numSamples))] <- rnbinom(numSamples,RP[i,1],1-RP[i,2])
      }else{
        Simulated_Data[i,((numSamples+1):(2*numSamples))] <- rnbinom(numSamples,RP[i,5],1-RP[i,6])
      }
    }
  }
}
return(list(Simulated_Data, f=round(f)))
}
