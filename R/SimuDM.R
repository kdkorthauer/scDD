#' simuDM
#'
#' Simulation for Differential Modalities Case
#'
#' @inheritParams singleCellSimu
#' @inheritParams simuDE
#'
#' @references Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, Kendziorski C. A statistical approach for identifying differential distributions
#' in single-cell RNA-seq experiments. Genome Biology. 2016 Oct 25;17(1):222. \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1077-y}
#'
#' @return Simulated_Data Simulated dataset for DM


simuDM <- function(Dataset1, Simulated_Data, DEIndex, samplename, Zeropercent_Base, f, FC, coeff, RP, modeFC,
                   generateZero, constantZero, varInflation){

  numGenes <- nrow(Simulated_Data)
  numSamples <- ncol(Simulated_Data)/2
  numDE <- length(DEIndex)
  
  # calculate means and variances with fold changes
  MVC2 <- matrix(data=0,nrow=numGenes,ncol=2)
  if(generateZero %in% c("empirical", "constant")){
    MVC2[,1:2] <- t(sapply(1:length(samplename), function(x) calcMV(Dataset1[samplename[x],], FC=FC[x], FC.thresh=FC[x]^(1/2), 
                                                                    threshold = 3, include.zeroes=FALSE)))
  }else{
    MVC2[,1:2] <- t(sapply(1:length(samplename), function(x) calcMV(Dataset1[samplename[x],], FC=FC[x], FC.thresh=FC[x]^(1/2), 
                                                                    threshold = 3, include.zeroes=TRUE)))
  }
  
  # calculate r and p parameters of NB
  RPC2 <- t(apply(MVC2[,1:2,drop=FALSE], 1, function(x) calcRP(x[1], x[2])))
  
  # calculate r and p parameters of inflated variance NB and add to the RPC2 matrix
  if(!is.null(varInflation)){
    MVC2.Infl1 <- MVC2
    MVC2.Infl1[,1] <- MVC2[,1]*(1 + MVC2[,2]/(MVC2[,1]^2))^((varInflation[1]-1)/2) 
    MVC2.Infl1[,2] <- MVC2[,2]*( ((1 + MVC2[,2]/(MVC2[,1]^2))^(2*varInflation[1])-(1 + MVC2[,2]/(MVC2[,1]^2))^varInflation[1]) / 
                                   ((1 + MVC2[,2]/(MVC2[,1]^2))^2-(1 + MVC2[,2]/(MVC2[,1]^2))) )
    
    MVC2.Infl2 <- MVC2
    MVC2.Infl2[,1] <- MVC2[,1]*(1 + MVC2[,2]/(MVC2[,1]^2))^((varInflation[2]-1)/2) 
    MVC2.Infl2[,2] <- MVC2[,2]*( ((1 + MVC2[,2]/(MVC2[,1]^2))^(2*varInflation[2])-(1 + MVC2[,2]/(MVC2[,1]^2))^varInflation[2]) / 
                                   ((1 + MVC2[,2]/(MVC2[,1]^2))^2-(1 + MVC2[,2]/(MVC2[,1]^2))) )
    
    RPC2 <- cbind(RPC2, t(apply(MVC2.Infl1[,1:2,drop=FALSE], 1, function(x) calcRP(x[1], x[2]))))
    RPC2 <- cbind(RPC2, t(apply(MVC2.Infl2[,1:2,drop=FALSE], 1, function(x) calcRP(x[1], x[2]))))
  }
  
  dp <- c(1,0.5)
  
  
  ### Simulate data in condition1
  for(i in 1:numGenes){
    if(generateZero=="empirical"){
      p <- round(numSamples*(1-Zeropercent_Base[i,1]))
    } else if (generateZero=="constant"){
      p <- round((1-constantZero)*numSamples)
    }
    
    if(generateZero %in% c("empirical", "constant")){
      temp <- matrix(data=0,nrow=p,ncol=2)
      if(is.null(varInflation)){
        for(j in 1:p){
          while(temp[j,1]==0){temp[j,1] <- rnbinom(1,RP[i,1],1-RP[i,2])}
          while(temp[j,2]==0){temp[j,2] <- rnbinom(1,RPC2[i,1],1-RPC2[i,2])}
        }
      }else{
        for(j in 1:p){
          while(temp[j,1]==0){temp[j,1] <- rnbinom(1,RP[i,3],1-RP[i,4])}
          while(temp[j,2]==0){temp[j,2] <- rnbinom(1,RPC2[i,3],1-RPC2[i,4])}
        }
      }

      
      # place nonzeroes randomly
      n1 <- round(p*dp[1],digits=0)
      randp <- sample(1:numSamples, p, replace=FALSE)
      Simulated_Data[i,randp] <-  c(sample(temp[,1],n1,replace=FALSE), 
                                    sample(temp[,2],p-n1,replace=FALSE))
      if(p<numSamples){
        Simulated_Data[i,(1:numSamples)[-randp]] <-0
      }
    }else{
      temp <- matrix(data=0,nrow=numSamples,ncol=2)
      if(is.null(varInflation)){
        temp[,1] <- rnbinom(numSamples,RP[i,1],1-RP[i,2])
        temp[,2] <- rnbinom(numSamples,RPC2[i,1],1-RPC2[i,2])
      }else{
        temp[,1] <- rnbinom(numSamples,RP[i,3],1-RP[i,4])
        temp[,2] <- rnbinom(numSamples,RPC2[i,3],1-RPC2[i,4])
      }
      
      # sample from each mode
      n1 <- round(numSamples*dp[1],digits=0)
      Simulated_Data[i,1:numSamples] <-  c(sample(temp[,1],n1,replace=FALSE), 
                                           sample(temp[,2],numSamples-n1,replace=FALSE))
    }

    ### Simulate for condition2
    if(generateZero=="empirical"){
      p <- round(numSamples*(1-Zeropercent_Base[i,1]))
    } else if (generateZero=="constant"){
      p <- round((1-constantZero)*numSamples)
    }
    
    if(generateZero %in% c("empirical", "constant")){  
      if(i%in%DEIndex){
        temp <- matrix(data=0,nrow=p,ncol=2)
        if(is.null(varInflation)){
          for(j in 1:p){
            while(temp[j,1]==0){temp[j,1] <- rnbinom(1,RP[i,1],1-RP[i,2])}
            while(temp[j,2]==0){temp[j,2] <- rnbinom(1,RPC2[i,1],1-RPC2[i,2])}
          }
        }else{
          for(j in 1:p){
            while(temp[j,1]==0){temp[j,1] <- rnbinom(1,RP[i,5],1-RP[i,6])}
            while(temp[j,2]==0){temp[j,2] <- rnbinom(1,RPC2[i,5],1-RPC2[i,6])}
          }
        }
        n1 <- round(p*dp[2],digits=0)
        randp <- sample(1:numSamples, p, replace=FALSE)
        Simulated_Data[i,numSamples+randp] <-  c(sample(temp[,1],n1,replace=FALSE), 
                                                 sample(temp[,2],p-n1,replace=FALSE))
        if(p<numSamples){
          Simulated_Data[i,((numSamples+1):(2*numSamples))[-randp]] <-0
        }
      }else{
        temp <- matrix(data=0,nrow=p,ncol=2)
        if(is.null(varInflation)){
          for(j in 1:p){
            while(temp[j,1]==0){temp[j,1] <- rnbinom(1,RP[i,1],1-RP[i,2])}
            while(temp[j,2]==0){temp[j,2] <- rnbinom(1,RPC2[i,1],1-RPC2[i,2])}
          }
        }else{
          for(j in 1:p){
            while(temp[j,1]==0){temp[j,1] <- rnbinom(1,RP[i,5],1-RP[i,6])}
            while(temp[j,2]==0){temp[j,2] <- rnbinom(1,RPC2[i,5],1-RPC2[i,6])}
          }
        }
        n1 <- round(p*dp[1],digits=0)
        randp <- sample(1:numSamples, p, replace=FALSE)
        Simulated_Data[i,numSamples+randp] <-  c(sample(temp[,1],n1,replace=FALSE), 
                                                 sample(temp[,2],p-n1,replace=FALSE))
        if(p<numSamples){
          Simulated_Data[i,((numSamples+1):(2*numSamples))[-randp]] <-0
        }
      }
    }else{
      if(i%in%DEIndex){
        temp <- matrix(data=0,nrow=numSamples,ncol=2)
        if(is.null(varInflation)){
          temp[,1] <- rnbinom(numSamples,RP[i,1],1-RP[i,2])
          temp[,2] <- rnbinom(numSamples,RPC2[i,1],1-RPC2[i,2])
        }else{
          temp[,1] <- rnbinom(numSamples,RP[i,5],1-RP[i,6])
          temp[,2] <- rnbinom(numSamples,RPC2[i,5],1-RPC2[i,6])
        }
        n1 <- round(numSamples*dp[2],digits=0)
        Simulated_Data[i,((numSamples+1):(2*numSamples))] <-  c(sample(temp[,1],n1,replace=FALSE), 
                                                                sample(temp[,2],numSamples-n1,replace=FALSE))
      }else{
        temp <- matrix(data=0,nrow=numSamples,ncol=2)
        if(is.null(varInflation)){
          temp[,1] <- rnbinom(numSamples,RP[i,1],1-RP[i,2])
          temp[,2] <- rnbinom(numSamples,RPC2[i,1],1-RPC2[i,2])
        }else{
          temp[,1] <- rnbinom(numSamples,RP[i,5],1-RP[i,6])
          temp[,2] <- rnbinom(numSamples,RPC2[i,5],1-RPC2[i,6])
        }
        n1 <- round(numSamples*dp[1],digits=0)
        Simulated_Data[i,((numSamples+1):(2*numSamples))] <-  c(sample(temp[,1],n1,replace=FALSE), 
                                                                sample(temp[,2],numSamples-n1,replace=FALSE))
      }
    }
  }
  return(list(Simulated_Data, f=f))
}
