 #' simuDB
#'
#' Simulation for Differential "Both" Case - both Differential Modality and Differential Mean
#'
#' @inheritParams singleCellSimu 
#' @inheritParams simuDE
#' 
#' @references Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, Kendziorski C. A statistical approach for identifying differential distributions
#' in single-cell RNA-seq experiments. Genome Biology. 2016 Oct 25;17(1):222. \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1077-y}
#'
#' @return Simulated_Data Simulated dataset for DB




simuDB <- function(Dataset1, Simulated_Data, DEIndex, samplename, Zeropercent_Base, f, FC, coeff, RP, modeFC, DP,
                     generateZero, constantZero, varInflation){
  
  numGenes <- nrow(Simulated_Data)
  numSamples <- ncol(Simulated_Data)/2
  numDE <- length(DEIndex)
  
  ### Simulate data for EE and condition 1 in DE first
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
      randp <- sample(1:numSamples, p, replace=FALSE)
      Simulated_Data[i,randp] <- temp
      if(p<numSamples){
        Simulated_Data[i,(1:numSamples)[-randp]] <-0
      }
    }else{
      if(is.null(varInflation)){
        Simulated_Data[i,1:numSamples] <- rnbinom(numSamples,RP[i,1],1-RP[i,2])
      }else{
        Simulated_Data[i,1:numSamples] <- rnbinom(numSamples,RP[i,3],1-RP[i,4])
      }
    }
  }
  
  
  cutoff <- -Inf  # optional cutoff to specify which will influence direction of second mode

  # calculate means and variances with fold changes
  MVC2 <- matrix(data=0,nrow=numGenes,ncol=5)
  if(generateZero %in% c("empirical", "constant")){
    MVC2[,1:2] <- t(sapply(1:length(samplename), function(x) calcMV(Dataset1[samplename[x],], FC=FC[x]^(1/2), FC.thresh=FC[x]^(-1/2), 
                                                                  threshold = cutoff, include.zeroes=FALSE)))
    MVC2[,3:4] <- t(sapply(1:length(samplename), function(x) calcMV(Dataset1[samplename[x],], FC=FC[x], FC.thresh=FC[x]^(1/2), 
                                                                  threshold = cutoff, include.zeroes=FALSE)))
  }else{
    MVC2[,1:2] <- t(sapply(1:length(samplename), function(x) calcMV(Dataset1[samplename[x],], FC=FC[x]^(1/2), FC.thresh=FC[x]^(-1/2), 
                                                                    threshold = cutoff, include.zeroes=TRUE)))
    MVC2[,3:4] <- t(sapply(1:length(samplename), function(x) calcMV(Dataset1[samplename[x],], FC=FC[x], FC.thresh=FC[x]^(1/2), 
                                                                    threshold = cutoff, include.zeroes=TRUE)))
  }
  MVC2[,5] <- apply(Dataset1[samplename,,drop=FALSE], 1, function(x) sum(log(mean(x[x>0]))>=cutoff))
  
  
  # calculate r and p parameters of NB
  RPC2 <- cbind( t(apply(MVC2[,1:2,drop=FALSE], 1, function(x) calcRP(x[1], x[2]))),
                 t(apply(MVC2[,3:4,drop=FALSE], 1, function(x) calcRP(x[1], x[2]))))
  
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
    
    RPC2 <- cbind(RPC2, cbind( t(apply(MVC2[,1:2,drop=FALSE], 1, function(x) calcRP(x[1], x[2]))),
                               t(apply(MVC2[,3:4,drop=FALSE], 1, function(x) calcRP(x[1], x[2])))))
    RPC2 <- cbind(RPC2, cbind( t(apply(MVC2[,1:2,drop=FALSE], 1, function(x) calcRP(x[1], x[2]))),
                               t(apply(MVC2[,3:4,drop=FALSE], 1, function(x) calcRP(x[1], x[2])))))
  }
  
  ### Simulate data in condition 1
  for(i in 1:numGenes){
    if(generateZero=="empirical"){
      p <- round(numSamples*(1-Zeropercent_Base[i,1]))
    } else if (generateZero=="constant"){
      p <- round((1-constantZero)*numSamples)
    }
    
    if(generateZero %in% c("empirical", "constant")){
      temp <- matrix(data=0,nrow=p,ncol=1)
      if(MVC2[i,5]==0){  # gene is below cutoff
        if(is.null(varInflation)){
          for(j in 1:p){
            while(temp[j]==0){temp[j] <- rnbinom(1,RPC2[i,1],1-RPC2[i,2])}
          }
        }else{
          for(j in 1:p){
            while(temp[j]==0){temp[j] <- rnbinom(1,RPC2[i,5],1-RPC2[i,6])}
          } 
        }
        randp <- sample(1:numSamples, p, replace=FALSE)
        Simulated_Data[i,randp] <- temp
        if(p<numSamples){
          Simulated_Data[i,(1:numSamples)[-randp]] <-0
        }
      }else{  # gene is above cutoff
        if(is.null(varInflation)){
          for(j in 1:p){
            while(temp[j]==0){temp[j] <- rnbinom(1,RP[i,1],1-RP[i,2])}
          }
        }else{
          for(j in 1:p){
            while(temp[j]==0){temp[j] <- rnbinom(1,RP[i,3],1-RP[i,4])}
          }
        }
        randp <- sample(1:numSamples, p, replace=FALSE)
        Simulated_Data[i,randp] <- temp
        if(p<numSamples){
          Simulated_Data[i,(1:numSamples)[-randp]] <-0
        }
      }
    } else{
      if(MVC2[i,5]==0){  # gene is below cutoff 
        if(is.null(varInflation)){
          Simulated_Data[i,1:numSamples] <- rnbinom(numSamples,RPC2[i,1],1-RPC2[i,2])
        }else{
          Simulated_Data[i,1:numSamples] <- rnbinom(numSamples,RPC2[i,5],1-RPC2[i,6])
        }
      }else{  # gene is above cutoff
        if(is.null(varInflation)){
          Simulated_Data[i,1:numSamples] <- rnbinom(numSamples,RP[i,1],1-RP[i,2])
        }else{
          Simulated_Data[i,1:numSamples] <- rnbinom(numSamples,RP[i,3],1-RP[i,4])
        }
      }
    }
  }

  
  ### Simulate for condition2
  for(i in 1:numGenes){
    if(generateZero=="empirical"){
      p <- round(numSamples*(1-Zeropercent_Base[i,1]))
    } else if (generateZero=="constant"){
      p <- round((1-constantZero)*numSamples)
    }
    
    if(generateZero %in% c("empirical", "constant")){
      if(i%in%DEIndex){ 
        temp <- matrix(data=0,nrow=p,ncol=2)
        if(MVC2[i,5]==0){ # gene is below cutoff
          if(is.null(varInflation)){
            for(j in 1:p){
              while(temp[j,1]==0){temp[j,1] <- rnbinom(1,RP[i,1],1-RP[i,2])}
              while(temp[j,2]==0){temp[j,2] <- rnbinom(1,RPC2[i,3],1-RPC2[i,4])}
            }
          }else{
            for(j in 1:p){
              while(temp[j,1]==0){temp[j,1] <- rnbinom(1,RP[i,5],1-RP[i,6])}
              while(temp[j,2]==0){temp[j,2] <- rnbinom(1,RPC2[i,11],1-RPC2[i,12])}
            }
          }
        }else{  # gene is above cutoff
          if(is.null(varInflation)){
            for(j in 1:p){
              while(temp[j,1]==0){temp[j,1] <- rnbinom(1,RPC2[i,1],1-RPC2[i,2])}
              while(temp[j,2]==0){temp[j,2] <- rnbinom(1,RPC2[i,3],1-RPC2[i,4])}
            }
          }else{
            for(j in 1:p){
              while(temp[j,1]==0){temp[j,1] <- rnbinom(1,RPC2[i,9],1-RPC2[i,10])}
              while(temp[j,2]==0){temp[j,2] <- rnbinom(1,RPC2[i,11],1-RPC2[i,12])}
            }
          }
        }
        
        n1 <- round(p*DP[1],digits=0)
        randp <- sample(1:numSamples, p, replace=FALSE)
        Simulated_Data[i,numSamples+randp] <-  c(sample(temp[,1],n1,replace=FALSE), 
                                                 sample(temp[,2],p-n1,replace=FALSE))
        if(p<numSamples)
        {Simulated_Data[i,((numSamples+1):(2*numSamples))[-randp]] <-0}
        
      }else{  #EE genes (simulate just like in C1)
        temp <- matrix(data=0,nrow=p,ncol=1)
        if(MVC2[i,5]==0){ # gene is below cutoff
          if(is.null(varInflation)){
            for(j in 1:p){
              while(temp[j,1]==0){temp[j,1] <- rnbinom(1,RPC2[i,1],1-RPC2[i,2])}
            }
          }else{
            for(j in 1:p){
              while(temp[j,1]==0){temp[j,1] <- rnbinom(1,RPC2[i,9],1-RPC2[i,10])}
            }
          }
        }else{ # gene is above
          if(is.null(varInflation)){
            for(j in 1:p){
              while(temp[j,1]==0){temp[j,1] <- rnbinom(1,RP[i,1],1-RP[i,2])}
            }
          }else{
            for(j in 1:p){
              while(temp[j,1]==0){temp[j,1] <- rnbinom(1,RP[i,5],1-RP[i,6])}
            }
          }
        }
        
        randp <- sample(1:numSamples, p, replace=FALSE)
        Simulated_Data[i,randp+numSamples] <- temp
        if(p<numSamples){
          Simulated_Data[i,((numSamples+1):(2*numSamples))[-randp]] <-0
        }
      }
    }else{
      if(i%in%DEIndex){ 
        temp <- matrix(data=0,nrow=numSamples,ncol=2)
        if(MVC2[i,5]==0){ # gene is below cutoff
          if(is.null(varInflation)){
            temp[,1] <- rnbinom(numSamples,RP[i,1],1-RP[i,2])
            temp[,2] <- rnbinom(numSamples,RPC2[i,3],1-RPC2[i,4])
          }else{
            temp[,1] <- rnbinom(numSamples,RP[i,5],1-RP[i,6])
            temp[,2] <- rnbinom(numSamples,RPC2[i,11],1-RPC2[i,12])
          }
        }else{  # gene is above cutoff
          if(is.null(varInflation)){
            temp[,1] <- rnbinom(numSamples,RPC2[i,1],1-RPC2[i,2])
            temp[,2] <- rnbinom(numSamples,RPC2[i,3],1-RPC2[i,4])
          }else{
            temp[,1] <- rnbinom(numSamples,RPC2[i,9],1-RPC2[i,10])
            temp[,2] <- rnbinom(numSamples,RPC2[i,11],1-RPC2[i,12])
          }
        }
        
        n1 <- round(numSamples*DP[1],digits=0)
        Simulated_Data[i,((numSamples+1):(2*numSamples))] <-  c(sample(temp[,1],n1,replace=FALSE), 
                                                                sample(temp[,2],numSamples-n1,replace=FALSE))
        
      }else{  #EE genes (simulate just like in C1)
        if(MVC2[i,5]==0){ # gene is below cutoff
          if(is.null(varInflation)){
            temp <- rnbinom(numSamples,RPC2[i,1],1-RPC2[i,2])
          }else{
            temp <- rnbinom(numSamples,RPC2[i,9],1-RPC2[i,10])
          }
        }else{ # gene is above cutoff
          if(is.null(varInflation)){
            temp <- rnbinom(numSamples,RP[i,1],1-RP[i,2])
          }else{
            temp <- rnbinom(numSamples,RP[i,5],1-RP[i,6])
          }
        }
      }
    }
  }
  return(list(Simulated_Data, f=f))
}
