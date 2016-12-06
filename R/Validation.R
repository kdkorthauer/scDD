#' validation
#'
#' Draw validation plots to show that the simulated dataset emulates characteristics of observed dataset.
#'
#' @inheritParams singleCellSimu
#'
#' @param MV Mean and Variance matrix for observed data
#'
#' @param DEIndex Index for genes chosen to be DE (can be NULL)
#'
#' @param Zeropercent_Base Zero percentage for corresponding gene expression values 
#'
#' @param Simulated_Data Simulated dataset
#' 
#' @references Korthauer KD, Chu LF, Newton MA, Li Y, Thomson J, Stewart R, Kendziorski C. A statistical approach for identifying differential distributions
#' in single-cell RNA-seq experiments. Genome Biology. 2016 Oct 25;17(1):222. \url{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1077-y}
#'
#' @return Validation plots


validation <- function(MV, DEIndex, Zeropercent_Base, Simulated_Data, numGenes){
MV_Simu <- matrix(data=0,nrow=numGenes,ncol=4)

MV_Simu[,1]=apply(Simulated_Data, 1, function(x) mean(x[x!=0]))
MV_Simu[,2]=apply(Simulated_Data, 1, function(x) var(x[x!=0]))
MV_Simu[,3]=apply(Simulated_Data, 1, mean)
MV_Simu[,4]=apply(Simulated_Data, 1, var)

par(mfrow=c(1,1))
plot(MV[,1], MV[,2], log="xy", main="Mean-Variance Relationship", xlab="Nonzero Mean", ylab="Nonzero Variance",
    pch=20, cex=0.75, col="darkgrey")


par(mfrow=c(2,2))
plot(MV[,1],MV_Simu[,1],log="xy",main="Nonzero Mean", xlab="Empirical", ylab="Simulated",
     pch=20, cex=0.25, col="darkgrey")
abline(0,1,lty=2)
if(length(DEIndex>0)){ 
  points(MV[DEIndex,1],MV_Simu[DEIndex,1],col="red") 
}

plot(MV[,2],MV_Simu[,2],log="xy",main="Nonzero Variance", xlab="Empirical", ylab="Simulated",
     pch=20, cex=0.25, col="darkgrey")
abline(0,1,lty=2)
if(length(DEIndex>0)){ 
  points(MV[DEIndex,2],MV_Simu[DEIndex,2],col="red")
}

plot(MV[,3],MV_Simu[,3],log="xy",main="Overall Mean", xlab="Empirical", ylab="Simulated",
     pch=20, cex=0.25, col="darkgrey")
abline(0,1,lty=2)
if(length(DEIndex>0)){ 
  points(MV[DEIndex,3],MV_Simu[DEIndex,2],col="red") 
}

plot(MV[,4],MV_Simu[,4],log="xy",main="Overall Variance", xlab="Empirical", ylab="Simulated",
     pch=20, cex=0.25, col="darkgrey")
abline(0,1,lty=2)
if(length(DEIndex>0)){ 
  points(MV[DEIndex,4],MV_Simu[DEIndex,4],col="red")
}

}
