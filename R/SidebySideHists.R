#' sideHist
#'
#' Plots two histograms side by side with smoothed density overlay
#'
#' 
#' @param x First numeric vector of data to plot.
#'
#' @param y Second numeric vector of data to plot.
#'
#' @param logT Logical that indicates whether to take the log(x+1) transformation.
#' 
#' @param title.gene Character vector that contains the gene name that you are plotting
#'
#'

sideHist <- function(x, y, logT=TRUE, title.gene=""){

# grab some pretty ggplot colors
ggplotColours <- function(n=6, h=c(0, 360) +15){
  if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

valsA <- y
valsB <- x

if (logT){
  valsA <- log(valsA+1)
  valsB <- log(valsB+1)
}

binwidth <- 0.35
bin.points <- seq(range(c(valsA,valsB))[1]-0.5, range(c(valsA,valsB))[2]+0.5, by=binwidth)

A <- hist(valsA, breaks=bin.points, plot=FALSE)
B <- hist(valsB, breaks=bin.points, plot=FALSE)

gcol <- ggplotColours(2)

densrange <-  c(-max(c(B$density, A$density)), max(c(A$density, B$density)))

plot(NULL, type = "n", xlim = densrange,
     ylim = c(range(c(A$breaks, B$breaks))), ylab="log(EC + 1)", xlab="density", xaxt="n",
     main = paste0("Gene ", title.gene))
axis(1, at=seq(-2,2,by=0.1), labels=abs(seq(-2,2,by=0.1)))
rect(0, A$breaks[1:(length(A$breaks) - 1)], A$density, A$breaks[2:length(A$breaks)], col=gcol[1])
rect(0, B$breaks[1:(length(B$breaks) - 1)], -B$density, B$breaks[2:length(B$breaks)], col=gcol[2])
lines(density(valsA, adjust=1.25)$y, density(valsA, adjust=1.25)$x, col="darkgrey", lwd=2)
lines(-density(valsB, adjust=1.25)$y, density(valsB, adjust=1.25)$x, col="darkgrey", lwd=2)

}
