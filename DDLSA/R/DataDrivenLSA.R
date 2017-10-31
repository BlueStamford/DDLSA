################################
# DataDrivenLSA.R  #
################################
#
#
# The following functions are implemented:
#
#	. approSignificance(x,y,maxDelay)
#
#	. DataDrivenLSA(x,y,maxDelay)
#
# Last updated: Jan 30, 2017, Fang Zhang.  @All rights reserved.
#############################################################################
################################
# approSignificance(x,y,maxDelay):
#
#	This function computes the approximated statistical significancelocal of
# local similarity score with maximal delay maxDelay.
#
################################
theoDistribution<-function(localScore,Length,Variance,maxDelay){
  normalizedLocalScore <- localScore/sqrt(Variance*Length)
  if (normalizedLocalScore == 0) {
    theoAppro <- 1
    } else{
      partialSum <- rep(0,1000)
      partialSum[1] <- (1/normalizedLocalScore^2 + 1/(pi^2)) * exp(-pi^2/(2 * normalizedLocalScore^2))
      i <- 2
      threshold <- 1
      while (threshold > 1e-5){
        threshold <- (1/normalizedLocalScore^2 + 1/((2*i - 1)^2*pi^2)) * exp(-(2*i-1)^2 * pi^2/(2*normalizedLocalScore^2))
        partialSum[i]<-partialSum[i-1] + threshold
        i<-i+1
      }
      theoAppro<-1-8^(2*maxDelay+1)*max(partialSum)^(2*maxDelay+1)
    }
  return(theoAppro)
}
##################################################################
approSignificance<-function(x,y,maxDelay){
  localScore <- LocalSimilarity(x,y,maxDelay,rankScale= TRUE)[1]
  theoAppro <- theoDistribution(localScore, length(x),1,maxDelay)
  return(theoAppro)
}
##################################################################
# DataDrivenLSA(x,y,maxDelay)
#
#  This function computes the data-driven statistical significance of
#  local similarity for two sequences.
#
# INPUT:
# ======
#
#	x, y	: sequences to copute LS score
#	maxDelay	: maximum time shift allowed in computing LS score.
#	scale		: If TRUE, perform normalization first; False, otherwise.
#
# RETURN:
# ======
#
#  A seven element vector contains: c(Maxscore, approximated p-value,
#                                     startX, startY, delay, length, PosOrNeg)
#
Omega <- function(x,y){
  timepoints <- length(x)
  z <- x * y
  residual <- z - mean(z)
  rhohat <- sum(residual[2:timepoints]*residual[1:(timepoints-1)])/sum((residual[1:(timepoints-1)])^2)
  alphahat <- 4*rhohat^2/(1-rhohat^2)^2
  bandwidth <- ceiling(1.1447*(alphahat*timepoints)^(1/3))
  covx <- acf(x,type="covariance",plot=FALSE,lag.max=timepoints-1)[[1]][,1,1]
  covy <- acf(y,type="covariance",plot=FALSE,lag.max=timepoints-1)[[1]][,1,1]
  covz <- covx*covy
  if (!is.na(bandwidth)){
    if (bandwidth > 1) omega <- sum(c(covz[bandwidth:1],covz[2:bandwidth])*(1-abs(seq(-bandwidth+1,bandwidth-1)/bandwidth)))
    else omega <- covz[1]
      }
  else omega <- covz[1]
  return(omega)
}

DataDrivenLSA<-function(x,y,maxDelay,scale = TRUE){

  if(scale == TRUE)
  {
    x<-scale(x)
    y<-scale(y)
  }

  timepoints <- length(x)
  originalSimilarity<- LocalSimilarity(x, y, maxDelay)
  localscore <- originalSimilarity[1]
  xOmegay <- Omega(x,y)
  variance<-var(x)*var(y)
  if (is.na(xOmegay)) {
     approximation <- theoDistribution(localscore,timepoints,variance,maxDelay)
     }else {
       approximation<-theoDistribution(localscore,timepoints,xOmegay,maxDelay)
     }
  value <- t(c(localscore, approximation, originalSimilarity[2:6]))
  colnames(value)<-c('Maxscore','approximated p-value', 'startX', 'startY', 'delay', 'length', 'PosOrNeg')
  return(value)
}