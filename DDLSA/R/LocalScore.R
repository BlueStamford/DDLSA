################################
# LocalScore.R  #
################################
#
#
# The following functions are implemented:
#
#	. rankNormalization(x)
#
#	. LocalSimilarity(x, y, maxDelay=3, rankscale=FALSE)
#
# Last updated: Jan 27, 2017, Fang Zhang.  @All rights reserved.
#############################################################################
################################
# rankNormalization(x):
#
#	This function performs rank normalization of given vector
#
###############################
rankNormalization <- function(x)
{
  xLength <- length(x)
  rankScore <- rep(0,xLength)
  rankScore <- scale(qnorm(rank(x)/(xLength+1)))
  return(rankScore)
}

##################################################################
# LocalSimilarity<- function(x, y, maxDelay=3, rankscale=FALSE)
#
#  This function computes the local similarity score for two sequences.
#
# INPUT:
# ======
#
#	x, y	: sequences to copute LS score
#	maxDelay	: maximum time shift allowed in computing LS score.
#	rankscale		: If TRUE, perform rankNormalization first; False, otherwise.
#
# RETURN:
# ======
#
#  A six element vector contains: c(scoreMax, startX, startY, delay, length, PosOrNeg)
#

LocalSimilarity <- function(x, y, maxDelay=3, rankScale = FALSE){

  if (rankScale == TRUE)
  {
    x <- rankNormalization(x)
    y <- rankNormalization(y)
  }

  scoreMax <- 0;
  PosOrNeg <- 0;
  startX <- 0;
  startY <- 0;
  numTimepoints <- length(x)
  posMatrix <- matrix(0, numTimepoints+1, numTimepoints+1)
  negMatrix <- matrix(0, numTimepoints+1, numTimepoints+1)
  for (i in 1:numTimepoints){
    for (j in max(1, i-maxDelay):min(i+maxDelay, numTimepoints)){
      posMatrix[i+1,j+1] <- max(0, posMatrix[i,j] + x[i] * y[j])
      negMatrix[i+1,j+1] <- max(0, negMatrix[i,j] - x[i] * y[j])
    }
  }
  posMatrix <- posMatrix[-1,]
  posMatrix <- posMatrix[,-1]
  negMatrix <- negMatrix[-1,]
  negMatrix <- negMatrix[,-1]
  scorePosmax <- max(posMatrix)
  scoreNegmax <- max(negMatrix)
  scoreMax <- max(scorePosmax, scoreNegmax)
  if( scorePosmax > scoreNegmax) PosOrNeg = 1 else PosOrNeg = -1
  if (scorePosmax > scoreNegmax) {
    Maxposition <- which(posMatrix == scorePosmax, arr.ind = TRUE)
    }else  Maxposition <- which(negMatrix == scoreNegmax, arr.ind = TRUE)
  delay <- Maxposition[1] - Maxposition[2]
  subinternal<-NULL
  for(i in max(1, delay):min(numTimepoints,numTimepoints + delay)){
    if (scorePosmax > scoreNegmax) {
      subinternal <- c(subinternal,posMatrix[i, i-delay])
    }else  subinternal <- c(subinternal,negMatrix[i,i-delay])
  }
  if (delay>0) {
    startX <- max(1,which(subinternal[1:Maxposition[1]] == 0) + 1) + delay
  }else{
    startX <- max(1,which(subinternal[1:Maxposition[1]] == 0) + 1)
  }
  startY <- startX - delay
  lengths <- Maxposition[1] - startX + 1
  value <- t(c(scoreMax, startX, startY, delay, lengths, PosOrNeg))
  colnames(value)<-c('Maxscore', 'startX', 'startY', 'delay', 'length', 'PosOrNeg')
  return(value)
}

######################
# permutationTest(x, y, numPermu=1000, maxDelay=3,scale=TRUE)
#
# This function computes the significance level of the LS score by permutation test
#
######################
permutationTest <- function(x, y, numPermu=1000, maxDelay=3,scale=TRUE){
	scoreArray <- rep(0.0, numPermu+1);

	if (scale == TRUE){
		x <- rankNormalization(x)
		y <- rankNormalization(y)
	}
  	numTimePoints <- length(x)
	scoreMax1 <- LocalSimilarity(x, y, maxDelay)[1];
	scoreArray[1] <- scoreMax1;
	highScoreCnt <- 0;

	for(idx in 1:numPermu)
	{
    		dtt1 <- x[sample(numTimePoints)]
		scoreTmp <- LocalSimilarity(dtt1, y, maxDelay)[1];

		scoreArray[idx+1] <- scoreTmp;
		highScoreCnt <- highScoreCnt + (scoreTmp >= scoreMax1);
	}

	pValue <- 1.0 * highScoreCnt / numPermu;

	return(pValue);
}
