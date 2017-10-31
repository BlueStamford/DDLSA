\name{LocalSimilarity}
\alias{LocalSimilarity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate Local Similarity Score of Two Time Series
}
\description{
  This function computes the local similarity score for two sequences.
}
\usage{
  LocalSimilarity(x, y, maxDelay=3, rankScale = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{
    two random variables with equal length.
  }
  \item{maxDelay}{
    maximum time shift allowed in computing LS score.
  }
  \item{rankScale}{
    If TRUE, perform rankNormalization first; False, otherwise.
  }
}
\value{
  \item{scoreMax}{
    Maximum local score between x and y.
  }
  \item{startX}{
    start position in the most similar subintervel of x.
  }
  \item{startY}{
    start position in the most similar subintervel of y.
  }
  \item{delay}{
    time delay between x and y when the maximum local similarity score is obtained. If delay is zero, no one leads anothe; if delay is positive, y leads x; otherwise, x leads y.
  }
  \item{length}{
    length of the most similar subinterval of x and y.
  }
  \item{PosOrNeg}{
    positive of negtive relationship of x and y when the max local similarity score is obtained. If PosOrNeg is 1, positive; otherwise negative.
  }
}
\examples{

  # Generate two random variables
  library(MASS)
  data <- mvrnorm(100, c(0, 0), Sigma = matrix(c(1, 0.6, 0.6, 1), 2, 2))

  # Rank normalize data
  x <- rankNormalization(data[,1])
  y <- rankNormalization(data[,2])

  # compute the local similarity score
  LocalSimilarity(x, y, maxDelay = 2)

  # example of MPHM datasets
  data(F4tongue)
  F4T <- as.matrix(F4tongue)
  length<-dim(F4T)[2]
  Count<-dim(F4T)[1]
  localScoreSignificance <- matrix(0,Count,Count)
  for (i in 1:Count){
    for (j in 1:Count){
      localScoreSignificance[i, j] <- approSignificance(F4T[i,], F4T[j,],3)
    }
  }
  # example of PML datasets
  data(PML)
  P <- as.matrix(PML)
  length_PML<-dim(P)[2]
  Count_PML<-dim(P)[1]
  localScoreSignificance_PML <- matrix(0,Count_PML,Count_PML)
  for (i in 1:Count_PML){
    for (j in 1:Count_PML){
      localScoreSignificance_PML[i, j] <- approSignificance(P[i,], P[j,],3)
    }
  }

  # R packages "qvalue" are needed to be installed before implementation of following code.
  # To install "qvalue", enter:
  # source("https://bioconductor.org/biocLite.R") (try http:// if https:// URLs are not supported)
  # biocLite("qvalue")
  # library(qvalue)
  # qv <- qvalue(localScoreSignificance[lower.tri(localScoreSignificance)])
  # summary(qv)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

