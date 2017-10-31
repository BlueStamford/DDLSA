\name{DataDrivenLSA}
\alias{DataDrivenLSA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate statistical significance of local similarity score of two time series
}
\description{
  This function computes the data-driven statistical significance of local similarity for two sequences.
}
\usage{
  DataDrivenLSA(x,y,maxDelay,scale = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{
    two random variables with equal length.
  }
  \item{maxDelay}{
    maximum time shift allowed in computing LS score.
  }
  \item{scale}{
    If TRUE, perform normalization first; False, otherwise.
  }
}
\value{
  \item{scoreMax}{
    Maximum local score between x and y.
  }
  \item{approximated p-value}{
    approximated significance of local score between x and y.
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

  # compute the significance of local similarity score
  DataDrivenLSA(x, y, maxDelay = 2)

  # Another example of real data
  data(F4tongue)
  F4T <- as.matrix(F4tongue)
  length<-dim(F4T)[2]
  Count<-dim(F4T)[1]
  PairWiseSignificance <- matrix(0,Count,Count)
  for (i in 1:Count){
    for (j in 1:Count){
      PairWiseSignificance[i, j] <- DataDrivenLSA(F4T[i,], F4T[j,], 3)[2]
    }
  }

  # example of PML datasets
  data(PML)
  P <- as.matrix(PML)
  length_PML<-dim(P)[2]
  Count_PML<-dim(P)[1]
  PairWiseSignificance_PML <- matrix(0,Count_PML,Count_PML)
  for (i in 1:Count_PML){
    for (j in 1:Count_PML){
      PairWiseSignificance_PML[i, j] <- DataDrivenLSA(P[i,], P[j,], 3)[2]
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

