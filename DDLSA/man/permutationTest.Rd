\name{permutationTest}
\alias{permutationTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   permutation test of local similarity score of two time series
}
\description{
  This function computes the significance level of the LS score by permutation test.
}
\usage{
  permutationTest(x, y, numPermu=1000, maxDelay,scale=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{
    two random variables with equal length.
  }
  \item{numPermu}{
    the number of permutation test.
  }
  \item{maxDelay}{
    maximum time shift allowed in computing LS score.
  }
  \item{scale}{
    If TRUE, perform normalization first; False, otherwise.
  }
}
\value{
  This function returns the significance level(p-value) of the local similarity score between two sequence by permutation test.
}
\examples{

  # Generate two random variables
  library(MASS)
  data <- mvrnorm(100, c(0, 0), Sigma = matrix(c(1, 0.6, 0.6, 1), 2, 2))

  # Rank normalize data
  x <- rankNormalization(data[,1])
  y <- rankNormalization(data[,2])

  # compute the significance of local similarity score by permutation test
  permutationTest(x, y, numPermu = 1000, maxDelay = 3,scale = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

