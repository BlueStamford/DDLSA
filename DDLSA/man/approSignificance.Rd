\name{approSignificance}
\alias{approSignificance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Approximate statistical significance of local similarity score of two time series
}
\description{
  This function computes the approximated statistical significance of local similarity for two sequences.
}
\usage{
  approSignificance(x, y, maxDelay)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{
    two random variables with equal length.
  }
  \item{maxDelay}{
    maximum time shift allowed in computing LS score.
  }
}
\value{
  This function returns the approximated statistical significance(p-value) of local similarity of two time series x, y with maximal delay maxDelay.
}
\examples{

  # Generate two random variables
  library(MASS)
  data <- mvrnorm(100, c(0, 0), Sigma = matrix(c(1, 0.6, 0.6, 1), 2, 2))

  # Rank normalize data
  x <- rankNormalization(data[,1])
  y <- rankNormalization(data[,2])

  # compute the approximated significance of local similarity score
  approSignificance(x, y, maxDelay=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

