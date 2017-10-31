\name{rankNormalization}
\alias{rankNormalization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Rank Normalization of a Vector
}
\description{
  This function performs rank normalization of a given vector.
}
\usage{
  rankNormalization(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Input variable vector
  }
}
\value{
  A normalized vector
}
\examples{

  # Generate one random variable
  x <- rnorm(100, 1, 2)

  # Rank normalize data
  rankNormalization(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

