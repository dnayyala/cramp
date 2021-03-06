% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chenzz_statistic.R
\name{chenzz.statistic}
\alias{chenzz.statistic}
\title{Chen et al. test statistic}
\usage{
chenzz.statistic(x, test = "identity")
}
\arguments{
\item{x}{data matrix with rows representing samples and columns representing variables}

\item{test}{The type of test being performed - "identity"(default) or "sphericity".}
}
\value{
A list with two values
\describe{
\item{test.statistic}{The test statistic value}
\item{p.value}{The p-value}
}
}
\description{
This function computes the test statistic for testing equality of covariance matrices as described in the Chen et al. paper (need citation)
}
\examples{
library(mvtnorm)
# Test for identity
x = rmvnorm(n = 20, mean = numeric(100), sigma = diag(100))
chenzz.statistic(x, test = "identity")
# Test for sphericity
y = rmvnorm(n = 10, mean = runif(100, 1, 3), sigma = diag(rgamma(100, 10, 3)))
chenzz.statistic(y, test = "sphericity")
}
\references{
S. X. Chen, L. X. Zhang, and P. S. Zhong. Tests for high-dimensional covariance matrices. Journal of the American Statistical Association, 105(490):810–819, 2010.
}
