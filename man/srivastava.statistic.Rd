% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/srivastava_tests.R
\name{srivastava.statistic}
\alias{srivastava.statistic}
\title{Srivastava et al. test statistics}
\usage{
srivastava.statistic(x, y = NULL, test = "identity")
}
\arguments{
\item{x}{data matrix of first group with rows representing samples and columns representing variables.}

\item{y}{data matrix of second group with same number of columns as \code{x}. If not provided (default value is NULL), then the one-sample tests are performed as specified.}

\item{test}{Type of test to be performed. When \code{y = NULL}, you can specify "identity" (Default) or "sphericity". If \code{y} is provided as a non-null data matrix, any specified value for \code{test} will be over-written to \code{test = "two.sample"}.}
}
\value{
A list with two values
\describe{
\item{test.statistic}{The test statistic value}
\item{p.value}{The p-value}
}
}
\description{
This function computes the test statistics from Srivastava et al. (2014) for tests involving the covariance matrices. If two sample matrices are provided, the function returns the result of the comparison of two covariance matrices. When only one data matrix is provided, tests for \code{sphericity} and \code{identity} can be performed.
}
\examples{
library(mvtnorm)
x = rmvnorm(n = 100, mean = numeric(10), sigma = diag(10))
y = rmvnorm(n = 100, mean = numeric(10), sigma = diag(10))
srivastava.statistic(x, y)
srivastava.statistic(x)
}
\references{
M. S. Srivastava, H. Yanagihara, and T. Kubokawa. Tests for covariance matrices in high dimension with less sample size. Journal of Multivariate Analysis, 130:289–309, 2014.
}
