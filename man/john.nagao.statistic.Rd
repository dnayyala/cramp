% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/john_nagao_tests.R
\name{john.nagao.statistic}
\alias{john.nagao.statistic}
\title{John and Nagao test statistics}
\usage{
john.nagao.statistic(x, test = "identity")
}
\arguments{
\item{x}{data matrix with rows representing samples and columns representing variables}

\item{test}{The type of test being performed - "identity"(default) or "sphericity". If \code{test == "identity"}, the Nagao test statistic is computed and when \code{test == "sphericity"}, the John test statistic is computed.}
}
\value{
A list with two values
\describe{
\item{test.statistic}{The test statistic value}
\item{p.value}{The p-value}
}
}
\description{
This function computes the test statistic for testing equality of covariance matrices as described by John (1972) and Nagao (1973)
}
\examples{
library(mvtnorm)
# Nagao's test for identity
x = rmvnorm(n = 20, mean = numeric(100), sigma = diag(100))
john.nagao.statistic(x, test = "identity")
# John's test for sphericity
y = rmvnorm(n = 10, mean = runif(100, 1, 3), sigma = diag(rgamma(100, 10, 3)))
john.nagao.statistic(y, test = "sphericity")
}
\references{
S. John. The Distribution of a Statistic Used for Testing Sphericity of Normal Distributions. Biometrika, 59(1):169–173, 1972.
H. Nagao. On some test criteria for covariance matrix. The Annals of Statistics, 1(4):700–709, 1973.
}
