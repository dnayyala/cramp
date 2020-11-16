cai.statistic <- function(x, y, alpha = 0.05){
    #' Cai test statistic
    #'
    #' This function computes the test statistic for testing equality of covariance matrices as described in the Cai et al. paper (need citation)
    #' @param x  data matrix of first group with rows representing the samples and columns representing variables.
    #' @param y  data matrix of second group with same number of columns as \code{x}
    #' @param alpha Significance level of the test (default value is 0.05)

    #' @return A list with two values
    #' \describe{
        #' \item{test.statistic}{The test statistic value}
        #' \item{decision}{A binary response of whether the null hypothesis is rejected or not rejected}
    #' }
    #' @export

    #' @examples
    #' library(mvtnorm)
    #' x = rmvnorm(n = 20, mean = numeric(100), sigma = diag(100))
    #' y = rmvnorm(n = 10, mean = numeric(100), sigma = diag(100))
    #' cai.statistic(x, y)
    #' @references T. Cai, W. Liu, and Y. Xia. Two-sample covariance matrix testing and support recovery in high-dimensional and sparse settings. Journal of the American Statistical Association, 108(501):265–277, 2013.

    if (ncol(x) != ncol(y)) stop("Dimensions do not match")
    n <- nrow(x)
    m <- nrow(y)
    p <- ncol(x)


    ## Compute the test statistic
    x.center <- t(x) - colMeans(x);
    y.center <- t(y) - colMeans(y);
    M <- matrix(0, nrow = p, ncol = p);
    for (i in 1:(p - 1)){
        for (j in (i + 1):p){
            theta1.ij <- mean( (x.center[i,]*x.center[j,] - mean(x.center[i,]*x.center[j,]))^2 )
            theta2.ij <- mean( (y.center[i,]*y.center[j,] - mean(y.center[i,]*y.center[j,]))^2 )
            M[i,j] <- ( mean(x.center[i,]*x.center[j,]) - mean(y.center[i,]*y.center[j,]) )^2/(theta1.ij/n + theta2.ij/m)
        }
    }
    test.statistic <- max(M)
    reject.crit <- -log(8*pi) - 2*log(log(1/(1 - alpha)))  +4*log(p) - log(log(p))

    return(list(test.statistic = test.statistic, decision = c("Do not reject H0", "Reject H0")[1 + (test.statistic > reject.crit)]) )
}
