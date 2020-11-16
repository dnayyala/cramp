lrt.statistic <- function(x, test = "identity"){
    #' One-sample likelihood ratio tests
    #'
    #' The traditional likelihood ratio test statistic for the one-sample hypothesis of covariance matrix.
    #' @param x  data matrix with rows representing samples and columns representing variables
    #' @param test The type of test being performed - "identity"(default) or "sphericity".
    #' @import mvtnorm stats
    NULL
    #' @export

    #' @return A list with two values
    #' \describe{
        #' \item{test.statistic}{The test statistic value}
        #' \item{p.value}{The p-value}
    #' }

    #' @examples
    #' library(mvtnorm)
    #' # Test for identity
    #' x = rmvnorm(n = 20, mean = numeric(100), sigma = diag(100))
    #' lrt.statistic(x, test = "identity")

    #' # Test for sphericity
    #' y = rmvnorm(n = 10, mean = runif(100, 1, 3), sigma = diag(rgamma(100, 10, 3)))
    #' lrt.statistic(y, test = "sphericity")
    #' @references T. W. Anderson. An introduction to multivariate statistical analysis. Wiley Series in Probability and Statistics, 3rd edition, 2003.

    n <- nrow(x)
    p <- ncol(x)

    lambda.vec <- eigen(var(x)*(n - 1)/n)$values
    nu = n - 1
    if (test == "identity"){
        test.statistic <- (1 - (2*p + 1 - 2/(p + 1))/(6*nu - 1))*nu*( -sum(log(lambda.vec)) + sum(lambda.vec) - p )
        p.value <- pchisq(test.statistic, df = p*(p + 1)/2, lower.tail = FALSE)
    } else if (test == "sphericity"){
        test.statistic <- -(nu - (2*p^2 + p + 2)/(6*p))*(p*log(p) + sum(log(lambda.vec)) - p*log(sum(lambda.vec)))
        p.value <- pchisq(test.statistic, df = p*(p + 1)/2 - 1, lower.tail = FALSE)
    }
    return(list(test.statistic = test.statistic, p.value = p.value))
}
