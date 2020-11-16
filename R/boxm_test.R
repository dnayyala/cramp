boxm.test <- function(x, y, type = "chi.squared"){
    #' Box-M test statistic
    #'
    #' This function computes the Box M test statistic for testing equality of covariance matrices. The p-value is computed using either a chi-squared or F approximation.
    #' @import mvtnorm stats
    NULL
    #' @export

    #' @param x  data matrix of first group with rows representing the samples and columns representing variables.
    #' @param y  data matrix of second group with same number of columns as \code{x}
    #' @param type Distribution to use for calculating the p-value. Possible values are "chi.squared"(Default) and "F"

    #' @return A list with two values
    #' \describe{
        #' \item{test.statistic}{The test statistic value}
        #' \item{p.value}{The p-value computed using the distribution specified in \code{type}}
    #' }

    #' @examples
    #' library(mvtnorm)
    #' x = rmvnorm(n = 100, mean = numeric(10), sigma = diag(10))
    #' y = rmvnorm(n = 100, mean = numeric(10), sigma = diag(10))
    #' boxm.test(x, y)
    #' @references T. W. Anderson. An introduction to multivariate statistical analysis. Wiley Series in Probability and Statistics, 3rd edition, 2003.

    if (ncol(x) != ncol(y)) stop("Dimensions do not match")
    if (ncol(x) >= nrow(x) || ncol(y) >= nrow(y)) stop("This is not a high dimensional test")
    n <- nrow(x)
    m <- nrow(y)
    p <- ncol(x)

    ## Compute the covariance matrices of the two groups and the pooled covariance matrix
    s1 <- cov(x)
    s2 <- cov(y)
    s.pooled <- ( (n - 1)*s1 + (m - 1)*s2 )/(n + m - 2)

    ## Log-M
    log.M <- ( (n - 1)*log(det(s1)) + (m - 1)*log(det(s2)) - (n + m - 2)*log(det(s.pooled)) )/2

    ## Different test statistics and p-values for the chi-squared and F tests
    if (type == "chi.squared"){
        c1 = (1/(n - 1) + 1/(m - 1) - 1/(n + m - 2))*(2*p^2 + 3*p - 1)/(6*(p + 1))
        test.statistic <- -2*(1 - c1)*log.M
        p.value <- pchisq(q = test.statistic, df = p*(p + 1)/2, lower.tail = FALSE)
    } else if (type == "F"){
        c1 = (1/(n - 1) + 1/(m - 1) - 1/(n + m - 2))*(2*p^2 + 3*p - 1)/(6*(p + 1))
        c2 <- (p - 1)*(p + 2)/6*( 1/(n - 1)^2 + 1/(m - 1)^2 - 1/(n + m - 2)^2   )
        a1 <- p*(p + 1)/2
        a2 <- (a1 + 2)/abs(c2 - c1^2)
        b1 <- (1 - c1 - a1/a2)/a1;
        b2 <- (1 - c1 - 2/a2)/a2;
        if (c2 > c1^2){
            test.statistic <- -2*b1*log.M
        } else {
            test.statistic <- -(a2*b2*log.M)/(a1*(1 + 2*b2*log.M))
        }
        p.value <- pf(q = test.statistic, df1 = a1, df2 = a2, lower.tail = FALSE)
    }

    return(list(test.statistic = test.statistic, p.value = p.value))
}
