schott.test <- function(x, y){
    #' Schott's test
    #'
    #' This function computes the Schott's two sample test statistic for testing equality of covariance matrices as described in Schott (need citation)
    #' @param x  data matrix with rows representing samples and columns representing variables
    #' @param y  data matrix of second group with same number of columns as \code{x}
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
    #' x = rmvnorm(n = 100, mean = numeric(10), sigma = diag(10))
    #' y = rmvnorm(n = 100, mean = numeric(10), sigma = diag(10))
    #' schott.test(x, y)

    #' @references J. R. Schott. A test for the equality of covariance matrices when the dimension is large relative to the sample sizes. Computational Statistics and Data Analysis, 51(12):6535–6542, 2007.
    
    if (ncol(x) != ncol(y)) stop("Dimensions do not match")

    n <- nrow(x)
    m <- nrow(y)
    p <- ncol(x)

    ## Compute the covariance matrices of the two groups and the pooled covariance matrix
    s1 <- cov(x)
    s2 <- cov(y)
    s.pooled <- ( (n - 1)*s1 + (m - 1)*s2 )/(n + m - 2)

    ## Compute the test statistic
    T <- (sum((s1 - s2)^2)
                      - ( (n - 1)*(n - 3)*sum(s1^2) + (n - 1)^2*sum(diag(s1))^2 )/((n-1)*(n-2)*(n+1))
                      - ( (m - 1)*(m - 3)*sum(s2^2) + (m - 1)^2*sum(diag(s2))^2 )/((m-1)*(m-2)*(m+1)) )

    theta.hat <- (4*( (n + m - 2)/((n - 1)*(m - 1)) )^2*(n + m - 2)^2/( (n+m)*(n+m-3) )
                  *(sum(s.pooled^2) - sum(diag(s.pooled))^2/(n + m - 2)) )

    test.statistic <- T/sqrt(theta.hat)
    p.value <- pnorm(q = test.statistic, lower.tail = FALSE)

    return(list(test.statistic = test.statistic, p.value = p.value))
}
