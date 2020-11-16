wald.statistic <- function(x, y){
    #' Two-sample Wald test
    #'
    #' Wald-type test statistic for comparing the covariance matrices of two independent groups.
    #' @param x  data matrix of first group with rows representing samples and columns representing variables.
    #' @param y  data matrix of second group.
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
    #' wald.statistic(x, y)

    #' @references J. R. Schott. A test for the equality of covariance matrices when the dimension is large relative to the sample sizes. Computational Statistics and Data Analysis, 51(12):6535–6542, 2007.

    if (ncol(x) != ncol(y)) stop("Dimensions do not match")
    if (ncol(x) >= nrow(x) || ncol(y) >= nrow(y)) stop("This is not a high dimensional test")
    n <- nrow(x)
    m <- nrow(y)
    p <- ncol(x)

    ## Compute the covariance matrices of the two groups and the pooled covariance matrix
    s1 <- cov(x)
    s2 <- cov(y)
    s.pooled <- ( (n - 1)*s1 + (m - 1)*s2 )/(n + m - 2)

    ## Compute the test statistics
    s.pooled.inv <- solve(s.pooled)
    test.statistic <- (n + m - 2)/2*( ((n-1)/(n+m-2) - (n-1)^2/(n+m-2)^2)*sum(diag( (s1%*%s.pooled.inv)%*%(s1%*%s.pooled.inv)  ))
                                    + ((m-1)/(n+m-2) - (m-1)^2/(n+m-2)^2)*sum(diag( (s2%*%s.pooled.inv)%*%(s2%*%s.pooled.inv)  ))
                                    - 2*(n-1)*(m-1)/(n+m-2)^2* sum(diag( (s1%*%s.pooled.inv)%*%(s2%*%s.pooled.inv)  ))  )
    p.value <- pchisq(q = test.statistic, df = p*(p + 1)/2, lower.tail = FALSE)
    return( list( test.statistic = test.statistic, p.value = p.value) )
}
