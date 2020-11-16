srivastava.statistic <- function(x, y = NULL, test = "identity"){
    #' Srivastava et al. test statistics
    #'
    #' This function computes the test statistics from Srivastava et al. (2014) for tests involving the covariance matrices. If two sample matrices are provided, the function returns the result of the comparison of two covariance matrices. When only one data matrix is provided, tests for \code{sphericity} and \code{identity} can be performed.
    #' @import mvtnorm stats
    NULL
    #' @export

    #' @param x  data matrix of first group with rows representing samples and columns representing variables.
    #' @param y  data matrix of second group with same number of columns as \code{x}. If not provided (default value is NULL), then the one-sample tests are performed as specified.
    #' @param test Type of test to be performed. When \code{y = NULL}, you can specify "identity" (Default) or "sphericity". If \code{y} is provided as a non-null data matrix, any specified value for \code{test} will be over-written to \code{test = "two.sample"}.

    #' @return A list with two values
    #' \describe{
        #' \item{test.statistic}{The test statistic value}
        #' \item{p.value}{The p-value}
    #' }

    #' @examples
    #' library(mvtnorm)
    #' x = rmvnorm(n = 100, mean = numeric(10), sigma = diag(10))
    #' y = rmvnorm(n = 100, mean = numeric(10), sigma = diag(10))
    #' srivastava.statistic(x, y)
    #' srivastava.statistic(x)

    #' @references M. S. Srivastava, H. Yanagihara, and T. Kubokawa. Tests for covariance matrices in high dimension with less sample size. Journal of Multivariate Analysis, 130:289–309, 2014.
    n <- nrow(x)
    p <- ncol(x)
    if (!is.null(y)){
        if (ncol(x) != ncol(y)) stop("Dimensions do not match")
        m <- nrow(y)
        test = "two.sample"
    }

    ## Center the data and transpose
    x.center <- t(x) - colMeans(x)
    ## Term a1
    a1 = sum( x.center^2  )/((n - 1)*p)
    ### Matrix forms for Term a2
    a2.term1 <- sum( (t(x.center)%*%x.center)^2 )
    a2.term2 <- sum( (rowSums(x.center^2))^2 )
    a2.term3 <- sum( rowSums(x.center^2) )^2
    a2 <- ((n-2)*(n-1)*a2.term1 - n*(n-1)*a2.term2 + a2.term3)/(p*n*(n-1)*(n-2)*(n-3))


    if (test == "sphericity"){
        test.statistic <- ((n - 1)/2)*(a2/a1^2 - 1)
        p.value <- pnorm(test.statistic, lower.tail = FALSE)
    } else if (test == "identity"){
        test.statistic <- ((n - 1)/2)*(a2 - 2*a1 + 1)
        p.value <- pnorm(test.statistic, lower.tail = FALSE)
    } else if (test == "two.sample"){
        ## Center the data and transpose
        y.center <- t(y) - colMeans(y)
        ## Term b1 -- Equal to term a1 for the second population
        b1 = sum( y.center^2  )/((m - 1)*p)
        # ## Term b2 -- Equal to term a2 for the second population
        b2.term1 <- sum( (t(y.center)%*%y.center)^2 )
        b2.term2 <- sum( (rowSums(y.center^2))^2 )
        b2.term3 <- sum( rowSums(y.center^2) )^2
        b2 <- ((m-2)*(m-1)*b2.term1 - m*(m-1)*b2.term2 + b2.term3)/(p*m*(m-1)*(m-2)*(m-3))
        ## Third term of the numerator
        c = sum(diag(  (x.center%*%t(x.center))%*%(y.center%*%t(y.center)) ) )/( (n - 1)*(m - 1)*p  )

        ## Compute the test statistic and the p-value
        test.statistic <- ( a2 + b2 - 2*c )/(2*(1/(n-1) + 1/(m-1))*((n - 1)*a2 + (m - 1)*b2)/(n + m - 2) )
        p.value <- pnorm(test.statistic, lower.tail = FALSE)
    }
    return(list(test.statistic = test.statistic, p.value = p.value))

}
