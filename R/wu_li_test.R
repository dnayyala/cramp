wu.li.statistic <- function(x, y = NULL, nrand = 1e3, test = "identity"){
    #' Wu-Li test statistic
    #'
    #' This function computes the test statistics for one and two sample hypothesis of covariance matrices as described in Wu and Li (need citation)
    #' @import mvtnorm stats
    NULL
    #' @export

    #' @param x  data matrix of first group with rows representing samples and columns representing variables.
    #' @param y  data matrix of second group with same number of columns as \code{x}. If not provided (default value is NULL), then the one-sample tests are performed as specified.
    #' @param nrand number of random projections to be generated for the one-dimension projections.
    #' @param test Type of test to be performed. When \code{y = NULL}, you can specify "identity" (Default) or "sphericity". If \code{y} is provided as a non-null data matrix, any specified value for \code{test} will be over-written to \code{test = "two.sample"}.

    #' @return A list with two values
    #' \describe{
        #' \item{test.statistic}{The test statistic value}
        #' \item{crit.value}{The critical value to accept or reject the null hypothesis based on the \code{test.statistic} returned}
    #' }
    #' @examples
    #' library(mvtnorm)
    #' x = rmvnorm(n = 100, mean = numeric(10), sigma = diag(10))
    #' y = rmvnorm(n = 100, mean = numeric(10), sigma = diag(10))
    #' wu.li.statistic(x, y, nrand = 1e2)

    #' @seealso \code{gev.critvalue}
    #' @references T.-L. Wu and P. Li. Projected tests for high-dimensional covariance matrices. Journal of Statistical Planning and Inference, 207:73 – 85, 2020.
    n <- nrow(x)
    p <- ncol(x)
    if (!is.null(y)){
        if (ncol(x) != ncol(y)) stop("Dimensions do not match")
        m <- nrow(y)
        test = "two.sample"
    }

    if (test == "identity"){
        ## Generate and normalize the random vectors
        T.n <- numeric(nrand)
        for (rep in 1:nrand){
            R <- rnorm(n = p, mean = 0, sd = 1); R <- R/sqrt(sum(R^2))
            X <- x%*%R
            T.n[rep] <- sum( (X - mean(X))^2 )
        }
        test.statistic <- max(T.n)

        ## Computing the critical value
        if (!exists("sig.level")){
            sig.level <- 0.05
        }
        crit.value <- gev.critvalue(T.n, sig.level)
        return(list(test.statistic = test.statistic, crit.value = mean(crit.value)))
    } else if (test == "two.sample"){
        ## Generate and normalize the random vectors
        T.nm <- numeric(nrand)
        for (rep in 1:nrand){
            R <- rnorm(n = p, mean = 0, sd = 1); R <- R/sqrt(sum(R^2))
            T.nm[rep] <- (var(x%*%R)/var(y%*%R))*((n - 1)/n)*(m/(m-1))
        }
        test.statistic <- max(T.nm)

        ## Computing the critical value
        if (!exists("sig.level")){
            sig.level <- 0.05
        }
        crit.value <- gev.critvalue(T.nm, nrand, sig.level)
        return(list(test.statistic = test.statistic, crit.value = mean(crit.value)))
    }
}
