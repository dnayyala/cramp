li.chen.statistic <- function(x, y){
    #'Li-Chen test statistic
    #'
    #' This function computes the test statistic for testing equality of covariance matrices as described by Li and Chen (need citation)
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
    #' x = rmvnorm(n = 20, mean = numeric(100), sigma = diag(100))
    #' y = rmvnorm(n = 20, mean = numeric(100), sigma = diag(100))
    #' li.chen.statistic(x, y)

    #' @references J. Li and S. X. Chen. Two sample tests for high-dimensional covariance matrices. Annals of Statistics, 40(2):908–940, 2012.
    ## Error message if the dimensions do not match
    if (ncol(x) != ncol(y)) stop("Dimensions do not match")

    n <- nrow(x)
    m <- nrow(y)
    p <- ncol(x)

    ## Numerator terms
    A.x <- A.y <- numeric(3)
    C <- numeric(4)

    ## Terms from the first group.
    ## Also include the quantities for the third term
    for (i in 1:n){
        for (j in 1:m){
            C[1] <- C[1] + sum(x[i,]*y[j,])^2
        }
        for (j in setdiff(1:n, i)){
            A.x[1] <- A.x[1] + sum(x[i,]*x[j,])^2
            for (k in 1:m){
                C[2] <- C[2] + sum(x[i,]*y[k,])*sum(x[j,]*y[k,])
                for (l in setdiff(1:m, k)){
                    C[4] <- C[4] + sum(x[i,]*y[k,])*sum(x[j,]*y[l,])
                }
            }
            for (k in setdiff(1:n, c(i,j))){
                A.x[2] <- A.x[2] + sum(x[i,]*x[j,])*sum(x[i,]*x[k,])
                for (l in setdiff(1:n, c(i,j,k))){
                    A.x[3] <- A.x[3] + sum(x[i,]*x[j,])*sum(x[k,]*x[l,])
                }
            }
        }
    }

    ## Terms from the second group
    for (i in 1:m){
        for (j in setdiff(1:m, i)){
            A.y[1] <- A.y[1] + sum(y[i,]*y[j,])^2
            for (k in 1:n){
                C[3] <- C[3] + sum(y[i,]*x[k,])*sum(y[j,]*x[k,])
            }
            for (k in setdiff(1:m, c(i,j))){
                A.y[2] <- A.y[2] + sum(y[i,]*y[j,])*sum(y[i,]*y[k,])
                for (l in setdiff(1:m, c(i,j,k))){
                    A.y[3] <- A.y[3] + sum(y[i,]*y[j,])*sum(y[k,]*y[l,])
                }
            }
        }
    }

    T.term1 <- A.x[1]/(n*(n - 1)) - 2*A.x[2]/(n*(n - 1)*(n - 2)) + A.x[3]/(n*(n - 1)*(n - 2)*(n - 3))
    T.term2 <- A.y[1]/(m*(m - 1)) - 2*A.y[2]/(m*(m - 1)*(m - 2)) + A.y[3]/(m*(m - 1)*(m - 2)*(m - 3))
    T.term3 <- C[1]/(n*m) - C[2]/(n*(n-1)*m) - C[3]/(m*(m-1)*n) + C[4]/(n*(n-1)*m*(m-1))

    test.statistic <- (T.term1 + T.term2 - 2*T.term3)/sqrt(2*T.term1/n + 2*T.term2/m)
    p.value <- pnorm(test.statistic, lower.tail = FALSE)

    return(list(test.statistic = test.statistic, p.value = p.value))

}
