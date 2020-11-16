chenzz.statistic <- function(x, test = "identity"){
    #' Chen et al. test statistic
    #'
    #' This function computes the test statistic for testing equality of covariance matrices as described in the Chen et al. paper (need citation)
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
    #' chenzz.statistic(x, test = "identity")

    #' # Test for sphericity
    #' y = rmvnorm(n = 10, mean = runif(100, 1, 3), sigma = diag(rgamma(100, 10, 3)))
    #' chenzz.statistic(y, test = "sphericity")
    #' @references S. X. Chen, L. X. Zhang, and P. S. Zhong. Tests for high-dimensional covariance matrices. Journal of the American Statistical Association, 105(490):810–819, 2010.

    n <- nrow(x)
    p <- ncol(x)

    ## Five individual terms

    y1n <- y2n <- y3n <- y4n <- y5n <- 0;
    for (i in 1:n){
        y1n <- y1n + sum(x[i,]^2)
        for (j in setdiff(1:n, i)){
            y2n <- y2n + sum(x[i,]*x[j,])^2
            y3n <- y3n + sum(x[i,]*x[j,])
            for (k in setdiff(1:n, c(i,j))){
                y4n <- y4n + sum(x[i,]*x[j,])*sum(x[j,]*x[k,])
                for (l in setdiff(1:n, c(i,j,k))){
                    y5n <- y5n + sum(x[i,]*x[j,])*sum(x[k,]*x[l,])
                }
            }
        }
    }
    y1n <- y1n/n
    y2n <- y2n/(choose(n,2)*factorial(2))
    y3n <- y3n/(choose(n,2)*factorial(2))
    y4n <- y4n/(choose(n,3)*factorial(3))
    y5n <- y5n/(choose(n,4)*factorial(4))

    ## Combine terms to calculate T1 and T2
    T1n <- y1n - y3n
    T2n <- y2n - 2*y4n + y5n
    if (test == "identity"){
        test.statistic <- T2n/p - 2*T1n/p + 1
    } else if (test == "sphericity"){
        test.statistic <- p*(T2n/T1n^2) - 1
    }

    return(list(test.statistic = n*test.statistic/2, p.value = pnorm(n*test.statistic/2, lower.tail = FALSE )))
}
