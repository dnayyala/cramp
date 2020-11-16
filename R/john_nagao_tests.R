john.nagao.statistic <- function(x, test = "identity"){
    #'John and Nagao test statistics
    #'
    #' This function computes the test statistic for testing equality of covariance matrices as described by John (1972) and Nagao (1973)
    #' @param x  data matrix with rows representing samples and columns representing variables
    #' @param test The type of test being performed - "identity"(default) or "sphericity". If \code{test == "identity"}, the Nagao test statistic is computed and when \code{test == "sphericity"}, the John test statistic is computed.
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
    #' # Nagao's test for identity
    #' x = rmvnorm(n = 20, mean = numeric(100), sigma = diag(100))
    #' john.nagao.statistic(x, test = "identity")

    #' # John's test for sphericity
    #' y = rmvnorm(n = 10, mean = runif(100, 1, 3), sigma = diag(rgamma(100, 10, 3)))
    #' john.nagao.statistic(y, test = "sphericity")

    #' @references S. John. The Distribution of a Statistic Used for Testing Sphericity of Normal Distributions. Biometrika, 59(1):169–173, 1972.
    #' H. Nagao. On some test criteria for covariance matrix. The Annals of Statistics, 1(4):700–709, 1973.

    n <- nrow(x)
    p <- ncol(x)

    S <- var(x)
    if (test == "sphericity"){
        ## Calculate the John test statistic
        test.statistic <- (n*p*(p - 1)/2)*( p*sum(S^2)/(sum(diag(S))^2)  - 1)/(p - 1)
        p.value <- pchisq(test.statistic, df = p*(p + 1)/2 - 1, lower.tail = FALSE)
    } else if (test == "identity"){
        ## Calculate the Nagao test statistic
        # test.statistic <- (1/p)*sum( (S - diag(p))^2 )
        test.statistic <- (n/2)*sum( (S - diag(p))^2 )
        p.value <- pchisq(test.statistic, df = p*(p + 1)/2 - 1, lower.tail = FALSE)
    }
    return(list(test.statistic = test.statistic, p.value = p.value))
}
