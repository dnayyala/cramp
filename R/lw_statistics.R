lw.statistic <- function(x, test = "identity"){
    #' Ledoit-Wolf test statistic
    #'
    #' This function computes the Ledoit-Wolf test statistic for testing the shape of covariance matrix in the one-sample case
    #' @param x  data matrix with rows representing samples and columns representing variables
    #' @param test The type of test being performed - "identity"(default) or "sphericity". The test statistics are modified versions of the \code{\link{john.nagao.statistic}}. If \code{test == "identity"}, the test statistic is equal to the Nagao test statistic and when \code{test == "sphericity"}, the Ledoit-Wolf test statistic computed is a modified version of the John test statistic.
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

    #' @seealso \code{\link{john.nagao.statistic}}

    #' @references O. Ledoit and M. Wolf. Some hypothesis tests for the covariance matrix when the dimension is large compared to the sample size. The Annals of Statistics, 30(4):1081–1102, 2002.
    ## Test statistics from Ledoit-Wolf
    n <- nrow(x)
    p <- ncol(x)

    S <- var(x)
    if (test == "sphericity"){
        mu <- 1 - (1 - 2/p)/n;
        U <- (sum(S^2))/(sum(diag(S)))^2
        p.value = pchisq(n*p*mu*(p*U - 1)/2, df = p*(p + 1)/2 - 1, lower.tail = FALSE)
        return(list(test.statistic = n*p*mu*(p*U - 1)/2, p.value = p.value))
    } else if (test == "identity"){
        W <- (1/p)*sum( (S - diag(p))^2 ) - (p/n)*(sum(diag(S))/p)^2 + (p/n)
        p.value = pchisq(n*p*W/2, df = p*(p + 1)/2, lower.tail = FALSE)
        return(list(test.statistic = W, p.value = p.value))
    }
}
