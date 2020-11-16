gev.critvalue <- function(T.n, nrand, sig.level){
    #' Wu-Li critical value
    #'
    #' Function to compute critical values for the \code{\link{wu.li.statistic}}
    #' @import ismev stats
    NULL

    #' @param T.n,nrand,sig.level Intermediary quantities computed in the \code{\link{wu.li.statistic}}
    #' @return The critical value for rejecting/accepting the null hypothesis as defined in Wu and Li (2020)
    #' @seealso \code{\link{wu.li.statistic}}
    #' @references T.-L. Wu and P. Li. Projected tests for high-dimensional covariance matrices. Journal of Statistical Planning and Inference, 207:73 – 85, 2020. ISSN 0378-3758.
    ## Computing the critical value
    cutoff <- quantile(T.n, 0.95)
    T.n <- sort(T.n)
    crit.value <- numeric(10)
    for (rep in 1:10){
        index.set <- sort(sample(nrand, size = nrand, replace = TRUE))
        T.n.tilde <-T.n[index.set]
        T.n.tilde.sub <- T.n.tilde[T.n.tilde >= cutoff]
        index.set <- index.set[T.n.tilde >= cutoff]
        m.u <- length(index.set)
        Wk.vec <- index.set[-1] - index.set[-m.u]
        m.c <- sum(Wk.vec - 1 != 0)

        ## Estimate of theta
        theta.hat <- 1 - (m.u - m.c - 1)/(2*m.c - (1 - m.u/nrand)*sum(Wk.vec - 1))
        ## Parameter estimates
        param.estims <- gev.fit(T.n.tilde.sub, show = FALSE);
        mu.hat <- param.estims$mle[1]; sigma.hat <- param.estims$mle[2]; eta.hat <- param.estims$mle[3];

        ## Compute the critical value for the current iteration
        if (eta.hat != 0){
            crit.value[rep] <- mu.hat + sigma.hat/eta.hat*( (-log(1 - sig.level)/theta.hat)^(-eta.hat) - 1 )
        } else {
            crit.value[rep] <- mu.hat - sigma.hat*log( -log(1 - sig.level)/theta.hat )
        }
    }
    return(crit.value)
}
