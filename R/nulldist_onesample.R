nulldist.onesample <- function(n, p, m = 1, nnull = 1e2, nrep = 1e2, mean = NULL, parallel = TRUE, ncores = NULL){
    #' One sample null distribution generator
    #'
    #' This function generates random samples and evaluates the test statistics for the identity hypothesis and reports a bootstrap random sample of average p-values. The generated values can be used to calculate the critical value for rejecting/accepting the null hypothesis

    #' @import foreach parallel doParallel mvtnorm
    NULL
    #' @export

    #' @param n sample size of your data
    #' @param p original dimension of your data
    #' @param m projected dimension of your data (Default value is 1)
    #' @param nnull number of repetitions of averages of pvalues
    #' @param nrep number of repetitions for each average of pvalues
    #' @param mean (Optional) Mean vector for the null distribution
    #' @param parallel Whether to parallelize or not (default is \code{TRUE})
    #' @param ncores Number of cores to use (If NULL and parallel is TRUE, use 75% of all cores). If \code{parallel = FALSE}, this argument will be ignored.

    #' @return A list with three elements corresponding to different test statistics
    #' \describe{
    #'      \item{LRT}{likelihood ratio test statistic}
    #'      \item{LW}{Ledoit-Wolf test statistic}
    #'      \item{John}{John-Nagao test statistic}
    #'}

    #' @examples
    #' \dontrun{
    #'   nulldist.onesample(n = 20, p = 100, m = 5)
    #' }

    if (is.null(mean)){
        mean <- numeric(p)
    }

    ## Parallelize if required
    if (parallel){
        ## Detemine number of cores if NULL. Use 75% cores if not specified
        if (is.null(ncores)){
            ncores <- ceiling(0.75*detectCores())
        }

        ## Create the cluster
        parclust <- makeCluster(ncores)
        registerDoParallel(parclust)
        # clusterEvalQ(parclust,
        #     sapply(c("ortho_randproj.R", "lrt_tests.R", "lw_statistics.R", "john_nagao_tests.R"), source))
        clusterEvalQ(parclust, library(cramp))

        ## Run the loops in parallel and record the average p-values
        pvalue.null <- foreach(i=1:nnull, .packages = c("mvtnorm", "cramp"), .combine = "rbind") %dopar% {
            # X <- rmvnorm(n = n, mean = mean, sigma = diag(p));
            X <- matrix(rnorm(n*p), nrow = n, ncol = p)
            lrt <- lw <- john <- numeric(nrep)

            ## Loop for the multiple random projections
            for (j in 1:nrep){
                Rstar <- ortho.randproj(nrow = m, ncol = p)
                X.proj <- t(Rstar%*%t(X))
                lrt[j] <- lrt.statistic(x = X.proj, test = "identity")$p.value
                lw[j] <- lw.statistic(x = X.proj, test = "identity")$p.value
                john[j] <- john.nagao.statistic(x = X.proj, test = "identity")$p.value
            }
            c(mean(lrt), mean(lw), mean(john))
        }
        stopCluster(parclust)
    } else {
        sapply(c("ortho_randproj.R", "lrt_tests.R", "lw_statistics.R", "john_nagao_tests.R"), source)
        pvalue.null <- matrix(0, nrow = nnull, ncol = 3)
        for (i in 1:nnull){
            # X <- rmvnorm(n = n, mean = mean, sigma = diag(p));
            X <- matrix(rnorm(n*p), nrow = n, ncol = p)
            lrt <- lw <- john <- numeric(nrep)

            ## Loop for the multiple random projections
            for (j in 1:nrep){
                Rstar <- ortho.randproj(nrow = m, ncol = p)
                X.proj <- t(Rstar%*%t(X))
                lrt[j] <- lrt.statistic(x = X.proj, test = "identity")$p.value
                lw[j] <- lw.statistic(x = X.proj, test = "identity")$p.value
                john[j] <- john.nagao.statistic(x = X.proj, test = "identity")$p.value
            }
            pvalue.null[i,] <- c(mean(lrt), mean(lw), mean(john))
        }
    }

    colnames(pvalue.null) <- c("LRT", "LW", "John")
    return(pvalue.null)

}
