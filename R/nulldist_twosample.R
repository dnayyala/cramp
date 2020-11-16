nulldist.twosample <- function(n1, n2 = NULL, p, m = 1,
                                nnull = 1e2, nrep = 1e2, mean1 = NULL, mean2 = NULL, parallel = TRUE, ncores = NULL){
    #' Two sample null distribution generator
    #'
    #' This function generates random samples and evaluates the test statistics for the two-sample covariance matrix test and reports a bootstrap random sample of average p-values. The generated values can be used to calculate the critical value for rejecting/accepting the null hypothesis
    #' @import foreach parallel doParallel mvtnorm
    NULL
    #' @export

    #' @param n1 sample size of the data from the first group
    #' @param n2 sample size of the data from the second group. If \code{n2 = NULL} (Default value), then n2 is set equal to n1
    #' @param p original dimension of your data
    #' @param m projected dimension of your data (Default value is 1)
    #' @param nnull number of repetitions of averages of pvalues
    #' @param nrep number of repetitions for each average of pvalues
    #' @param mean1 (Optional) Mean vector of the first group for the null distribution
    #' @param mean2 (Optional) Mean vector of the second group for the null distribution
    #' @param parallel Whether to parallelize or not (default is \code{TRUE})
    #' @param ncores Number of cores to use (If NULL and parallel is TRUE, use 75% of all cores). If \code{parallel = FALSE}, this argument will be ignored.

    #' @return A list with three elements corresponding to different test statistics
    #' \describe{
    #'      \item{Box-M}{Box's M test statistic}
    #'      \item{Wald}{Wald test statistic}
    #'}

    #' @examples
    #' \dontrun{
    #'  nulldist.twosample(n = 20, p = 100, m = 5)
    #' }

    ## If only one sample size is provided, use that for the second group
    if (is.null(n2)){
        n2 <- n1
    }
    ## If mean vectors  NULL, set them to the zero vector
    if (is.null(mean1)){
        mean1 <- numeric(p)
    }
    if (is.null(mean2)){
        mean2 <- numeric(p)
    }

    if (parallel){
        ## Detemine number of cores if NULL - Use 75% cores if not specified
        if (is.null(ncores)){
            ncores <- ceiling(0.75*detectCores())
        }

        ## Create the cluster
        parclust <- makeCluster(ncores)
        registerDoParallel(parclust)
        # clusterEvalQ(parclust,
        #     sapply(c("ortho_randproj.R", "boxm_test.R", "wald_statistic.R"), source))
        clusterEvalQ(parclust, library(cramp))

        ## Run the loops in parallel and record the average p-pvalues
        pvalue.null <- foreach(i = 1:nnull, .packages = c("mvtnorm", "ismev", "cramp"), .combine = "rbind") %dopar% {
            # X <- rmvnorm(n = n1, mean = numeric(p), sigma = diag(p))
            # Y <- rmvnorm(n = n2, mean = numeric(p), sigma = diag(p))
            X <- matrix(rnorm(n1*p), nrow = n1, ncol = p)
            Y <- matrix(rnorm(n1*p), nrow = n1, ncol = p)
            boxm <- wald <- numeric(nrep)

            for (j in 1:nrep){
                Rstar <- ortho.randproj(nrow = m, ncol = p)
                X.proj <- t(Rstar%*%t(X))
                Y.proj <- t(Rstar%*%t(Y))
                boxm[j] <- boxm.test(x = X.proj, y = Y.proj)$p.value
                wald[j] <- wald.statistic(x = X.proj, y = Y.proj)$p.value
            }
            c(mean(boxm), mean(wald))
        }
        stopCluster(parclust)

    } else {
        sapply(c("ortho_randproj.R", "boxm_test.R", "wald_statistic.R"), source)
        pvalue.null <- matrix(0, nrow = nnull, ncol = 2)
        for (i in 1:nnull){
            # X <- rmvnorm(n = n1, mean = numeric(p), sigma = diag(p))
            # Y <- rmvnorm(n = n2, mean = numeric(p), sigma = diag(p))
            X <- matrix(rnorm(n1*p), nrow = n1, ncol = p)
            Y <- matrix(rnorm(n1*p), nrow = n1, ncol = p)
            boxm <- wald <- numeric(nrep)

            ## Loop for the multiple random projections
            for (j in 1:nrep){
                Rstar <- ortho.randproj(nrow = m, ncol = p)
                X.proj <- t(Rstar%*%t(X))
                Y.proj <- t(Rstar%*%t(Y))
                boxm[j] <- boxm.test(x = X.proj, y = Y.proj)$p.value
                wald[j] <- wald.statistic(x = X.proj, y = Y.proj)$p.value
            }
            pvalue.null[i,] <- c(mean(boxm), mean(wald))
        }
    }

    colnames(pvalue.null) <- c("Box-M", "Wald")
    return(pvalue.null)

}
