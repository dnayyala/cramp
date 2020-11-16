cramp.statistic <- function(x, y = NULL, nrep = 1e3, m = 1, parallel = TRUE, ncores = NULL){
    #' CRAMP test statistic
    #'
    #' This function computes the test statistic based on random projections and the Hotelling's \eqn{T^2} test statistic.
    #' @import foreach parallel doParallel mvtnorm
    NULL
    #' @export

    #' @param x  data matrix of first group with rows representing samples and columns representing variables.
    #' @param y  data matrix of second group with same number of columns as \code{x}. If not provided (default value is NULL), then the one-sample tests are performed as specified.
    #' @param nrep number of repetitions for each average of pvalues
    #' @param m Dimension to randomly project the data (Default value is 1)
    #' @param parallel Logical variable - whether to parallelize or not (default is \code{TRUE})
    #' @param ncores Number of cores to use (If NULL and parallel is TRUE, use 75% of all cores). If \code{parallel = FALSE}, this argument will be ignored.

    #' @return The average p-value based on \code{nrep} bootstrap random projections.

    #' @examples
    #' \dontrun{
    #'      x = rmvnorm(n = 20, mean = numeric(1000), sigma = diag(runif(1000)))
    #'      y = rmvnorm(n = 20, mean = numeric(1000), sigma = diag(runif(1000)))
    #'      cramp.statistic(x, y, m = 5)
    #' }


    p <- ncol(x)
    if (!is.null(y)){
        if (ncol(x) != ncol(y)) stop("Dimensions do not match")
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

        if (is.null(y)){
            pvalue.null <- foreach(i=1:nrep, .packages = c("mvtnorm", "cramp"), .combine = "rbind") %dopar% {
                Rstar <- ortho.randproj(nrow = m, ncol = p)
                X.proj <- t(Rstar%*%t(x))
                lrt <- lrt.statistic(x = X.proj, test = "identity")$p.value
                lw <- lw.statistic(x = X.proj, test = "identity")$p.value
                john <- john.nagao.statistic(x = X.proj, test = "identity")$p.value
                c(lrt, lw , john)
            }
            colnames(pvalue.null) <- c("LRT", "LW", "John-Nagao")
            stopCluster(parclust)
        } else {
            pvalue.null <- foreach(i = 1:nrep, .packages = c("mvtnorm", "ismev", "cramp"), .combine = "rbind") %dopar% {
                Rstar <- ortho.randproj(nrow = m, ncol = p)
                X.proj <- t(Rstar%*%t(x))
                Y.proj <- t(Rstar%*%t(y))
                boxm<- boxm.test(x = X.proj, y = Y.proj)$p.value
                wald <- wald.statistic(x = X.proj, y = Y.proj)$p.value
                c(mean(boxm), mean(wald))
            }
            colnames(pvalue.null) <- c("Box-M", "Wald")
            stopCluster(parclust)
        }
    } else {
        ## Run it in serial
        if (is.null(y)){
            pvalue.null <- matrix(0, nrow = nrep, ncol = 3)
            colnames(pvalue.null) <- c("LRT", "LW", "John-Nagao")
            for (i in 1:nrep){
                Rstar <- ortho.randproj(nrow = m, ncol = p)
                X.proj <- t(Rstar%*%t(x))
                pvalue.null[i,1] <- lrt.statistic(x = X.proj, test = "identity")$p.value
                pvalue.null[i,2] <- lw.statistic(x = X.proj, test = "identity")$p.value
                pvalue.null[i,3] <- john.nagao.statistic(x = X.proj, test = "identity")$p.value
            }
        } else {
            pvalue.null <-  matrix(0, nrow = nrep, ncol = 3)
            colnames(pvalue.null) <- c("Box-M", "Wald")
            for (i in 1:nrep){
                Rstar <- ortho.randproj(nrow = m, ncol = p)
                X.proj <- t(Rstar%*%t(x))
                Y.proj <- t(Rstar%*%t(y))
                pvalue.null[i,1] <- boxm.test(x = X.proj, y = Y.proj)$p.value
                pvalue.null[i,2] <- wald.statistic(x = X.proj, y = Y.proj)$p.value
            }
        }
    }
    return(colMeans(pvalue.null))

}
