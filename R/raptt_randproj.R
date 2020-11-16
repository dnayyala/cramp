raptt.randproj <- function(nrow, ncol, seed = NULL){
    #' RAPTT projection matrix generator
    #'
    #' Generate random matrices of given dimension which have full row and are orthogonal as defined in Srivastava et al. (2016)
    #' @import mvtnorm stats
    NULL
    #' @export

    #' @param nrow Number of rows in the random matrix to be generated (dimension of projected space)
    #' @param ncol Number of colums (dimension of original space)
    #' @param seed Set the seed to replicate the random matrix generated. Default is \code{NULL} and no seed is used.

    #' @return An orthogonal matrix of with \code{nrow} rows and \code{ncol} columns.
    #' @examples
    #' a = raptt.randproj(nrow = 5, ncol = 10)

    #' @references Srivastava, R., Li, P. and Ruppert, D. RAPTT: An Exact Two-Sample Test in High Dimensions Using Random Projections. Journal of Computational and Graphical Statistics, 25 (3), 2016.
    if (!is.null(seed)){
        set.seed(seed)
    }

    if (nrow > ncol){
        stop("Cannot project into a higher dimensional space. The number of rows should be smaller than the number of columns")
    }
    block.length <- ceiling(ncol/nrow);

    random.samples <- rnorm(n = ncol)
    permutation <- sample(ncol, ncol, replace = FALSE)

    random.matrix <- matrix(0, nrow = nrow, ncol = ncol)
    ## Fill the matrix blocks with random samples
    for (i in 1:(nrow - 1)){
        range <- (i-1)*block.length + (1:block.length)
        random.matrix[i, range] <- random.samples[range]
        random.matrix[i,] <- random.matrix[i,]/sqrt(sum(random.matrix[i,]^2))
    }
    range <- ((nrow - 1)*block.length + 1):ncol
    random.matrix[nrow, range] <- random.samples[range]
    random.matrix[nrow,] <- random.matrix[nrow,]/sqrt(sum(random.matrix[nrow,]^2))

    ## Permute the blocks
    random.matrix <- random.matrix[, permutation]

    return(random.matrix)
}
