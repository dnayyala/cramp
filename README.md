# CRAMP
**C**ovariance matrix testing using **RA**ndom **M**atrix **P**rojection

Estimation and hypothesis tests for the covariance matrix in high dimensions is a challenging problem as the traditional multivariate asymptotic theory is no longer valid. When the dimension is larger than or increasing with the sample size, standard likelihood based tests for the covariance matrix have poor performance. Existing high dimensional tests are either computationally expensive or have very weak control of type I error. This package provides the R functions for testing hypotheses involving one or more covariance matrices using random projections. Projecting the high dimensional data randomly into lower dimensional subspaces  alleviates of the *curse* of dimensionality, allowing for the use of traditional multivariate tests.

For the one-sample case, tests for uniformity and sphericity are provided. In the two sample case, the test is for equality of two covariance matrices. Functions to compute several traditional Gaussian likelihood-based and high-dimensional test statistics are provided in this package.

## Installation

To install the package, run the following commands:

```library(devtools)```

```install_github("dnayyala/cramp")```
## References

Ayyala, D. N., Ghosh, S. and Linder, D. F., *Covariance matrix testing in high dimension using random projects*, **Computational Statistics**, *37, 1111–1141 (2022)*. [Download](https://doi.org/10.1007/s00180-021-01166-4)
