
# SSDL

<!-- badges: start -->
<!-- badges: end -->
The SSDL package implements the Sketched Stochastic Dictionary Learning (SSDL) algorithm.
The goal of this algorithm is to extract a set of the most repetitive patterns (a.k.a dictionary) from a large-scale dataset. SSDL method operates on a compressed data version referred to as a data sketch. The [chickn](https://CRAN.R-project.org/package=chickn) package is used to carry out the data compression. SSDL operates with Filebacked Big Matrices (FBM) of [bigstatsr](https://github.com/privefl/bigstatsr) package to store and to access to matrix-like data that cannot be loaded in memory. The package provides the dictionary initialization routines based on the TV norm criterion and on the Compressive Orthogonal Matching Pursuit method. 

## Installation

You can install the released version of SSDL from [Gitlab](https://gitlab.com/Olga.Permiakova/ssdl/) with:

``` r
devtools::install_gitlab(repo = "Olga.Permiakova/ssdl")
```

## Example

This is a basic example which shows how to learn a dictionary using SSDL package:

``` r
library(SSDL)
library(chickn)
# Convert data matrix into a Filebacked Big Matrix
X = matrix(abs(rnorm(n = 1000)), ncol = 100, nrow = 10)
X_fbm = bigstatsr::as_FBM(X, backingfile = file.path(tempdir(), 'X_fbm'))$save()
# Compute the data sketch
## Generate frequency matrix
m = 256 # number of the frequency vectors
out_freq = chickn::GenerateFrequencies(Data = X_fbm, m = m, N0 = ncol(X_fbm),
                                       ncores = 1, niter= 5, nblocks = 2, sigma_start = 0.001)
W = out_freq$W # the obtained frequency matrix
## Compute the data sketch
SK= chickn::Sketch(X_fbm, W)
# Initialize a dictionary 
D0 = X_fbm[,sample(ncol(X_fbm),20)]
# Apply SSDL
result = SSDL::CDL(SK_Data = SK, Data = X_fbm, Frequencies = W, K = 20,
                   D =D0, pos.dic = TRUE, maxEpoch = 4, batch_size = 20,
                   lambda = 0, learn_rate = 0.2, alpha = 0.9, gamma = 0.1, ncores = 2,
                   DIR_tmp = tempdir())
# Obtained dictionary                   
D = result$D
#
```

