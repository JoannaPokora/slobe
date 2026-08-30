# slobe

`slobe` is an R package for fitting Adaptive Bayesian SLOPE (ABSLOPE) model using the SLOBE algorithm.

The package supports Gaussian and logistic regression and provides an R interface for the updated implementation of ABSLOPE. The method combines the SLOPE penalty with an adaptive spike-and-slab prior to perform variable selection and coefficient estimation.

## Installation

The package can be installed directly from GitHub:

```r
# install.packages("remotes")
remotes::install_github("JoannaPokora/slobe")
