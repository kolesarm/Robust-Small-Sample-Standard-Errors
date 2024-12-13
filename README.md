[![R-CMD-check](https://github.com/kolesarm/Robust-Small-Sample-Standard-Errors/workflows/R-CMD-check/badge.svg)](https://github.com/kolesarm/Robust-Small-Sample-Standard-Errors/actions) [![Coverage status](https://codecov.io/gh/kolesarm/Robust-Small-Sample-Standard-Errors/branch/master/graph/badge.svg)](https://app.codecov.io/github/kolesarm/Robust-Small-Sample-Standard-Errors?branch=master) [![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/dfadjust)](https://cran.r-project.org/package=dfadjust) [![Download statistics](https://cranlogs.r-pkg.org/badges/grand-total/dfadjust)](https://cran.r-project.org/package=dfadjust)

# dfadjust

This package implements the small-sample degrees of freedom adjustments for
robust and cluster-robust standard errors in linear regression described in
[Imbens and Kolesár (2016)](https://doi.org/10.1162/REST_a_00552).

SAS version by Nicolas Moreau is available
[here](http://cemoi.univ-reunion.fr/econometrie-avec-r-et-sas).

See vignette [dfadjust](doc/dfadjust.pdf) for description of the package
(available through `vignette("dfadjust")` once package is installed), and the
package [manual](doc/manual.pdf) for documentation of the package functions.


## Example

No clustering:
``` r
set.seed(42)
x <- sin(1:10)
y <- rnorm(10)
fm <- lm(y~x)
dfadjustSE(fm)
```
Clustering:
``` r
clustervar <- as.factor(c(rep(1, 6), rep(2, 2), rep(3, 2)))
dfadjustSE(fm, clustervar)
```
Here we defined the first six observations to be in cluster 1, the next two in
cluster 2, and the last three in cluster three.

The package handles cluster fixed effects, and large clusters. Computing the
adjustment with one million observations and 50 clusters, takes about 5 seconds:

``` r
N <- 10^6
x <- sin(1:N)
y <- seq(N)
clustervar <- as.factor(rep(1:50, each=N/50))
fm <- lm(y~x+clustervar)
## Inference on x, i.e. second coefficient
dfadjustSE(fm, ell=2, clustervar=clustervar)
```

## Installation

You can install the released version of `dfadjust` from
[CRAN](https://CRAN.R-project.org/package=dfadjust) with:

``` r
install.packages("dfadjust")
```

Alternatively, you can get the current development version from GitHub:
``` r
if (!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("kolesarm/Robust-Small-Sample-Standard-Errors")
```
