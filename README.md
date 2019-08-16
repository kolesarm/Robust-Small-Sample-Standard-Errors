[![Travis build status](https://travis-ci.org/kolesarm/Robust-Small-Sample-Standard-Errors.svg?branch=master)](https://travis-ci.org/kolesarm/Robust-Small-Sample-Standard-Errors) [![Coverage status](https://codecov.io/gh/kolesarm/Robust-Small-Sample-Standard-Errors/branch/master/graph/badge.svg)](https://codecov.io/github/kolesarm/Robust-Small-Sample-Standard-Errors?branch=master)

# dfadjust

This package implements the small-sample degrees of freedom adjustments for
robust and cluster-robust standard errors in linear regression described in
[Imbens and Kolesár (2016)](http://www.mitpressjournals.org/doi/abs/10.1162/REST_a_00552).

SAS version by Nicolas Moreau is available [here](http://cemoi.univ-reunion.fr/econometrie-avec-r-et-sas).

See vignette [dfadjust](doc/dfadjust.pdf) for description of the package
(available through `vignette("dfadjust")` once package is installed), and the
package [manual](doc/manual.pdf) for documentation of the package functions.


## Installation

<!-- You can install the released version of dfadjust from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->
<!-- install.packages("dfadjust") -->
<!-- ``` -->

You can install the package manually by downloading the source code here, or
using the function `install_github()` from the `devtools` package:

``` r
install.packages("devtools") ## if devtools package not installed
devtools::install_github("https://github.com/kolesarm/Robust-Small-Sample-Standard-Errors")
```

## Example

No clustering:
``` r
set.seed(42)
x <- sin(1:10)
y <- rnorm(10)
fm <- lm(y~x)
BMlmSE(fm)
```
Clustering:
``` r
clustervar <- as.factor(c(rep(1,6),rep(2,2),rep(3,2)))
dfadjustSE(fm, clustervar)
```
Here we defined the first six observations to be in cluster 1, the next two in
cluster 2, and the last three in cluster three.
