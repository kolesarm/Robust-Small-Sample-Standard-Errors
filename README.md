# dfadjust

<!-- badges: start -->
<!-- badges: end -->

## Description

The =R= function ~BMlmSE(model, clustervar=NULL, ell=NULL, IK=TRUE)~ implements
the Bell-McCaffrey standard error degrees-of-freedom adjustment and associated
standard errors for regression parameters in the linear regression
Y_i=X_i'\beta+\epsilon_i as described in [[http://www.mitpressjournals.org/doi/abs/10.1162/REST_a_00552][Imbens and Kolesár (2016), ReStat]]

SAS version by Nicolas Moreau is available [[http://cemoi.univ-reunion.fr/econometrie-avec-r-et-sas][here]].



## Dependencies
The code needs the =sandwich= package to be installed

## Arguments

- The first argument, =model=, takes a model fitted using the =lm()= function
- =clustervar= is a factor variable that defines clusters. The command computes
  cluster-robust standard errors if the variable is supplied, otherwise it
  computes heteroscedasticity-robust standard errors.
- The vector =ell=, of the same length as the dimension of the covariates,
  specifies which linear combination of coefficients to compute. If supplied,
  the command computes the standard error for =ell='\beta. Otherwise the
  function computes standard errors for each element of \beta.
- The flag =IK=, only relevant if cluster-robust standard errors are being
  computed, specifies whether to compute the degrees-of-freedom adjustment using
  the Imbens-Kolesár method (if =TRUE=), or the Bell-McCaffrey method (if =FALSE=)

## Value
The function returns a list with the following components:
- =vcov= :: Variance-covariance matrix estimator. For the case without
            clustering, it corresponds to the HC2 estimator (see MacKinnon and
            White, 1985 and the reference manual for the =sandwich= package).
            For the case with clustering, it corresponds to a generalization of
            the HC2 estimator, called LZ2 in Imbens and Kolesár.
- =dof= :: Degrees-of-freedom adjustment
- =adj.se= :: Adjusted standard errors. For \beta_j, they are defined as
              : adj.se[j]=sqrt(vcov[j,j]se*qt(0.975,df=dof)
              so that the Bell-McCaffrey confidence intervals are given as
              =coefficients(fm)[j] +- 1.96* adj.se=
- =se.Stata= :: Square root of the cluster-robust variance estimator used in
                =STATA=.



## Installation

You can install the released version of dfadjust from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("dfadjust")
```

## Example

** Examples
No clustering:
``` r
set.seed(42)
x <- sin(1:10)
y <- rnorm(10)
fm <- lm(y~x)
BMlmSE(fm)
#+end_src
Clustering:
#+begin_src R
clustervar <- as.factor(c(rep(1,6),rep(2,2),rep(3,2)))
BMlmSE(fm, clustervar)
```
Here we defined the first six observations to be in cluster 1, the next two in
cluster 2, and the last three in cluster three.
