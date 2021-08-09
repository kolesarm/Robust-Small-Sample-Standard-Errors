# dfadjust 1.0.3

## Minor improvements and fixes

- Fix incorrect computation of p-values in the `print.dfadjustSE` method

# dfadjust 1.0.2

## Minor improvements and fixes

- Fix incorrect computation of CR2 variance estimator and degrees of freedom
  adjustment if data not sorted by cluster

# dfadjust 1.0.1

## Minor improvements and fixes

- Fix problem with failing tests when platform didn't use long double

# dfadjust 1.0.0

## New Features

- The function `dfadjustSE` implements small-sample degrees of freedom
  adjustment discussed in [Imbens and Kolesár
  (2016)](https://www.doi.org/10.1162/REST_a_00552), using
  both heteroskedasticity-robust and clustered standard errors. For clustered
  standard errors, the package implements both the Imbens and Kolesár (2016) and
  the Bell and McCaffrey (2002, Survey Methodology) degrees of freedom
  adjustments.
- This implementation can handle models with fixed effects, as well as datasets
  with a large number of observations (for heteroskedasticity-robust standard
  errors) or datasets with large clusters (for clustered standard errors)
