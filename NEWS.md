# dfadjust 1.0.0

## New Features

- The function `dfadjustSE` implements small-sample degrees of freedom
  adjustment discussed in [Imbens and Kolesár
  (2016)](http://www.mitpressjournals.org/doi/abs/10.1162/REST_a_00552), using
  both heteroskedasticity-robust and clustered standard errors. For clustered
  standard errors, the package implements both the Imbens and Kolesár (2016) and
  the Bell and McCaffrey (2002, Survey Methodology) degrees of freedom
  adjustments.
- This implementation can handle models with fixed effects, as well as datasets
  with a large number of observations (for heteroskedasticity-robust standard
  errors) or datasets with large clusters (for clustered standard errors)
