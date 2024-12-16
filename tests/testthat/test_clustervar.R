context("Test clustering")

test_that("Test clustering if clusters out of order", {

    dt <- data.frame(y=c(cos(0:9), tan(1:10)),
                     x=c(sin(1:10), sin(2:11)),
                     clusters=factor(c(2, 1, 3, 2, 2, 3, 3, 1, 1, 2)))
    ## foreign::write.dta(dt, "dt.dta")
    ## Stata: foreign::write.dta(dt, "dt.dta")
    ## eststo r0: reg y x, cluster(clusters)
    ## eststo r1: reg y x, r
    ## eststo r2: reg y x, vce(hc2)
    ## esttab r0 r1 r2, b(%11.10f) se(%11.10f)

    r0 <- lm(y~x, data=dt)
    r1 <- dfadjustSE(r0, cluster=dt$clusters)
    ## Check against stata
    expect_equal(as.vector(r1$coefficients[, 1:2]),
                 c(-0.4456202905, 0.3241183110, 0.4827602183, 0.1814501834))
    expect_lt(max(abs(r1$coefficients[, 4]-c(1.053441197, 0.533126248))), 1e-9)
    expect_lt(max(abs(r1$coefficients[, 5]-c(2.001430933, 1.783242385))), 1e-9)

    ## Check same result as when ordered
    dt2 <- dt[order(dt$clusters), ]
    r2 <- dfadjustSE(lm(y~x, data=dt2), cluster=dt2$clusters)
    expect_equal(r1, r2)

    ## check non-clustered data
    r1 <- dfadjustSE(lm(y~x, data=dt), cluster=as.factor(1:20)[c(10:1, 11:20)],
                     IK=FALSE)
    r2 <- dfadjustSE(lm(y~x, data=dt), cluster=as.factor(1:20), IK=FALSE)
    r3 <- dfadjustSE(lm(y~x, data=dt))
    r4 <- dfadjustSE(lm(y~x, data=dt), cluster=as.factor(1:20)[c(10:1, 11:20)],
                     IK=TRUE)
    r5 <- dfadjustSE(lm(y~x, data=dt), cluster=as.factor(1:20), IK=TRUE)

    expect_equal(r1, r2)
    expect_equal(r1, r3)
    expect_equal(r4, r5)
    expect_equal(r1$coefficients, r5$coefficients)
    ## Check against Stata
    expect_equal(unname(r5$coefficients[, "HC2 se"]),
                 c(0.4227339762, 0.4421048122))
    expect_equal(unname(r5$coefficients[, "HC1 se"]),
                 c(0.4278077455, 0.4392155956))
    expect_equal(unname(r5$coefficients[, 5]), c(17.6142349, 12.10389102))

})

test_that("ell specifications", {
    x <- sin(1:100)
    y <- 1:100
    fm <- lm(y ~ x + I(x^2))
    clustervar <- as.factor(c(rep(1, 40), rep(5, 20),
                              rep(2, 20), rep(3, 10), rep(4, 10)))
    r0 <- dfadjustSE(fm, clustervar, ell=2)
    r1 <- dfadjustSE(fm, clustervar, ell=c(0, 1, 0))
    expect_equal(drop(unname(r1$coefficients-r0$coefficients)), rep(0, 5))

    r2 <- dfadjustSE(fm, clustervar, ell=c(1, 3))
    r3 <- dfadjustSE(fm, clustervar)
    expect_equal(r3$coefficients[c(1, 3), ], r2$coefficients)
})

test_that("collinear specifications", {
    x <- sin(1:100)
    y <- 1:100

    fm1 <- lm(y ~ x)
    fm2 <- lm(y ~ x +I(x))
    r1 <- dfadjustSE(fm1, clustervar)
    r2 <- dfadjustSE(fm2, clustervar)
    expect_equal(r1, r2)
    fm3 <- lm(y ~ x +I(x^2) + I(x) + I(x^2+0)+I(sin(x)))
    fm4 <- lm(y ~ x +I(x^2) + I(sin(x)))
    r3 <- dfadjustSE(fm3, clustervar)
    r4 <- dfadjustSE(fm4, clustervar)
    expect_equal(r4, r3)
    r5 <- dfadjustSE(fm3, clustervar, ell=6)
    r6 <- dfadjustSE(fm4, clustervar, ell=4)
    expect_equal(r5, r6)
})
