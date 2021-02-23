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
})
