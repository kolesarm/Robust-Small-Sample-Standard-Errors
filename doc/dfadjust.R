## ----include=FALSE, cache=FALSE-----------------------------------------------
library("knitr")
knitr::opts_knit$set(self.contained = FALSE)
knitr::opts_chunk$set(tidy = TRUE, collapse=TRUE, comment = "#>",
                      tidy.opts=list(blank=FALSE, width.cutoff=55))
oldoptions <- options(digits=3)

## ----setup--------------------------------------------------------------------
library(dfadjust)

## -----------------------------------------------------------------------------
set.seed(7)
d1 <- data.frame(y=rnorm(1000), x1=c(rep(1, 3), rep(0, 997)),
                 x2=c(rep(1, 150), rep(0, 850)),
                 x3=rnorm(1000),
                 cl=as.factor(c(rep(1:10, each=50), rep(11, 500))))

## -----------------------------------------------------------------------------
r1 <- lm(y~x1, data=d1)
## No clustering
dfadjustSE(r1)

## -----------------------------------------------------------------------------
r1 <- lm(y~x2, data=d1)
# Default Imbens-Kolesár method
dfadjustSE(r1, clustervar=d1$cl)
# Bell-McCaffrey method
dfadjustSE(r1, clustervar=d1$cl, IK=FALSE)

## -----------------------------------------------------------------------------
r1 <- lm(y~x3+cl, data=d1)
dfadjustSE(r1, clustervar=d1$cl, ell=c(0, 1, rep(0, r1$rank-2)))
dfadjustSE(r1, clustervar=d1$cl, ell=c(0, 1, rep(0, r1$rank-2)), IK=FALSE)

## -----------------------------------------------------------------------------
d2 <- do.call("rbind", replicate(500, d1, simplify = FALSE))
d2$y <- rnorm(length(d2$y))
r2 <- lm(y~x2, data=d2)
summary(r2)
# Default Imbens-Kolesár method
dfadjustSE(r2, clustervar=d2$cl)
# Bell-McCaffrey method
dfadjustSE(r2, clustervar=d2$cl, IK=FALSE)

## ----cleanup, include=FALSE---------------------------------------------------
options(oldoptions)

