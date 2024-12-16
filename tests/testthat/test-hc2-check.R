context("Test formulas")

elt <- function(a, ep) testthat::expect_lt(max(abs(a)), ep)

test_that("HC1 and HC2 formulas match sandwich", {
    y <- 1:10
    x <- cbind(sin(1:10), sin(2:11))
    fm <- lm(y~x)

    HC1 <- sandwich::vcovHC(fm, type="HC1")
    HC2 <- sandwich::vcovHC(fm, type="HC2")
    r <- dfadjustSE(fm)

    elt(r$vcov-HC2, 10^3*.Machine$double.eps)
    elt(r$coefficients[, "HC1 se"]-sqrt(diag(HC1)), 10^3*.Machine$double.eps)
    elt(r$coefficients[, "df"]-c(4.94624731, 3.90754737, 4.77265484),
        1e-8)
})

test_that("clustervar must be a factor", {
    y <- 1:10
    x <- cbind(sin(1:10), sin(2:11))
    fm <- lm(y~x)
    expect_error(dfadjustSE(fm, clustervar=1:10))
})


## Compute the inverse square root of a symmetric matrix
MatSqrtInverse <- function(A) {
    ei <- eigen(A, symmetric=TRUE)
    d <- pmax(ei$values, 0)
    d2 <- 1/sqrt(d)
    d2[d == 0] <- 0
    ## diag(d2) is d2 x d2 identity if d2 is scalar, instead we want 1x1 matrix
    ei$vectors %*% (if (length(d2)==1) d2 else diag(d2)) %*% t(ei$vectors)
}

BMlmSE <- function(model, clustervar=NULL, ell=NULL, IK=TRUE) {
    X <- model.matrix(model)
    sum.model <- summary.lm(model)
    n <- sum(sum.model$df[1:2])
    K <- model$rank
    XXinv <- sum.model$cov.unscaled # XX^{-1}
    u <- residuals(model)

    ## Compute DoF given G'*Omega*G without calling eigen as suggested by
    ## Winston Lin
    df <- function(GG) sum(diag(GG))^2 / sum(GG * GG)
    ## Previously:
    ## lam <- eigen(GG, only.values=TRUE)$values
    ## sum(lam)^2/sum(lam^2)

    ## no clustering
    if (is.null(clustervar)) {
        Vhat <- sandwich::vcovHC(model, type="HC2")
        VhatStata <- Vhat*NA

        M <- diag(n)-X %*% XXinv %*% t(X)       # annihilator matrix
        ## G'*Omega*G
        GOG <- function(ell) {
            Xtilde <- drop(X %*% XXinv %*% ell / sqrt(diag(M)))
            crossprod(M * Xtilde)
        }
    } else {
        if (!is.factor(clustervar)) stop("'clustervar' must be a factor")

        ## Stata
        S <- length(levels(clustervar)) # number clusters
        uj <- apply(u*X, 2, function(x) tapply(x, clustervar, sum))
        VhatStata <- S / (S-1) * (n-1) / (n-K) *
            sandwich::sandwich(model, meat = crossprod(uj)/n)

        ## HC2
        tXs <- function(s) {
            Xs <- X[clustervar==s, , drop=FALSE] #nolint
            MatSqrtInverse(diag(NROW(Xs))-Xs%*% XXinv %*% t(Xs)) %*% Xs
        }
        tX <- lapply(levels(clustervar), tXs) # list of matrices

        tu <- split(u, clustervar)
        tutX <- vapply(seq_along(tu), function(i) crossprod(tu[[i]], tX[[i]]),
                       numeric(K))
        Vhat <- sandwich::sandwich(model, meat = tcrossprod(tutX)/n)

        ## DOF adjustment
        tHs <- function(s) {
            Xs <- X[clustervar==s, , drop=FALSE] #nolint
            index <- which(clustervar==s)
            ss <- outer(rep(0, n), index)     # n x ns matrix of 0
            ss[cbind(index, seq_along(index))] <- 1
            ss-X %*% XXinv %*% t(Xs)
        }
        tH <- lapply(levels(clustervar), tHs) # list of matrices

        Moulton <- function() {
            ## Moulton estimates
            ns <- tapply(u, clustervar, length)
            ssr <- sum(u^2)
            s1 <- vapply(seq_along(tu), function(i) sum(tu[[i]] %o% tu[[i]]),
                         numeric(1))
            rho <- max((sum(s1)-ssr) / (sum(ns^2)-n), 0)
            c(sig.eps=max(ssr/n - rho, 0), rho=rho)
        }

        GOG <- function(ell) {
            G <- sapply(seq_along(tX),
                        function(i)  tH[[i]] %*% tX[[i]] %*% XXinv %*% ell)
            GG <- crossprod(G)

            ## IK method
            if (IK==TRUE) {
                Gsums <- apply(G, 2,
                               function(x) tapply(x, clustervar, sum)) # Z'*G
                GG <- Moulton()[1]*GG + Moulton()[2]*crossprod(Gsums)
            }
            GG
        }
    }

    if (!is.null(ell)) {
        se <- drop(sqrt(crossprod(ell, Vhat) %*% ell))
        dof <- df(GOG(ell))
        seStata <- drop(sqrt(crossprod(ell, VhatStata) %*% ell))
    } else {
        se <- sqrt(diag(Vhat))
        dof <- vapply(seq_len(K), function(k) df(GOG(diag(K)[, k])), numeric(1))
        seStata <- sqrt(diag(VhatStata))
    }
    names(dof) <- names(se)

    list(vcov=Vhat, dof=dof, adj.se=se*qt(0.975, df=dof)/qnorm(0.975), se=se,
         seStata=seStata)
}

test_that("New implementation matches old", {
    cl3 <- c(rep(1, 60), rep(2, 20), rep(3, 20))
    d0 <- list(data.frame(y=tan(1:10),
                          x=cbind(sin(1:10), cos(2:11)),
                          cl=as.factor(c(rep(1, 6), rep(2, 2), rep(3, 2)))),
               data.frame(y=1:100,
                          x=cbind(sin(1:100), (1:100)^3),
                          cl=as.factor(c(rep(1, 60), rep(2, 20), rep(3, 20)))),
               data.frame(y = (1:100)^2+cl3,
                          x=cbind(sin(1:100), 1:100),
                          cl=as.factor(cl3)))
    ep <- 100*.Machine$double.eps
    ## Increase numerical tolerance if platform doesn't have long double.
    ## BLAS/LAPACK 3.10.1 in R-devel also seemes to be creating issues on M1 Mac
    bigep <- if (capabilities("long.double")) 10*ep else 10^4*ep

    for (j in seq_along(d0)) {
        fm <- lm(y~x.1+x.2, data=d0[[j]])

        ## No clustering
        r <- dfadjustSE(fm)
        ro <- BMlmSE(fm)
        elt((r$vcov-ro$vcov)/r$vcov, ep)
        elt(r$coefficients[, "df"]-ro$dof, bigep)
        elt(r$coefficients[, "Adj. se"]-ro$adj.se, bigep)
        elt((r$coefficients[, "HC2 se"]-ro$se)/ro$se, ep)
        elt(r$coefficients[, 5]-ro$dof, bigep)
        ## If each observation in its cluster, we get HC2
        cl <- as.factor(seq_along(fm$model$y))
        elt(dfadjustSE(fm, cl, IK=FALSE)$coefficients-r$coefficients, bigep)

        ## Clustering
        r <- dfadjustSE(fm, d0[[j]]$cl, IK=FALSE)
        ro <- BMlmSE(fm, d0[[j]]$cl, IK=FALSE)
        elt((r$vcov-ro$vcov)/r$vcov, bigep)

        elt((r$coefficients[, "HC2 se"]-ro$se)/ro$se, bigep)
        elt((r$coefficients[, "df"]-ro$dof)/ro$dof, bigep)
        elt((r$coefficients[, "Adj. se"]-ro$adj.se) /
                r$coefficients[, "Adj. se"], 10*bigep)
        elt((r$coefficients[, "HC1 se"]-ro$seStata)/ro$seStata, ep)

        r <- dfadjustSE(fm, d0[[j]]$cl, IK=TRUE, rho0=TRUE)
        ro <- BMlmSE(fm, d0[[j]]$cl, IK=TRUE)
        elt((r$vcov-ro$vcov)/r$vcov, bigep)
        elt((r$coefficients[, "HC2 se"]-ro$se)/ro$se, bigep)
        elt(r$coefficients[, "df"]-ro$dof, ep)
        elt((r$coefficients[, "Adj. se"]-ro$adj.se)/ro$adj.se, bigep)
        elt((r$coefficients[, "HC1 se"]-ro$seStata)/ro$seStata, ep)

        r <- dfadjustSE(fm, d0[[j]]$cl, IK=TRUE, ell=c(1, 1, 0), rho0=TRUE)
        ro <- BMlmSE(fm, d0[[j]]$cl, IK=TRUE, ell=c(1, 1, 0))
        elt((r$vcov-ro$vcov)/r$vcov, bigep)
        elt((r$coefficients[, "HC2 se"]-ro$se)/ro$se, bigep)
        elt(r$coefficients[, "df"]-ro$dof, ep)
        elt((r$coefficients[, "Adj. se"]-ro$adj.se)/ro$adj.se, bigep)
        elt((r$coefficients[, "HC1 se"]-ro$seStata)/ro$seStata, ep)
    }

    r <- dfadjustSE(lm(d0[[3]]$y~1))
    ro <- BMlmSE(lm(d0[[3]]$y~1))
    elt((r$vcov-ro$vcov)/ro$vcov, bigep)
    elt(r$coefficients[, "df"]-ro$dof, bigep)
    elt((r$coefficients[, "Adj. se"]-ro$adj.se)/ro$adj.se, ep)
    elt((r$coefficients[, "HC2 se"]-ro$se)/ro$se, ep)
    ## previously incorrectly 9.05e-37 for p-value
    expect_equal(capture.output(print(r, digits=3))[4],
                 "(Intercept)     3385    303    303     306 99 2.83e-19")

    ## Noninvertible cases
    d1 <- d0[[3]]
    d1$x.3 <- FALSE
    d1$x.3[1] <- TRUE
    fm1 <- lm(y~x.1+x.2+x.3, data=d1)
    rq <- dfadjustSE(fm1)
    r0 <- dfadjustSE(lm(y~x.1+x.2, data=d1))
    elt(rq$coefficients[,  4]-c(188.97371028, 107.8865573,
                                3.36194461, 211.67556014), 1e-8)
    elt(r0$coefficients[, 4]-c(188.164883779, 109.636001038,
                               3.371563406), 1e-9)
    expect_warning(rr <- BMlmSE(fm1)$se)
    expect_true(all(is.nan(rr)))
    ## Fixed effects: here minimum eigenvalue is only numerically positive
    d2 <- d0[[1]]
    d2$x <- c(rep(1, 4), rep(0, 5), 1)
    fm2 <- lm(y~x+cl, data=d2)
    p1 <- dfadjustSE(fm2, d2$cl)
    p2 <- BMlmSE(fm2, d2$cl)
    elt(p2$adj.se - p1$coefficients[, "Adj. se"], 1e-6)

    ## P-values
    expect_equal(capture.output(print(p1, digits=3))[5],
                 "x               1.71 0.3501 0.3705     2.4  1  0.1362")

    ## Test scaling
    p3 <- dfadjustSE(lm.fit(x=model.matrix(fm2),
                            y=1e-8*fm2$model$y), d2$cl)
    p4 <- dfadjustSE(lm.fit(x=model.matrix(fm2),
                            y=1e8*fm2$model$y), d2$cl)
    expect_lt(max(abs(p3$vcov*1e16-p1$vcov)), 1e-12)
    expect_lt(max(abs(p4$vcov*1e-16-p1$vcov)), 1e-12)
    p5 <- dfadjustSE(lm.fit(x=model.matrix(fm2)*1e8,
                            y=fm2$model$y), d2$cl)
    p6 <- dfadjustSE(lm.fit(x=model.matrix(fm2)*1e-8,
                            y=fm2$model$y), d2$cl)
    expect_lt(max(abs(p5$vcov*1e16-p1$vcov)), 1e-12)
    expect_lt(max(abs(p6$vcov*1e-16-p1$vcov)), 1e-12)
})
