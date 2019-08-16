context("Test formulas")

test_that("HC1 and HC2 formulas match sandwich", {
    set.seed(42)
    y <- rnorm(10)
    x <- cbind(sin(1:10), sin(2:11))
    fm <- lm(y~x)

    HC1 <- sandwich::vcovHC(fm, type="HC1")
    HC2 <- sandwich::vcovHC(fm, type="HC2")
    r <- dfadjustSE(fm)

    expect_lt(max(abs(r$vcov-HC2)), 100*.Machine$double.eps)
    expect_lt(max(abs(r$coefficients[, "HC1 se"]-sqrt(diag(HC1)))),
              100*.Machine$double.eps)
})


test_that("New implementation matches old", {

## Compute the inverse square root of a symmetric matrix
MatSqrtInverse <- function(A) {

    ei <- eigen(A, symmetric=TRUE)
    if (min(ei$values) <= 0)
        warning("Gram matrix doesn't appear to be positive definite")

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
    df <- function(GG)
        sum(diag(GG))^2 / sum(GG * GG)
    ## Previously:
    ## lam <- eigen(GG, only.values=TRUE)$values
    ## sum(lam)^2/sum(lam^2)

    ## no clustering
    if(is.null(clustervar)) {
        Vhat <- sandwich::vcovHC(model, type="HC2")
        Vhat.Stata <- Vhat*NA

        M <- diag(n)-X %*% XXinv %*% t(X)       # annihilator matrix
        ## G'*Omega*G
        GOG <- function(ell) {
            Xtilde <- drop(X %*% XXinv %*% ell / sqrt(diag(M)))
            crossprod(M * Xtilde)
        }
    } else {
        if(!is.factor(clustervar)) stop("'clustervar' must be a factor")

        ## Stata
        S <- length(levels(clustervar)) # number clusters
        uj <- apply(u*X, 2, function(x) tapply(x, clustervar, sum))
        Vhat.Stata <- S/(S-1) * (n-1)/(n-K) *
            sandwich::sandwich(model, meat = crossprod(uj)/n)

        ## HC2
        tXs <- function(s) {
            Xs <- X[clustervar==s, , drop=FALSE] #nolint
            MatSqrtInverse(diag(NROW(Xs))-Xs%*% XXinv %*% t(Xs)) %*% Xs
        }
        tX <- lapply(levels(clustervar), tXs) # list of matrices

        tu <- split(u, clustervar)
        tutX <- sapply(seq_along(tu), function(i) crossprod(tu[[i]], tX[[i]]))
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
            rho <- max((sum(sapply(seq_along(tu), function(i)
                sum(tu[[i]] %o% tu[[i]])))-ssr) / (sum(ns^2)-n), 0)
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
        se.Stata <- drop(sqrt(crossprod(ell, Vhat.Stata) %*% ell))
    } else {
        se <- sqrt(diag(Vhat))
        dof <- sapply(seq(K), function(k) df(GOG(diag(K)[, k])))
        se.Stata <- sqrt(diag(Vhat.Stata))
    }
    names(dof) <- names(se)

    list(vcov=Vhat, dof=dof, adj.se=se*qt(0.975, df=dof)/qnorm(0.975),
                se=se, se.Stata=se.Stata)
}


    set.seed(42)
    cl3 <- c(rep(1, 60), rep(2, 20), rep(3, 20))
    d0 <- list(data.frame(y=rnorm(10),
                     x=cbind(sin(1:10), cos(2:11)),
                     cl=as.factor(c(rep(1, 6), rep(2, 2), rep(3, 2)))),
               data.frame(y=rnorm(100),
                     x=cbind(sin(1:100), rnorm(100)),
                     cl=as.factor(c(rep(1, 60), rep(2, 20), rep(3, 20)))),
               data.frame(y=rnorm(100)+cl3,
                     x=cbind(sin(1:100), rnorm(100)),
                     cl=as.factor(cl3)))
    ep <- 50*.Machine$double.eps
    for (j in seq_along(d0)) {
        fm <- lm(y~x.1+x.2, data=d0[[j]])

        ## No clustering
        r <- dfadjustSE(fm)
        rold <- BMlmSE(fm)
        expect_lt(max(abs(r$vcov-rold$vcov)), ep)
        expect_lt(max(abs(r$coefficients[, "df"]-rold$dof)), 10*ep)
        expect_lt(max(abs(r$coefficients[, "Adj. se"]-rold$adj.se)), ep)
        expect_lt(max(abs(r$coefficients[, "HC2 se"]-rold$se)), ep)

        ## If each observation in its cluster, we get HC2
        cl <- as.factor(seq_along(fm$model$y))
        expect_lt(max(abs(dfadjustSE(fm, cl, IK=FALSE)$coefficients-
                                        r$coefficients)), 10*ep)

        ## Clustering
        r <- dfadjustSE(fm, d0[[j]]$cl, IK=FALSE)
        rold <- BMlmSE(fm, d0[[j]]$cl, IK=FALSE)
        expect_lt(max(abs(r$vcov-rold$vcov)), ep)
        expect_lt(max(abs(r$coefficients[, "HC2 se"]-rold$se)),
                  ep)
        expect_lt(max(abs(r$coefficients[, "df"]-rold$dof)), ep)
        expect_lt(max(abs(r$coefficients[, "Adj. se"]-rold$adj.se)), ep)
        expect_lt(max(abs(r$coefficients[, "HC1 se"]-rold$se.Stata)), ep)

        r <- dfadjustSE(fm, d0[[j]]$cl, IK=TRUE, rho0=TRUE)
        rold <- BMlmSE(fm, d0[[j]]$cl, IK=TRUE)
        expect_lt(max(abs(r$vcov-rold$vcov)), ep)
        expect_lt(max(abs(r$coefficients[, "HC2 se"]-rold$se)), ep)
        expect_lt(max(abs(r$coefficients[, "df"]-rold$dof)), ep)
        expect_lt(max(abs(r$coefficients[, "Adj. se"]-rold$adj.se)), ep)
        expect_lt(max(abs(r$coefficients[, "HC1 se"]-rold$se.Stata)), ep)

        r <- dfadjustSE(fm, d0[[j]]$cl, IK=TRUE, ell=c(1, 1, 0), rho0=TRUE)
        rold <- BMlmSE(fm, d0[[j]]$cl, IK=TRUE, ell=c(1, 1, 0))
        expect_lt(max(abs(r$vcov-rold$vcov)), ep)
        expect_lt(max(abs(r$coefficients[, "HC2 se"]-rold$se)), ep)
        expect_lt(max(abs(r$coefficients[, "df"]-rold$dof)), ep)
        expect_lt(max(abs(r$coefficients[, "Adj. se"]-rold$adj.se)), ep)
        expect_lt(max(abs(r$coefficients[, "HC1 se"]-rold$se.Stata)), ep)
    }

    r <- dfadjustSE(lm(d0[[3]]$y~1))
    rold <- BMlmSE(lm(d0[[3]]$y~1))
    expect_lt(max(abs(r$vcov-rold$vcov)), ep)
    expect_lt(max(abs(r$coefficients[, "df"]-rold$dof)), 10*ep)
    expect_lt(max(abs(r$coefficients[, "Adj. se"]-rold$adj.se)), ep)
    expect_lt(max(abs(r$coefficients[, "HC2 se"]-rold$se)), ep)
    expect_equal(capture.output(print(r, digits=3))[4],
                 "(Intercept)     1.64  0.128  0.128    0.13 99 9.05e-37")

    ## Noninvertible cases
    d0[[3]]$x.3 <- FALSE
    d0[[3]]$x.3[1] <- TRUE
    fm1 <- lm(y~x.1+x.2+x.3, data=d0[[3]])
    r <- dfadjustSE(fm1)
    r0 <- dfadjustSE(lm(y~x.1+x.2, data=d0[[3]]))
    expect_lt(max(abs(r0$coefficients[, 1:4]-r$coefficients[1:3, 1:4])),
              0.05)
    expect_true(all(is.nan(BMlmSE(fm1)$se)))
    ## Fixed effects: here minimum eigenvalue is only numerically positive
    d0[[1]]$x <- c(rep(1, 4), rep(0, 5), 1)
    fm2 <- lm(y~x+cl, data=d0[[1]])
    p <- dfadjustSE(fm2, d0[[1]]$cl)$coefficients[, "Adj. se"]
    expect_lt(max(abs(BMlmSE(fm2, d0[[1]]$cl)$adj.se-p)), 1e-6)
})
