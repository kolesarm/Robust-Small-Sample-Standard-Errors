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
    expect_lt(max(abs(r$se.Stata-sqrt(diag(HC1)))), 100*.Machine$double.eps)
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
            Xs <- X[clustervar==s, , drop=FALSE]
            MatSqrtInverse(diag(NROW(Xs))-Xs%*% XXinv %*% t(Xs)) %*% Xs
        }
        tX <- lapply(levels(clustervar), tXs) # list of matrices

        tu <- split(u, clustervar)
        tutX <- sapply(seq_along(tu), function(i) crossprod(tu[[i]], tX[[i]]))
        Vhat <- sandwich::sandwich(model, meat = tcrossprod(tutX)/n)

        ## DOF adjustment
        tHs <- function(s) {
            Xs <- X[clustervar==s, , drop=FALSE]
            index <- which(clustervar==s)
            ss <- outer(rep(0, n), index)     # n x ns matrix of 0
            ss[cbind(index, 1:length(index))] <- 1
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
    y <- rnorm(10)
    x <- cbind(sin(1:10), cos(2:11))
    clustervar <- as.factor(c(rep(1, 6), rep(2, 2), rep(3, 2)))
    fm <- lm(y~x)

    ## No clustering
    r <- dfadjustSE(fm)
    rold <- BMlmSE(fm)
    expect_lt(max(abs(r$vcov-rold$vcov)), 50*.Machine$double.eps)
    expect_lt(max(abs(r$dof-rold$dof)), 50*.Machine$double.eps)
    expect_lt(max(abs(r$adj.se-rold$adj.se)), 50*.Machine$double.eps)
    expect_lt(max(abs(r$se-rold$se)), 50*.Machine$double.eps)
    ## Clustering
    r <- dfadjustSE(fm, clustervar)
    rold <- BMlmSE(fm, clustervar)
    expect_lt(max(abs(r$vcov-rold$vcov)), 50*.Machine$double.eps)
    expect_lt(max(abs(r$se-rold$se)), 50*.Machine$double.eps)
    ## expect_lt(max(abs(r$adj.se-rold$adj.se)), 50*.Machine$double.eps)
    ## expect_lt(max(abs(r$dof-rold$dof)), 50*.Machine$double.eps)
    ## TODO: IK


})
