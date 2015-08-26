## Compute Bell-McCaffrey Standard Errors
require('sandwich')

MatSqrtInverse <- function(A) {
    ##  Compute the inverse square root of a matrix
    ei <- eigen(A)
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

    df <- function(GG) {                # Compute DoF given G'*Omega*G
        lam <- eigen(GG, only.values=TRUE)$values
        sum(lam)^2/sum(lam^2)
    }

    if(is.null(clustervar)) {           # no clustering
        Vhat <- vcovHC(model, type="HC2")
        Vhat.Stata <- Vhat*NA

        M <- diag(n)-X %*% XXinv %*% t(X)       # annihilator matrix
        GOG <- function(ell) {           # G'*Omega*G
            Xtilde <- drop(X %*% XXinv %*% ell / sqrt(diag(M)))
            crossprod(M * Xtilde)
        }
    } else {
        if(!is.factor(clustervar)) stop("'clustervar' must be a factor")

        ## Stata
        S <- length(levels(clustervar)) # number clusters
        uj <- apply(u*X, 2, function(x) tapply(x, clustervar, sum))
        Vhat.Stata <- S/(S-1) * (n-1)/(n-K) * sandwich(model, meat = crossprod(uj)/n)

        ## HC2
        tXs <- function(s) {
            Xs <- X[clustervar==s, , drop=FALSE]
            MatSqrtInverse(diag(NROW(Xs))-Xs%*% XXinv %*% t(Xs)) %*% Xs
        }
        tX <- lapply(levels(clustervar), tXs) # list of matrices

        tu <- split(u, clustervar)
        tutX <- sapply(seq_along(tu),function(i) crossprod(tu[[i]],tX[[i]]))
        Vhat <- sandwich(model, meat = tcrossprod(tutX)/n)

        ## DOF adjustment
        tHs <- function(s) {
            Xs <- X[clustervar==s, , drop=FALSE]
            index <- which(clustervar==s)
            ss <- outer(rep(0,n),index)     # n x ns matrix of 0
            ss[cbind(index,1:length(index))] <- 1
            ss-X %*% XXinv %*% t(Xs)
        }
        tH <- lapply(levels(clustervar), tHs) # list of matrices

        Moulton <- function() {
            ## Moulton estimates
            ns <- tapply(u, clustervar,length)
            ssr <- sum(u^2)
            rho <- max((sum(sapply(seq_along(tu), function(i)
                sum(tu[[i]] %o% tu[[i]])))-ssr) / (sum(ns^2)-n), 0)
            c(sig.eps=max(ssr/n - rho, 0), rho=rho)
        }

        GOG <- function(ell) {
            G <- sapply(seq_along(tX),
                        function(i)  tH[[i]] %*% tX[[i]] %*% XXinv %*% ell)
            GG <- crossprod(G)

            if (IK==TRUE) {            # IK method
                Gsums <- apply(G, 2, function(x) tapply(x, clustervar, sum)) # Z'*G
                GG <- Moulton()[1]*GG + Moulton()[2]*crossprod(Gsums)
            }
            GG
        }
    }

    if (!is.null(ell)) {
        se <- drop(sqrt(crossprod(ell,Vhat) %*% ell))
        dof <- df(GOG(ell))
        se.Stata <- drop(sqrt(crossprod(ell,Vhat.Stata) %*% ell))
    } else {
        se <- sqrt(diag(Vhat))
        dof <- sapply(seq(K), function(k) df(GOG(diag(K)[,k])))
        se.Stata <- sqrt(diag(Vhat.Stata))
    }
    names(dof) <- names(se)

    return(list(vcov=Vhat, dof=dof, adj.se=se*qt(0.975,df=dof)/qnorm(0.975),
                se.Stata=se.Stata))
}
