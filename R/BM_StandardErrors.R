tol <- 1e-9

## TODO: HC2 when diaghat = 1, X must be full rank, X has one column, X is just the intercept
## TODO: If each observation in its cluster, we get HC2

## Compute DoF given G'*Omega*G without calling eigen as suggested by
## Winston Lin. sum(GG * GG) = sum(diag(GG %*% GG))
df <- function(GG)
    sum(diag(GG))^2 / sum(GG * GG)

#' Standard Errors with adjusted degrees of freedom
#' @param model Fitted model returned by the \code{\link{lm}} function
#' @param clustervar Factor variable that defines clusters. If \code{NULL} (or
#'     not supplied), the command computes heteroscedasticity-robust standard
#'     errors, rather than cluster-robust standard errors.
#' @param ell A vector of the same length as the dimension of covariates,
#'     specifying which linear combination \eqn{\ell'\beta} of coefficients
#'     \eqn{\beta} to compute. If \code{NULL}, compute standard errors for each
#'     regressor coefficient
#' @param IK Logical flag only relevant if cluster-robust standard errors are
#'     being computed. Specifies whether to compute the degrees-of-freedom
#'     adjustment using the Imbens-Kolesár method (if \code{TRUE}), or the
#'     Bell-McCaffrey method (if \code{FALSE})
#' @return Returns a list with the following components \describe{
#'
#' \item{vcov}{Variance-covariance matrix estimator. For the case without
#' clustering, it corresponds to the HC2 estimator (see MacKinnon and White,
#' 1985 and the reference manual for the \code{sandwich} package). For the case
#' with clustering, it corresponds to a generalization of the HC2 estimator,
#' called LZ2 in Imbens and Kolesár.}
#'
#' \item{dof}{Degrees-of-freedom adjustment}
#'
#' \item{se}{Standard error}
#'
#' \item{adj.se}{Adjusted standard errors. For \beta_j, they are defined as
#' \code{adj.se[j]=sqrt(vcov[j,j]se*qt(0.975,df=dof)} so that the Bell-McCaffrey
#' confidence intervals are given as \code{coefficients(fm)[j] +- 1.96* adj.se=}}
#'
#' \item{se.Stata}{Square root of the cluster-robust variance estimator used in
#' STATA}
#'
#' }
#' @examples
#' ## No clustering:
#' set.seed(42)
#' x <- sin(1:10)
#' y <- rnorm(10)
#' fm <- lm(y~x)
#' dfadjustSE(fm)
#' ## Clustering, defining the first six observations to be in cluster 1, the
#' #next two in cluster 2, and the last three in cluster three.
#' clustervar <- as.factor(c(rep(1, 6), rep(2, 2), rep(3, 2)))
#' dfadjustSE(fm, clustervar)
#' @export
dfadjustSE <- function(model, clustervar=NULL, ell=NULL, IK=TRUE) {
    Q <- qr.Q(model$qr)
    R <- qr.R(model$qr)
    n <- NROW(Q)
    u <- residuals(model)
    K <- model$rank

    sandwich <- function(meat) backsolve(R, t(backsolve(R, meat)))

    ## no clustering
    if(is.null(clustervar)) {
        ## Compute meat of HC1 and HC2
        diaghat <- try(hatvalues(model), silent = TRUE)
        AQ <- (1-diaghat >= tol) * (1/sqrt(pmax(1-diaghat, tol))) * Q
        HC2 <- crossprod(u * AQ)
        HC1 <- n/(n-K)*crossprod(u*Q)

        ## G'*G
        df0 <- function (ell) {
            lt <- backsolve(R, ell, transpose=TRUE)
            a <-  drop(AQ %*% lt)
            B <- a*Q
            num <- (sum(a^2)-sum(B^2))^2
            den <- sum(a^4)-2*sum((a*B)^2)+sum(crossprod(B)^2)
            num/den
        }
        ## GOG <- function(ell) {
        ##     lt <- backsolve(R, ell, transpose=TRUE)
        ##     a <-  drop(AQ %*% lt)
        ##     diag(a^2)-tcrossprod(a * Q)
        ## }
    } else {
        if(!is.factor(clustervar)) stop("'clustervar' must be a factor")

        sum.model <- summary.lm(model)
        XXinv <- sum.model$cov.unscaled # XX^{-1}

        ## Compute meat of HC1 and HC2
        S <- length(levels(clustervar)) # number of clusters
        uj <- apply(u*Q, 2, function(x) tapply(x, clustervar, sum))
        HC1 <- S/(S-1) * (n-1)/(n-K) * crossprod(uj)
        as <- function(s) {
            Qs <- Q[clustervar==s, , drop=FALSE]
            e <- eigen(crossprod(Qs))
            Ds <- e$vectors %*% ((1-e$values >= tol) *
                                 (1/sqrt(pmax(1-e$values, tol))) * t(e$vectors))
            crossprod(crossprod(u[clustervar==s],  Qs %*% Ds))
        }
        HC2 <- Reduce("+", lapply(levels(clustervar), as)) # list of matrices

        ## ## DOF adjustment
        ## tHs <- function(s) {
        ##     Xs <- X[clustervar==s, , drop=FALSE]
        ##     index <- which(clustervar==s)
        ##     ss <- outer(rep(0, n), index)     # n x ns matrix of 0
        ##     ss[cbind(index, 1:length(index))] <- 1
        ##     ss-X %*% XXinv %*% t(Xs)
        ## }
        ## tH <- lapply(levels(clustervar), tHs) # list of matrices

        ## Moulton <- function() {
        ##     ## Moulton estimates
        ##     ns <- tapply(u, clustervar, length)
        ##     ssr <- sum(u^2)
        ##     rho <- max((sum(sapply(seq_along(tu), function(i)
        ##         sum(tu[[i]] %o% tu[[i]])))-ssr) / (sum(ns^2)-n), 0)
        ##     c(sig.eps=max(ssr/n - rho, 0), rho=rho)
        ## }

        ## GOG <- function(ell) {
        ##     G <- sapply(seq_along(tX),
        ##                 function(i)  tH[[i]] %*% tX[[i]] %*% XXinv %*% ell)
        ##     GG <- crossprod(G)

        ##     ## IK method
        ##     if (IK==TRUE) {
        ##         Gsums <- apply(G, 2,
        ##                        function(x) tapply(x, clustervar, sum)) # Z'*G
        ##         GG <- Moulton()[1]*GG + Moulton()[2]*crossprod(Gsums)
        ##     }
        ##     GG
        ## }
        df0 <- function(ell) Inf
    }
    Vhat <- sandwich(HC2)
    Vhat.Stata <- sandwich(HC1)


    if (!is.null(ell)) {
        se <- drop(sqrt(crossprod(ell, Vhat) %*% ell))
        dof <- df0(ell)
        se.Stata <- drop(sqrt(crossprod(ell, Vhat.Stata) %*% ell))
    } else {
        se <- sqrt(diag(Vhat))
        dof <- sapply(seq(K), function(k) df0(diag(K)[, k]))
        se.Stata <- sqrt(diag(Vhat.Stata))
    }

    names(dof) <- names(se) <- names(se.Stata) <- colnames(Vhat) <-
        rownames(Vhat) <- names(model$coefficients)

    list(vcov=Vhat, dof=dof, adj.se=se*qt(0.975, df=dof)/qnorm(0.975),
                se=se, se.Stata=se.Stata)
}
