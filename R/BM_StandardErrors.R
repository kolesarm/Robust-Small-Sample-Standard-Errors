## TODO: HC2 when diaghat = 1, X must be full rank, X has one column, X is just
## the intercept TODO: If each observation in its cluster, we get HC2, TODO:
## negative rho TODO: example with few treated, with few treated clusters, and
## with fixed effects. Check partialling out doesn't affect anything.

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
#' @param tol Numerical tolerance for determining whether an eigenvalue equals zero.
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
#' \item{adj.se}{Adjusted standard errors. For \eqn{\beta_j}, they are defined as
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
dfadjustSE <- function(model, clustervar=NULL, ell=NULL, IK=TRUE, tol=1e-9) {
    Q <- qr.Q(model$qr)
    R <- qr.R(model$qr)
    n <- NROW(Q)
    u <- stats::residuals(model)
    K <- model$rank
    ## Moulton estimates
    rho <- sig <- NA

    sandwich <- function(meat) backsolve(R, t(backsolve(R, meat)))

    ## no clustering
    if(is.null(clustervar)) {
        ## Compute meat of HC1 and HC2
        diaghat <- try(stats::hatvalues(model), silent = TRUE)
        AQ <- (1-diaghat >= tol) * (1/sqrt(pmax(1-diaghat, tol))) * Q
        HC2 <- crossprod(u * AQ)
        HC1 <- n/(n-K)*crossprod(u*Q)

        ## G'*G
        df0 <- function (ell) {
            a <-  drop(AQ %*% backsolve(R, ell, transpose=TRUE))
            B <- a*Q
             (sum(a^2)-sum(B^2))^2 /
                 (sum(a^4)-2*sum((a*B)^2)+sum(crossprod(B)^2))
        }
    } else {
        if(!is.factor(clustervar)) stop("'clustervar' must be a factor")

        ## Compute meat of HC1 and HC2
        S <- length(levels(clustervar)) # number of clusters
        uj <- apply(u*Q, 2, function(x) tapply(x, clustervar, sum))
        HC1 <- S/(S-1) * (n-1)/(n-K) * crossprod(uj)
        AQf <- function(s) {
            Qs <- Q[clustervar==s, , drop=FALSE]
            e <- eigen(crossprod(Qs))
            Ds <- e$vectors %*% ((1-e$values >= tol) *
                                 (1/sqrt(pmax(1-e$values, tol))) * t(e$vectors))
            Qs %*% Ds
        }
        AQ <- lapply(levels(clustervar), AQf) # list of matrices
        AQ <- do.call(rbind, AQ)
        uj <- apply(u*AQ, 2, function(x) tapply(x, clustervar, sum))
        HC2 <- crossprod(uj)

        if (IK) {
            ## Moulton estimates
            ssr <- sum(u^2)
            rho <- (sum(tapply(u, clustervar, sum)^2)-ssr) /
                (sum(tapply(u, clustervar, length)^2)-n)
            ## Don't allow for negative correlation
            rho <- max(rho, 0)
            sig <- max(ssr/n - rho, 0)
        }

        df0 <- function(ell) {
            a <-  drop(AQ %*% backsolve(R, ell, transpose=TRUE))
            as <- as.vector(tapply(a^2, clustervar, sum))
            B  <- apply(a*Q, 2, function(x) tapply(x, clustervar, sum))
            if (!IK) {
                (sum(as)-sum(B^2))^2 /
                    (sum(as^2)-2*sum(as*B^2)+sum(crossprod(B)^2))
            } else {
                D <- as.vector(tapply(a, clustervar, sum))
                F  <- apply(Q, 2, function(x) tapply(x, clustervar, sum))
                GG <- sig*(diag(as)-tcrossprod(B)) +
                    rho*tcrossprod(diag(D)-tcrossprod(B, F))
                sum(diag(GG))^2 / sum(GG^2)
            }
        }
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

    list(vcov=Vhat, dof=dof,
         adj.se=se*stats::qt(0.975, df=dof)/stats::qnorm(0.975),
         se=se, se.Stata=se.Stata, rho=rho, sig=sig)
}
