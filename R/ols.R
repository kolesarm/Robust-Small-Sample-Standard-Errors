#' Standard Errors with adjusted degrees of freedom
#' @param model Fitted model returned by the \code{lm} function
#' @param clustervar Factor variable that defines clusters. If \code{NULL} (or
#'     not supplied), the command computes heteroscedasticity-robust standard
#'     errors, rather than cluster-robust standard errors.
#' @param ell A vector of the same length as the dimension of covariates,
#'     specifying which linear combination \eqn{\ell'\beta}{ell'beta} of
#'     coefficients \eqn{\beta}{beta} to compute. If \code{NULL}, compute
#'     standard errors for each regressor coefficient.
#' @param IK Only relevant for cluster-robust standard errors. Specifies whether
#'     to compute the degrees-of-freedom adjustment using the Imbens-Kolesár
#'     (2016) method (if \code{TRUE}), or the Bell-McCaffrey (2002) method (if
#'     \code{FALSE}).
#' @param tol Numerical tolerance for determining whether an eigenvalue equals
#'     zero.
#' @param rho0 Impose positive \eqn{\rho}{rho} when estimating the Moulton
#'     (1986) model when implementing the \code{IK} method?
#' @return Returns a list with the following components \describe{
#'
#' \item{vcov}{Variance-covariance matrix estimator. For independent errors, it
#'  corresponds to the HC2 estimator (see MacKinnon and White, 1985, or the
#'  reference manual for the \code{sandwich} package). For clustered errors, it
#'  corresponds to a version the generalization of the HC2 estimator, called LZ2
#'  in Imbens and Kolesár.}
#'
#' \item{coefficients}{Matrix of estimated coefficients, along with HC1, and HC2
#' standard errors, Adjusted standard errors, and effective degrees of freedom.
#' Adjusted standard error is HC2 standard error multiplied by \code{qt(0.975,
#' df=dof)/qnorm(0.975)} so that one can construct 95% confidence intervals by
#' adding and subtracting 1.96 times the adjusted standard error.}
#'
#' \item{rho, sig}{Estimates of \eqn{\rho} and \eqn{\sigma} of the Moulton
#' (1986) model for the regression errors. Only computed if \code{IK} method is
#' used}
#'
#' }
#' @references{
#'
#' \cite{Robert M. Bell and Daniel F. McCaffrey. Bias reduction in standard
#' errors for linear regression with multi-stage samples. Survey Methodology,
#' 28(2):169–181, December 2002.}
#'
#' \cite{Guido W. Imbens and Michal Kolesár. Robust standard errors in small
#' samples: Some practical advice. Review of Economics and Statistics,
#' 98(4):701–712, October 2016. \doi{10.1162/REST_a_00552}}
#'
#' \cite{James G. MacKinnon and Halbert White. Some
#' Heteroskedasticity-Consistent Covariance Matrix Estimators with Improved
#' Finite Sample Properties. Journal of Econometrics, (29)3:305–325, September
#' 1985. \doi{10.1016/0304-4076(85)90158-7}}
#'
#' \cite{Brent R. Moulton. Random group effects and the precision of regression
#' estimates. Journal of Econometrics, 32(3):385–397, August 1986.
#' \doi{10.1016/0304-4076(86)90021-7}.}
#'
#' }
#' @examples
#' ## No clustering:
#' set.seed(42)
#' x <- sin(1:100)
#' y <- rnorm(100)
#' fm <- lm(y ~ x + I(x^2))
#' dfadjustSE(fm)
#' ## Clustering, with 5 clusters
#' clustervar <- as.factor(c(rep(1, 40), rep(1, 20),
#'                         rep(2, 20), rep(3, 10), rep(4, 10)))
#' dfadjustSE(fm, clustervar)
#' @export
dfadjustSE <- function(model, clustervar=NULL, ell=NULL, IK=TRUE, tol=1e-9,
                       rho0=FALSE) {
    Q <- qr.Q(model$qr)
    R <- qr.R(model$qr)
    n <- NROW(Q)
    u <- stats::residuals(model)
    K <- model$rank
    ## Moulton estimates
    rho <- sig <- NA

    sandwich <- function(meat) backsolve(R, t(backsolve(R, meat)))

    ## no clustering
    if (is.null(clustervar)) {
        ## Compute meat of HC1 and HC2. Hatvalues scale invariant
        diaghat <- try(stats::hatvalues(model), silent = TRUE)
        AQ <- (1-diaghat >= tol) * (1/sqrt(pmax(1-diaghat, tol))) * Q
        HC2 <- crossprod(u * AQ)
        HC1 <- n / (n-K) * crossprod(u*Q)

        ## G'*G
        df0 <- function(ell) {
            a <-  drop(AQ %*% backsolve(R, ell, transpose=TRUE))
            B <- a*Q
            (sum(a^2)-sum(B^2))^2 /
                (sum(a^4)-2*sum((a*B)^2)+sum(crossprod(B)^2))
        }
    } else {
        if (!is.factor(clustervar)) stop("'clustervar' must be a factor")
        ## Order u, Q by clusters
        ord <- order(clustervar)
        u <- u[ord]
        Q <- Q[ord, ]
        clustervar <- clustervar[ord]

        ## Compute meat of HC1 and HC2
        S <- nlevels(clustervar) # number of clusters
        uj <- apply(u*Q, 2, function(x) tapply(x, clustervar, sum))
        HC1 <- S / (S-1) * (n-1) / (n-K) * crossprod(uj)
        ## A_s * Q_s; Q matrix scale invariant
        AQf <- function(s) {
            Qs <- Q[clustervar==s, , drop=FALSE] # nolint
            e <- eigen(crossprod(Qs))
            Ds <-  (1-e$values >= tol) *
                (1/sqrt(pmax(1-e$values, tol))) * t(e$vectors)
            Qs %*% e$vectors %*% Ds
        }

        AQ <- lapply(levels(clustervar), AQf) # list of matrices
        AQ <- do.call(rbind, AQ)
        uj <- apply(u*AQ, 2, function(x) tapply(x, clustervar, sum))
        HC2 <- crossprod(uj)

        if (IK) {
            ## Moulton estimates, set rho=0 if no clustering
            ssr <- sum(u^2)
            den <- sum(tapply(u, clustervar, length)^2)-n
            rho <- 0
            if (den>0)
                rho <- (sum(tapply(u, clustervar, sum)^2)-ssr) / den

            ## Don't allow for negative correlation
            if (rho0) rho <- max(rho, 0)
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
                Fm  <- apply(Q, 2, function(x) tapply(x, clustervar, sum))
                GG <- sig * (diag(as)-tcrossprod(B)) +
                    rho*tcrossprod(diag(D)-tcrossprod(B, Fm))
                sum(diag(GG))^2 / sum(GG^2)
            }
        }
    }

    Vhat <- sandwich(HC2)
    VhatStata <- sandwich(HC1)

    if (!is.null(ell)) {
        se <- drop(sqrt(crossprod(ell, Vhat) %*% ell))
        dof <- df0(ell)
        seStata <- drop(sqrt(crossprod(ell, VhatStata) %*% ell))
        beta <- sum(ell*model$coefficients)
    } else {
        se <- sqrt(diag(Vhat))
        dof <- vapply(seq(K), function(k) df0(diag(K)[, k]), numeric(1))
        seStata <- sqrt(diag(VhatStata))
        beta <-  model$coefficients
    }

    r <- cbind("Estimate"=beta,
               "HC1 se"=seStata,
               "HC2 se"=se,
               "Adj. se"=se*stats::qt(0.975, df=dof)/stats::qnorm(0.975),
               "df"=dof)
    rownames(r) <- names(beta)
    colnames(Vhat) <- rownames(Vhat) <- names(model$coefficients)

    structure(list(vcov=Vhat, coefficients=r,
                   rho=rho, sig=sig), class="dfadjustSE")
}


#' @export
print.dfadjustSE <- function(x, digits = getOption("digits"), ...) {
    xx <- x$coefficients
    pv <- 2*stats::pt(-abs(xx[, "Estimate"] / xx[, "HC2 se"]), df=xx[, "df"])
    cat("\nCoefficients:\n")
    print(cbind(xx, "p-value"=pv), digits=digits)

    invisible(x)
}
