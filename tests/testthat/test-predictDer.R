## *****************************************************************************
## AUTHOR: Yves Deville <deville.yves@alpestat.com>
##
## Check the derivatives (Jacobian arrays) returned by the 'predict' method
## for the class "km" as refactored in dolka.
##
## CAUTION. The 'predict.km' function of dolka may be no longer exported in
## the future. Switch then to a simple 'predict'.
##
## The check is for the following stats trend, mean, variance, and covariance.
##
## *****************************************************************************

context("Predict with derivatives")

library(numDeriv)
library(dolka)
library(testthat)

set.seed(12345)
n <- 16
d <- 2
nNew <- 5
trueResp <- match.fun("branin")
types <- c("UK", "SK")
stats <- c("trend", "mean", "sd2", "cov")
dims <- list("trend" = c(nNew, nNew,  d),
             "mean" = c(nNew, nNew,  d),
             "sd2" = c(nNew, nNew,  d),
             "cov" = c(nNew, nNew, nNew,  d))

prec <- c("trend" = 1e-5, "mean" = 1e-5, "sd2" = 1e-5, "cov" = 1e-3)

krigeStat <- function(x, object, which = stats, type = "UK", ...) {
    which <- match.arg(which)
    XNew <- matrix(x, ncol = object@d,
                   dimnames = list(NULL, colnames(object@X)))
    p <- predict(object = object,
                 which = which,
                 cov = TRUE,
                 newdata = XNew, type = type)[[which]]
    p
}

## =============================================================================
## Generate design and response
## =============================================================================

X <- array(runif(n * d),
           dim = c(n, d),
           dimnames = list(NULL, paste0("x", 1:d)))
y <- apply(X, 1, trueResp)

myKm <- km(~.,
           design = data.frame(X),
           response = y,
           covtype = "gauss")

## =============================================================================
## New design for tests
## =============================================================================

XNew <- matrix(runif(nNew * d), nrow = nNew)
colnames(XNew) <- colnames(myKm@X)

## =============================================================================
## Compute the Jacobian arrays with 'numDeriv'
## =============================================================================
for (type in types) {
    predNew <- predict(myKm, newdata = XNew, type = type,
                       cov = TRUE, deriv = TRUE)
    for (stat in stats) {
        J <- jacobian(krigeStat, x = XNew,
                      which = stat,
                      object = myKm, type = type)
        dim(J) <- dims[[stat]]
        errDer <- predNew[[paste0(stat, ".deriv")]] - J
        test_that(desc = sprintf(paste0("Jacobian array of prediction",
                                        " stat =\"%s\", type = \"%s\""),
                                 stat, type),
                  code = expect_true(max(abs(errDer)) < prec[[stat]]))
    }
}

if (FALSE) {

    ## =========================================================================
    ## CAUTION. 'DiceOptim' uses a quite unusual rule for Jacobians. E.g., for
    ## the conditional covariance 'Sigma'...
    ##
    ## o DiceOptim
    ##
    ##          [i, j, k, ell] = (d Sigma[k, ell] / (d X[i, j])
    ##
    ## with dimension c(nNew, d, nNew, nNew)
    ##
    ## o Gradient: output first!!!
    ##
    ##          [i, j, k, ell] = (d Sigma[i, j] / (d X[k, ell])
    ##
    ## with dimension c(nNew, nNew, nNew, d)
    ## =========================================================================

    require(DiceOptim)
    DOKD <- DiceOptim:::krigingDeriv(x = XNew, model = myKm, type = "UK")
    DOJacCov <- KPKD[[4]]
    dimnames(DOJacCov) <- list(NULL, colnames(myKm@X), NULL, NULL)
    max(abs(aperm(DOJacCov, perm = c(3, 4, 1, 2)) - JacCov))

}
