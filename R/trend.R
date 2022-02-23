##' Evaluate the trend component of a \code{km} object on the
##' given design \code{X}.
##'
##' @title Trend part of \code{km} object
##'
##' @param object A \code{km} object.
##'
##' @param X A "design": a numeric matrix or data frame wih numeric
##'     columns compatible with \code{object}.
##' 
##' @param deriv Logical. If \code{TRUE} the derivative will be
##'     returned as the \code{"deriv"} attribute of the result.  This
##'     will be an array with dimension either \code{c(n, p, n, d)} or
##'     \code{c(n, p, d)} if \code{diag = FALSE}.
##'    
##' @param diagOnly Logical. If \code{TRUE} the structure of the
##'     result takes into account the fact that the row \eqn{i} of the
##'     trend matrix \eqn{\mathbf{F}}{F} depends only on the row
##'     \enw{i} of the design matrix \eqn{\mathbf{X}}{X}. Then the
##'     deivative is a 3-dimensional array with its element \eqn{i, j,
##'     k} equal to \eqn{\partial F_{i,j}/ \partial X_{i,k}}{dF[i,
##'     j]/dX[i, k]}.
##' 
##' @return
##' @export
##' @examples
##' library(DiceKriging)
##' # a 16-points factorial design, and the corresponding response
##' d <- 2; n <- 16
##' design.fact <- expand.grid(x1 = seq(0, 1, length = 4),
##'                            x2 = seq(0, 1, length = 4))
##' y <- apply(design.fact, 1, branin)
##' # kriging model 1 : matern5_2 covariance structure, no trend, no
##' # nugget effect
##' fit <- km(~ x1 + x2, design = design.fact, response = y)
##' X <- matrix(runif(n = 40), ncol = 2,
##'             dimnames = list(NULL, c("x1", "x2")))
##' fitTrend <- trend(fit, X = X, deriv = 1)
##' fitTrend <- trend(fit, X = X[1, ], deriv = 1)
trend <- function(object, X, deriv = 0, diagOnly = FALSE) {

    d <- object@d
    p <- object@p
    
    if (!is.numeric(X)) {
        stop("'X' must be numeric")
    }
    if (is.null(dX <- dim(X))) {
        if (length(X) != d) {
            stop("'X' must be of length 'object@d'")
        }
        X <- matrix(X, ncol = d)
        n <- 1
    } else {
        if (dX[2] != 2) {
            stop("when 'X' is a matrix, it must have ",
                 "'object@d' columns")
        }
        n <- dX[1]
    }
    
    if (p == 0) {
        FX <- matrix(0, nrow = n, ncol = 0)
        if (deriv) attr(FX, "deriv") <- F
        return(FX)
    }
    
    trendFun <- function(x, asVector = FALSE) {
        dim(x) <- c(n, d)
        colnames(x) <- colnames(object@X)
        FX <- model.matrix(object@trend.formula,
                         data = as.data.frame(x))
        if (asVector) FX <- as.vector(FX)
        FX
    }

    FX <- trendFun(X)
    if (deriv) {
        if (diagOnly) {
            deriv <- array(0.0, dim = c(n, p, d),
                           dimnames = list(NULL, colnames(FX),
                                           colnames(object@X)))
            for (i in 1:n) {
                res <- numDeriv::jacobian(trendFun,
                                          x = as.vector(X[i, ]),
                                          asVector = TRUE)
                deriv[i, , ] <- attr(res, "deriv")
            }
            
        } else {
            deriv <- numDeriv::jacobian(trendFun, x = as.vector(X),
                                        asVector = TRUE)
            dim(deriv) <- c(n, p, n, d)
            dimnames(deriv) <- list(NULL, colnames(FX), NULL,
                                    colnames(object@X))
        }
        attr(FX, "deriv") <- deriv
    }
    FX
}
