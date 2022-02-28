## *****************************************************************************

##' Return the diagonal of the cross-product
##' \eqn{\mathbf{X}^\top\mathbf{X}}{t(X) \%*\% X} for a given matrix
##' \eqn{\mathbf{X}}{X}. As is well-known the cross-product is
##' computed by the single argument of \code{crossprod}; however we
##' often have to compute \emph{only the diagonal} of the
##' cross-product. For instance, the variance of a kriging predition
##' may have to be computed without formin the covariance.
##'
##' Several methods can be used to compute the result in R (using
##' \code{apply}, ...). To our best knowledge the method used here is
##' the most efficient in general framework.
##' 
##' @title Diagonal of the Cross-Product of a Matrix
##' 
##' @param X A numeric matrix.
##'
##' @return The diagonal vector of \code{crossprod(X)}.
##'
##' @author Thanks to Clément Walter for comparisons and the suggested
##'     code.
##' 
dcrossprod <- function(X) {

    if (is.null(dX <- dim(X))) return(crossprod(X))
    
    drop(rep(1, dX[1]) %*% (X^2))
    
}
