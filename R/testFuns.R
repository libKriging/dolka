## *****************************************************************************
##' @description Branin function with its exact gradient returned into
##'     a list.
##'
##' @title Branin Function with Gradient
##' 
##' @param x A vector with length \eqn{d = 2}.
##' 
##' @return A list with two numeric elements named \code{objective}
##'     and \code{gradient} corresponding to the function value
##'     (length 1) and its gradient (length 2).
##'
##' @seealso \code{\link[DiceKriging]{branin}}.
##'
##' @export
##' 
##' @examples
##' X <- matrix(runif(100), ncol = 2, dimnames = list(NULL, c("x1", "x2")))
##' Obj <- apply(X, 1, function(x) braninGrad(x)$objective)
##' ObjDK <- apply(X, 1, DiceKriging::branin)
##' max(abs(Obj - ObjDK))
##' Grad <- t(apply(X, 1, function(x) braninGrad(x)$gradient))
##' GradNum <- t(apply(X, 1, function(x) grad(DiceKriging::branin, x)))
##' max(abs(Grad - GradNum))
##' 
braninGrad <- function(x) {
    
    x1 <- x[1] * 15 - 5   
    x2 <- x[2] * 15
    b <- 5 / (4 * pi * pi)
    c <- 5 / pi
    d <- 10 * (1 - 1 / 8 / pi)
    A <- x2 - b * x1 * x1  + c * x1 - 6
    obj <- A * A + d * cos(x1) + 10
    grad <- c("x1" = 15 * (2 * A * (-2 * b * x1 + c) - d * sin(x1)),
              "x2" = 15 * 2 * A )
    list(objective = obj,
         gradient = grad)
    
}
