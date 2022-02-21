trend <- function(object, X, deriv = 0) {

    if (!is.numeric(X)) {
        stop("'X' must be numeric")
    }
    if (is.null(dim(X))) {
        if (length(X) != object@d) {
            stop("'X' must be of length 'object@d'")
        }
    }
    X <- matrix(X, ncol = d)
    nNew <- nrow(X)
    colnames(X) <- colnames(object@X)
    rownames(x) <- NULL
    
    F <- function(x) {
        dim(x) <- c(nNew, d)
        model.matrix(model@trend.formula,
                     data = as.data.frame(x))
    }
    
    numDeriv::jacobian(F, x = as.vector(x))
    


}
