## *****************************************************************************
##' @description Contours of a kriging function or statistic related
##'     to a kriging model object with two inputs. This can be the
##'     kriging mean, the kriging (conditional) standard deviation and
##'     more.
##'
##' @details When several functions are given, they will be displayed
##'     on different facets of the plot. These functions should then
##'     have the same dimension because \emph{only one colour scale}
##'     is used.  For instance it makes no sense to give the mean and
##'     the variance which do not have the same dimension, nor even
##'     the mean and the standard deviation because these two
##'     functions may have different order of magnitude. We can
##'     usefully want the true function and the kriging mean, or the
##'     kriging mean and the kriging trend.
##' 
##' @title Contours of a Function Related to a Kriging Model object
##'     with Two Inputs
##' 
##' @param object An object with class \code{"km"} or equivalent
##'     with two inputs.
##' 
##' @param which Character vector. Name of the kriging
##'     function(s)/statistic(s) to be shown. Can contain
##'     \code{"mean"}, \code{"trend"}, \code{"sd"}, \code{"sd2"} or
##'     \code{"var"}. Can have length 0 as well, in order to display
##'     the contours of a function having or not a formal argument
##'     \code{object}. This function is then provided by using
##'     \code{other}.
##' 
##' @param other Character. Name of a function. This can be a function
##'     to be compared to the stat(s) given in \code{which} (usually
##'     the kriging mean). If \code{which} has length zero we can also
##'     specify a gradient function. The function must be a function
##'     of one vector argument. It optionally can have one more
##'     argument with name \code{"model"} or \code{"object"} in which
##'     case \code{object} will be passed will be passed to this
##'     argument as is required to diplay a one-point Bayesian
##'     optimisation criterion (EI, AEI, ...), see \bold{Examples}.
##'
##' @param lower,upper Numeric vectors with length 2 definign the
##'     bounds of the rectangular region for the contours.
##'
##' @param nGrid Integer. Number of grid points. Can be of length 2 if
##'     different grids are to be used for the two dimensions.
##'
##' @param grad Logical. If \code{TRUE} the gradient of the kriging
##'     function/statistic will be shown using arrows. This will only
##'     be possible when \emph{only one} function is used among
##'     \code{"mean"} \code{"trend"} \code{"sd"} or \code{"sd2"}, or
##'     when a function is given in \code{other}. In the second case
##'     the gradient function must be given in \code{otherGrad}.
##'   
##' @param otherGrad Character. Name of the gradient function of the
##'     function given in \code{other}. This will be ignored if
##'     \code{other} is not provided or if \code{grad} is
##'     \code{FALSE}. \bold{Not implemented yet}.
##'
##' @param whereGrad Either a character with length one or a data
##'     frame or matrix with two columns. This defines the points at
##'     which the gradient will be evaluated and plotted as a small
##'     arrow in a vector-field fashion. The default is
##'     \code{"levels"} indicating that the gradient is evaluated
##'     along some level curves (not those used to define the filled
##'     contours). Another possible value is \code{"grid"} Finally if
##'     a matrix or data frame is provided, the gradient will only be
##'     evaluated and plotted at the specified points.
##' 
##' @param ... Not used yet.
##'
##' @return A graphical object inheriting from \code{"ggplot"}.
##'
##' @export
##' @import ggplot2
##' 
##' @examples
##' set.seed(1521)
##' n <- 9
##' design <- matrix(runif(n * 2), ncol = 2,
##'                  dimnames = list(NULL, c("x1", "x2")))
##' design <- data.frame(design)
##' y <- apply(design, 1, branin)
##' model <- km(~1, design = design, response = y)
##' 
##' ## Contours only
##' ## =============
##' contours(model) + ggtitle("kriging mean (default in 'contours')")
##' contours(model, other = "branin") + ggtitle("'branin' and kriging mean")
##' contours(model, which = c("mean", "trend")) + ggtitle("Kriging mean and trend")
##' contours(model, which = "sd") + ggtitle("Kriging sd")
##' 
##' ## Display gradients
##' ## =================
##' contours(model, which = "mean", grad = TRUE) +
##'      ggtitle("Kriging mean, gradients at contours")
##' contours(model, which = "sd", grad = TRUE) +
##'     ggtitle("Kriging sd, gradients at contours")
##' contours(model, which = "sd", grad = TRUE, whereGrad = "grid") +
##'     ggtitle("Kriging sd, gradients at grid points")
##'
##' ## more involved examples
##' ## ======================
##' \dontrun{
##' braninG <- function(x) branin_with_grad(x)$gradient
##' contours(model, which = character(0), grad = TRUE,
##'          other = "branin", otherGrad = "braninG", whereGrad = "grid") +
##'     ggtitle("Branin function and its gradient")
##' contours(model, which = character(0), grad = TRUE,
##'          other = "EI", otherGrad = "EI.grad", whereGrad = "grid") +
##'     ggtitle("Expected Improvement and its gradient")
##' }
contours <- function(object,
                     which = "mean",
                     other = NULL,
                     lower = c(0, 0),
                     upper = c(1, 1),
                     nGrid = 50,
                     grad = FALSE,
                     otherGrad = NULL,
                     whereGrad = "levels",
                     ...) {
    
    x1 <- x2 <- stat <- statDer1 <- statDer2 <- z <- NULL
    
    ## XXX TODO : manage input names, here hard-coded as c("x1", "x2")
    
    if (!requireNamespace("tidyr", quietly = TRUE)) {
        stop("This function requires the 'tidyr' package")
    }
    if (!requireNamespace("grDevices", quietly = TRUE)) {
        stop("This function requires the 'grDevices' package")
    }
    colors <- grDevices::colorRampPalette(c("seagreen","orange",
                                            "red","brown"))(20)
    if (object@d != 2) {
        stop("'contours' only works when object@d == 2")
    }

    Which <- c("mean", "sd", "var")
    i <- match(which, Which)
    
    nGrid <- rep(nGrid, length.out = 2)
    xGrid <- list(x1 = seq(from = lower[1], to = upper[1], length.out = nGrid[1]),
                  x2 = seq(from = lower[2], to = upper[2], length.out = nGrid[2]))
    df <- expand.grid(x1 = xGrid[[1]], x2 = xGrid[[2]])
    
    colnames(df) <- inputNames <- colnames(object@X)
    inputNames <- c("x1", "x2")

    mat <- matrix(as.numeric(unlist(df)), ncol = 2)
    pred <- predict(object, newdata = mat, deriv = FALSE, type = "UK")
    
    Which <- c("mean", "sd", "var")
    i <- match(which, Which)
    o <- !is.null(other)
    
    ## =========================================================================
    ## The gradient is evaluated at as smaller number of points
    ## =========================================================================
    
    if (grad) {
        
        if (length(which) + o > 1 ||
            (length(which) > 0 &&
             !(which %in% c("trend", "mean", "sd", "var")))) {
            stop("When 'grad' is TRUE only one function can be ",
                 "used along with its gradient")   
        }

        if (o && is.null(otherGrad)) {
            stop("When 'other' is given and 'grad' is TRUE, ",
                 "otherGrad must be provided")
        }

        if (o) {
            nmo <- as.character(other)
            other <- match.fun(other)
            if ("model" %in% names(formals(other))) {
                stat <- apply(df, MARGIN = 1, FUN = other, model = object)
            } else if ("object" %in% names(formals(other))) {
                stat <- apply(df, MARGIN = 1, FUN = other, object = object)
            } else {
                stat <- apply(df, MARGIN = 1, FUN = other)
            }
            df <- data.frame(df, stat = stat)
            
        } else {
            if (which == "var") which <- "s2"
            df <- data.frame(df, stat = pred[[which]])
        }

        ## =====================================================================
        ## Define the data frame where the gradient will be evaluated
        ## =====================================================================
        
        if (whereGrad == "levels") {
           
            cL <- grDevices::contourLines(x = xGrid$x1, y = xGrid$x2,
                                          z = matrix(df$stat, nrow = nGrid[1],
                                                     ncol = nGrid[2]),
                                          nlevels = 10)
            
            dfDer <- data.frame(x1 = unlist(sapply(cL, function(x) x$x)),
                                x2 = unlist(sapply(cL, function(x) x$y)))
            
        } else if (whereGrad == "grid") {
            hDer1 <- upper[1] - lower[1]
            hDer2 <- upper[2] - lower[2]
            xDer1 <- seq(lower[1] + hDer1 / 20, by = hDer1 / 10, to = upper[1])
            xDer2 <- seq(lower[2] + hDer2 / 20, by = hDer2 / 10, to = upper[2])
            dfDer <- expand.grid(x1 = xDer1, x2 = xDer2)
        } else {
            if (ncol(whereGrad) != 2) stop("'whereGrad' must have 2 columns")
            dfDer <- whereData
            colnames(whereData) <- c("x1", "x2")
        }
        
        if (o) {
            ## nmo <- as.character(other)
            ## other <- match.fun(other)
            ## otherGrad <- match.fun(otherGrad)
            ## stat  <- apply(df, MARGIN = 1, FUN = other)
            ## statDer <- t(apply(df, MARGIN = 1, FUN = otherGrad))
            ## dfDer <- data.frame(dfDer, stat = stat,
            ##                     statDer1 = statDer[ , 1],
            ##                     statDer2 = statDer[ , 2])
            if ("model" %in% names(formals(other))) {
                stat <- apply(dfDer, MARGIN = 1, FUN = other, model = object)
            } else if ("object" %in% names(formals(other))) {
                stat <- apply(dfDer, MARGIN = 1, FUN = other, object = object)
            } else {
                stat <- apply(dfDer, MARGIN = 1, FUN = other)
            }
            
            predDer <- data.frame(stat = stat)
            which <- "stat"
            
            if (!is.null(otherGrad)) {
                otherGrad <- match.fun(otherGrad)
                if ("model" %in% names(formals(other))) {
                    statDer <- t(apply(dfDer, MARGIN = 1, FUN = otherGrad, model = object))
                } else if ("object" %in% names(formals(other))) {
                    statDer <- t(apply(dfDer, MARGIN = 1, FUN = otherGrad, object = object))
                } else {
                    statDer <- t(apply(dfDer, MARGIN = 1, FUN = otherGrad))
                }
                colnames(statDer) <- c("x1", "x2")
            }
            
        } else {
            predDer <- predict(object, newdata = dfDer,
                               deriv = TRUE, type = "UK")    
            statDer <- predDer[[paste0(which, ".deriv")]]
            dDer <- dim(statDer)
            statDer <- statDer[slice.index(statDer, 1) == slice.index(statDer, 2)]
            statDer <- matrix(statDer, nrow = dDer[1])
        }
        
        ratio <- 30 * max(abs(statDer[ , 1]) / (upper[1] - lower[1]),
                          abs(statDer[ , 2]) / (upper[2] - lower[2])) 
            
        dfDer <- data.frame(dfDer, stat = predDer[[which]],
                            statDer1 = statDer[ , 1],
                            statDer2 = statDer[ , 2])
        v <- ggplot() +
            geom_contour_filled(data = df, 
                                mapping = aes(x = x1, y = x2, z = stat),
                                bins = 20, alpha = 0.7) +
            scale_fill_manual(values = colors,
                              name ="Value", drop = FALSE) +
                geom_point(data = as.data.frame(object@X),
                           mapping = aes(x = x1, y = x2),
                           size = 3) +
            geom_segment(data = dfDer,
                         aes(x = x1,
                             y = x2,
                             xend = x1 + statDer1 / ratio,
                             yend = x2 + statDer2 / ratio),
                         colour = "royalBlue",
                         arrow = arrow(length = unit(0.1, "cm")), size = 0.25)
        v
        return(v)
        
    }
    
    if ("mean" %in% which) df <- data.frame(df, mean = pred$mean)
    if ("trend" %in% which) df <- data.frame(df, trend = pred$trend)
    if ("sd" %in% which) df <- data.frame(df, sd = pred$sd)
    if ("var" %in% which || "sd2" %in% which) df <- data.frame(df, sd2 = pred$sd2)
    
    if (!is.null(other)) {
        nmo <- as.character(other)
        other <- match.fun(other)
        if ("model" %in% names(formals(other))) {
            df[[nmo]] <- apply(df[ , inputNames], MARGIN = 1, FUN = other, model = object)
        } else if ("object" %in% names(formals(other))) {
            df[[nmo]] <- apply(df[ , inputNames], MARGIN = 1, FUN = other, object = object)
        } else {
            df[[nmo]] <- apply(df[ , inputNames], MARGIN = 1, FUN = other)
        }
        ## df[[nmo]] <- apply(df, MARGIN = 1, FUN = other)
    }
    
    df2 <- tidyr::gather(df, key = "which", value = "z", -inputNames)
    
    colors <- grDevices::colorRampPalette(c("seagreen","orange",
                                            "red","brown"))(20)
    v <- ggplot() +
        geom_contour_filled(data = df2, 
                            mapping = aes(x = x1, y = x2, z = z),
                            bins = 20, alpha = 0.7) +
        ##geom_contour(data = df2, 
        ##             mapping = aes(x = x1, y = x2, z = z),
        ##             bins = 20) +
        scale_fill_manual(values = colors,
                          name ="Value", drop = FALSE) +
        geom_point(data = as.data.frame(object@X),
                   mapping = aes(x = x1, y = x2),
                   size = 3) +
        facet_grid(which ~ ., scales = "free")  
    v
}
               
