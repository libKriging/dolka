##' @title Analytical expression of the Expected Improvement criterion
##' 
##' @description Computes the Expected Improvement at current
##'     location. The current minimum of the observations can be
##'     replaced by an arbitrary value (plugin), which is useful in
##'     particular in noisy frameworks.
##' 
##' 
##' @param x A vector representing the input for which one wishes to
##'     calculate EI,
##'
##' @param model An object of class \code{\link[DiceKriging]{km}}.
##'
##' @param plugin Optional scalar: if provided, it replaces the
##'     minimum of the current observations.
##'
##' @param type \code{"UK"} (default) or \code{"SK"}, depending
##'     whether uncertainty related to trend estimation has to be
##'     taken into account.
##'
##' @param minimization Logical specifying if EI is used in
##'     minimiziation or in maximization.
##'
##' @param envir An optional environment specifying where to assign
##'     intermediate values for future gradient calculations. Default
##'     is \code{NULL}.
##'
##' @param proxy Optional logical. If \code{TRUE} EI is replaced by
##'     the kriging mean (to minimize).
##' 
##' @return The expected improvement, defined as \deqn{EI(x) := E[(
##'     min(Y(X)) - Y(x))^{+} | Y(X)=y(X)],} where X is the current
##'     design of experiments and Y is the random process assumed to
##'     have generated the objective function y.  If a plugin is
##'     specified, it replaces \deqn{min(Y(X))} in the previous
##'     formula.
##'
##' @export EI
##'
##' @author David Ginsbourger, Olivier Roustant and Victor Picheny.
##'
##' @seealso \code{\link{max_EI}}, \code{\link{EGO.nsteps}}, \code{\link{qEI}}
##'
##' @references
##' 
##' D.R. Jones, M. Schonlau, and W.J. Welch (1998), Efficient global
##' optimization of expensive black-box functions, \emph{Journal of Global
##' Optimization}, 13, 455-492.
##' 
##' J. Mockus (1988), \emph{Bayesian Approach to Global
##' Optimization}. Kluwer academic publishers.
##' 
##' T.J. Santner, B.J. Williams, and W.J. Notz (2003), \emph{The
##' design and analysis of computer experiments}, Springer.
##' 
##' M. Schonlau (1997), \emph{Computer experiments and global
##' optimization}, Ph.D. thesis, University of Waterloo.
##' 
##' @keywords models
##' 
##' @examples
##' set.seed(123)
##' ## =========================================================================
##' ## 	EI Surface Associated with an Ordinary Kriging Model for the Branin    
##' ##  Function Known at a 9-Points Factorial Design  
##' ## ========================================================================
##' 
##' ## a 9-points factorial design, and the corresponding response
##' ## ===========================================================
##' d <- 2; n <- 9
##' design.fact <-
##'     expand.grid(x1 = seq(0, 1, length = 3), x2 = seq(0, 1, length = 3))
##' y <- apply(design.fact, 1, branin) 
##' 
##' ## model fit
##' ## =========
##' fit1 <- km(~1, design = design.fact, response = y, covtype = "gauss",
##'            control = list(pop.size = 50, trace = FALSE), parinit = c(0.5, 0.5))
##' 
##' ## graphics
##' ## =========
##' n.grid <- 12
##' x.grid <- y.grid <- seq(0, 1, length = n.grid)
##' design.grid <- expand.grid(x.grid, y.grid)
##' #response.grid <- apply(design.grid, 1, branin)
##' EI.grid <- apply(design.grid, MARGIN = 1, FUN = EI, model = fit1)
##' z.grid <- matrix(EI.grid, nrow = n.grid, ncol = n.grid)
##' contour(x.grid, y.grid, z.grid, nlevels = 25)
##' title("Expected Improvement for the Branin function known at 9 points")
##' points(design.fact[ , 1], design.fact[ , 2], pch = 17, col = "blue")
##' 
EI <- function (x, model, plugin = NULL, type = c("UK", "SK"),
                minimization = TRUE, envir = NULL,
                proxy = FALSE) {

    type <- match.arg(type)
    
    if (is.null(plugin)){ 
        if (minimization) {
            plugin <- min(model@y)
        } else {
            plugin <- -max(model@y)
        }
    }
    m <- plugin
    
    ## =========================================================================
    ## Convert x in proper format(s)
    ## =========================================================================
    
    if (is.data.frame(x)) {
        d <- length(x)
        if (d != model@d) stop("'x' does not have the right size") 
        newdata.num <- as.numeric(x)
        newdata <- data.frame(t(newdata.num))
    } else {
        if (is.null(dim(x))) {
            d <- length(x)
            if (d != model@d){ stop("x does not have the right size") }
            newdata <- data.frame(t(as.numeric(x)))
        } else {
            d <- ncol(x)
            if (d != model@d){ stop("x does not have the right size") }
            newdata <- data.frame(x)
        }
    }
    colnames(newdata) <- colnames(model@X)
    
    ## =========================================================================
    ## The derivatives are computed only if they are to be stored in
    ## the environment. Note that since 'newdata' has only one row
    ## computing the covariance or not does not really matter.
    ## =========================================================================
    
    if (!is.null(envir)) {
        pred <- predict(object = model, newdata = newdata, type = type,
                        checkNames = FALSE, deriv = TRUE)
    } else {
        pred <- predict(object = model, newdata = newdata, type = type,
                        checkNames = FALSE, light.return = TRUE)
    }
  
    kriging.mean <- pred$mean
    if(!minimization) {
        kriging.mean <- -kriging.mean
    }
    kriging.sd   <- pred$sd
    
    if (proxy) {
        xcr <- xcr.prob <- xcr.dens <- NULL
        res <- m - kriging.mean
    } else {
        xcr <- (m - kriging.mean) / kriging.sd
        xcr.prob <- pnorm(xcr)
        xcr.dens <- dnorm(xcr)	        
        res <- (m - kriging.mean) * xcr.prob + kriging.sd * xcr.dens
    }
    
    ## =========================================================================
    ## Store the prediction and the auxiliary variables required (not
    ## very expansive to compute however).
    ## =========================================================================
    
    if (!is.null(envir))  {
        assign("pred", pred, envir = envir)
        assign("xcr", xcr, envir = envir)
        assign("xcr.prob", xcr.prob, envir = envir)
        assign("xcr.dens", xcr.dens, envir = envir)
    }
    
    return(res)
    
}


## =============================================================================
##'
##' @title Analytical gradient of the Expected Improvement criterion
##' 
##' @description Computes the gradient of the Expected Improvement at
##'     the current location.  The current minimum of the observations
##'     can be replaced by an arbitrary value (plugin), which is
##'     useful in particular in noisy frameworks.
##' 
##' @param x A vector representing the input for which one wishes to
##'     calculate \code{\link{EI}}.
##'
##' @param model An object of class \code{\link[DiceKriging]{km}}.
##'
##' @param plugin Optional scalar: if provided, it replaces the
##'     minimum of the current observations.
##'
##' @param type Character \code{"UK"} (default) or \code{"SK"},
##'     depending whether uncertainty related to trend estimation has
##'     to be taken into account.
##' 
##' @param minimization Logical specifying if EI is used in
##'     minimiziation or in maximization.
##'
##' @param envir Optional environment specifying where to get
##'     intermediate values calculated in \code{\link{EI}}.
##'
##' @param proxy Optional logical. If \code{TRUE}, EI is replaced by the
##'     kriging mean (to minimize).
##'
##' @return The gradient of the expected improvement criterion with
##'     respect to \eqn{x}.  Returns 0 at design points (where the
##'     gradient does not exist).
##'
##' @author David Ginsbourger, Olivier Roustant and Victor Picheny.
##' 
##' @seealso \code{\link{EI}}
##'
##' @section Caution: XXX the code uses \code{model@covariance@s2}
##'     which may not exist if some other class of object or even of
##'     covariance is used.
##' 
##' @references
##' 
##' D. Ginsbourger (2009), \emph{Multiples metamodeles pour
##' l'approximation et l'optimisation de fonctions numeriques
##' multivariables}, Ph.D. thesis, Ecole Nationale Superieure des
##' Mines de Saint-Etienne, 2009.
##' 
##' J. Mockus (1988), \emph{Bayesian Approach to Global
##' Optimization}. Kluwer academic publishers.
##' 
##' T.J. Santner, B.J. Williams, and W.J. Notz (2003), \emph{The
##' design and analysis of computer experiments}, Springer.
##' 
##' M. Schonlau (1997), \emph{Computer experiments and global
##' optimization}, Ph.D. thesis, University of Waterloo.
##'
##' @keywords models optimize
##'
##' @export EI.grad
##' @importFrom stats pnorm
##' 
##' @examples
##' set.seed(123)
##' ## =========================================================================
##' ## 	EI Surface Associated with an Ordinary Kriging Model for the Branin    
##' ##  Function Known at a 9-Points Factorial Design  
##' ## ========================================================================
##' 
##' ## a 9-points factorial design, and the corresponding response
##' ## ===========================================================
##' d <- 2; n <- 9
##' design.fact <- expand.grid(x1 = seq(0, 1, length=3), x2 = seq(0, 1, length=3))
##' y <- apply(design.fact, 1, branin) 
##' 
##' ## model fit
##' ## =========
##' fit1 <- km(~1, design = design.fact, response = y, covtype = "gauss",
##'            control = list(pop.size = 50, trace = FALSE), parinit = c(0.5, 0.5))
##'
##' ## graphics
##' ## ========
##' n.grid <- 9  # Increase to 50 for a nicer picture
##' x.grid <- y.grid <- seq(0, 1, length = n.grid)
##' design.grid <- expand.grid(x1 = x.grid, x2 = y.grid)
##' EI.grid <- apply(design.grid, MARGIN = 1, FUN = EI, model = fit1)
##' 
##' z.grid <- matrix(EI.grid, n.grid, n.grid)
##' 
##' contour(x.grid, y.grid, z.grid, 20)
##' title("Expected Improvement for the Branin function known at 9 points")
##' points(design.fact[ , 1], design.fact[ , 2], pch = 17, col = "blue")
##' 
##' # graphics
##' n.gridx <- 5  # increase to 15 for nicer picture
##' n.gridy <- 5  # increase to 15 for nicer picture
##' x.grid2 <- seq(0, 1, length = n.gridx) 
##' y.grid2 <- seq(0, 1, length = n.gridy) 
##' design.grid2 <- expand.grid(x1 = x.grid2, x2 = y.grid2)
##' 
##' EI.envir <- new.env()	
##' environment(EI) <- environment(EI.grad) <- EI.envir 
##'
##' for(i in 1:nrow(design.grid2)) {
##' 	x <- design.grid2[i,  ]
##' 	ei <- EI(x, model = fit1, envir = EI.envir)
##' 	eigrad <- EI.grad(x , model = fit1, envir = EI.envir)
##' 	if (!(is.null(ei))) {
##' 	suppressWarnings(arrows(x$x1, x$x2,
##' 	   x$x1 + eigrad[1] * 2.2 * 1e-4,
##'        x$x2 + eigrad[2] * 2.2 * 1e-4, 
##' 	   length = 0.04, code = 2, col = "orange", lwd = 2))
##' 	}
##' }
##' 
EI.grad <- function(x, model, plugin = NULL, type = c("UK", "SK"),
                    minimization = TRUE, envir = NULL,
                    proxy = FALSE){ 

    type <- match.arg(type)
    
    ## =========================================================================
    if (is.null(plugin)){ 
        if (minimization) plugin <- min(model@y)
        else plugin <- -max(model@y)
    }
    m <- plugin
    
    ## =========================================================================
    ## Convert x in proper format(s)
    ## =========================================================================

    d <- length(x)
    if (d != model@d) stop("'x' does not have the right length") 
    newdata.num <- as.numeric(x)
    newdata <- data.frame(t(newdata.num))
    colnames(newdata) <- colnames(model@X)
    
    ## Get quantities related to the prediction
    if (is.null(envir)) {
        
        pred <- predict(object = model, newdata = newdata, type = type,
                        checkNames = FALSE, se.compute = TRUE,
                        cov.compute = FALSE, deriv = TRUE)
        
        kriging.mean <- pred$mean
        if (!minimization) kriging.mean <- -kriging.mean

        kriging.sd <- pred$sd
        xcr <- (m - kriging.mean) / kriging.sd
        xcr.prob <- pnorm(xcr)
        xcr.dens <- dnorm(xcr)    

    } else {
        
        ## Note that we do not need kriging.mean
        pred <- envir$pred
        xcr <- envir$xcr
        xcr.prob   <- envir$xcr.prob
        xcr.dens   <- envir$xcr.dens
        kriging.sd <- pred$sd

    }
    
    ## =========================================================================
    ## Pursue calculation only if standard deviation is non-zero
    ## =========================================================================
    
    if (kriging.sd / sqrt(model@covariance@sd2) < 1e-06)  {
        EI.grad <- rep(0, d)
    } else  { 
    
        kriging.mean.grad <- as.vector(pred$mean.deriv)
        if (!minimization) kriging.mean.grad <- -kriging.mean.grad
        
        if (proxy) {
            EI.grad <- - kriging.mean.grad
        } else {
            kriging.sd.grad <- pred$sd.deriv
            EI.grad <- - kriging.mean.grad * xcr.prob + kriging.sd.grad * xcr.dens
        }
    }

    return(EI.grad)
    
}
