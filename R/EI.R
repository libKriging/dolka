##' @title Expected Improvement Criterion and its Gradient
##' 
##' @description The function \code{EI} computes the Expected
##'     Improvement at current location \code{x} while \code{EI.grad}
##'     compute the gradient at \code{x}. The current minimum of the
##'     observations in \code{model} can be replaced by an arbitrary
##'     value (plugin), which is useful in particular in noisy
##'     frameworks.
##'
##' @details The Expected Improvement (EI) is defined as \deqn{EI(x)
##'     := E\left[\{\min Y(X) - Y(x)\}_+ | Y(X) = y(X) \right],}{ E[{
##'     min Y(X) - Y(x) }_+ | Y(X) = y(X)]} where \eqn{X} is the
##'     current design of experiments and \eqn{Y} is the random
##'     process assumed to have generated the objective function
##'     \eqn{y} and \eqn{z_{+} := \max\{z, \, 0)}{z_+ = max(z, 0)}
##'     denotes the positive part of a real number \eqn{z}. The value
##'     of EI is non-negative but can be numerically zero close to the
##'     inputs used in \code{model}. The EI and its gradient are
##'     computed using their closed forms.
##' 
##' 
##' @param x A numeric vector representing the input for which one
##'     wishes to calculate EI. The length \eqn{d} of this vector must
##'     be equal to \eqn{d}, the dimension of the input space used for
##'     the kriging results in \code{model}.
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
##'     minimization or in maximization.
##'
##' @param envir An optional environment specifying where to assign
##'     intermediate values for future gradient calculations. Default
##'     is \code{NULL}.
##'
##' @param proxy Optional logical. If \code{TRUE}, EI is replaced by
##'     the kriging mean, to be \emph{minimized}.
##' 
##' @return The expected improvement as defined in \bold{Details}
##'     (for \code{EI}) or its gradient (for \code{EI.grad}).  If
##'     \code{plugin} is specified, its provided value will replace
##'     \eqn{\min Y(X)}{min Y(X)} in the formula. The EI and its
##'     gradient are numeric vectors with length \eqn{1} and \eqn{d}.
##'
##' @export EI
##'
##' @author David Ginsbourger, Olivier Roustant and Victor Picheny.
##'
##' @seealso \code{\link{max_EI}}, \code{\link{EGO.nsteps}},
##'     \code{\link{qEI}}
##'
##' @references
##' D. Ginsbourger (2009), \emph{Multiples métamodèles pour
##' l'approximation et l'optimisation de fonctions numériques
##' multivariables}, Ph.D. thesis, Ècole Nationale Supérieure des
##' Mines de Saint-Ètienne.
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
##' ## =========================================================================
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
##' ## computing the EI
##' ## ================
##' x <- c(x1 = 0.2, x2 = 0.4)
##' EI(x, model = fit1)
##' EI.grad(x, model = fit1)
##' 
##' ## graphics
##' ## ========
##' contours(object = fit1, which = character(0), grad = TRUE,
##'          other = "EI", otherGrad = "EI.grad",
##'          whereGrad = "grid", nGrid = 30) +
##'     ggtitle("Expected Improvement and its gradient")
##' 
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
        d <- ncol(x)
        if (d != model@d) stop("'x' does not have the right number of columns") 
        newdata.num <- as.numeric(x)
        newdata <- matrix(newdata.num, ncol = d)
    } else {
        if (is.null(dim(x))) {
            d <- length(x)
            if (d != model@d){ stop("x does not have the right size") }
            newdata <- matrix(as.numeric(x), ncol = d)
        } else {
            d <- ncol(x)
            if (d != model@d){
                stop("x does not have the right number of colums")
            }
            newdata <- matrix(as.numeric(x), ncol = d)
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
##' @keywords models optimize
##'
##' @export EI.grad
##' @importFrom stats pnorm
##' @rdname EI
##' 
##' 
EI.grad <- function(x, model, plugin = NULL, type = c("UK", "SK"),
                    minimization = TRUE, envir = NULL,
                    proxy = FALSE){ 

    type <- match.arg(type)
    
    ## =========================================================================
    if (is.null(plugin)) { 
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
        xcr.prob <- envir$xcr.prob
        xcr.dens <- envir$xcr.dens
        kriging.sd <- pred$sd

    }
    
    ## =========================================================================
    ## Pursue calculation only if standard deviation is non-zero
    ## =========================================================================
    
    if (kriging.sd / sqrt(model@covariance@sd2) < 1e-06)  {
        EI.grad <- rep(0, d)
    } else  { 
    
        kriging.mean.grad <- as.vector(pred$mean.deriv)
        if (!minimization) kriging.mean.grad <- - kriging.mean.grad
        
        if (proxy) {
            EI.grad <- - kriging.mean.grad
        } else {
            kriging.sd.grad <- as.vector(pred$sd.deriv)
            EI.grad <- - kriging.mean.grad * xcr.prob + kriging.sd.grad * xcr.dens
        }
    }

    return(EI.grad)
    
}
