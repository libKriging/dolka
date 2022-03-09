##' @title Expected Improvement Criterion With Gradient
##' 
##' @description The function \code{EI_with_grad} computes the
##'     Expected Improvement at current location \code{x} and its
##'     gradient if wanted. The current minimum of the observations in
##'     \code{model} can be replaced by an arbitrary value (plugin),
##'     which is useful in particular in noisy frameworks.
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
##' @param proxy Optional logical. If \code{TRUE}, EI is replaced by
##'     the kriging mean, to be \emph{minimized}.
##' 
##' @param deriv Logical. If \code{TRUE} the result is a list with the
##'     two elements \code{objective} and \code{gradient}. Else
##'     the result is the value of the objective.
##' 
##' @param out_list Logical When \code{out_list} is \code{TRUE} the
##'     result is a \emph{list} with one element \code{objective}. If
##'     \code{deriv} is \code{TRUE} the list has a second element
##'     \code{gradient}. When \code{out_list} is \code{FALSE} the
##'     result is the numerical value of the objective, possibly
##'     having an attribute named \code{"gradient"}.
##'
##' @return The expected improvement as defined in \bold{Details}
##'     (for \code{EI}) or its gradient (for \code{EI.grad}).  If
##'     \code{plugin} is specified, its provided value will replace
##'     \eqn{\min Y(X)}{min Y(X)} in the formula. The EI and its
##'     gradient are numeric vectors with length \eqn{1} and \eqn{d}.
##'
##' @export EI_with_grad
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
##' 
EI_with_grad <- function (x, model, plugin = NULL, type = c("UK", "SK"),
                          minimization = TRUE, proxy = FALSE, deriv = TRUE,
                          out_list = deriv) {
    
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
        if (d != model@d) stop("'x' does not have the right lenght") 
        newdata.num <- as.numeric(x)
        newdata <- data.frame(t(newdata.num))
    } else {
        if (is.null(dim(x))) {
            d <- length(x)
            if (d != model@d){ stop("'x' does not have the right lenght") }
            newdata <- data.frame(t(as.numeric(x)))
        } else {
            d <- ncol(x)
            if (d != model@d){ stop("'x' does not have the right length") }
            newdata <- data.frame(x)
        }
    }
    colnames(newdata) <- colnames(model@X)
    
    ## =========================================================================
    ## The derivatives are computed only if they are to be stored in
    ## the environment. Note that since 'newdata' has only one row
    ## computing the covariance or not does not really matter.
    ## =========================================================================
    
    pred <- predict(object = model, newdata = newdata, type = type,
                    checkNames = FALSE, deriv = deriv)
    
    kriging.mean <- pred$mean
    
    if(!minimization) {
        kriging.mean <- -kriging.mean
    }
    kriging.sd   <- pred$sd
    
    if (proxy) {
        xcr <- xcr.prob <- xcr.dens <- NULL
        EI <- m - kriging.mean
    } else {
        xcr <- (m - kriging.mean) / kriging.sd
        xcr.prob <- pnorm(xcr)
        xcr.dens <- dnorm(xcr)	        
        EI <- (m - kriging.mean) * xcr.prob + kriging.sd * xcr.dens
    }

    if (!deriv) {
        if (out_list) EI <- list(objective = EI)
        return(EI)
    }
        
    ## =========================================================================
    ## Store the prediction and the auxiliary variables required (not
    ## very expansive to compute however).
    ## =========================================================================

    kriging.sd <- pred$sd
    xcr <- (m - kriging.mean) / kriging.sd
    xcr.prob <- pnorm(xcr)
    xcr.dens <- dnorm(xcr)    
 
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

    if (out_list) {
        return(list(objective = EI, gradient = EI.grad))
    } else {
        attr(EI, "gradient") <- EI.grad
        return(EI)
    }
}
