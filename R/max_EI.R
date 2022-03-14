## *****************************************************************************

##' @title Maximization of the Expected Improvement Criterion
##' 
##' @description Given an object in
##'     \code{\link[DiceKriging]{km-class}} and a set of tuning parameters
##'     (\code{lower}, \code{upper}, \code{parinit}, and
##'     \code{genoud_args}), \code{max_EI} performs the maximization of
##'     the Expected Improvement criterion and delivers the next point
##'     to be visited in an EGO-like procedure.
##' 
##' @details The latter maximization relies on a genetic algorithm
##'     using derivatives, \code{\link[rgenoud]{genoud}}. This
##'     function plays a central role in the package since it is in
##'     constant use in the proposed algorithms. It is important to
##'     remark that the information needed about the objective
##'     function reduces here to the vector of response values
##'     embedded in \code{model} (no call to the objective function or
##'     simulator).
##' 
##' The current minimum of the observations can be replaced by an
##' arbitrary value given in \code{plugin}, which is useful in
##' particular in noisy frameworks.
##' 
##' @param model An object of class \code{"km"} as created by using
##'     \code{\link[DiceKriging]{km}}.
##' 
##' @param plugin Optional scalar: if provided, it replaces the
##'     minimum of the current observations.
##' 
##' @param type Character \code{"UK"} (default) or \code{"SK"} giving
##'     the kriging type.
##'
##' @param lower,upper Numeric vectors of length \eqn{d} giving the
##'     lower and upper bounds for the variables to be optimized over.
##' 
##' @param parinit Optional numeric vector of initial values for the
##'     variables to be optimized over.
##'
##' @param minimization Logical specifying if EI is used in
##'     minimization or in maximization. Of course, this concerns
##'     the objective function to be optimized, not the EI which
##'     is always maximized.
##' 
##' @param genoud_args Optional named list of arguments for the
##'     \code{\link[rgenoud]{genoud}} optimization.
##'      Some arguments can not be used (such as \code{fn} or \code{gr})
##'      and an error will occur if they are given.
##'      Also note that for some arguments the default values of \code{genoud}
##'      have been changed:
##'      \itemize{
##'         \item{\code{pop.size}  }{default
##'             \code{ifelse(d < 6, 3 * 2^d, 32 * d)} where \code{d} is
##'              the number of inputs,}
##'          \item{\code{max.generations}  }{default: 12,}
##'          \item{\code{wait.generations}  }{default: 2,}
##'          \item{\code{BFGSburnin}  }{ default: 2.}
##'      }
##'      The values given in \code{genoud_args} are passed to the corresponding
##'      arguments of the function \code{\link[rgenoud]{genoud}}. 
##' 
##' @return A list with components:
##' \itemize{
##'     \item{\code{par} }{A numeric matrix with 1 row and \code{d} columns 
##'         the best set of parameters found.}
##'     \item{\code{value} }{A numeric vector with length one,
##'        giving the value of expected improvement at \code{par}.}
##' }
##' The matrix \code{par} will often have to be coerced into a data frame
##' (see \bold{Examples} or into a numeric vector.
##' @author David Ginsbourger, Olivier Roustant and Victor Picheny.
##' 
##' @references
##' 
##' D. Ginsbourger (2009), \emph{Multiples métamodèles pour
##' l'approximation et l'optimisation de fonctions numériques
##' multivariables}, Ph.D. thesis, Ècole Nationale Supérieure des
##' Mines de Saint-Étienne, 2009.
##' 
##' D.R. Jones, M. Schonlau, and W.J. Welch (1998), "Efficient global
##' optimization of expensive black-box functions", \emph{Journal of
##' Global Optimization}, \bold{13}, 455-492.
##' 
##' W.R. Mebane, Jr. and J.S. Sekhon
##' (2011). "Genetic Optimization Using Derivatives: The rgenoud Package for R".
##' \emph{Journal of Statistical Software}, \bold{42}(11), 1-26. URL
##' \url{http://www.jstatsoft.org/v42/i11/}.
##' 
##' @keywords optimize
##'
##' @export max_EI
##' 
##' @examples
##' library(rgenoud)
##' set.seed(123)
##' 
##' ## =========================================================================
##' ##  "ONE-SHOT" EI-MAXIMIZATION OF THE BRANIN FUNCTION 
##' ## 	KNOWN AT A 9-POINTS FACTORIAL DESIGN         
##' ## =========================================================================
##' 
##' ## a 9-points factorial design, and the corresponding response
##' ## ===========================================================
##' d <- 2; n <- 9
##' design.fact <- expand.grid(x1 = seq(0, 1, length = 3),
##'                            x2 = seq(0, 1, length = 3)) 
##' y.branin <- apply(design.fact, 1, branin)
##' 
##' ## model identification
##' ## ====================
##' fit.branin <- km(~1, design = design.fact, response = y.branin, 
##'                  covtype="gauss", control = list(pop.size = 50, trace = FALSE),
##'                  parinit = c(0.5, 0.5))
##' 
##' ## EGO one step
##' ## ============
##' lower <- rep(0.0, d); upper <- rep(1.0, d)    
##' oEGO <- max_EI(fit.branin, lower = lower, upper = upper, 
##'                genoud_args = list(pop.size = 20, BFGSburnin = 2))
##' print(oEGO)
##' 
##' ## graphics
##' ## =========
##' g <- contours(fit.branin, which = "", other = "branin") +
##'          ggtitle("Branin Function. EI maximized at the point shown with a lozenge.") +
##'               geom_point(data = data.frame(oEGO$par),
##'                          mapping = aes(x = x1, y = x2), shape = 23,
##'                          col = "purple", fill = "lavender", size = 4)
##' g
##' 
##' ## =========================================================================
##' ## "ONE-SHOT" EI-MAXIMIZATION OF THE CAMELBACK FUNCTION
##' ##	KNOWN AT A 16-POINTS FACTORIAL DESIGN            
##' ## =========================================================================
##' \dontrun{
##' ## a 16-points factorial design, and the corresponding response
##' d <- 2; n <- 16
##' design.fact <- expand.grid(x1 = seq(0, 1, length = 4),
##'                            x2 = seq(0, 1, length = 4)) 
##' yCB <- apply(design.fact, 1, camelback)
##' 
##' ## model fit
##' ## =========
##' fitCB <- km(~1, design = design.fact, response = yCB,
##'             covtype = "gauss",
##'             control = list(pop.size = 50, trace = FALSE),
##'             parinit = c(0.5, 0.5))
##' 
##' ## EI maximization
##' ## ================
##' lower <- rep(0.0, d); upper <- rep(1.0, d)   
##' oCB <- max_EI(fitCB, lower = lower, upper = upper, 
##'               genoud_args = list(pop.size = 20, BFGSburnin = 2))
##' print(oCB)
##' 
##' ## graphics
##' ## ========
##' g <- contours(fitCB, which = "", other = "camelback") +
##'          ggtitle("Camel back Function. EI maximized at the lozenge.") +
##'               geom_point(data = data.frame(oCB$par),
##'                          mapping = aes(x = x1, y = x2), shape = 23,
##'                          col = "purple", fill = "lavender", size = 4)
##' g
##' }
##' 
##' ## =========================================================================
##' ##  "ONE-SHOT" EI-MAXIMIZATION OF THE GOLDSTEIN-PRICE FUNCTION 
##' ## 	     KNOWN AT A 9-POINTS FACTORIAL DESIGN              
##' ## =========================================================================
##' 
##' \dontrun{
##' ## a 9-points factorial design, and the corresponding response
##' d <- 2; n <- 9
##' design.fact <- expand.grid(x1 = seq(0, 1, length = 3),
##'                            x2 = seq(0, 1, length = 3))
##' 
##' y.gP <- apply(design.fact, 1, goldsteinPrice)
##' 
##' ## model identification
##' fit.gP <- km(~1, design = design.fact, response = y.gP, 
##'               covtype="gauss",
##'               control = list(pop.size = 50, max.generations = 50, 
##'                              wait.generations = 5, BFGSburnin = 10,
##'                              trace = FALSE),
##'               parinit = c(0.5, 0.5), optim.method = "gen")
##' 
##' ## EI maximization
##' ## ===============
##' lower <- rep(0.0, d); upper <- rep(1.0, d);    
##' oEGO.gP <- max_EI(model = fit.gP, lower = lower, upper = upper,
##'                   genoud_args = list(pop.size = 50, max.generations = 50,
##'                                      wait.generations = 5, BFGSburnin = 10))
##' print(oEGO.gP)
##' 
##' ## graphics
##' ## ========
##' g <- contours(fit.gP, which = "", other = "goldsteinPrice") +
##'          ggtitle("GoldsteinPrice Function. EI maximized at the lozenge") +
##'               geom_point(data = data.frame(oEGO.gP$par),
##'                          mapping = aes(x = x1, y = x2), shape = 23,
##'                          col = "purple", fill = "lavender", size = 4)
##' g
##' 
##' }
##' 
max_EI <-function(model, plugin = NULL, type = "UK",
                  lower, upper, parinit = NULL,
                  minimization = TRUE,
                  genoud_args = NULL) {
    
    type <- match.arg(type)
    if (is.null(plugin)) plugin <- min(model@y)
    
    d <- ncol(model@X)
    if (is.null(genoud_args)) genoud_args <- list()
    
    args <- formals(genoud)
    
    args$fn <- EI_with_grad
    args$nvars <- d
    args$max <- TRUE

    ## =========================================================================
    ##  XXXY make a function for that? In 'DiceOptim', the default
    ##  values of 'genoud' are overloaded for all the criteria, but this is
    ##  not very well documented.
    ##  =========================================================================

    ind <- match(genoud_args, c("fn", "max", "gr")) 
    ind1 <- !is.na(ind)
    if (any(ind1)) { 
        stop("\"fn\", \"max\" and \"gr\" ncan not be given in ",
             "'genoud_args'")
    }
    
    if (is.null(genoud_args$print.level)) args$print.level <- 1
    N <- ifelse(d <= 6, 3 * 2^d, 32 * d) 
    ## provide value (not provided by default)
    if (is.null(genoud_args$control$maxit)) {
        args$control$maxit <- N
    } 
    ## replace default
    if (is.null(genoud_args$pop.size)) {
        args$pop.size <- N
    } else {
        args$pop.size <- genoud_args$pop.size
    }
    ## replace defaut (0.001)
    if (is.null(genoud_args$solution.tolerance)) {
        args$solution.tolerance <- 1e-21
    } else {
        args$solution.tolerance <- genoud_args$solution.tolerance
    }
    ## replace default (100)
    if (is.null(genoud_args$max.generations)) {
        args$max.generations <- 12
    } else {
        args$max.generations <- genoud_args$max.generations
    }
    ## replace default (10)
    if (is.null(genoud_args$wait.generations)) {
        args$wait.generations <- 2
    } else {
        args$wait.generations <- genoud_args$wait.generations
    }
    ## replace default (0)
    if (is.null(genoud_args$BFGSburnin)) {
        args$BFGSburnin <- 2
    } else {
        args$BFGSburnin <- genoud_args$BFGSburnin
    }
    ## doublon in args
    if (!is.null(genoud_args$starting.values)) {
        stop("Please use 'parinit' instead of 'genoud_args$starting.values'")
    }
    ## If 'parinit' is not provide, 'starting.values' will be NULL (default)
    if (is.null(parinit)) {
        args$starting.values <- lower + runif(d) * (upper - lower)
    }
    
    args$Domains <- cbind(lower, upper)

    ## CAUTION it is required to evaluate this line here
    if (is.null(genoud_args$boundary.enforcement)) {
        args$boundary.enforcement <- 0 
    } else {
        args$boundary.enforcement <- genoud_args$boundary.enforcement
    }
        
    if (is.null(genoud_args$optim.method)) {
        args$optim.method <- ifelse(args$boundary.enforcement < 2, "BFGS",  "L-BFGS-B")
    } else {
        args$optim.method <- genoud_args$optim.method
    }
          
    args$genoud_args <- NULL
    
    args <- c(args,                ## formals of 'genoud'
              model = model,       ## further arguments for 'EI', passed by dots
              plugin = plugin,
              type = type,
              minimization = minimization)
    args[["..."]] <- NULL    
    
    o <- do.call(genoud_cache, args)
    
    o$par <- t(as.matrix(o$par))
    colnames(o$par) <- colnames(model@X)
    o$value <- as.matrix(o$value)
    colnames(o$value) <- "EI"  

    return(list(par = o$par,
                value = o$value))
    
}

 
