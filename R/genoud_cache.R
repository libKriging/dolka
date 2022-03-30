## *****************************************************************************

##' @description Run the \code{\link[rgenoud]{genoud}} function from
##'     the \pkg{rgenoud} package with the objective and the gradient
##'     defined in \emph{one} function. The result is global
##'     minimization of a real-valued function of a numeric vector
##'     with length \eqn{d} to be provided by \code{fn}. The R
##'     function \code{fn} must provide the values of the objective
##'     and its gradient.
##'
##'     The formal arguments of this function and their default values
##'     are those of gthe \code{\link[rgenoud]{genoud}} function at
##'     the time when this function is written, based on
##'     version \code{5.8.3.0} of \pkg{rgenoud}.
##' 
##' @title Run the \code{genoud} Optimization Function with a Cached
##'     Gradient
##' 
##' @param fn A function to be minimized. This function should return
##'     a list with two named numeric elements \code{objective} and
##'     \code{gradient} representing the value of the objective
##'     (length 1) and the gradient (length \eqn{d}).
##' 
##' @param nvars,max,pop.size,max.generations See \code{\link[rgenoud]{genoud}}.
##' @param wait.generations,hard.generation.limit See \code{\link[rgenoud]{genoud}}.
##' @param starting.values See \code{\link[rgenoud]{genoud}}.
##' @param MemoryMatrix See \code{\link[rgenoud]{genoud}}.
##' @param Domains See \code{\link[rgenoud]{genoud}}.
##' @param default.domains See \code{\link[rgenoud]{genoud}}.
##' @param solution.tolerance See \code{\link[rgenoud]{genoud}}.
##' @param gr \bold{Can not be used}. This argument is for compatibility 
##' @param boundary.enforcement See \code{\link[rgenoud]{genoud}}.
##' @param lexical See \code{\link[rgenoud]{genoud}}.
##' @param gradient.check See \code{\link[rgenoud]{genoud}}.
##' @param BFGS See \code{\link[rgenoud]{genoud}}.
##' @param data.type.int  See \code{\link[rgenoud]{genoud}}.
##' @param hessian See \code{\link[rgenoud]{genoud}}.
##' @param unif.seed See \code{\link[rgenoud]{genoud}}.
##' @param int.seed See \code{\link[rgenoud]{genoud}}.
##' @param print.level See \code{\link[rgenoud]{genoud}}. 
##' @param share.type See \code{\link[rgenoud]{genoud}}. 
##' @param instance.number See \code{\link[rgenoud]{genoud}}.
##' @param output.path See \code{\link[rgenoud]{genoud}}.
##' @param output.append See \code{\link[rgenoud]{genoud}}.
##' @param project.path See \code{\link[rgenoud]{genoud}}.
##' @param P1,P2,P3,P4,P5,P6 See \code{\link[rgenoud]{genoud}}.
##' @param P7,P8,P9,P9mix See \code{\link[rgenoud]{genoud}}.
##' @param BFGSburnin See \code{\link[rgenoud]{genoud}}.
##' @param BFGSfn See \code{\link[rgenoud]{genoud}}.
##' @param BFGShelp See \code{\link[rgenoud]{genoud}}.
##' @param control See \code{\link[rgenoud]{genoud}}.
##' @param optim.method See \code{\link[rgenoud]{genoud}}.
##' @param transform See \code{\link[rgenoud]{genoud}}.
##' @param debug See \code{\link[rgenoud]{genoud}}.
##' @param cluster See \code{\link[rgenoud]{genoud}}.
##' @param balance See \code{\link[rgenoud]{genoud}}.
##' @param ... Further arguments to be passed to the function given in
##'     \code{fn}.
##' 
##' @return The result of a call to \code{\link[rgenoud]{genoud}}.
##'
##' @note The optimization is expected to be faster than
##'     \code{\link[rgenoud]{genoud}} if the function and its gradient
##'     are costly to evaluate separately. However this may not be the
##'     case for quite simple functions, due to the cost of setting
##'     the "cache" mechanism.
##'
##' @section Caution: This function is hightly experimental and has
##'     not been tested enough yet. It should be refactored on the
##'     basis of \code{\link{optim_cache}} which is more achieved.
##' 
##' @importFrom stats runif
##' @importFrom rgenoud genoud
##' @export
##' 
##' @examples
##'
##' \dontrun{
##' ## Note that in this example, gradient caching would not be worth it.
##' 
##' dom <- cbind(lower = rep(0.0, 2), upper = rep(1.0, 2))
##' library(rgenoud)
##' 
##' ## emulate a costly-to-evaluate-alone gradient
##' ## ===========================================
##' braninDer <- function(x) {
##'    Sys.sleep(0.01)
##'    branin_with_grad(x)$gradient
##' }
##' 
##' ## separate objective and gradient functions
##' ## =========================================
##' te <- system.time(res <- genoud(fn = branin, nvars = 2,
##'                                 unif.seed = 123,
##'                                 int.seed = 456,
##'                                 Domains = dom,
##'                                 gr = braninDer))
##'
##' ## gradient "cached"
##' ## ================
##' teCache <- system.time(resCache <- genoud_cache(fn = branin_with_grad, nvars = 2,
##'                                                 unif.seed = 123,
##'                                                 int.seed = 456,
##'                                                 Domains = dom))
##' rbind("genoud" = te, "genoud_cache" = teCache)
##' c("genoud" = res$value, "genoud_cache" = resCache$value)
##'
##' ## Check the use of ...
##' ## ====================
##' braninShift <- function(x, shift = 0) {
##'     res <- branin_with_grad(x)
##'     res$objective <- res$objective + shift
##'     res
##'  }
##'  resShifted <- genoud_cache(fn = braninShift, nvars = 2,
##'                             unif.seed = 123,
##'                             int.seed = 456,
##'                             Domains = dom,
##'                             shift = 100)
##' c("genoud_cache" = resCache$value, "genoudShifted" = resShifted$value - 100)
##' }
genoud_cache <- function(fn, nvars, max = FALSE, pop.size = 1000,
                         max.generations = 100, 
                         wait.generations = 10,
                         hard.generation.limit = TRUE,
                         starting.values = NULL, 
                         MemoryMatrix = TRUE,
                         Domains = NULL,
                         default.domains = 10, 
                         solution.tolerance = 0.001,
                         gr = NULL,
                         boundary.enforcement = 0, 
                         lexical = FALSE,
                         gradient.check = TRUE,
                         BFGS = TRUE,
                         data.type.int = FALSE, 
                         hessian = FALSE,
                         unif.seed = round(runif(1, 1, 2147483647L)), 
                         int.seed = round(runif(1, 1, 2147483647L)),
                         print.level = 2, 
                         share.type = 0,
                         instance.number = 0,
                         output.path = "stdout", 
                         output.append = FALSE,
                         project.path = NULL,
                         P1 = 50, P2 = 50, P3 = 50,
                         P4 = 50, P5 = 50, P6 = 50,
                         P7 = 50, P8 = 50, P9 = 0, 
                         P9mix = NULL,
                         BFGSburnin = 0,
                         BFGSfn = NULL,
                         BFGShelp = NULL, 
                         control = list(),
                         optim.method = ifelse(boundary.enforcement < 2,
                                               "BFGS", "L-BFGS-B"),
                         transform = FALSE,
                         debug = FALSE, 
                         cluster = FALSE,
                         balance = FALSE,
                         ...) {
    
    if (!requireNamespace("rgenoud", quietly = TRUE)) {
        stop("This function requires the 'rgenoud' package")
    }

    if (!is.null(gr)) {
        stop("'gr' must not be given in this function!")
    }
    
    mc <- as.list(match.call())

    ## =========================================================================
    ## 'nmArgs' contains the names of the formals in the call that are
    ## not formals of `genoud`. These should be passed to `fn` hence
    ## must be formals of `fn`. 
    ##
    ## Caution the first name in the call is "" (to mean the present
    ## function). Note that we do not pass all the arguments ofd 'fn'
    ## to 'fnCache' and 'grCache', but only the first one. The other
    ## arguments are found in 'envf'. We use names with dots, since 'res' or
    ## 'gradient' could be formal arguments of 'fn'.
    ## =========================================================================
    
    nmArgs <- setdiff(names(mc), names(formals(genoud_cache)))
    nmF <- names(formals(fn))
    nmArgs[1] <- nmF[1]
    if (!all(nmArgs %in% nmF)) {
        Pb <- setdiff(nmArgs, nmF)
        stop("Arguments ", Pb, " can not be passed to 'fn'") 
    }

    L <- list(NA)
    names(L) <- nmF[1]
    
    if (length(nmArgs) > 1) {
        for (i in 2:length(nmArgs)) {
            nm <- nmArgs[i]
            L[[nm]] <- eval(mc[[nm]], envir = parent.frame())
        }
    }

    envf <- new.env()
    
    fnCache <- function(x, envir) {
        L[[1]] <- x
        .res <- do.call(fn, L)
        assign(".gradient", value = .res$gradient, envir = envir)
        .res$objective
    }
    
    environment(fnCache) <- envf

    grCache <- function(x, envir) {
        get(".gradient", envir = envir)
    }
    
    environment(grCache) <- envf

    ## =========================================================================
    ## The arguments of the call to 'genoud' include the arguments
    ## given in the call which are not arguments of 'fn', plus 'fn' and 'gr'.
    ## =========================================================================

    args <- list()
    nms2 <- setdiff(intersect(names(mc), names(formals(genoud))), c("fn", "gr"))
    if (length(nms2) > 0) {
        for (i in 1:length(nms2)) {
            nm <- nms2[i]
            args[[nm]] <- eval(mc[[nm]], envir = parent.frame())
        }
    }
    args[["fn"]] <- fnCache
    args[["gr"]] <- grCache
    args[["envir"]] <- envf
    
    res <- do.call(genoud, args = args)

    res
    
}


