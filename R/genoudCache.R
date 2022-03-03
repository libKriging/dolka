## *****************************************************************************

##' @description Run the \code{\link[rgenoud]{genoud}} function from
##'     the \pkg{rgenoud} package with the objective and the gradient
##'     defined in \emph{one} function. The result is global
##'     minimization of a real-valued function of a numeric vector
##'     with length \eqn{d} to be provided by \code{fn}. The R
##'     function \code{fn} must provide the values of the objective
##'     and its gradient.
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
##'     not been tested enough yet.
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
##' dom <- cbind(lower = rep(0, 2), upper = rep(0, 2))
##' library(rgenoud)
##' 
##' ## emulate a costly-to-evaluate-alone gradient
##' ## ===========================================
##' braninDer <- function(x) {
##'    Sys.sleep(0.01)
##'    braninGrad(x)$gradient
##' }
##' 
##' ## separate objective and gradient functions
##' ## =========================================
##' te <- system.time(res <- genoud(fn = branin, nvars = 2, Domains = dom,
##'                                 gr = braninDer))
##'
##' ## gradient "cached"
##' ## ================
##' teCache <- system.time(resCache <- genoudCache(fn = braninGrad, nvars = 2,
##'                                                Domains = dom))
##' rbind("genoud" = te, "genoudCache" = teCache)
##' c("genoud" = res$value, "genoudCache" = resCache$value)
##' }
genoudCache <- function(fn, nvars, max = FALSE, pop.size = 1000,
                        max.generations = 100, 
                        wait.generations = 10,
                        hard.generation.limit = TRUE,
                        starting.values = NULL, 
                        MemoryMatrix = TRUE,
                        Domains = NULL,
                        default.domains = 10, 
                        solution.tolerance = 0.001,
                        ## gr = NULL,
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
        stop("This function requires the '' package")
    }

    mc <- as.list(match.call())
    Ldots <- list(...)

    if (!is.na(match("gr", names(Ldots)))) {
        stop("'gr' is not a formal argument of the 'genoudCache' function ",
             "and it can not either be passed in dots '...'")
    }

    envir <- new.env()
    
    fnCache <- function(x, ...) {
        res <- fn(x, ...)
        assign("gradient", value = res$gradient, envir = envir)
        res$objective
    }
    
    grCache <- function(x, ...) {
        get("gradient", envir = envir)
    }

    ## is this really useful?
    environment(fnCache) <- environment(grCache) <- envir

    args <- list()
    ## copy/replace the default values
    if (length(mc) > 1) {
        nms <- names(mc)
        for (i in 2:length(mc)) {
            nm <- nms[i]
            args[[nm]] <- mc[[nm]]
        }
    }
    args[["fn"]] <- fnCache
    args[["gr"]] <- grCache
    
    res <- do.call(genoud, args = args)

    res
    
}
