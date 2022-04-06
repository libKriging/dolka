
## *****************************************************************************
##
##' @title Maximization of Multipoint Expected Improvement Criterion (qEI)
##'
##' @description Maximization of the \code{\link{qEI}} criterion. This
##'     is a refactoring of the function from the \pkg{DiceOptim}
##'     package. However only the \code{"genoud"} method can be used
##'     for now. This function calls the function \code{\link{genoud}}
##'     using analytical formulae of \code{\link{qEI}} and its
##'     gradient \code{\link{qEI.grad}}.
##'
##' @section Caution: this function may well be renamed
##'     \code{max_qEI_genoud} in a future version.
##' 
##' @details
##' The parameters of list \code{optimcontrol} include
##' \code{optimcontrol$method}, with character value \code{"BFGS"}
##' (default) or \code{"genoud"}, specifying the method used to
##' maximize the criterion. This is irrelevant when \code{crit} is
##' \code{"CL"} because this method always uses \code{genoud}.
##'
##' \itemize{ 
##'
##'     \item{\bold{When \code{crit == "CL"}}} {
##'
##'          The elements of \code{optimcontrol} can have the following
##'          names and values.
##' 
##'          \itemize{
##'
##'              \item{\bold{\code{parinit}} }{Optional numeric
##'                   matrix of initial values. Must have
##'                   \code{model@d} columns, the number of rows is not
##'                   constrained.
##'               }
##'  
##'               \item{\bold{\code{L}} }{ The "liar". Either a
##'                   character in \code{c("max", "min", "mean")}, or
##'                   a scalar value specifying the liar. For the
##'                   character values: \code{"min"} sets the liar to
##'                   \code{min(model@y)}, \code{"max"} set the liar to
##'                   \code{max(model@y)}, and "mean" set it to the
##'                   prediction of the model. When \code{L} is
##'                   \code{NULL}, it is replaced by \code{"min"} if
##'                   \code{minimization} is \code{TRUE}, or by \code{"max"}
##'                   otherwise.  }
##'
##'               \item{\bold{Formal arguments of \code{\link{genoud}}}.}{
##'                   These include \code{pop.size} with
##'                   default : \code{3 * 2^d} when \code{d <
##'                   6} and \code{32 * d} otherwise, where \code{d}
##'                   is the number of inputs.  One can also use \code{max.generations}
##'                   (default: 12), \code{wait.generations} (default:
##'                   2) and \code{BFGSburnin} (default: 2).}  } }
##'
##'     \item{\bold{When \code{optimcontrol$method == "BFGS"}} } {
##'
##'          The elements of \code{optimcontrol} can have the
##'          following names and values.
##' 
##'          \itemize{
##' 
##'               \item{\bold{\code{nStarts}} } {(default: 4).}
##' 
##'               \item{\bold{\code{fastCompute}} } {Logical. If \code{TRUE}
##'                    (default), a fast
##'                    approximation method based on a semi-analytic
##'                    formula is used, see Marmin (2014) for details.}
##' 
##'              \item{\bold{\code{samplingFun}} } {A function which
##'                   samples a batch of starting point
##'                   Default : \code{\link{sampleFromEI}}.}
##'  
##'              \item{\bold{\code{parinit}} }{Optional 3d-array of
##'                   initial (or candidate)
##'                   batches. For each \code{k}, the slice \code{parinit[ , , k]}
##'                   is a matrix of size \code{c(npoints, d)}
##'                   representing one batch. The number of initial
##'                   batches \code{dim(parinit)[3]} is not
##'                   contrained and does not have to be equal to
##'                   \code{nStarts}. If there is too few initial
##'                   batches for \code{nStarts}, missing batches are
##'                   drawn with \code{samplingFun} (default
##'                   \code{NULL}).}
##'         }
##'     } 
##'
##'     \item{\bold{When \code{optimcontrol$method = "genoud"}} }{
##'        \itemize{
##' 
##'            \item{\bold{\code{optimcontrol$fastCompute}} }{Logical. If
##'                \code{TRUE} (default), a fast approximation method
##'                based on a semi-analytic formula is used, see
##'                Marmin(2014) for details.}
##' 
##'            \item{\bold{\code{parinit}} }{
##'                Optional numeric matrix of candidate starting
##'                points (one row corresponds to one point).}
##
##'            \item{\bold{The formal arguments of the \code{\link{genoud}} function}} {
##'                 Main parameters are \code{pop.size} (default:
##'                 \code{50 * d * npoints}),
##'                 \code{max.generations} (default: 5),
##'                 \code{wait.generations} (default: 2),
##'                 \code{BFGSburnin} (default: 2).
##'            }
##'        }
##'    }
##' }
##' 
##' @param model An object in \code{\link[DiceKriging]{km-class}} as created
##'    by using \code{\link[DiceKriging]{km}}.
##'
##' @param npoints Integer representing the desired number of
##'     iterations.
##' 
##' @param lower,upper Numeric vectors of lower and upper bounds.
##'
##' @param crit Character \code{"exact"} or \code{"CL"} specifying the
##'     criterion used. \code{"exact"} triggers the maximization of
##'     the multipoint expected improvement at each iteration (see
##'     \code{\link{max_qEI}}), \code{"CL"} applies the \emph{Constant
##'     Liar} heuristic.
##'
##' @param minimization Logical specifying if the qEI to be maximized
##'     is used in minimization or in maximization.
##' 
##' @param optimcontrol an optional list of control parameters for
##'     optimization. See details.
##' 
##' @param trace Integer. Level of verbosity.
##' 
##'
##' @return A list with components:
##' \itemize{
##' \item{par}{ A matrix with its rows containing the \code{npoints}
##'             input vectors found.}
##' \item{value}{ A value giving the qEI computed in \code{par}.}
##' }
##' @author Sébastien Marmin, Clément Chevalier and David Ginsbourger.
##'
##' @seealso \code{\link{qEI_with_grad}}.
##'
##' @references
##' 
##' C. Chevalier and D. Ginsbourger (2014) Learning and Intelligent
##' Optimization - 7th International Conference, Lion 7, Catania, Italy,
##' January 7-11, 2013, Revised Selected Papers, chapter Fast computation of
##' the multipoint Expected Improvement with applications in batch selection,
##' pages 59-69, Springer.
##' 
##' D. Ginsbourger, R. Le Riche, L. Carraro (2007), A Multipoint Criterion for
##' Deterministic Parallel Global Optimization based on Kriging. The
##' International Conference on Non Convex Programming, 2007.
##' 
##' D. Ginsbourger, R. Le Riche, and L. Carraro. Kriging is well-suited to
##' parallelize optimization (2010), In Lim Meng Hiot, Yew Soon Ong, Yoel
##' Tenne, and Chi-Keong Goh, editors, \emph{Computational Intelligence in
##' Expensive Optimization Problems}, Adaptation Learning and Optimization,
##' pages 131-162. Springer Berlin Heidelberg.
##' 
##' J. Mockus (1988), \emph{Bayesian Approach to Global Optimization}. Kluwer
##' academic publishers.
##' 
##' M. Schonlau (1997), \emph{Computer experiments and global optimization},
##' Ph.D. thesis, University of Waterloo.
##' 
##' @keywords optimize
##'
##' @importFrom stats optim
##' @importFrom DiceOptim sampleFromEI
##' 
##' @export max_qEI
##' 
##' @examples
##' set.seed(0)
##' ## 3-points EI maximization.
##' ## 9-points factorial design, and the corresponding response
##' d <- 2; n <- 9
##' design.fact <- expand.grid(x1 = seq(0, 1, length = 3),
##'                            x2 = seq(0, 1, length = 3)) 
##' response.branin <- apply(design.fact, 1, branin)
##' lower <- c(0, 0); upper <- c(1, 1)
##' 
##' ## number of point in the batch
##' batchSize <- 3
##' 
##' ## model fit
##' fitted.model <- km(~1, design = design.fact, response = response.branin, 
##'                    covtype = "gauss",
##'                    control = list(pop.size = 50, trace = FALSE),
##'                    parinit = c(0.5, 0.5))
##' 
##' # maximization of qEI
##' 
##' ## With a multistarted BFGS algorithm
##' ## ==================================
##' \dontrun{
##'     maxBFGS <- max_qEI(model = fitted.model, npoints = batchSize,
##'                        lower = lower, upper = upper, 
##'                        crit = "exact",
##'                        optimcontrol = list(nStarts = 3, method = "BFGS"))
##' 
##'     print(maxBFGS$value)
##' }
##' 
##' ## With a genetic algorithm using derivatives
##' ## ==========================================
##' maxGen  <- max_qEI(model = fitted.model, npoints = batchSize,
##'                    lower = lower, upper = upper, 
##'                    crit = "exact",
##'                    optimcontrol = list(nStarts = 3, method = "genoud",
##'                                        pop.size = 100, max.generations = 15))
##' print(maxGen$value)
##' 
##' ## With the constant liar heuristic
##' ## ================================
##' \dontrun{
##'     maxCL   <- max_qEI(model = fitted.model, npoints = batchSize,
##'                        lower = lower, upper = upper,
##'                        crit = "CL",
##'                        optimcontrol = list(pop.size = 20))
##'     print(maxCL$value)
##' }
max_qEI <- function(model,
                    npoints,
                    lower,
                    upper,
                    crit = c("exact", "CL"),
                    minimization = TRUE,
                    optimcontrol = NULL,
                    trace = 1) {

    message("This function is still not well-tested. ",
            "Use 'DiceOptim::max_qEI' in case of doubt")
    
    crit <- match.arg(crit)
    
    if (crit != "exact") {
        stop("for now, one can only use 'crit = \"exact\"'")
    }
    
    if (is.null(optimcontrol$method)) optimcontrol$method <- "BFGS"
    optim.method <- optimcontrol$method
    d <- model@d 
    parinit <- optimcontrol$parinit
    
    if (crit == "CL") {
        if (trace) {
            cat("o Using the 'max_qEI.CL' function\n")
            cat("=================================\n")
        }

        ## XXXY Not available for now (requires several imports from
        ## DiceOptim)
        ## res <- max_qEI.CL(model, npoints, optimcontrol$L,
        ##                   lower, upper,
        ##                   parinit = parinit,
        ##                   minimization = minimization,
        ##                   control = optimcontrol)
        ## res <- list(par = res$par,
        ##             value = matrix(qEI(x = res$par, model = model)))
        
    } else if (crit == "exact") {
        
        ## EI.envir <- new.env()
        ## environment(qEI) <- environment(qEI.grad) <- EI.envir
        
        LOWER <- c(apply(matrix(lower, d, 1), 1, rep, npoints))
        UPPER <- c(apply(matrix(upper, d, 1), 1, rep, npoints))
        
        if (optim.method == "BFGS") {
            if (is.null(parinit))
                parinit <- array(NaN, dim = c(npoints, model@d, 0))
            if (is.null(optimcontrol$nStarts))
                optimcontrol$nStarts <- 4
            if (is.null(optimcontrol$fastCompute))
                optimcontrol$fastCompute <- TRUE
            if (is.null(optimcontrol$sampleFun))
                optimcontrol$sampleFun <- sampleFromEI
            if (is.null(optimcontrol$gradNum))
                optimcontrol$gradNum <- FALSE
            if (is.null(optimcontrol$maxit))
                optimcontrol$maxit <- 100
            maxi <- 0
            startPoints <-
                optimcontrol$sampleFun(model = model,
                                       minimization = TRUE,
                                       n = npoints * optimcontrol$nStarts,
                                       lower = lower, upper = upper)
            for (i in 1:(optimcontrol$nStarts)) {
                if (i > length(parinit[1, 1, ])){
                    x <- startPoints[((i - 1) * npoints + 1):(i * npoints),]
                } else {
                    x <- parinit[ , , i]
                }
                if (!(optimcontrol$gradNum)) {
                    if (trace) {
                        cat("Optimizing 'qEI_with_grad' using 'optim_cache'\n")
                        cat("==============================================\n")
                    }
                    o <- optim_cache(par = c(x), fn = qEI_with_grad,
                                     control = list(trace = 0,
                                                    REPORT = 1,
                                                    fnscale = -1,
                                                    maxit = optimcontrol$maxit),
                                     method = "L-BFGS-B",
                                     lower = LOWER, upper = UPPER,
                                     model = model,
                                     ## envir = EI.envir,
                                     fastCompute = optimcontrol$fastCompute,
                                     minimization = minimization)
                } else {
                    cat("Optimizing 'qEI' using 'optim'\n")
                    cat("==============================\n")
                    o <- optim(par = c(x), fn = qEI, gr = NULL,
                               control = list(trace = 0,
                                              REPORT = 1,
                                              fnscale = -1,
                                              maxit = optimcontrol$maxit),
                               method = "L-BFGS-B",
                               lower = LOWER, upper = UPPER,
                               model = model,
                               ## envir = EI.envir,
                               fastCompute = optimcontrol$fastCompute,
                               minimization = minimization)
                }
                par <- matrix(o$par, ncol = d)
                value <- as.matrix(o$value)
                if (value >= maxi) {
                    maxi <- value
                    parmax <- par
                }
            }
            colnames(parmax) <- colnames(model@X)
            colnames(maxi) <- "EI"
            res <- list(par = parmax, value = maxi)
            ## res <- max_qEI.BFGS(model,npoints,lower,upper,
            ## parinit = NULL,minimization = minimization, control = control)
            
        } else if (optim.method == "genoud") {
            
            args <- formals(genoud)
            
            args$fn <- qEI_with_grad
            args$nvars <- npoints *d 
            args$max <- TRUE
            ind <- match(optimcontrol, c("fn", "max", "gr")) 
            ind1 <- !is.na(ind)
            if (any(ind1)) { 
                stop("\"fn\", \"max\" and \"gr\" ncan not be given in ",
                     "'optimcontrol'")
            }
            
            if (is.null(optimcontrol$pop.size)) {
                args$pop.size <- 50 * d
            } 
            if (is.null(optimcontrol$max.generations)) {
                args$max.generations <- 5
            } 
            if (is.null(optimcontrol$wait.generations)) {
                args$wait.generations <- 2
            } 
            if (is.null(optimcontrol$BFGSburnin)) {
                args$BFGSburnin <- 2
            } 
            if (is.null(optimcontrol$fastCompute)) {
                fastCompute <- TRUE
            } 
            if (is.null(optimcontrol$print.level)) {
                args$print.level <- 0
            }

            ## Added to cope with the dependence between default
            ## values
            
            ## if (is.null(optimcontrol$boundary.enforcement)) {
            ##     args$boundary.enforcement <- 0 
            ## } else {
            ##     args$boundary.enforcement <- optimcontrol$boundary.enforcement
            ## }
            args$boundary.enforcement <- 2
            
            if (is.null(optimcontrol$optim.method)) {
                args$optim.method <-
                    ifelse(args$boundary.enforcement < 2, "BFGS",  "L-BFGS-B")
            } else {
                args$optim.method <- optimcontrol$optim.method
            }
          
            args$Domains <- cbind(LOWER, UPPER)
            
            nmsGen <- names(formals(genoud))
            for(nm in names(optimcontrol)) {
                if (nm %in% nmsGen) {
                    args[[nm]] <- optimcontrol[[nm]]
                }
            }   

            args <- c(args,                ## formals of 'genoud'
                      model = model,       ## further arguments for 'EI', passed by dots
                      fastCompute = fastCompute,
                      minimization = minimization)

            args[["..."]] <- NULL            
            
            cat("Optimizing 'qEI_with_max' using 'genoud'\n")
            cat("========================================\n")
            
            o <- do.call(genoud_cache, args)
            
            ## o <- genoud_cache(qEI_with_grad,
            ##                   nvars = (d * npoints),
            ##                   max = TRUE,
            ##                   ## gr = qEI.grad,
            ##                   pop.size = optimcontrol$pop.size, 
            ##                   max.generations = optimcontrol$max.generations, 
            ##                   wait.generations = optimcontrol$wait.generations, 
            ##                   hard.generation.limit = TRUE,
            ##                   starting.values = c(optimcontrol$parinit), 
            ##                   MemoryMatrix = TRUE,
            ##                   Domains = domaine,
            ##                   default.domains = 10, 
            ##                   solution.tolerance = 1e-09,
            ##                   boundary.enforcement = 2, 
            ##                   lexical = FALSE,
            ##                   gradient.check = FALSE,
            ##                   BFGS = TRUE, 
            ##                   data.type.int = FALSE,
            ##                   hessian = FALSE, 
            ##                   unif.seed = floor(runif(1, max = 10000)), 
            ##                   int.seed = floor(runif(1, max = 10000)),
            ##                   print.level = 0,
            ##                   share.type = 0,
            ##                   instance.number = 0, 
            ##                   output.path = "stdout",
            ##                   output.append = FALSE,
            ##                   project.path = NULL, 
            ##                   P1 = 50, P2 = 50, P3 = 50,
            ##                   P4 = 50, P5 = 50, P6 = 50, 
            ##                   P7 = 50, P8 = 50, P9 = 0,
            ##                   P9mix = NULL,
            ##                   BFGSburnin = optimcontrol$BFGSburnin, 
            ##                   BFGSfn = NULL,
            ##                   BFGShelp = NULL,
            ##                   cluster = FALSE,
            ##                   balance = FALSE, 
            ##                   debug = FALSE,
            ##                   model = model,
            ##                   fastCompute = optimcontrol$fastCompute,
            ##                   minimization = minimization)
            
            o$par <- matrix(o$par, ncol = d)
            colnames(o$par) <- colnames(model@X)
            o$value <- as.matrix(o$value)
            colnames(o$value) <- "EI"
            res <- list(par = o$par, value = o$value)

            ## res <- max_qEI.gen(model,npoints,lower,upper,
            ## parinit = NULL,minimization = minimization, control = control)
            
        } else {
            stop(paste(paste("\'", optim.method, "\'", sep = ""),"is unknown."))
        }
    }
    
    return(res)
}


## max_qEI.CL <- function(model, npoints, L,
##                        lower, upper,
##                        parinit = NULL,
##                        minimization = TRUE,
##                        control = NULL) {
##     KB <- FALSE
##     n1 <- nrow(model@X)
##     lengthParinit <- length(parinit[1,])
##     if (is.null(L)) {
##         liar <- minimization * min(model@y) + (!minimization) * max(model@y)
##     } else if (L == "max") {
##         liar <- max(model@y)
##     } else if (L == "min") {
##         liar <- min(model@y)
##     } else if (L == "mean") {
##         KB <- TRUE
##     } else {
##         if ((!is.numeric(L)) || (length(L) != 1) || is.na(as.numeric(L)))
##             stop("control$L must be NULL, \"max\", \"min\", \"mean\" or ",
##                  "a scalar specifying the plugin.")
##         liar <- L
##     }
##     for (s in 1:npoints) {
##         if (s > lengthParinit){
##             startPoint <- NULL
##         } else {
##             startPoint <- parinit[ , s]
##         }

##         oEGO <- max_EI(model = model,
##                        lower = lower, upper = upper,
##                        parinit = startPoint,
##                        minimization = minimization,
##                        control = c(control, list(print.level = 0)))


##         ## XXXY Is'nt this something like 'update.km'? This could be a
##         ## cheap version with no re-estimation.
##         model@X <- rbind(model@X, oEGO$par)
        
##         if (KB) liar <- predict(object = model,newdata = oEGO$par,
##                                 se.compute = FALSE, cov.compute = FALSE,
##                                 light.return = TRUE, checkNames = FALSE)$mean

##         model@y <- rbind(model@y, liar, deparse.level = 0)   		
##         model@F <- trendMatrix.update(model, Xnew = data.frame(oEGO$par))	
        
##         if (model@noise.flag) {
##             ## heterogenous case : use 0 nugget for new points
##             model@noise.var = c(model@covariance@nugget, 0)
##         }	
##         model <- computeAuxVariables(model)
##     }
##     return(list(par = model@X[(n1 + 1):(n1 + npoints), , drop = FALSE],
##                 value = model@y[(n1 + 1):(n1 + npoints), , drop = FALSE]))
## }
