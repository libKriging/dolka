## =============================================================================
##
##' @description Analytical expression of the multipoint expected
##'     improvement criterion, also known as the \eqn{q}-\emph{point
##'     expected improvement} and denoted by \eqn{q}-EI or \code{qEI}.
##' 
##' @title Computes the Multipoint Expected Improvement (qEI) Criterion
##'  
##' @param x A numeric matrix representing the set of input points
##'     (one row corresponds to one point) where to evaluate the qEI
##'     criterion.
##'
##' @param model An object of class \code{km}.
##'
##' @param plugin Optional scalar: if provided, it replaces the
##'     minimum of the current observations.
##'
##' @param type \code{"SK"} or \code{"UK"} (by default), depending
##'     whether uncertainty related to trend estimation has to be
##'     taken into account.
##'
##' @param minimization Logical specifying if EI is used in
##'     minimiziation or in maximization.
##'
##' @param fastCompute Logical. If \code{TRUE}, a fast approximation
##'     method based on a semi-analytic formula is used. See the
##'     reference Marmin (2014) for details.
##'
##' @param eps Numeric value of \eqn{epsilon} in the fast computation
##'     trick.  Relevant only if \code{fastComputation} is
##'     \code{TRUE}.
##'
##' @param deriv Logical. When \code{TRUE} the dervivatives of the
##'     kriging mean vector and of the kriging covariance (Jacobian
##'     arrays) are computed and stored in the environment given in
##'     \code{envir}.
##' 
##' @param out_list Logical. Logical When \code{out_list} is
##'     \code{TRUE} the result is a \emph{list} with one element
##'     \code{objective}. If \code{deriv} is \code{TRUE} the list has
##'     a second element \code{gradient}. When \code{out_list} is
##'     \code{FALSE} the result is the numerical value of the
##'     objective, possibly having an attribute named
##'     \code{"gradient"}.
##'
##' @param trace Integer level of verbosity.
##'
##' @return The multipoint Expected Improvement, defined as
##'     \deqn{qEI(X_{new}) := E\left[\{ \min Y(X)  - \min Y(X_{new}) \}_{+}
##'     \vert Y(X) = y(X) \right],}{qEI(Xnew) := E[{ min Y(X)  - min Y(Xnew) }_{+}
##'     \vert Y(X) = y(X)],}
##' where \eqn{X} is the current design of experiments,
##' \eqn{X_new}{Xnew} is a new candidate design, and
##' \eqn{Y} is a random process assumed to have generated the
##' objective function \eqn{y}.
##'
##' @author Sebastien Marmin, Clement Chevalier and David Ginsbourger.
##' 
##' @seealso \code{\link{EI}}
##'
##' @references
##' 
##' C. Chevalier and D. Ginsbourger (2014) Learning and Intelligent
##' Optimization - 7th International Conference, Lion 7, Catania,
##' Italy, January 7-11, 2013, Revised Selected Papers, chapter Fast
##' computation of the multipoint Expected Improvement with
##' applications in batch selection, pages 59-69, Springer.
##' 
##' D. Ginsbourger, R. Le Riche, L. Carraro (2007), A Multipoint
##' Criterion for Deterministic Parallel Global Optimization based on
##' Kriging. The International Conference on Non Convex Programming,
##' 2007.
##' 
##' S. Marmin (2014). Developpements pour l'evaluation et la
##' maximisation du critere d'amelioration esperee multipoint en
##' optimisation globale. Master's thesis, Mines Saint-Etienne
##' (France) and University of Bern (Switzerland).
##' 
##' D. Ginsbourger, R. Le Riche, and L. Carraro. Kriging is
##' well-suited to parallelize optimization (2010), In Lim Meng Hiot,
##' Yew Soon Ong, Yoel Tenne, and Chi-Keong Goh, editors,
##' \emph{Computational Intelligence in Expensive Optimization
##' Problems}, Adaptation Learning and Optimization, pages
##' 131-162. Springer Berlin Heidelberg.
##' 
##' J. Mockus (1988), \emph{Bayesian Approach to Global
##' Optimization}. Kluwer academic publishers.
##' 
##' M. Schonlau (1997), \emph{Computer experiments and global
##' optimization}, Ph.D. thesis, University of Waterloo.
##' 
##' @keywords models parallel optimization
##'
##' @importFrom mnormt dmnorm pmnorm
##' @export
##'
##' @examples
##' \donttest{
##' set.seed(007)
##' ## Monte-Carlo validation
##' 
##' ## a 4-d, 81-points grid design, and the corresponding response
##' ## ============================================================
##' d <- 4; n <- 3^d
##' design <- expand.grid(rep(list(seq(0, 1, length = 3)), d))
##' names(design) <- paste0("x", 1:d)
##' y <- apply(design, 1, hartman4)
##' 
##' ## learning
##' ## ========
##' model <- km(~1, design = design, response = y, control = list(trace = FALSE))
##' 
##' ## pick up 10 points sampled from the 1-point expected improvement
##' ## ===============================================================
##' q <- 10
##' X <- sampleFromEI(model, n = q)
##' 
##' ## simulation of the minimum of the kriging random vector at X
##' ## ===========================================================
##' t1 <- proc.time()
##' newdata <- as.data.frame(X)
##' colnames(newdata) <- colnames(model@X)
##' 
##' krig  <- predict(object = model, newdata = newdata, type = "UK",
##'                  se.compute = TRUE, cov.compute = TRUE)
##' mk <- krig$mean
##' Sigma.q <- krig$cov
##' mychol <- chol(Sigma.q)
##' nsim <- 300000
##' white.noise <- rnorm(n = nsim * q)
##' minYsim <- apply(crossprod(mychol, matrix(white.noise, nrow = q)) + mk,
##'                  MARGIN = 2, FUN = min)
##' 
##' ## simulation of the improvement (minimization)
##' ## ============================================
##' qImprovement <- min(model@y) - minYsim
##' qImprovement <- qImprovement * (qImprovement > 0)
##' 
##' ## empirical expectation of the improvement and confidence interval (95%)
##' ## ======================================================================
##' EIMC <- mean(qImprovement)
##' seq <- sd(qImprovement) / sqrt(nsim)
##' confInterv <- c(EIMC - 1.96 * seq, EIMC + 1.96 * seq)
##'
##' ## evaluate time
##' ## =============
##' tMC <- proc.time() - t1
##' 
##' ## MC estimation of the qEI
##' ## ========================
##' print(EIMC) 
##' 
##' ## qEI with analytical formula and with fast computation trick
##' ## ===========================================================
##' tForm <- system.time(qEI(X, model, fastCompute = FALSE))
##' tFast <- system.time(qEI(X, model))
##' 
##' rbind("MC" = tMC, "form" = tForm, "fast" = tFast)
##' 
##' }
##' 
qEI_with_grad <- function(x, model, plugin = NULL,
                          type = c("UK", "SK"),
                          minimization = TRUE,
                          fastCompute = TRUE,
                          eps = 1e-5,
                          deriv = TRUE,
                          out_list = TRUE,
                          trace = 1) {

    type <- match.arg(type)
    d <- model@d

    if (!is.matrix(x)) x <- matrix(x, ncol = d)
    nx <- length(x)

    if (trace) cat("length(x) = ", nx, "\n")
    
    xb <- rbind(model@X, x)

    ## =========================================================================
    ## XXXY Remove the rows of 'x' that are already present in
    ## 'model@X'. If no row remains, return 0 with derivative 0.
    ##
    ## Note that 'model@X' can have duplicated rows (possible with
    ## DiceKriging::km).
    ##
    ## CAUTION This check is not the same as the one done in
    ## 'qEI.grad' which only checks the presence of duplicated rows in
    ## rbind(model@X, x).
    ## =========================================================================
    
    nExp <- model@n
    x <- unique(round(xb, digits = 8))
    
    if ((nExp + 1) > length(x[ , 1])) {
        cat("AAAAAAAA\n")
        this.qEI <- 0
        if (deriv) {
            if (out_list) {
                this.qEI <- list("objective" = this.qEI,
                                 "gradient" = rep(0, nx))
            } else attr(this.qEI, "gradient") <- rep(0, nx)
        }
        return(this.qEI)
    }
    
    x <- matrix(x[(nExp + 1):length(x[ , 1]), ], ncol = d)
    q <- nrow(x)
    ## end 
        
    if (is.null(plugin) &&  minimization) plugin <- min(model@y)
    if (is.null(plugin) && !minimization) plugin <- max(model@y)

    if (nrow(x) == 1) {
        ## XXXY in the original code 'minimization' is not
        ## passed to 'EI'. BUG?
        if (trace) {
            cat("'x' has one row. Using 'EI_with_grad'\n")
            print(x)
        }
        return(EI_with_grad(x = x,
                            model = model,
                            plugin = plugin,
                            type = type,
                            minimization = minimization,
                            deriv = deriv,
                            out_list = out_list))
    }
    
    ## Compute the kriging prediction with derivatives
    pred  <- predict(object = model, newdata = x, type = type,
                     se.compute = TRUE,
                     cov.compute = TRUE,
                     deriv = deriv,
                     checkNames = FALSE)
    
    sigma <- pred$cov
    mu <- pred$mean
    
    pk <- first_term <- second_term <- rep(0, times = q)
    
    if (!fastCompute) {
        symetric_term <- matrix(0, q, q)
        non_symetric_term <- matrix(0, q, q)
    }
    
    for (k in 1:q) {
        ## covariance matrix of the vector Z^(k)
        Sigma_k <- covZk(sigma = sigma, index = k )
        
        ## mean of the vector Z^(k)
        mu_k <- mu - mu[k] 
        mu_k[k] <- -mu[k]
        if (minimization) mu_k <- -mu_k
        
        b_k <- rep(0, times = q)
        b_k[k] <- -plugin
        if(minimization) b_k <- -b_k
        pk[k] <- pmnorm(x = b_k - mu_k, varcov = Sigma_k, maxpts = q * 200)[1]
     
        first_term[k] <- (mu[k] - plugin) * pk[k]
        if (minimization) first_term[k] <- -first_term[k]
        
        if (fastCompute) {
	    second_term[k] <- (pmnorm(x = b_k + eps * Sigma_k[ , k] - mu_k,
                                      varcov = Sigma_k,
                                      maxpts = q * 200)[1] - pk[k]) / eps
        } else {
            for(i in 1:q){
        	non_symetric_term[k, i] <- Sigma_k[i, k]
                if (i >= k) {
                    mik <- mu_k[i]
                    sigma_ii_k <- Sigma_k[i,i]
	            bik <- b_k[i]
	            phi_ik <- dnorm(x = bik, mean = mik, sd = sqrt(sigma_ii_k))
                    
                    ##need c.i^(k) and Sigma.i^(k)
	            cik <- get_cik(b = b_k , m = mu_k , sigma = Sigma_k , i = i)
	            sigmaik <- get_sigmaik(sigma = Sigma_k , i = i)
	            Phi_ik <- pmnorm(x = cik, varcov = sigmaik,
                                     maxpts = (q - 1) * 200)[1]
	            symetric_term[k, i] <- phi_ik * Phi_ik
                }
            }
        }
    }
    
    if (!fastCompute) {
        symetric_term <- symetric_term + t(symetric_term)
        diag(symetric_term) <- 0.5 * diag(symetric_term) 
        second_term <- sum(symetric_term * non_symetric_term)
    }
    
    thisqEI <- sum(first_term, second_term)

    if (!deriv) {
        if (out_list) return(list("objective" = thisqEI))
        else return(thisqEI)
    }

    ## *************************************************************************
    ## Derivation part. Some code has been removed or is kept commented out.
    ## *************************************************************************

    if (FALSE) {
        if (!minimization && is.null(plugin)) {
            plugin <- -max(model@y)
        } 
        
        xb <- rbind(model@X, x)
        ux <- unique(round(xb, digits = 8))
        if(length(xb[ , 1]) != length(ux[ , 1])) {
            return (matrix(0, q, d))
        }

        ## XXXY BUG: if !minimization, this should not be done because it cancels
        ## the code about 10 lines above.
        if (is.null(plugin)) plugin <- min(model@y)
        
        if (q == 1) {
            if (!minimization) {
                stop("'qEI.grad' doesn't work in \'minimization = FALSE\' when ",
                     "dim = 1 (in progress).")
            }
            return(EI.grad(x, model, plugin, type, envir))
        }
    }
    
    ## ==========================================================================
    ## Kriging mean vector and covariance matrix
    ##
    ## o kriging.mean: vector with length q
    ##
    ## o kriging.cov:  2-dimensional array with dim c(q, q)
    ##
    ## These are 'mu' and 'sigma' above
    ## ==========================================================================

    kriging.mean <- pred[["mean"]]
    kriging.cov <- pred[["cov"]]
    
    ## ==========================================================================
    ## Jacobian arrays of the mean vector and the covariance matrix
    ##
    ## o kriging.mean.jacob: 3-dimensional array with dim c(q, d, q)
    ##
    ##          kriging.mean.jacob[l, j, k] = d(m_k) / d(x_lj)
    ##
    ## o kriging.cov.jacob: 4-dimensional array with dim  c(q, d, q, q)
    ## 
    ##          kriging.cov.jacob[l, j, k, i]= d(Sigma_ki) / (d x_lj)
    ##
    ## CAUTION Since the indexation of this code taken from DiceOptim
    ## differs from that in the predict method, 'aperm' is used.
    ## ==========================================================================
    
    kriging.mean.jacob  <- aperm(pred[["mean.deriv"]], perm = c(1, 3, 2))
    kriging.cov.jacob <- aperm(pred[["cov.deriv"]], perm = c(3, 4, 1, 2))
    
    if (!minimization) {
        kriging.mean <- -kriging.mean;
        kriging.mean.jacob <- -kriging.mean.jacob
    }
    
    ## Initialisation
    EI.grad <- matrix(0, q, d) # result
    b <- matrix(rep(0, q), q, 1)
    L <- -diag(rep(1, q))
    Dpk <- matrix(0, q, d)
    
    if (fastCompute) termB <- matrix(0, nrow = q, ncol = d)

    Sigk_dx <- array(0, dim = c(q, d, q, q))
    mk_dx <- array(0, dim = c(q, d, q))

    ## =========================================================================
    ## XXXY Notes YD:
    ##
    ## o Since 'bk' and 'mk' are used through indexed forms, they both
    ## could be vectors.
    ## 
    ## o Why not use bk - mk?
    ## =========================================================================
    
    ## First sum of the formula
    for (k in 1:q) {
        
        bk <- b
        bk[k, 1] <- plugin                   # creation of vector b^(k)
        Lk <- L
        Lk[ , k] <- rep(1, q)                # linear transform: Y to Zk
        tLk <- t(Lk)
        mk <- Lk %*% kriging.mean            # mean of Zk (m^(k) in the formula)
        Sigk <- Lk %*% kriging.cov %*% tLk   # covariance of Zk
        Sigk <- 0.5 * (Sigk + t(Sigk))       # numerical symetrization
        
        ## term 'A1' XXXY comment out
        ## ==========================
        
        ## if (is.null(envir)) {
        ##     EI.grad[k, ] <-  EI.grad[k, ] - kriging.mean.jacob[k, , k] *
        ##         pmnorm(x = t(bk), mean = t(mk), varcov = Sigk, maxpts = q * 200)
        ## } else {

        EI.grad[k, ] <-  EI.grad[k, ] - kriging.mean.jacob[k, , k] * pk[k]
        
        ## }
        
        ## compute gradient ans hessian matrix of the CDF term 'pk'
        gradpk <- GPhi(x = bk - mk, mu = b, Sigma = Sigk)
        hesspk <- HPhi(x = bk, mu = mk, Sigma = Sigk, gradient = gradpk)
        
        ## term 'A2'
        ## ========
        for (l in 1:q) {
            for (j in 1:d) {
                Sigk_dx[l, j, , ] <- Lk %*% kriging.cov.jacob[l, j, , ] %*% tLk
                mk_dx[l, j, ] <- Lk %*% kriging.mean.jacob[l, j, ]
                Dpk[l, j] <- 0.5 * sum(hesspk * Sigk_dx[l, j, , ]) -
                    crossprod(gradpk, mk_dx[l, j, ])
            }
        }
        EI.grad <- EI.grad + (plugin - kriging.mean[k]) * Dpk
        
        ## term 'B'
        ## =======
        if (fastCompute) {
            
            gradpk1 <- GPhi(x = bk - mk + Sigk[ , k] * eps, mu = b, Sigma = Sigk)
            hesspk1 <- HPhi(x = bk + Sigk[ , k] * eps,  mu = mk, Sigma = Sigk,
                           gradient = gradpk1)
            for (l in 1:q) {
                for (j in 1:d) {
                    f1 <- -crossprod(mk_dx[l, j, ], gradpk1) +
                        eps * crossprod(Sigk_dx[l, j, , k], gradpk1) +
                        0.5 * sum(Sigk_dx[l, j, , ] * hesspk1)
                    f  <- -crossprod(mk_dx[l, j, ], gradpk) +
                        0.5 * sum(Sigk_dx[l, j, , ] * hesspk)
                    termB[l, j] <- (f1 - f) / eps
                }
            }
        } else {
            B1 <- B2 <- B3 <- matrix(0, q, d)
            for (i in 1:q) {
                Sigk_ik <- Sigk[i, k]
                Sigk_ii <- Sigk[i, i]
                mk_i <- mk[i]
                mk_dx_i <- mk_dx[ , , i]
                bk_i <- bk[i, 1]
                ck_pi <- bk[-i, 1] - mk[-i, ] - (bk[i, 1] - mk[i, 1]) / Sigk_ii * Sigk[-i, i]
                Sigk_pi <- 0.5 * (Sigk[-i, -i] - 1 / Sigk_ii * Sigk[-i, i] %*% t(Sigk[-i, i]) +
                                  t(Sigk[-i, -i] - 1 / Sigk_ii * Sigk[-i, i] %*% t(Sigk[-i, i])))
                Sigk_dx_ii <- Sigk_dx[ , , i, i]
                Sigk_dx_ik <- Sigk_dx[ , , i, k]
                phi_ik <- dnorm(bk[i, 1], mk_i, sqrt(Sigk_ii))
                dphi_ik_dSig <- ((bk_i - mk_i)^2 / (2 * Sigk_ii^2) - 0.5 * 1 / Sigk_ii) * phi_ik
                dphi_ik_dm <- (bk_i - mk_i) / Sigk_ii * phi_ik

                ## XXXY comment out
                ## if (is.null(envir)) {
                ##     Phi_ik <- pmnorm(x = ck_pi, mean = rep(0, q - 1), varcov = Sigk_pi,
                ##                     maxpts = (q - 1) * 200)
                ## } else {
                Phi_ik <- symetric_term[k, i] / phi_ik
                ## }
                
                GPhi_ik <- GPhi(x = matrix(ck_pi, q - 1, 1), mu = matrix(0, q - 1, 1),
                                Sigma = Sigk_pi)
                HPhi_ik <- HPhi(x = matrix(ck_pi, q - 1, 1), mu = matrix(0, q - 1, 1),
                                Sigma = Sigk_pi, gradient = GPhi_ik)
                Sigk_mi <- Sigk[-i, i]
                
                for (l in 1:q) {
                    for (j in 1:d) {
                        ## B1
                        B1[l, j] <- B1[l, j] + Sigk_dx_ik[l, j] * phi_ik * Phi_ik
                        ## B2
                        B2[l, j] <- B2[l, j] + Sigk_ik *
                            (mk_dx_i[l,j] * dphi_ik_dm +
                             dphi_ik_dSig * Sigk_dx_ii[l, j]) * Phi_ik
                        ## B3
                        dck_pi <- -mk_dx[l, j, -i] +
                            (mk_dx_i[l, j] * Sigk_ii + (bk_i - mk_i) * Sigk_dx_ii[l, j]) /
                            Sigk_ii^2 * Sigk_mi -
                            (bk[i, 1] - mk_i) / Sigk_ii * Sigk_dx[l, j, -i, i]
                        SigtCross <- tcrossprod(Sigk_dx[l, j, -i, i], Sigk_mi)
                        dSigk_pi <- Sigk_dx[l, j, -i, -i] + Sigk_dx_ii[l, j] / Sigk_ii^2 *
                            tcrossprod(Sigk_mi, Sigk_mi) - Sigk_ii^-1 *
                            (SigtCross + t(SigtCross))
                        B3[l, j] <- B3[l, j] +
                            Sigk_ik * phi_ik * (crossprod(GPhi_ik, dck_pi) +
                                                0.5 * sum(HPhi_ik * dSigk_pi))
                    }
                }
            }
        }
        if (fastCompute) {
            EI.grad <- EI.grad + termB
        } else {
            EI.grad <- EI.grad + B1 + B2 + B3
        }
    }

    if (is.nan(sum(EI.grad))) EI.grad <- matrix(0, q, d)

    EI.grad <- as.vector(EI.grad)
    
    if (out_list){ 
        res <- list("objective" = thisqEI, "gradient" = EI.grad)
    } else {
        res <- thisqEI
        attr(res, "gradient") <- EI.grad
    }
    
    res 
    
}




























