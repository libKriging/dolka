
## *****************************************************************************
##' Overload the \code{predict} method of the class \code{"km"} of the
##' \pkg{DiceKriging} package in order make possible the computation
##' of derivatives.
##'
##' When \code{deriv} is \code{TRUE} "Jacobian" arrays are returned
##' with the following rule. For a function
##' \eqn{\mathbf{F}(\mathbf{X})}{F(X)} with a
##' \eqn{q \times d}{c(q, d)}
##' matrix argument and a \eqn{n \times m}{c(n, m)} matrix value,
##' the Jacobian array has dimension
##' \eqn{(n \times m) \times (q \times d)}{c(n, m, q, d)} and element
##' \deqn{D \mathbf{F}(\mathbf{X})[i, j, k, \ell] = \frac{\partial F_{i,j}}{\partial X_{k,\ell}}}{DF(X)[i, j, k, ell] = dF[i, j] / dX[k, ell]}. This rule is compatible with the
##' R arrays rule: if the function is considered as a function
##' of a vector argument \code{as.vector(X)} with the vector value
##' \code{as.vector(F(X))}, then the 
##' 
##' @title Predicted values, confidence intervals
##'
##' @param object,newdata,type see
##'     \code{\link[DiceKriging]{predict.km}}.
##'
##' @param se.compute,cov.compute see
##'     \code{\link[DiceKriging]{predict.km}}.
##'
##' @param light.return,bias.correct see
##'     \code{\link[DiceKriging]{predict.km}}.
##' 
##' @param checkNames see \code{\link[DiceKriging]{predict.km}}.
##'
##' @param deriv Logical. If \code{TRUE} further elements are added to
##'     the returned list, all concerning the the derivatives.
##' 
##' @param ... Not used yet.
##'
##' @import DiceKriging
##' @export
##'
##' @section Caution: XXXY remettre "method predict km" dans le roxygen
##' 
##' @return A list with the elements of \code{\link[DiceKriging]{predict.km}}
##' plus the following elements that relate to the derivatives w.r.t. the input 
##' \itemize{
##' \item{\code{trend.deriv} }{
##' Derivative of the trend component. This is an array with dimension
##'      \eqn{[n, n, d]}{c(n, n, d)}.
##' }
##' \item{\code{mean.deriv}, \code{s2.deriv} }{
##' Derivatives of the kriging mean and kriging variance. These are
##'      arrays with dimension \eqn{[n, n, d]}{c(n, n, d)}.
##' }
##' \item{\code{cov.deriv} }{
##'    Derivative of the kriging covariance. This is a
##'     four-dimensional array with dimension \eqn{[n, n, n, d]}{c(n, n, n, d)}.
##'    }
##' }
##' @examples
##' ## a 16-points factorial design, and the corresponding response
##' d <- 2; n <- 16
##' X  <- expand.grid(x1 = seq(0, 1, length = 4), x2 = seq(0, 1, length = 4))
##' y <- apply(X, MARGIN = 1, FUN = branin)
##' ## kriging model 1 : gaussian covariance structure, no trend,
##' ##                   no nugget effect
##' myKm <- km(~1 + x1 + x2, design = X, response = y, covtype = "gauss")
##' ## predicting at new points
##' XNew <- expand.grid(x1 = s <- seq(0, 1, length = 15), x2 = s)
##' pred <- predict(myKm, newdata = XNew[10, ], type = "UK", deriv = TRUE)
##' newdata <- XNew[10, ]
##' c.newdata <- covMat1Mat2(object = myKm@covariance,
##'                          X1 = myKm@X, X2 = matrix(newdata, nrow = 1),
##'                          nugget.flag = myKm@covariance@nugget.flag)
##' covVector.dx(x = newdata, X = myKm@X,
##'              object = myKm@covariance,
##'              c = c.newdata)
##' trend.deltax(x = newdata, model = myKm)
##' 
predict.km <- function(object, newdata, type,
                       se.compute = TRUE,
                       cov.compute = FALSE,
                       light.return = FALSE,
                       bias.correct = FALSE,
                       checkNames = TRUE,
                       deriv = FALSE,
                       ...) {
    
    ## newdata : n x d
    
    ## extract some slots and rename
    nugget.flag <- object@covariance@nugget.flag 
    X <- object@X
    y <- object@y
    L <- t(object@T)
    zStar <- object@z
    FStar <- object@M
    betaHat <- object@trend.coef
    
    ## =========================================================================
    ## Checks and coercion
    ## =========================================================================

    if (checkNames) {
        XNew <- checkNames(X1 = X, X2 = newdata, X1.name = "the design",
                           X2.name = "newdata")
    } else {
        XNew <- as.matrix(newdata)
        dNew <- ncol(XNew)
        if (!identical(dNew, object@d)) {
            stop("'newdata' must have the same numbers of columns as the ",
                 "experimental design used in 'object'")
        }
        if (!identical(colnames(XNew), colnames(X))) {
            colnames(XNew) <- colnames(X)
        }
    }
    
    nNew <- nrow(XNew)
    FNew <- model.matrix(object@trend.formula, data = data.frame(XNew))
    muNewHat <- FNew %*% betaHat
  
    ## =========================================================================
    ## compute 'KOldNew := K(XOld, XNew)'  with dim c(n, nNew) and then
    ## KOldNewStar := L^{-1} %*% KOldNew.
    ## =========================================================================
    
    KOldNew <- covMat1Mat2(object@covariance, X1 = X, X2 = XNew,
                           nugget.flag = object@covariance@nugget.flag)
    KOldNewStar <- backsolve(L, KOldNew, upper.tri = FALSE)
    
    meanNewHat <- muNewHat + t(KOldNewStar) %*% zStar
    meanNewHat <- as.numeric(meanNewHat)
    
    ## ========================================================================= 
    ## Prepare an output list 'OL' stands for "output.list"
    ## =========================================================================
    
    OL <- list()
    OL$trend <- muNewHat
    OL$mean <- meanNewHat
    
    if (!light.return) {
        OL$c <- as.vector(KOldNew)
        OL$Tinv.c <- as.vector(KOldNewStar)
    } 

    
    ## =========================================================================
    ## 'sd2Total' is the variance
    ## =========================================================================
    
    if ((se.compute) || (cov.compute)) {
        if (!is(object@covariance, "covUser") ) {
            sd2Total <- object@covariance@sd2
        } else {
            s2Total <- rep(NA, nNew)
            for(i in 1:nNew) {
                sd2Total[i] <- object@covariance@kernel(XNew[i, ], XNew[i, ])
            }
        }
        if (object@covariance@nugget.flag) {
            sd2Total <- sd2Total + object@covariance@nugget
        }

        if (type == "UK") {
            ## =============================================================
            ## Compute 'RStar := qr.R(qr(FStar))',
            ## Then the transpose of 'ENewDag := ENew %*% RStar^{-1}'
            ## with 'ENew := FNew - FNewHat'
            ## =============================================================
        
            RStar <- chol(crossprod(FStar)) 
            ENewDagT <- backsolve(t(RStar),
                                  t(FNew - t(KOldNewStar) %*% FStar),
                                  upper.tri = FALSE)
            
            OL$ENewDagT <- ENewDagT
            
        }
    }
    
    if (se.compute) {

        ## =====================================================================
        ## Compute K(XNew, XOld) %*% K(XOld, XOld)^(-1) %*% K(XOld, XNew)
        ## This is diagCrossprod(A) with A := KOldNewStar with
        ## dimension c(n, nNew)
        ## =====================================================================
        
        ## s2Part1 <- apply(KOldNewStar, 2, crossprod)         ## less efficient
        s2Pred1 <- dcrossprod(KOldNewStar)
        
        if (type == "SK") {

            s2Pred <- pmax(sd2Total - s2Pred1, 0)
            s2Pred <- as.numeric(s2Pred)
            q95 <- qnorm(0.975)
            
        } else if (type == "UK") {

            ## =================================================================
            ## Compute 'RStar := qr.R(qr(FStar))',
            ## Then the transpose of 'ENewDag := ENew %*% RStar^{-1}' with
            ## 'ENew := FNew - FNewHat'
            ## =================================================================
            
            ## RStar <- chol(crossprod(FStar)) 
            ## ENewDagT <- backsolve(t(RStar),
            ##                      t(FNew - t(KOldNewStar) %*% FStar),
            ##                      upper.tri = FALSE)
            
            ## =================================================================
            ## Compute 'diagCrossprod(A) = diag(t(A) %*% A)' with
            ## 'A := ENewDagT with dimension c(p, nNew) 
            ## =================================================================
            
            s2Pred2 <- dcrossprod(ENewDagT)
            ## s2Pred2 <- apply(ENewDagT, 2, crossprod)    ## less efficient
            
            s2Pred <- pmax(sd2Total - s2Pred1 + s2Pred2, 0)
            s2Pred <- as.numeric(s2Pred)

            if (bias.correct) {
                s2Pred <- s2Pred * object@n / (object@n - object@p)
            }
            q95 <- qt(0.975, object@n - object@p)
        }
        
        sPred <- sqrt(s2Pred)
        lower95 <- meanNewHat - q95 * sPred
        upper95 <- meanNewHat + q95 * sPred
        
        OL$sd2 <- s2Pred
        OL$sd <- sPred
        OL$lower95 <- lower95
        OL$upper95 <- upper95
    }
    
    if (cov.compute) {		

        ## =====================================================================
        ## Note that 'Sigma' may not be numerically positive definite
        ## Also remind that the parameters of 'object@covariance' are
        ## (constrained) ML estimates, including for the variance. So
        ## the variance uses the denominator 'n', not 'n - p'.
        ## =====================================================================
        
        KNewNew <- covMatrix(object@covariance, XNew)[["C"]]
        Sigma <- KNewNew - crossprod(KOldNewStar)

        ## XXX
        OL$Sigma1 <- Sigma
        
        if (type == "UK") {	
            Sigma <- Sigma + crossprod(ENewDagT)
            if (bias.correct) {
                Sigma <- Sigma * object@n / (object@n - object@p)
            }
        }
        ## XXX
        OL$Sigma2 <- Sigma
        OL$bias.correct <- bias.correct
        OL$cov <- Sigma        
    }

    ## =========================================================================
    ## Compute the derivative of the kriging mean the kriging variance and
    ## the kriging standard deviation.
    ## =========================================================================
    
    if (deriv) {
        
        ## =====================================================================
        ## Prepare some 'empty' arrays for output
        ## =====================================================================
        
        muNewHatDer <- meanNewHatDer <- s2Der <-
            array(0.0, dim = c(nNew, nNew, object@d),
                  dimnames = list(NULL, NULL, colnames(object@X)))
        
        ## =====================================================================
        ## Prepare more auxiliary variables.
        ##  ________________________________________________________________
        ## |     name       |         formula         |      dim            |
        ## | -------------------------------------------------------------- |
        ## | KOldNewStarDer | L^{-1} %*% D(KOldNew)   | c(n, nNew, nNew, d) |
        ## | ENewDagDer     | D(ENew) %*% RStar^{-1}  | c(nNew, nNew, p, d) |
        ## |________________________________________________________________|
        ##
        ## where 'D' stands for the derivative. Remind that 'KOldNew' and
        ## 'ENewDagT' have been stored previously. Storing the derivatives
        ## is only useful when the derivative of the covariance is needed,
        ## but it costs nothing but space. 
        ## =====================================================================
        
        KOldNewStarDer <-
            array(0.0, dim = c(object@n, nNew, nNew, object@d),
                  dimnames = list(NULL, NULL, NULL, colnames(object@X)))
        ENewDagDer <-
            array(0.0, dim = c(nNew, nNew, object@p, object@d),
                  dimnames = list(NULL, NULL, NULL, colnames(object@X)))
        
        for (i in 1:nNew) {
            
            ## =================================================================
            ## 'FNewiDer' and 'KOldNewiDer' are matrices with
            ## dimension c(nNew, d)
            ## =================================================================
            
            FNewiDer <- trend.deltax(x = XNew[i, ], model = object)
            KOldNewi <- as.vector(KOldNew[ , i])
            KOldNewiDer <- covVector.dx(x = as.vector(XNew[i, ]),
                                        X = X,
                                        object = object@covariance,
                                        c = KOldNewi)
            
            KOldNewStarDer[ , i, i, ] <-
                backsolve(L, KOldNewiDer, upper.tri = FALSE)
            
            ## Gradient of the kriging trend and mean
            muNewHatDer[i, i, ] <- crossprod(FNewiDer, betaHat)
            ## dim in product c(d, n) and NULL(length d)
            meanNewHatDer[i, i, ] <- muNewHatDer[i, i, ] +
                crossprod(KOldNewStarDer[ , i, i,  ], zStar) 

            ## dim in product c(d, n) and NULL(length n)
            s2Der[i, i,  ] <-
                - 2 * crossprod(KOldNewStarDer[ , i, i, ],
                                drop(KOldNewStar[ , i]))
            
            ## dim in product c(d, n) and c(n, p)
            
            if (type == "UK") {
                ENewDagDer[i, i, , ] <-
                    backsolve(t(RStar),
                              FNewiDer -
                              t(crossprod(KOldNewStarDer[ , i, i, ], FStar)),
                              upper.tri = FALSE)
                ## dim in product NULL (length p) and c(p, d) because of 'drop'
                s2Der[i, i, ] <- s2Der[i, i, ] + 2 * drop(ENewDagT[ , i]) %*%
                    drop(ENewDagDer[i, i, , ])
            }
            
        }
        
        OL$trend.deriv <- muNewHatDer
        OL$mean.deriv <- meanNewHatDer
        
        OL$sd2.deriv <- s2Der
        OL$sd.deriv <- s2Der / (2 * OL$sd)
        
        if (cov.compute) {
            
            covDeriv <-
                array(0.0, dim = c(nNew, nNew, nNew, object@d),
                      dimnames = list(NULL, NULL, NULL, colnames(object@X)))
            
            for (i in 1:nNew) {
                for (j in 1:nNew) {
                    covDeriv[i, j, i, ] <-
                        covVector.dx(x = as.vector(XNew[i, ]),
                                     X = XNew[j, , drop = FALSE],
                                     object = object@covariance,
                                     c = KNewNew[i, j])
                    
                    ## dim in product c(d, n) and NULL(length n)
                    covDeriv[i, j, i, ]  <- covDeriv[i, j, i, ] -
                        crossprod(drop(KOldNewStarDer[ , i, i, ]),
                                  drop(KOldNewStar[ , j]))
                    
                    if (type == "UK") {
                        ## dim in product NULL (length p) and c(p, d)
                        covDeriv[i, j, i, ] <- covDeriv[i, j, i, ] +
                            drop(ENewDagT[ , j]) %*% drop(ENewDagDer[i, i, , ])
                    }
                    
                }
            }
            for (i in 1:nNew) {
                for (ell in 1:object@d) {
                    covDeriv[ , , i, ell] <-  covDeriv[ , , i, ell] +
                        t(covDeriv[ , , i, ell])
                }
            }
            OL$cov.deriv <- covDeriv
        }
    }
    
    return(OL)
    
}

##' @export
##' @method predict km
if (FALSE) {
    setMethod("predict", "km", 
              function(object, newdata, type, se.compute = TRUE,
                       cov.compute = FALSE, light.return = FALSE, bias.correct = FALSE,
                       checkNames = TRUE, ...) {
                  predict.km(object = object, newdata = newdata, type = type,
                             se.compute = se.compute, cov.compute = cov.compute,
                       light.return = light.return,
                       bias.correct = bias.correct, checkNames = checkNames, ...)
              }
              )
}

