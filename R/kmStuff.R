##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Predicted values, confidence intervals
##' @param object,newdata,type see \code{\link[DiceKriging]{predict.km}}. 
##' @param se.compute,cov.compute see \code{\link[DiceKriging]{predict.km}}. 
##' @param light.return,bias.correct see \code{\link[DiceKriging]{predict.km}}. 
##' @param checkNames see \code{\link[DiceKriging]{predict.km}}. 
##'
##' @param deriv Logical. If \code{TRUE} further elements are added to
##'     the returned list, all concerning the the derivatives.
##' 
##' @param ... Not used yet.
##'
##' @import DiceKriging
##' 
##' @method predict km
##' 
##' @return A list with the elements of \code{\link[DiceKriging]{predict.km}}
##' plus the following elements that relate to the derivatives w.r.t. the input 
##' \itemize{
##' \item{\code{trend.deriv} }{Derivative of the trend component.}
##' \item{\code{mean.deriv} }{Derivative of the kriging mean.}
##' \item{\code{s2.deriv} }{Derivative of the kriging variance.}
##' }
##' @examples
##' ## a 16-points factorial design, and the corresponding response
##' d <- 2; n <- 16
##' fact.design  <- expand.grid(x1 = seq(0, 1, length = 4), x2 = seq(0, 1, length = 4))
##' branin.resp <- apply(fact.design, 1, branin)
##' ## kriging model 1 : gaussian covariance structure, no trend,
##' ##                   no nugget effect
##' myKm <- km(~1 + x1 + x2, design=fact.design, response=branin.resp, covtype="gauss")
##' ## predicting at testdata points
##' testdata <- expand.grid(x1 = s <- seq(0, 1, length = 15), x2 = s)
##' pred <- predict(myKm, newdata = testdata[10, ], type = "UK", deriv = TRUE)
##' newdata <- testdata[10, ]
##' c.newdata <- covMat1Mat2(object = myKm@covariance,
##'                          X1 = myKm@X, X2 = matrix(newdata, nrow = 1),
##'                          nugget.flag = myKm@covariance@nugget.flag)
##' covVector.dx(x = newdata, X = myKm@X, object = myKm@covariance, c = c.newdata)
##' trend.deltax(x = newdata, model = myKm)
##' 
predict.km <- function(object, newdata, type,
                       se.compute = TRUE, cov.compute = FALSE,
                       light.return = FALSE,
                       bias.correct = FALSE, checkNames = TRUE,
                       deriv = FALSE, ...) {
    ## newdata : n x d
    
    nugget.flag <- object@covariance@nugget.flag 
    
    X <- object@X
    y <- object@y
    T <- object@T
    z <- object@z
    M <- object@M
    beta <- object@trend.coef
    
    if (checkNames) {
        newdata <- checkNames(X1 = X, X2 = newdata, X1.name = "the design",
                              X2.name = "newdata")
    } else {
        newdata <- as.matrix(newdata)
        d.newdata <- ncol(newdata)
        if (!identical(d.newdata, object@d)) {
            stop("newdata must have the same numbers of columns than the experimental design")
        }
        if (!identical(colnames(newdata), colnames(X))) {
            ##  warning("column names mismatch between 'newdata' and
            ## the experimental design - the columns of 'newdata' are
            ## interpreted in the same order as the experimental
            ## design names")
            colnames(newdata) <- colnames(X)
        }
    }
    
    F.newdata <- model.matrix(object@trend.formula, data = data.frame(newdata))
    y.predict.trend <- F.newdata %*% beta
    
    c.newdata <- covMat1Mat2(object@covariance, X1 = X, X2 = newdata,
                             nugget.flag = object@covariance@nugget.flag)
    ## compute c(x) for x = newdata ; remark that for prediction (or
    ## filtering), cov(Yi, Yj)=0 even if Yi and Yj are the outputs
    ## related to the equal points xi and xj.
    
    Tinv.c.newdata <- backsolve(t(T), c.newdata, upper.tri = FALSE)
    y.predict.complement <- t(Tinv.c.newdata) %*% z
    y.predict <- y.predict.trend + y.predict.complement
    y.predict <- as.numeric(y.predict)
    
    output.list <- list()
    output.list$trend <- y.predict.trend
    output.list$mean <- y.predict
    
    if (!light.return) {
        output.list$c <- c.newdata
        output.list$Tinv.c <- Tinv.c.newdata
    } 
    
    ## A FAIRE : 
    ## REMPLACER total.sd2 par cov(Z(x),Z(x)) ou x = newdata
    ## partout dans les formules ci-dessous
    ## c'est utile dans le cas non stationnaire
    
    if ((se.compute) || (cov.compute)) {
        if (!is(object@covariance, "covUser") ) {
            total.sd2 <- object@covariance@sd2
        } else {
            m <- nrow(newdata)
            total.sd2 <- rep(NA,m)
            for(i in 1:m) {
                total.sd2[i] <- object@covariance@kernel(newdata[i, ], newdata[i, ])
            }
        }
        if (object@covariance@nugget.flag) {
            total.sd2 <- total.sd2 + object@covariance@nugget
        }
    }
    
    
    if (se.compute) {
        
        ## compute c(x)'*C^(-1)*c(x)   for x = newdata
        s2.predict.1 <- apply(Tinv.c.newdata, 2, crossprod)
        
        if (type == "SK") {
            s2.predict <- pmax(total.sd2 - s2.predict.1, 0)
            s2.predict <- as.numeric(s2.predict)
            q95 <- qnorm(0.975)
        }
        else if (type == "UK") {
            
            T.M <- chol(t(M) %*% M)   # equivalently : qrR <- qr.R(qr(M))
            ## 's2.predict.mat' is a matrix with n rows and p columns 
            s2.predict.mat <- backsolve(t(T.M),
                                        t(F.newdata - t(Tinv.c.newdata) %*% M),
                                        upper.tri = FALSE)
            
            s2.predict.2 <- apply(s2.predict.mat, 2, crossprod)
            s2.predict <- pmax(total.sd2 - s2.predict.1 + s2.predict.2, 0)
            s2.predict <- as.numeric(s2.predict)
            if (bias.correct) {
                s2.predict <- s2.predict * object@n/(object@n - object@p)
            }
            q95 <- qt(0.975, object@n - object@p)
        }
        
        s.predict <- sqrt(s2.predict)
        lower95 <- y.predict - q95 * s.predict
        upper95 <- y.predict + q95 * s.predict
        
        output.list$sd <- s.predict
        output.list$lower95 <- lower95
        output.list$upper95 <- upper95
    }
    
    if (cov.compute) {		
        
        C.newdata <- covMatrix(object@covariance, newdata)[[1]]
        cond.cov <- C.newdata - crossprod(Tinv.c.newdata)
        
        if (type == "UK") {	
            T.M <- chol(t(M)%*%M)   # equivalently : qrR <- qr.R(qr(M))
            s2.predict.mat <- backsolve(t(T.M),
                                        t(F.newdata - t(Tinv.c.newdata) %*% M),
                                        upper.tri = FALSE)
            cond.cov <- cond.cov + crossprod(s2.predict.mat)
            if (bias.correct) cond.cov <- cond.cov * object@n / (object@n - object@p)
        }
        
        output.list$cov <- cond.cov
        
    }

    ## =========================================================================
    ## Compute the derivative of the kriging mean the kriging variance and
    ## the kriging standard deviation.
    ## =========================================================================
    
    if (deriv) {
        
        newdata <- as.numeric(newdata)
        
        if (!is.null(dnd <- dim(newdata))) {
            if ((dnd[1] == 1) && (dnd[2] == object@d)) {
                newdata <- drop(newdata)
            } else {
                stop("for now 'deriv = TRUE' can only be used when ",
                     "'newdata' is a numeric vector with lenght 'd' ",
                     "or a matrix with dimension c(1, d)")
            }
        }

        ## 'derFnew' is a vector with length d
        F.newdata.deriv <- trend.deltax(x = newdata, model = object)
        
        ## 'c.deriv' and 'cStar.deriv' are  matrices with dimension c(n, d)
        c.deriv <- covVector.dx(x = newdata, X = X,
                             object = object@covariance, c = c.newdata)
        cStar.deriv <- backsolve(t(T), c.deriv, upper.tri = FALSE)

        ## Gradient of the kriging mean
        output.list$trend.deriv <- crossprod(F.newdata.deriv, beta)
        output.list$mean.deriv <- output.list$trend.deriv +
            crossprod(cStar.deriv, z) 

        ## 'Tinv.c.newdata' is something like 'cStar'
        sd2.deriv <-  -2 * crossprod(Tinv.c.newdata, cStar.deriv)
        
        if (type == "UK") {
            
            ## remind of the defintiion of 's2.predict.mat'
            ## s2.predict.mat <- backsolve(t(T.M),
            ##                            t(F.newdata - t(Tinv.c.newdata) %*% M),
            ##                            upper.tri = FALSE)
            
            s2.predict.der <- backsolve(t(T.M),
                                        F.newdata.deriv - t(M) %*% cStar.deriv,
                                        upper.tri = FALSE)
            
            sd2.deriv <-  sd2.deriv + 2 * crossprod(s2.predict.mat, s2.predict.der)
        }
        
        output.list$sd2.deriv <- sd2.deriv
        output.list$sd.deriv <- sd2.deriv / (2 * output.list$sd)

    }
    
    return(output.list)
    
}

##' @export
##' @method predict km
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

if (FALSE) {
    
    krigeMean <- function(x, model, type = "UK", ...) {
        newdata <- data.frame(matrix(as.numeric(x), nrow = 1))
        colnames(newdata) <- colnames(model@X)
        
        p <- predict(object = model,
                     newdata = newdata, type = type)$mean
        p
    }
    
    krigeVar <- function(x, model, type = "UK", ...) {
        newdata <- data.frame(matrix(as.numeric(x), nrow = 1))
        colnames(newdata) <- colnames(model@X)
        
        p <- predict(object = model,
                     newdata = newdata, type = type)$sd^2
        p
    }
    krigeSd <- function(x, model, type = "UK", ...) {
        newdata <- data.frame(matrix(as.numeric(x), nrow = 1))
        colnames(newdata) <- colnames(model@X)
        
        p <- predict(object = model,
                     newdata = newdata, type = type)$sd
        p
    }
    
    library(numDeriv)
    x <- runif(2)
    type <- "SK"
    krigeMean(x, mod = myKm, deriv = TRUE)
    gMean <- grad(krigeMean, x = x, model = myKm, type = type)
    gVar <- grad(krigeVar, x = x, model = myKm, type = type)
    gSd <- grad(krigeSd, x = x, model = myKm, type = type)
    p <- predict(myKm,
                 newdata = matrix(x, ncol = 2, dimnames = list(NULL, c("x1", "x2"))),
                 type = type, deriv = TRUE)

}
