# *****************************************************************************
## AUTHOR: Yves Deville <deville.yves@alpestat.com>
##
## Test the function 'qEI_with_grad'.
##
## o The precision used here is quite loose. The relative precision
## seems to be much better but it is not used here because the
## gradient hase often elements close to zero.
##
## o When 'nNew' a.k.a 'q' increases, the precision on the gradient
## becomes poor. With 'nNew = 6' the test fails.
## 
## *****************************************************************************

library(testthat)
# context("qEI_with_grad")

library(dolka)
library(numDeriv)
set.seed(123)
prec <- 1e-1

Funs <- c("branin", "goldsteinPrice", "camelback")
Types <- c("UK", "SK")

d <- 2; n <- 9; nNew <- 2
X <- expand.grid(x1 = seq(0, 1, length = 3),
                 x2 = seq(0, 1, length = 3))

XNew <- matrix(runif(nNew * d), nrow = nNew, ncol = d)

for (i in seq_along(Funs)) {
    fun <- get(Funs[i], mode = "function")
    y <- apply(X, 1, fun)
    fit <- km(~1, design = X, response = y, 
              covtype = "gauss", control = list(pop.size = 50, trace = FALSE),
              parinit = c(0.5, 0.5))
    
    ## function of a vector to check the gradient
    qEIFun <- function(xNew, type) {
        dim(xNew) <- c(nNew, d)
        qEI_with_grad(xNew, model = fit,
                      deriv = FALSE,
                      type = type,
                      out_list = FALSE)
    }
    for (type in Types) {

        gradNum <- grad(func = qEIFun, x = as.vector(XNew), type = type)
        dim(gradNum) <- c(nNew, d)

        res <- qEI_with_grad(XNew, model = fit, deriv = TRUE,
                             type = type,
                             out_list = TRUE)
        
        ## res <- DiceOptim::qEI.grad(XNew, model = fit, type = type)
        
        errDer <- abs(gradNum - res$gradient)
        ## print(errDer)
        Den <- (abs(gradNum) + abs(res$gradient)) / 2.0
        ## errDer <- errDer / Den
        test_that(desc = sprintf(paste0("Derivative of qEI",
                                        " fun =\"%s\", type = \"%s\""),
                                 Funs[i], type),
                  code = expect_true(max(errDer) < prec))
    }
}
