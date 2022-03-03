# *****************************************************************************
## AUTHOR: Yves Deville <deville.yves@alpestat.com>
##
## Check the value of EI and its gradient
##
## *****************************************************************************

## context("EI and its derivatives")

library(numDeriv)
library(dolka)
library(testthat)

set.seed(12345)
d <- 2
X <- expand.grid(x1 = seq(from = 0, to = 1, length.out = 3),
                 x2 = seq(from = 0, to = 1, length.out = 3))
y <- apply(X, 1, branin)  

## model fit
## =========
fit1 <- km(~1, design = X, response = y, covtype = "gauss",
           control = list(pop.size = 50, trace = FALSE), parinit = c(0.5, 0.5))

if (FALSE) {
    nNew <- 100
    XNew <- matrix(runif(nNew * d), ncol = d,
                   dimnames = list(NULL, paste0("x", 1:d)))
} else {
    XNew <- expand.grid(x1 = seq(from = 0.01, to = 0.99, length.out = 20),
                        x2 = seq(from = 0.01, to = 0.99, length.out = 20))  
    nNew <- nrow(XNew)
}

## Compare with DiceOptim
## ======================
gradEnv <- new.env()
environment(EI) <- environment(EI.grad) <- gradEnv

DOEIs <- EIs <- rep(0, nNew)
DOEIGrads <- EIGrads <- NumEIGrads <- array(NaN, dim = c(nNew, d))

for (i in 1:nNew) {
    DOEIs[i] <- DiceOptim::EI(x = XNew[i, ], model = fit1, envir = gradEnv)
    DOEIGrads[i, ] <- DiceOptim::EI.grad(x = XNew[i , ], model = fit1, envir = gradEnv)
    EIs[i] <- EI(x = XNew[i, ], model = fit1, envir = gradEnv)
    EIGrads[i, ] <- EI.grad(x = XNew[i , ], model = fit1, envir = gradEnv)
    NumEIGrads[i, ] <- grad(func = EI, x = XNew[i, ], model = fit1, envir = NULL,
                            method = "simple")
}

test_that(desc = "'EI' and 'DiceOptim::EI' are consistent",
          code = expect_true(max(abs(EIs - DOEIs)) < 1e-5))

test_that(desc = "'EI.grad' and 'DiceOptim::EI' are consistent",
          code = expect_true(max(abs(EIGrads - DOEIGrads)) < 1e-5))

if (FALSE) {
    max(abs(EIGrads - NumEIGrads))
    error <- apply(EIGrads - NumEIGrads, 1, function(x) max(abs(x)))
    library(scatterplot3d)
    scatterplot3d(x = XNew[ , 1], y = XNew[ , 2], z = error,
                  type = "h", xlab = "x1", ylab = "x2", zlab = "error")
    scatterplot3d(x = XNew[ , 1], y = XNew[ , 2], z = EIs,
                  type = "h", xlab = "x1", ylab = "x2", zlab = "EI")
}
