## *****************************************************************************
## AUTHOR: Yves Deville <deville.yves@alpestat.com>
##
## Test that the function 'max_EI' in the dolka package returns the
## same result as 'DiceOptim::max_EI' on a example with several
## iterations. Of course, the test could be limited to one iteration.
##
## NOTE that 'DiceOptim::max_EI' generates warnings hence so does
## 'dolka::max_EI' which refactors the code striving to keep the
## results unchanged.
## 
## *****************************************************************************

library(testthat)
context("'max_EI' consistency between dolka and DiceOptim")
library(DiceKriging)
## library(rlibkriging)

## =============================================================================
## CAUTION dolka must be attached AFTER rlibkriging or the
## predict method for 'km' will not work as expected.
## =============================================================================
library(dolka)

f <- DiceKriging::branin
d <- 2
lower <- rep(0.0, d)
upper <- rep(1.0, d)

set.seed(123)
n0 <- 9
X0 <- as.matrix(expand.grid(x1 = seq(0.0, 1.0, len = 3),
                            x2 = seq(0.0, 1.0, len = 3)))
## only for the plot...
colnames(X0) <- c("x1", "x2")
y0 <- apply(X0, MARGIN = 1, f)

model <- opt <- mEI <- newX <- newy <- list()
nIter <- 10

for (pkg in c("dolka", "DiceOptim")) {
    
    set.seed(456)

    if (FALSE) {
    ## get the 'max_EI' function of the current package 
        mEI[[pkg]] <- get("max_EI_genoud", mode = "function",
                          envir = asNamespace("DiceOptim"))
    }
    
    ## Initial 'km' model
    model[[pkg]] <- km(design = X0, response = y0)

    ## EGO iterations
    for (iter in 1:nIter) {
        if (pkg == "dolka") {
            opt[[pkg]] <-
                suppressWarnings(dolka::max_EI_genoud(model = model[[pkg]],
                                              lower = lower, upper = upper))
            ## genoud_args = list(print.level = 0))
        } else {
            opt[[pkg]] <-
                suppressWarnings(DiceOptim::max_EI(model = model[[pkg]],
                                                   lower = lower, upper = upper,
                                                   control = list(print.level = 0)))
        }
        ## opt[[pkg]] <- mEI[[pkg]](model = model[[pkg]],
        ##     lower = lower, upper = upper)
        newX[[pkg]] <- opt[[pkg]]$par
        newy[[pkg]] <- f(newX[[pkg]]) 
        model[[pkg]] <- update(model[[pkg]],
                               newX = newX[[pkg]], newy = newy[[pkg]],
                               newnoise.var = 1e-8)
    }
}

test_that(desc = paste0("'dolka::max_EI_genoud' leads to the same results as\n",
                        "'DiceOptim::max_EI_genoud' on several EGO iterations\n",
                        "to minimize the branin function"),
          code = expect_true(all.equal(model[["DiceOptim"]]@X,
                                       model[["dolka"]]@X)))
