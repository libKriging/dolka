
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package dolka

**dolka** is a transitional R package aimed at helping to use
**rlibkriging** as an alternative to **DiceKriging** when doing Bayesian
optimization with **DiceOptim** or when performing kriging inversion
with **KrigInv**.

The name is built as an acronym for *DiceOptim with the LibKriging
Alternative*.

The idea is to rely on R methods, be they S3 or S4, to allow an easy
switch from an object with class `"km"` to an other object with a
different class for instance `"KM"`in the **rlibriging** R package. The
key methods to be enhanced are the `predict` and `update` methods since
most of the tasks in Bayes optimization and kriging inversion rely on
these methods. The method `simulate` is also important. In a first step,
these methods will be re-implemented in **dolka** based on existing
source codes from the packages cited, in view of overloading the
existing methods for the `"km"` class. Then, we will suggest some light
refactoring to the authors of **DiceOptim**, **KrigInv** and maybe to
those of other packages using **DiceKriging**, if they want their
package to allow the use of both **DiceKriging** and **libKriging** via
its R interface **rlibkriging**.

Of major importance in Bayesian optimisation, the `predict` method
should be able to provide the derivative of the kriging mean and the
kriging covariance w.r.t. one or several “new” design points. This is
made possible by using the `deriv` argument in the `predict`method for
the class `"km"`provided by **dolka**, which overloads the existing
method when **dolka** is attached.

# News

-   The `genoud_cache` and `optim_cache` functions are wrapers for the
    optimization function `genoud` from the **rgenoud** package and
    `optim` from **stats**. In both cases the function to be optimized
    and its gradient are provided as *one single function* returning a
    named list with the elements `objective` and `gradient` as is the
    rule in the **nloptr** package. Compared to using two separated
    functions, this usually allows some savings due to the use of
    auxiliary variables. To stick to the use of `genoud` and `optim` the
    two required functions `fn`and `gr` are created “on the fly” in the
    wrapers, the gradient being “cached” in an environment.

-   The Bayesian optimization criteria to be optimized at each step
    (Expected Improvement, …) are and will be provided in versions
    `_with_grad` returning the list with objective and the gradient as
    described above. The computation of the gradient is made optional
    thanks to the `deriv` logical argument, and the functions can return
    simply the numeric value of the objective, possibly with a
    `gradient` attribute. This versatility allows the use of different
    optimization functions e.g. taken from **nloptr**. For now the
    Expected Improvement `EI` and the multi-points a.k.a *q-point*
    Expected Improvement `qEI`.
