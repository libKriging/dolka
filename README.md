<!-- README.md is generated from README.Rmd. Please edit that file -->

dolka package
=============

**dolka** is a transitional package aimed at helping to use
**rlibkriging** as an alternative to **DiceKriging** when doing Bayesian
optimization with **DiceOptim** or when performing kriging inversion
with **KrigeInv**.

The name is built as an acronym for *DiceOptim with the LibKriging
Alternative*.

The idea is to rely on R methods, be they S3 or S4, to allow an easy
switch from an object with class `"km"` to an other object with a
different class for instance `"KM"`in the **rlibriging** R package. The
key methods to be enhanced are the `predict` and `update` methods since
most of the tasks in Bayes optimization and kriging inversion rely on
these methods. In a first step these methods will be re-implemented in
**dolka** (based on existing source codes from the packages) in view of
overloading the existing methods for the `"km"` class. Then, we will
suggest some light refactoring to the authors of **DiceOptim**,
**KrigeInv** and maybe to those of other packages using **DiceKriging**.

Of major importance in Bayesian optimisation, the `predict` method
should be able to provide the derivative of the kriging mean and the
kriging covariance w.r.t. one or several “new” design points. This is
made possible by using the `deriv` argument.