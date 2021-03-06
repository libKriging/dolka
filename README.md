
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

# INSTALLATION

For now, we do not provide precompiled binary versions (`.zip` for
Windows). These may be provided soon in the
[**Releases**](https://github.com/libKriging/dolka/releases) tab.

## Simple installation using **remotes** or **devtools**

Since **dolka** does not for now need compilation (it doe not embed C or
Cpp program files), the **remotes** package can be used.

    remotes::install_github("libKriging/dolka")

If you have the
[**Rtools**](https://cran.r-project.org/bin/windows/Rtools/) installed,
you can also rely on **devtools** package

    devtools::install_github("libKriging/dolka")

In both cases many options are available e.g. to select a specific
branch or commit. See the manuals of **remotes** or **devtools** for
more information.

## Clone, build and install

### Cloning the repository

If you do not have yet a local `dolka` repository, use `git clone` to
clone the `dolka` repository

``` bash
git clone https://github.com/liKriging/dolka
```

This will create a `dolka` sub-directory of the current directory,
i.e. the directory from which the git command was issued.

### Installation on Unix and MacOs systems

With these systems you can install a package from its source. Move to
the parent directory of your cloned repository and use the following
command from a terminal to create a tarball source file

``` bash
R CMD build dolka
```

This will produce a so-called *source tarball* file `dolka_x.y.z` where
`x`, `y` and `z` stand for the major, minor and patch version numbers.
Then you can install from a command line

``` bash
R CMD INSTALL dolka_x.y.z.tar.gz
```

Of course if you are using the **RStudio** IDE you can install the
source tarball using the menu **Tools/Install Packages…** Note that a
*source tarball* file is called a *Package Archive File* in **RStudio**.

### Installation on Windows systems

Provided that you have the
[**Rtools**](https://cran.r-project.org/bin/windows/Rtools/) installed,
you can use the same commands as Linux and MacOS sytems.

Moreover you can create a so-called **precompiled binary** file using

``` bash
R CMD INSTALL --build dolka_x.y.z.tar.gz
```

This will create a `dolka.zip` file that can be used on a Windows
plateform which may not be equipped with the **Rtools**. For instance,
with **RStudio** you can use the menu `Tools/Install Packages` and
select `Install from:` to install a precompiled `.zip` file.
