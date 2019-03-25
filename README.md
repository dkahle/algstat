<!-- README.md is generated from README.Rmd. Please edit that file -->

algstat
=======

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/algstat)](https://cran.r-project.org/package=algstat)
[![Travis build
status](https://travis-ci.org/dkahle/algstat.svg?branch=master)](https://travis-ci.org/dkahle/algstat)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/dkahle/algstat?branch=master&svg=true)](https://ci.appveyor.com/project/dkahle/algstat)
<!-- badges: end -->

**algstat** is a collection of tools to help you use algebraic
statistical methods in R. It is intended to be an end user package for
that purpose, and consequently depends on other packages that make
connections to math software to do key computations. Currently, [the
**latte** package](https://github.com/dkahle/latte.git) is used to
connect R to [LattE](https://www.math.ucdavis.edu/~latte/) with
[4ti2](http://www.4ti2.de) for lattice problems and the computation of
Markov bases, [the **m2r**
package](https://github.com/coneill-math/m2r.git) connects R to
[Macaulay2](http://www.math.uiuc.edu/Macaulay2/) for algebraic
computations, and [the **bertini**
package](https://github.com/dkahle/bertini) connects R to
[Bertini](https://bertini.nd.edu), which is used to numerically solve
systems of polynomial equations. These are at varying stages of
development. **m2r** can be used for that purpose to some extent
currently by using Macaulay2’s connection to
[PHCPack](http://homepages.math.uic.edu/~jan/download.html), so be sure
to look there for now if that’s what you’re trying to do.

If you have something you’re wanting implemented here, please file an
issue!

Exact inference with log-linear models
======================================

*Note: this section assumes you have
[**latte**](https://github.com/dkahle/latte.git) installed and working
on your machine.*

One of the most well-developed parts of the package allows users to
perform (conditional) exact tests for log-linear models. There are
several great references on the math behind this, such as [Diaconis and
Sturmfels’ original
paper](http://projecteuclid.org/euclid.aos/1030563990), the [Lectures on
Algebraic
Statistics](http://smile.amazon.com/Lectures-Algebraic-Statistics-Oberwolfach-Seminars/dp/3764389044/ref=sr_1_1?ie=UTF8&qid=1430536908&sr=8-1&keywords=lectures+on+algebraic+statistics),
and [Markov Bases in Algebraic
Statistics](http://smile.amazon.com/Markov-Bases-Algebraic-Statistics-Springer/dp/1461437180/ref=sr_1_fkmr0_1?ie=UTF8&qid=1430536933&sr=8-1-fkmr0&keywords=aoki%2C+hada%2C+and+takemura),
so we’ll keep the technical discussion to a minimum.

### Fisher’s exact test

We’ll begin by doing Fisher’s exact test on a built-in dataset called
politics.

``` r
library("algstat")
# Loading required package: mpoly
# Loading required package: latte
#   Please cite latte! See citation("latte") for details.
# Loading required package: bertini
#   Please cite bertini! See citation("bertini") for details.
# Loading required package: m2r
#   Please cite m2r! See citation("m2r") for details.
#   M2 found in /Applications/Macaulay2-1.10/bin
# Please cite algstat! See citation("algstat") for details.
data(politics)
politics
#            Party
# Personality Democrat Republican
#   Introvert        3          7
#   Extrovert        6          4
```

Here’s how you typically do Fisher’s exact test in R:

``` r
fisher.test(politics)
# 
#   Fisher's Exact Test for Count Data
# 
# data:  politics
# p-value = 0.3698
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.03005364 2.46429183
# sample estimates:
# odds ratio 
#   0.305415
```

Since the independence model is log-linear, this exact same procedure
can be done with **algstat**. The go-to function here is `loglinear()`:

``` r
loglinear(~ Personality + Party, data = politics)
# Computing Markov moves (4ti2)... done.
# Running chain (C++)... done.
# Call:
# loglinear(model = ~Personality + Party, data = politics)
# 
# Fitting method:
# Iterative proportional fitting (with stats::loglin)
# 
# MCMC details:
# N = 10000 samples (after thinning), burn in = 1000, thinning = 10
# 
#                  stat p.value     se mid.p.value
#       P(table)         0.3739 0.0048      0.2219
#    Pearson X^2 1.8182  0.3739 0.0048      0.2219
# Likelihood G^2 1.848   0.3739 0.0048      0.2219
#  Freeman-Tukey 1.8749  0.3739 0.0048      0.2219
#   Cressie-Read 1.8247  0.3739 0.0048      0.2219
#     Neyman X^2 2.0089  0.3739 0.0048      0.2962
```

Exact inference in algebraic statistics is done using MCMC to sample
from the conditional distribution of the data given its sufficient
statistics under the model. Consequently, the p-values estimated are
only determined up to Monte Carlo error. The standard p-value is given
under the column `p.value` in the row labeled `P(samp)`.

The asymptotic test of independence analogous to Fisher’s exact test
(which does not condition on the marginals being known) can be done in
either of two ways.

The first way uses the `loglin()` function from the **stats** package.
It outputs the likelihood ratio statistic (`Likelihood G^2` in the
output above) and Pearson’s chi-squared statistic (`Pearson X^2` above),
but you have to calculate the p-value yourself.

``` r
(loglinMod <- stats::loglin(politics, list(1, 2)))
# 2 iterations: deviation 0
# $lrt
# [1] 1.848033
# 
# $pearson
# [1] 1.818182
# 
# $df
# [1] 1
# 
# $margin
# $margin[[1]]
# [1] "Personality"
# 
# $margin[[2]]
# [1] "Party"
pchisq(loglinMod$pearson, df = 1, lower.tail = FALSE)
# [1] 0.1775299
```

The second way is the `loglm()` function in the **MASS** package, which
is a nice wrapper of `loglin()` (in fact, **algstat**’s `loglinear()`
function uses the IPF implementation from `loglin()`, although it
doesn’t need to). It’s syntax looks identical to `loglinear()`’s above:

``` r
MASS::loglm(~ Personality + Party, data = politics)
# Call:
# MASS::loglm(formula = ~Personality + Party, data = politics)
# 
# Statistics:
#                       X^2 df  P(> X^2)
# Likelihood Ratio 1.848033  1 0.1740123
# Pearson          1.818182  1 0.1775299
```

### Fisher’s exact test on RxC tables

Doing Fisher’s exact test on larger problems is a significantly more
complicated problem. The documentation for `fisher.test()` illustrates
how it can be used on RxC tables in general, not just on 2x2 tables.
Here’s an example from its documentation drawn from Agresti (2002,
p.57):

``` r
Job <- matrix(
  c(1,2,1,0, 3,3,6,1, 10,10,14,9, 6,7,12,11), nrow = 4, ncol = 4,
  dimnames = list(
    "income" = c("< 15k", "15-25k", "25-40k", "> 40k"),
    "satisfaction" = c("VeryD", "LittleD", "ModerateS", "VeryS")
  )
)

Job
#         satisfaction
# income   VeryD LittleD ModerateS VeryS
#   < 15k      1       3        10     6
#   15-25k     2       3        10     7
#   25-40k     1       6        14    12
#   > 40k      0       1         9    11

fisher.test(Job)
# 
#   Fisher's Exact Test for Count Data
# 
# data:  Job
# p-value = 0.7827
# alternative hypothesis: two.sided
```

Here’s the **algstat** counterpart:

``` r
loglinear(~ income + satisfaction, data = Job)
# Care ought be taken with tables with sampling zeros to ensure the MLE exists.
# Computing Markov moves (4ti2)... done.
# Running chain (C++)... done.
# Call:
# loglinear(model = ~income + satisfaction, data = Job)
# 
# Fitting method:
# Iterative proportional fitting (with stats::loglin)
# 
# MCMC details:
# N = 10000 samples (after thinning), burn in = 1000, thinning = 10
# 
#                  stat p.value     se mid.p.value
#       P(table)         0.789  0.0041      0.7884
#    Pearson X^2 5.9655  0.782  0.0041      0.782 
# Likelihood G^2 6.7641  0.7849 0.0041      0.7849
#  Freeman-Tukey 8.6189  0.7813 0.0041      0.7813
#   Cressie-Read 6.0752  0.7837 0.0041      0.7837
#     Neyman X^2 6.2442  0.6087 0.0049      0.6087
```

The asymptotic test can be performed as well. The chi-square
approximation is actually very good here:

``` r
MASS::loglm(~ income + satisfaction, data = Job)
# Call:
# MASS::loglm(formula = ~income + satisfaction, data = Job)
# 
# Statistics:
#                       X^2 df  P(> X^2)
# Likelihood Ratio 6.764053  9 0.6616696
# Pearson          5.965515  9 0.7433647
```

### Fisher’s exact test on multi-way tables

`fisher.test()` does not work with multi-way tables and is prone to
crashing even in large-celled two-way tables (see `?loglinear` for an
example). Thus, the only way to do exact inference in multi-way tables
is to use `loglinear()`. We’ll illustrate this using the drugs dataset
from `loglinear()`’s documentation, taken from Agresti (2002, p.322), on
which we’ll test the no-three-way interaction model:

``` r
data(drugs)
ftable(drugs)
#                     Alcohol Yes  No
# Cigarette Marijuana                
# Yes       Yes               991   3
#           No                538  43
# No        Yes                44   2
#           No                456 279

loglinear(subsets(1:3, 2), data = drugs)
# Computing Markov moves (4ti2)... done.
# Running chain (C++)... done.
# Call:
# loglinear(model = subsets(1:3, 2), data = drugs)
# 
# Fitting method:
# Iterative proportional fitting (with stats::loglin)
# 
# MCMC details:
# N = 10000 samples (after thinning), burn in = 1000, thinning = 10
# 
#                  stat p.value     se mid.p.value
#       P(table)         0.6065 0.0049      0.4667
#    Pearson X^2 0.5279  0.6065 0.0049      0.4667
# Likelihood G^2 0.4845  0.6065 0.0049      0.4667
#  Freeman-Tukey 0.4672  0.6065 0.0049      0.4667
#   Cressie-Read 0.512   0.6065 0.0049      0.4667
#     Neyman X^2 0.4294  0.6065 0.0049      0.4667
```

Note that here we’ve used the more concise syntax of facet
specification; if you want to understand the model specification better,
read the documentation in `?loglinear`. You can perform the asymptotic
test with `loglm()` like this:

``` r
MASS::loglm(~ 1*2 + 2*3 + 1*3, data = drugs)
# Call:
# MASS::loglm(formula = ~1 * 2 + 2 * 3 + 1 * 3, data = drugs)
# 
# Statistics:
#                        X^2 df  P(> X^2)
# Likelihood Ratio 0.4845145  1 0.4863845
# Pearson          0.5279994  1 0.4674492
```

Other statistical applications of LattE and 4ti2
================================================

*Note: this section assumes you have
[**latte**](https://github.com/dkahle/latte.git) installed and working
on your machine.*

Most [LattE](https://www.math.ucdavis.edu/~latte/) programs are
available as functions in **latte**, which is imported by **algstat**.
Checkout the readme for **latte**
[here](https://github.com/dkahle/latte).

There are many statistical applications and potential applications of
LattE in R. One example is found in the `count` program, implemented in
`latte::latte_count()`. `latte::latte_count()` counts the number of
integer points in a [convex
polytope](https://en.wikipedia.org/wiki/Convex_polytope). This can be
useful for counting the number of contingency tables with fixed
marginals. **algstat** uses `latte::latte_count()` in the
`count_tables()` function, which determines the number of contingency
tables in the [fiber (isostatistical
region)](http://en.wikipedia.org/wiki/Fiber_(mathematics)) of a table
given an [exponential family
model](http://en.wikipedia.org/wiki/Exponential_family).

``` r
count_tables(politics) # the independence model is the default
# [1] 10
```

For example, we can determine the number of tables with the same row
sums of `politics` as follows:

``` r
(A <- hmat(varlvls = c(2, 2), facets = 1:2)[1:2,])
#    11 12 21 22
# 1+  1  1  0  0
# 2+  0  0  1  1
count_tables(politics, A)
# [1] 121
```

Numerically solving systems of polynomial equations
===================================================

*Note: this section assumes you have [Bertini](https://bertini.nd.edu)
installed and algstat has registered it.*

**algstat** also provides back-end connections to
[Bertini](https://bertini.nd.edu) to solve systems of polynomial
equations. While this work is still being implemented, here’s a peak at
what it can currently do.

First, **algstat** can run raw Bertini programs using `bertini()`. It
also has a nice print method to display the results. For example, here’s
how you would find the intersection of the line f(x) = x and the unit
circle using Bertini:

``` r
code <- "
INPUT

variable_group x, y;
function f, g;

f = x^2 + y^2 - 1;
g = y - x;

END;
"
bertini(code)
# 2 solutions (x,y) found. (2 real, 0 complex)
#     (-0.707,-0.707) (R)
#     ( 0.707, 0.707) (R)
```

Even better, **algstat** can team up with
[**mpoly**](http://github.com/dkahle/mpoly) (working under the hood) to
solve systems of polynomial equations using `poly_solve()`:

``` r
ggvariety(mp("(y - x^2) (y - (2 - x^2))"), xlim = c(-2,2), ylim = c(0,2), n = 351) +
  ggplot2::theme_classic()
```

![](tools/poly-solve-1.png)

``` r
poly_solve(c("y = x^2", "y = 2 - x^2"), varorder = c("x", "y"))
# 2 solutions (x,y) found. (2 real, 0 complex)
#     (-1,1) (R)
#     ( 1,1) (R)
```

Acknowledgements
================

This material is based upon work supported by the National Science
Foundation under Grant Nos.
[1622449](https://nsf.gov/awardsearch/showAward?AWD_ID=1622449) and
[1622369](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1622369).

Installation
============

Installing algstat
------------------

There are currently two ways to get **algstat**. Here’s the preferred
way until we push a new release to CRAN:

-   From Github:

``` r
if(!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("dkahle/mpoly")
devtools::install_github("coneill-math/m2r")
devtools::install_github("dkahle/latte")
devtools::install_github("dkahle/bertini")
devtools::install_github("dkahle/algstat")
```

-   From CRAN: `install.packages("algstat")` (don’t get this now; it’s
    currently out of date)

Installing supporting software
------------------------------

Coming soon! See the links above for direct information.
