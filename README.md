<!-- README.md is generated from README.Rmd. Please edit that file -->
algstat
=======

**algstat** is a collection of tools to help do algebraic statistics in R. Many (but not all) of the tools make use of back-end connections to software used in various math communities, such as [Macaulay2](http://www.math.uiuc.edu/Macaulay2/) (for algebraic computations), [Bertini](https://bertini.nd.edu) (for solving systems of polynomial equations), and [LattE](https://www.math.ucdavis.edu/~latte/) with [4ti2](http://www.4ti2.de) (for lattice problems and, in particular, the computation of Markov bases).

This brief intro is currently under construction.

Exact inference with log-linear models
======================================

One of the most well-developed parts of the package allows users to perform (conditional) exact tests for log-linear models. There are several great references on the math behind this, such as [Diaconis and Sturmfels' original paper](http://projecteuclid.org/euclid.aos/1030563990), the [Lectures on Algebraic Statistics](http://smile.amazon.com/Lectures-Algebraic-Statistics-Oberwolfach-Seminars/dp/3764389044/ref=sr_1_1?ie=UTF8&qid=1430536908&sr=8-1&keywords=lectures+on+algebraic+statistics), and [Markov Bases in Algebraic Statistics](http://smile.amazon.com/Markov-Bases-Algebraic-Statistics-Springer/dp/1461437180/ref=sr_1_fkmr0_1?ie=UTF8&qid=1430536933&sr=8-1-fkmr0&keywords=aoki%2C+hada%2C+and+takemura), so we'll keep the technical discussion to a minimum.

We'll begin by doing Fisher's exact test on a built-in dataset called politics.

``` r
library(algstat)
#> Loading required package: mpoly
#> Loading required package: stringr
#> 
#> Attaching package: 'algstat'
#> 
#> The following object is masked _by_ '.GlobalEnv':
#> 
#>     tabFill
data(politics)
politics
#>            Party
#> Personality Democrat Republican
#>   Introvert        3          7
#>   Extrovert        6          4

fisher.test(politics)
#> 
#>  Fisher's Exact Test for Count Data
#> 
#> data:  politics
#> p-value = 0.3698
#> alternative hypothesis: true odds ratio is not equal to 1
#> 95 percent confidence interval:
#>  0.03005364 2.46429183
#> sample estimates:
#> odds ratio 
#>   0.305415
```

Since the independence model is log-linear, this exact same procedure can be done with **algstat**. The go-to function here is `loglinear()` (formerly `hierarchical()`):

``` r
loglinear(~ Personality + Party, data = politics)
#> Computing Markov moves... done.
#> Running chain... done.
#> Call:
#> loglinear(model = ~Personality + Party, data = politics)
#> 
#> Fitting method:
#> Iterative proportional fitting (with stats::loglin)
#> 
#> MCMC details:
#> N = 10000 samples (after thinning), burn in = 1000, thinning = 10
#> 
#>       Distance   Stat     SE p.value     SE mid.p.value
#>        P(samp)                0.3645 0.0048      0.2178
#>    Pearson X^2 1.8182 0.0149  0.3645 0.0048      0.2178
#> Likelihood G^2 1.848  0.0159  0.3645 0.0048      0.2178
#>  Freeman-Tukey 1.8749 0.0172  0.3645 0.0048      0.2178
#>   Cressie-Read 1.8247 0.0151  0.3645 0.0048      0.2178
#>     Neyman X^2 2.0089 0.0234  0.3645 0.0048      0.2938
```

Exact inference in algebraic statistics is done using MCMC to sample from the conditional distribution of the data given its sufficient statistics under the model. Consequently, the p-values estimated are only determined up to Monte Carlo error. The standard p-value is given under the column `p.value` in the row labeled `P(samp)`. The analogous asymptotic test can be done in either of two ways.

The first way uses the `loglin()` function from the **stats** package. It outputs the likelihood ratio statistic (`Likelihood G^2` in the output above) and Pearson's chi-squared statistic (`Pearson X^2` above), and you have to calculate the p-value yourself.

``` r
(loglinMod <- stats::loglin(politics, list(1, 2)))
#> 2 iterations: deviation 0
#> $lrt
#> [1] 1.848033
#> 
#> $pearson
#> [1] 1.818182
#> 
#> $df
#> [1] 1
#> 
#> $margin
#> $margin[[1]]
#> [1] "Personality"
#> 
#> $margin[[2]]
#> [1] "Party"
pchisq(loglinMod$pearson, df = 1, lower.tail = FALSE)
#> [1] 0.1775299
```

The second way is `loglm()` in the **MASS** package, which is a nice wrapper of `loglin()` (in fact, **algstat**'s `loglinear()` function uses the IPF implementation from `loglin()`, although it doesn't need to). It's syntax looks identical to `loglinear()`'s above, but `loglinear()` can actually take a wide range of model specification syntaxes. Here's `loglm()`:

``` r
MASS::loglm(~ Personality + Party, data = politics)
#> Call:
#> MASS::loglm(formula = ~Personality + Party, data = politics)
#> 
#> Statistics:
#>                       X^2 df  P(> X^2)
#> Likelihood Ratio 1.848033  1 0.1740123
#> Pearson          1.818182  1 0.1775299
```

Installation
============

Installing algstat
------------------

-   From CRAN: `install.packages("algstat")` (this is not up-to-date)

-   From Github (dev version):

    ``` r
    # install.packages("devtools")
    # install.packages("mpoly")
    devtools::install_github("dkahle/algstat")
    ```

Installing supporting software
------------------------------

Coming soon! See the links above for direct information.
