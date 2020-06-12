# oglm

<!-- badges: start -->
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Travis build status](https://travis-ci.com/linogaliana/oglm.svg?branch=master)](https://travis-ci.com/linogaliana/oglm)
[![Codecov test coverage](https://codecov.io/gh/linogaliana/oglm/branch/master/graph/badge.svg)](https://codecov.io/gh/linogaliana/oglm?branch=master)
<!-- badges: end -->

A package to work with interval regression and ordered logit and probit models

:warning: `oglm` is a personal adaptation of `oglmx` :package: available on [CRAN](https://cran.r-project.org/web/packages/oglmx/index.html). Because the
:package: does not look actively maintened and I needed a few enhancement, I 
started to work on the source code

## Installation

You can install the github version of `oglm` with:

``` r
devtools::install_github("linogaliana/oglm")
```

## Difference with `oglmx` :package:

This :package: has a few difference with `oglmx` :package: that are listed here:

* The :package: presents some tests and is built using continuous integration (based on `travis`) ! 
* `predict` method
* ...
