# plotness

<!-- badges: start -->
[![R-CMD-check](https://github.com/zaynesember/plotness/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zaynesember/plotness/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/zaynesember/plotness/graph/badge.svg)](https://app.codecov.io/gh/zaynesember/plotness)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

An R package improving upon the `vcd` package's `distplot()` function. It lets
you produce **Poissonness**, **binomialness**, and **negative binomialness**
plots with `ggplot2` — or return just the data needed for the plot so you can
draw it with your preferred graphics library.

Ever need to make a Poissonness, binomialness, or negative binomialness plot but
want it to look prettier than the base R output of `vcd::distplot()`? Look no
further.

## Installation

```r
# install.packages("remotes")
remotes::install_github("zaynesember/plotness")
```

## Usage

`plotness()` takes a vector of counts, a 1-way frequency table, or a two-column
`(freq, count)` data frame / matrix, and a distribution `type`:

```r
library(plotness)

# Poissonness plot
plotness(rpois(15, 10), type = "poisson")

# Binomialness plot
plotness(rbinom(15, 10, prob = 0.5), type = "binomial")

# Negative binomialness plot
plotness(rnbinom(15, 10, prob = 0.5), type = "nbinomial")
```

Set `plot = FALSE` to get the underlying data frame instead of a ggplot, so you
can render it however you like:

```r
plotness(rpois(15, 10), type = "poisson", plot = FALSE)
#>   Counts Freq Metameter CI.center CI.width CI.lower CI.upper   y_line
#> ...
```

See `?plotness` for the full list of arguments (confidence intervals, leveled
Poissonness plots via `lambda`, custom titles/labels, and more).

### Getting the fitted estimates

`plotness_fit()` runs the same statistics and returns a `plotness_fit` object
that exposes the fitted parameter estimates (rather than only printing them on
the plot). `plotness()` is a thin wrapper around it plus `autoplot()`:

```r
fit <- plotness_fit(rpois(200, 6), type = "poisson")
fit                       # prints slope/intercept and the parameter estimates
coef(fit)                 # c(intercept, slope)
fit$estimates$par_estim   # exp(slope): the implied rate
ggplot2::autoplot(fit)    # the same figure plotness() would draw
```

## Acknowledgements

The statistical machinery is adapted from `vcd::distplot()` (Meyer, Zeileis &
Hornik). `plotness` reskins it with `ggplot2` and adds a data-only return mode.

## License

GPL (>= 2) © Zayne Sember. `plotness` adapts code from `vcd::distplot()`, which
is GPL-2, so the package is distributed under a compatible GPL license.
