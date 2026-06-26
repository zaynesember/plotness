# plotness

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

## Acknowledgements

The statistical machinery is adapted from `vcd::distplot()` (Meyer, Zeileis &
Hornik). `plotness` reskins it with `ggplot2` and adds a data-only return mode.

## License

MIT © Zayne Sember
