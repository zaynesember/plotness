# plotness 0.1.0.9000

## New features

* New `plotness_fit()` computes the fit and returns a `plotness_fit` object that
  exposes the fitted parameter estimates programmatically (via `fit$estimates`
  or `coef(fit)`) instead of only as text on the plot. It has `print()`,
  `coef()`, `as.data.frame()`, and an `ggplot2::autoplot()` method.
* `plotness()` is now a thin wrapper over `plotness_fit()` + `autoplot()`; its
  signature and dual return type are unchanged.
* Inputs are validated with clear, class-tagged error messages (negative or
  non-integer counts, counts exceeding `size` for a binomialness plot,
  `conf_level` outside `(0, 1)`, missing values, and degenerate fits).

## Plot changes

* The palette is now the colourblind-safe Okabe-Ito scheme, mapped by a named
  vector so colours no longer depend on factor ordering. Confidence-interval
  bars are grey and visually separate from the observed points.
* The fit summary is drawn as a subtitle rather than positioned inside the panel
  with quantile heuristics, so it can no longer overlap the data.

## Behaviour change

* Input with only one distinct count now errors (class
  `plotness_error_degenerate`) instead of returning a data frame with an
  undefined slope; exactly two distinct counts warn but still plot.
* Supplying `lambda` with `type = "binomial"` or `type = "nbinomial"` now warns
  (class `plotness_warning_lambda`) and ignores it; previously it was silently
  applied to the metameter.

## Licensing

* Relicensed to `GPL (>= 2)`. `plotness` adapts code from `vcd::distplot()`,
  which is GPL-2, so a compatible GPL license is required.

## Bug fixes

* `plotness()` no longer errors on two-column `(freq, count)` input that
  contains a zero frequency. The fitted "perfect distribution" line is now
  evaluated at every count so it always aligns with the returned rows.
* `xlim` now controls the plotted x-axis range instead of corrupting the fitted
  line. Previously passing `xlim` produced a degenerate two-point line.
* Fixed the negative-binomial example, which incorrectly used
  `type = "binomial"`.

## ggplot2 compatibility

* Use `linewidth` instead of the deprecated `size` aesthetic for the fitted
  line (ggplot2 >= 3.4.0).
* Removed an invalid `position` argument passed to `labs()`.
* Reference the `Counts` column directly instead of relying on partial matching.

## Documentation & packaging

* Corrected the spelling "Poissonness" throughout (title, plot titles, docs).
* Moved to `Authors@R`, declared dependencies in `Imports`, dropped the spurious
  `LazyData` field, and added `URL`/`BugReports`.
* Added a `testthat` test suite (including `vdiffr` visual-regression tests), a
  vignette, a `pkgdown` site, GitHub Actions CI, spell checking, `.Rbuildignore`,
  and `.gitignore`.
