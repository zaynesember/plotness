# plotness 0.1.0.9000

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
* Moved to `Authors@R`, declared `stats` in `Imports`, dropped the spurious
  `LazyData` field, and added `URL`/`BugReports`.
* Switched the license metadata to the CRAN `MIT + file LICENSE` form.
* Added a `testthat` test suite, `.Rbuildignore`, and `.gitignore`.
