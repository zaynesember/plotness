# plotness — Audit & TODO

Audit of the original grad-school version of `plotness`, performed 2026-06-25 on
the `dev` branch. Every runtime finding was reproduced against R 4.5.2 /
ggplot2 4.0.1, with the logic diffed against the upstream `vcd::distplot` it was
adapted from.

Status legend: `[ ]` open · `[x]` done

**As of this pass: all items resolved; `R CMD check` is clean (0/0/0).**

> A subsequent enhancement round (relicense to GPL (>= 2) + vcd attribution,
> input validation, the `plotness_fit()`/`autoplot()` API with exposed
> estimates, Okabe-Ito plot polish, vdiffr/value-level tests, CI + pkgdown +
> spell-check, and a styler/lintr pass) is recorded in `NEWS.md`. It passed an
> adversarial multi-agent review; `R CMD check --as-cran` stays 0/0/0.

---

## 🔴 Correctness bugs

- [x] **B1 — Fitted-line / `xlim` logic is broken; crashes on documented input.**
  `R/plotness.R:130-131`. The line was renamed from `xlim <-` to `x_lim <-`
  during the adaptation, leaving `xlim` `NULL` in the `predict()` call below it.
  `predict()` then silently falls back to `mycount` from the enclosing scope.
  Consequences (both reproduced):
  - Crash `replacement has N rows, data has M` on the documented 2-column
    `(freq, count)` input whenever any count has zero frequency.
  - Garbage line (2 recycled values) when a user actually passes `xlim`.
  - `x_lim` (line 130) is dead code.
  Fix: predict per row — `predict(fm, newdata = data.frame(mycount = RVAL$Counts))`
  — drop the dead `x_lim`, and route the user's `xlim` to the axis
  (`coord_cartesian(xlim = ...)`), not into `predict()`.

- [x] **B2 — Wrong `type` in an example.** `R/plotness.R:36`,
  `man/plotness.Rd:71`, and the README narrative. `plotness(rnbinom(...),
  type="binomial")` feeds negative-binomial data to the binomial path; should be
  `type="nbinomial"`.

## 🟡 ggplot2 modernization (verified warnings on ggplot2 4.0.1)

- [x] **B3 — Deprecated `size` aesthetic.** `R/plotness.R:162`. `geom_line(size=0.75)`
  deprecated since ggplot2 3.4.0; use `linewidth`.
- [x] **B4 — Invalid `labs()` argument.** `R/plotness.R:178`.
  `labs(color="", position="bottom")` → `Ignoring unknown labels: position`.
  Remove `position`.
- [x] **B5 — Partial-match `RVAL$Count`.** `R/plotness.R:139,144,149`. Column is
  `Counts`; emits a partial-match warning. Use `RVAL$Counts`.

## 🟡 Packaging / CRAN-readiness

- [x] **P1 — `LazyData: true` with no `data/` dir** → R CMD check WARNING. Remove it.
  `DESCRIPTION:15`.
- [x] **P2 — License.** `License: MIT` needs `MIT + file LICENSE`, and `LICENSE`
    must be CRAN's 2-line MIT template (a `LICENSE.md` can hold the full text).
    `DESCRIPTION:13`.
- [x] **P3 — Authorship.** Replace `Author:`/`Maintainer:` with `Authors@R`;
    update the stale `zsember@ucsd.edu` email. `DESCRIPTION:5-6`.
- [x] **P4 — `stats` not declared in `Imports:`** though imported in NAMESPACE.
    `DESCRIPTION:10`.
- [x] **P5 — No `.Rbuildignore`** → `plotness.Rproj`/`.Rproj.user` get bundled.
- [x] **P6 — No `.gitignore`.**
- [x] **P7 — Stale `RoxygenNote: 7.1.1`** → regenerate docs with current roxygen2.
- [x] **P8 — Missing `URL` / `BugReports`** in DESCRIPTION.

## 🟡 Spelling / docs

- [x] **D1 — "Poisonness" misspelled** throughout (Title, Description, the plot
    title at `R/plotness.R:138`, README). Correct statistical term is
    **"Poissonness"** (Hoaglin–Tukey).
- [x] **D2 — README typo** "displot" → "distplot" (`README.md:2`); README is a
    single thin paragraph — expand with install + usage.
- [x] **D3 — Duplicated `@description`** between roxygen and DESCRIPTION (cosmetic).

## 🟢 Tests & infrastructure (the "proper package" treatment)

- [x] **T1 — No tests.** Add `tests/testthat/` covering: each distribution type,
    `plot=FALSE` data-frame shape/names, the 2-column zero-freq input (regression
    test for B1), and `xlim` handling.
- [x] **T2 — No `NEWS.md`.**
- [x] **T3 — Vignette** demonstrating the three plot types, misfit diagnosis,
    and data-only mode (`vignettes/plotness.Rmd`).

---

### Note on what is *correct*

The core statistics are faithful to `vcd::distplot` and should not be touched:
metameter math, ML parameter estimates per type, the CI width formula, and the
CI centering on the smoothed `yhat` (not the observed point).
