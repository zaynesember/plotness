# Visual-regression tests. Baselines live in tests/testthat/_snaps/.
# These are skipped on CRAN and only fail locally/CI when a figure changes.

test_that("autoplot produces stable figures", {
  skip_if_not_installed("vdiffr")

  set.seed(1)
  pois <- plotness_fit(rpois(200, 6), type = "poisson")
  vdiffr::expect_doppelganger("poissonness", ggplot2::autoplot(pois))
  vdiffr::expect_doppelganger(
    "poissonness-bare",
    ggplot2::autoplot(pois, conf_int = FALSE, label = FALSE, legend = FALSE)
  )

  set.seed(2)
  binom <- plotness_fit(rbinom(200, 10, prob = 0.4),
    type = "binomial",
    size = 10
  )
  vdiffr::expect_doppelganger("binomialness", ggplot2::autoplot(binom))

  set.seed(3)
  nbin <- plotness_fit(rnbinom(200, size = 5, prob = 0.4), type = "nbinomial")
  vdiffr::expect_doppelganger("negative-binomialness", ggplot2::autoplot(nbin))
})
