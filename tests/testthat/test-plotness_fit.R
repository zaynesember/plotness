set.seed(42)

test_that("plotness_fit() returns a well-formed object", {
  fit <- plotness_fit(rpois(200, 6), type = "poisson")
  expect_s3_class(fit, "plotness_fit")
  expect_named(
    fit,
    c(
      "type", "call", "data", "tidy", "estimates", "model",
      "size", "conf_level", "n_distinct", "degenerate"
    )
  )
  expect_s3_class(fit$model, "lm")
  expect_identical(names(fit$data), names(plotness(rpois(50, 6), "poisson",
    plot = FALSE
  )))
})

test_that("estimates are exposed and consistent", {
  fit <- plotness_fit(rpois(500, 6), type = "poisson")
  e <- fit$estimates
  expect_named(
    e,
    c(
      "type", "par_name", "par_estim", "par_ml", "slope", "intercept",
      "size", "lambda", "conf_level"
    )
  )
  expect_type(e$par_estim, "double")
  expect_equal(e$par_name, "lambda")
  # exp(slope) recovers a rate near the simulated lambda = 6
  expect_equal(e$par_estim, 6, tolerance = 0.5)
  # coef() agrees with the stored estimates
  expect_named(coef(fit), c("intercept", "slope"))
  expect_equal(unname(coef(fit)["slope"]), e$slope)
})

test_that("binomial / nbinomial estimates use prob and agree with coef()", {
  bfit <- plotness_fit(rbinom(300, 12, 0.4), type = "binomial", size = 12)
  expect_equal(bfit$estimates$par_name, "prob")
  expect_equal(bfit$estimates$par_estim, 0.4, tolerance = 0.15)
  expect_true(is.finite(bfit$estimates$par_ml))
  expect_equal(unname(coef(bfit)["slope"]), bfit$estimates$slope)
  expect_equal(unname(coef(bfit)["intercept"]), bfit$estimates$intercept)
  # the regression and ML estimators are computed differently, so differ
  expect_false(isTRUE(all.equal(bfit$estimates$par_estim, bfit$estimates$par_ml)))

  # nbinomial: size estimated from the data
  nfit <- plotness_fit(rnbinom(300, size = 5, prob = 0.4), type = "nbinomial")
  expect_equal(nfit$estimates$par_name, "prob")
  expect_equal(unname(coef(nfit)["slope"]), nfit$estimates$slope)
})

test_that("nbinomial with explicit size uses the method-of-moments prob branch", {
  set.seed(99)
  y <- rnbinom(300, size = 5, prob = 0.4)
  nfit <- plotness_fit(y, type = "nbinomial", size = 5)
  expect_equal(nfit$size, 5)
  # par_ml = size / (size + weighted-mean count) on the positive-freq rows
  tab <- table(y)
  cnt <- as.numeric(names(tab))
  fr <- as.vector(tab)
  xbar <- weighted.mean(cnt[fr > 0], fr[fr > 0])
  expect_equal(nfit$estimates$par_ml, 5 / (5 + xbar))
})

test_that("as.data.frame(legacy=) toggles the column style but not the values", {
  fit <- plotness_fit(rpois(100, 6), type = "poisson")
  legacy <- as.data.frame(fit)
  tidy <- as.data.frame(fit, legacy = FALSE)
  expect_identical(names(legacy), names(fit$data))
  expect_identical(
    names(tidy),
    c(
      "count", "freq", "metameter", "ci_center", "ci_width",
      "ci_lower", "ci_upper", "y_line"
    )
  )
  # only the column names differ; the underlying values are identical
  expect_equal(unname(as.matrix(legacy)), unname(as.matrix(tidy)))
})

test_that("print() shows a summary and returns the object invisibly", {
  fit <- plotness_fit(rpois(100, 6), type = "poisson")
  out <- cli::cli_fmt(print(fit))
  expect_true(any(grepl("plotness_fit", out)))
  expect_true(any(grepl("slope", out)))
  expect_invisible(print(fit))
})

test_that("input validation rejects bad data with tagged conditions", {
  expect_error(plotness_fit(c(2, 5, 11), "binomial", size = 10),
    class = "plotness_error_count_gt_size"
  )
  expect_error(plotness_fit(rpois(50, 6), "poisson", conf_level = 1.2),
    class = "plotness_error_conf_level"
  )
  expect_error(plotness_fit(c(2, 2, 3, -1), "poisson"),
    class = "plotness_error_negative_count"
  )
  expect_error(plotness_fit(c(1.5, 2, 3, 3), "poisson"),
    class = "plotness_error_noninteger_count"
  )
  expect_error(plotness_fit(c(1, 2, NA, 3), "poisson"),
    class = "plotness_error_na"
  )
  expect_error(plotness_fit(rep(3, 20), "poisson"),
    class = "plotness_error_degenerate"
  )
  expect_error(plotness_fit(c(2, 5, 8), "binomial", size = -3),
    class = "plotness_error_size"
  )
  expect_error(plotness_fit(rbinom(50, 10, 0.4), "binomial", size = 2.5),
    class = "plotness_error_size"
  )
  expect_error(
    plotness_fit(data.frame(f = c(-1, 2, 3), c = c(0, 1, 2)), "poisson"),
    class = "plotness_error_negative_freq"
  )
})

test_that("the NA error reports the true number of missing observations", {
  expect_error(plotness_fit(c(1, 2, NA, NA, 3), "poisson"),
    class = "plotness_error_na", regexp = "2 missing"
  )
})

test_that("leveled Poissonness plot shifts the metameter (numeric)", {
  set.seed(7)
  x <- rpois(200, 6)
  lam <- 5
  base <- plotness_fit(x, type = "poisson")
  lev <- plotness_fit(x, type = "poisson", lambda = lam)

  expect_equal(lev$estimates$lambda, lam)
  expect_equal(rel_label(lev), "exp(slope) x lambda")

  # leveled metameter = unleveled + lambda - count * log(lambda) at observed rows
  obs <- !is.na(base$data$Metameter)
  cnt <- base$data$Counts[obs]
  expect_equal(
    lev$data$Metameter[obs],
    base$data$Metameter[obs] + lam - cnt * log(lam)
  )

  # exp(slope) * lambda is invariant to leveling, but the intercept changes
  expect_equal(lev$estimates$par_estim, base$estimates$par_estim)
  expect_false(
    isTRUE(all.equal(lev$estimates$intercept, base$estimates$intercept))
  )
})

test_that("two distinct counts warn but still fit", {
  expect_warning(plotness_fit(c(1, 1, 3, 3), "poisson"),
    class = "plotness_warning_degenerate"
  )
  fit <- suppressWarnings(plotness_fit(c(1, 1, 3, 3), "poisson"))
  expect_true(fit$degenerate)
})

test_that("lambda is rejected/ignored appropriately", {
  expect_error(plotness_fit(rpois(50, 6), "poisson", lambda = -1),
    class = "plotness_error_lambda"
  )
  expect_warning(
    plotness_fit(rbinom(80, 10, 0.5), "binomial", size = 10, lambda = 5),
    class = "plotness_warning_lambda"
  )
  # leveled Poissonness plot accepts a valid lambda
  expect_s3_class(
    plotness_fit(rpois(80, 6), "poisson", lambda = 6),
    "plotness_fit"
  )
})
