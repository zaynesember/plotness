set.seed(42)

test_that("each distribution type returns a ggplot", {
  expect_s3_class(plotness(rpois(50, 10), type = "poisson"), "ggplot")
  expect_s3_class(
    suppressWarnings(plotness(rbinom(50, 10, prob = 0.5), type = "binomial")),
    "ggplot"
  )
  expect_s3_class(
    plotness(rnbinom(50, 10, prob = 0.5), type = "nbinomial"),
    "ggplot"
  )
})

test_that("plot = FALSE returns the expected data frame", {
  counts <- rpois(50, 10)
  df <- plotness(counts, type = "poisson", plot = FALSE)
  expect_s3_class(df, "data.frame")
  expect_identical(
    names(df),
    c(
      "Counts", "Freq", "Metameter", "CI.center", "CI.width",
      "CI.lower", "CI.upper", "y_line"
    )
  )
  # one row per distinct observed count
  expect_equal(nrow(df), length(unique(counts)))
  expect_equal(sum(df$Freq), length(counts))
})

test_that("two-column input with a zero frequency aligns values to rows (regression: B1)", {
  # Counts with a gap: count 1 has zero frequency. Metameter/CI must be NA at
  # the gap and defined elsewhere, and the values must land on the right rows.
  gap <- data.frame(freq = c(5, 0, 8, 3), count = c(0, 1, 2, 3))
  expect_no_error(plotness(gap, type = "poisson", plot = TRUE))

  df <- plotness(gap, type = "poisson", plot = FALSE)
  expect_equal(nrow(df), 4)

  # the zero-frequency row (Counts == 1) is NA in every fitted column ...
  na_cols <- c("Metameter", "CI.center", "CI.width", "CI.lower", "CI.upper")
  expect_true(all(is.na(df[df$Counts == 1, na_cols])))
  # ... but the non-gap rows are fully populated, and y_line everywhere.
  expect_false(any(is.na(df[df$Counts != 1, na_cols])))
  expect_false(any(is.na(df$y_line)))

  # value anchor: metameter at Counts==2 is lgamma(2+1) + log(8/16) == 0
  expect_equal(df$Metameter[df$Counts == 2], lgamma(3) + log(8 / 16))
})

test_that("y_line equals the fitted line through every count", {
  fit <- plotness_fit(rpois(50, 10), type = "poisson")
  expect_equal(
    fit$data$y_line,
    fit$estimates$intercept + fit$estimates$slope * fit$data$Counts
  )
})

test_that("passing xlim does not corrupt the fitted line (regression: B1)", {
  df <- plotness(rpois(50, 10), type = "poisson", xlim = c(0, 20), plot = FALSE)
  # Before the fix, y_line collapsed to two recycled values.
  expect_equal(length(df$y_line), nrow(df))
  expect_gt(length(unique(round(df$y_line, 6))), 2)
})

test_that("xlim is applied to the plot via coord_cartesian", {
  p <- plotness(rpois(50, 10), type = "poisson", xlim = c(0, 20))
  expect_s3_class(p$coordinates, "CoordCartesian")
  expect_equal(p$coordinates$limits$x, c(0, 20))
  # without xlim the x limits are left unset
  p0 <- plotness(rpois(50, 10), type = "poisson")
  expect_null(p0$coordinates$limits$x)
})

test_that("plotness() facade matches autoplot() with options threaded through", {
  x <- rpois(80, 8)
  fit <- plotness_fit(x, type = "poisson")
  for (ci in c(TRUE, FALSE)) {
    for (lab in c(TRUE, FALSE)) {
      expect_equal(
        ggplot2::ggplot_build(
          plotness(x, "poisson", conf_int = ci, label = lab)
        )$data,
        ggplot2::ggplot_build(
          ggplot2::autoplot(fit, conf_int = ci, label = lab)
        )$data
      )
    }
  }
  p1 <- plotness(x, "poisson", title = "T", xlab = "XX", legend = FALSE)
  p2 <- ggplot2::autoplot(fit, title = "T", xlab = "XX", legend = FALSE)
  expect_identical(
    p1$labels[c("title", "x", "y", "subtitle")],
    p2$labels[c("title", "x", "y", "subtitle")]
  )
})

test_that("binomial without size warns and defaults to the maximum count", {
  expect_warning(
    plotness(rbinom(50, 10, prob = 0.5), type = "binomial"),
    "size was not given"
  )
})
