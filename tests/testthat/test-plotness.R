set.seed(42)

test_that("each distribution type returns a ggplot", {
  expect_s3_class(plotness(rpois(50, 10), type = "poisson"), "ggplot")
  expect_s3_class(
    suppressWarnings(plotness(rbinom(50, 10, prob = 0.5), type = "binomial")),
    "ggplot"
  )
  expect_s3_class(plotness(rnbinom(50, 10, prob = 0.5), type = "nbinomial"),
                  "ggplot")
})

test_that("plot = FALSE returns the expected data frame", {
  counts <- rpois(50, 10)
  df <- plotness(counts, type = "poisson", plot = FALSE)
  expect_s3_class(df, "data.frame")
  expect_identical(
    names(df),
    c("Counts", "Freq", "Metameter", "CI.center", "CI.width",
      "CI.lower", "CI.upper", "y_line")
  )
  # one row per distinct observed count
  expect_equal(nrow(df), length(unique(counts)))
  expect_equal(sum(df$Freq), length(counts))
})

test_that("two-column input with a zero frequency does not crash (regression: B1)", {
  # Counts with a gap: count 1 has zero frequency. The fitted-line column must
  # still align with every row of the result, plot or no plot.
  gap <- data.frame(freq = c(5, 0, 8, 3), count = c(0, 1, 2, 3))
  expect_no_error(plotness(gap, type = "poisson", plot = FALSE))
  expect_no_error(plotness(gap, type = "poisson", plot = TRUE))

  df <- plotness(gap, type = "poisson", plot = FALSE)
  expect_equal(nrow(df), 4)
  expect_false(any(is.na(df$y_line)))      # line is defined at every count
})

test_that("y_line is the fitted regression line, one value per count", {
  df <- plotness(rpois(50, 10), type = "poisson", plot = FALSE)
  expect_equal(length(df$y_line), nrow(df))
  # fitted line through a single lm() is monotone in Counts
  ord <- order(df$Counts)
  expect_true(!is.unsorted(df$y_line[ord]) || !is.unsorted(rev(df$y_line[ord])))
})

test_that("passing xlim does not corrupt the fitted line (regression: B1)", {
  df <- plotness(rpois(50, 10), type = "poisson", xlim = c(0, 20), plot = FALSE)
  # Before the fix, y_line collapsed to two recycled values.
  expect_equal(length(df$y_line), nrow(df))
  expect_gt(length(unique(round(df$y_line, 6))), 2)
})

test_that("binomial without size warns and defaults to the maximum count", {
  expect_warning(
    plotness(rbinom(50, 10, prob = 0.5), type = "binomial"),
    "size was not given"
  )
})
