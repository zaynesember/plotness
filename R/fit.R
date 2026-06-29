#' Fit a distribution diagnostic ("-ness") plot
#'
#' @description
#' `plotness_fit()` runs the statistics behind a Poissonness, binomialness, or
#' negative binomialness plot and returns a `plotness_fit` object. The object
#' carries the per-count plotting data, the fitted line, and the parameter
#' estimates, and has [autoplot()][autoplot.plotness_fit], `print()`, `coef()`,
#' and `as.data.frame()` methods.
#'
#' Most users can call [plotness()] directly; use `plotness_fit()` when you want
#' the fitted estimates programmatically (via `fit$estimates` or `coef(fit)`)
#' rather than only as text on the plot.
#'
#' @inheritParams plotness
#' @return An object of class `plotness_fit`: a list with elements
#'   \describe{
#'     \item{`type`}{the distribution, one of `"poisson"`, `"binomial"`,
#'       `"nbinomial"`.}
#'     \item{`data`}{the plotting data frame (one row per distinct count).}
#'     \item{`tidy`}{a snake_case copy of `data`.}
#'     \item{`estimates`}{a named list of fitted estimates (see below).}
#'     \item{`model`}{the underlying [stats::lm()] fit.}
#'     \item{`size`, `conf_level`}{the resolved size and confidence level.}
#'     \item{`n_distinct`, `degenerate`}{number of distinct observed counts and
#'       whether the fit is degenerate (exactly two distinct counts).}
#'   }
#'   `estimates` contains `par_name`, `par_estim` (the regression-implied
#'   estimate), `par_ml` (the [vcd::goodfit()] maximum-likelihood estimate),
#'   `slope`, `intercept`, `size`, `lambda`, and `conf_level`.
#' @seealso [plotness()], [autoplot.plotness_fit()]
#' @examples
#' fit <- plotness_fit(rpois(150, 8), type = "poisson")
#' fit
#' coef(fit)
#' fit$estimates$par_estim
#' @importFrom vcd goodfit
#' @importFrom rlang arg_match
#' @import stats
#' @export
plotness_fit <- function(x, type = c("poisson", "binomial", "nbinomial"),
                         size = NULL, lambda = NULL, conf_level = 0.95) {
  type <- arg_match(type)
  cl <- match.call()

  # ---- reshape input to (count, freq) --------------------------------------
  # Flag missing values rather than silently dropping them in table(); report
  # the true count of missing observations.
  if (is.vector(x)) {
    if (anyNA(x)) {
      cli::cli_abort(
        c(
          "{.arg x} contains {sum(is.na(x))} missing value{?s}.",
          "i" = "Remove or impute them before plotting."
        ),
        class = "plotness_error_na"
      )
    }
    x <- table(x)
  }
  if (is.table(x)) {
    if (length(dim(x)) > 1) {
      cli::cli_abort("{.arg x} must be a 1-way table.")
    }
    freq <- as.vector(x)
    count <- as.numeric(names(x))
  } else {
    if (!(!is.null(ncol(x)) && ncol(x) == 2)) {
      cli::cli_abort(
        "{.arg x} must be a 2-column matrix or data frame of (freq, count)."
      )
    }
    freq <- as.vector(x[, 1])
    count <- as.vector(x[, 2])
  }

  # ---- validation (shared by the plot and data paths) ----------------------
  check_plotness_input(count, freq, type, size, lambda, conf_level)

  # ---- binomial size fallback (counts already validated; max() is safe) ----
  if (type == "binomial" && is.null(size)) {
    size <- max(count)
    cli::cli_warn("size was not given, taken as maximum count")
  }

  # lambda only levels the Poissonness plot; ignore it elsewhere (already warned)
  if (type != "poisson") lambda <- NULL

  n_distinct <- length(unique(count[freq > 0]))
  degenerate <- n_distinct == 2

  myindex <- which(freq > 0)
  mycount <- count[myindex]
  myfreq <- freq[myindex]

  # ---- per-type metameter + fit (arithmetic preserved from vcd::distplot) ---
  switch(type,
    poisson = {
      par_ml <- suppressWarnings(goodfit(x, type = type)$par$lambda)
      phi <- function(nk, k, n_total, size = NULL) {
        ifelse(nk > 0, lgamma(k + 1) + log(nk / n_total), NA)
      }
      y <- phi(myfreq, mycount, sum(freq))
      if (!is.null(lambda)) y <- y + lambda - mycount * log(lambda)
      fm <- lm(y ~ mycount)
      par_estim <- exp(coef(fm)[2])
      if (!is.null(lambda)) par_estim <- par_estim * lambda
      par_estim <- unname(par_estim)
      par_name <- "lambda"
    },
    binomial = {
      par_ml <- suppressWarnings(
        goodfit(x, type = type, par = list(size = size))$par$prob
      )
      phi <- function(nk, k, n_total, size = NULL) {
        log(nk) - log(n_total * choose(size, k))
      }
      y <- phi(myfreq, mycount, sum(freq), size = size)
      fm <- lm(y ~ mycount)
      b <- exp(coef(fm)[2])
      par_estim <- unname(b / (1 + b))
      par_name <- "prob"
    },
    nbinomial = {
      if (is.null(size)) {
        gf <- suppressWarnings(goodfit(x, type = type)$par)
        size <- gf$size
        par_ml <- gf$prob
      } else {
        xbar <- weighted.mean(mycount, myfreq)
        par_ml <- size / (size + xbar)
      }
      phi <- function(nk, k, n_total, size = NULL) {
        log(nk) - log(n_total * choose(size + k - 1, k))
      }
      y <- phi(myfreq, mycount, sum(freq), size = size)
      fm <- lm(y ~ mycount)
      par_estim <- unname(1 - exp(coef(fm)[2]))
      par_name <- "prob"
    }
  )

  yhat <- ifelse(myfreq > 1.5, myfreq - 0.67, 1 / exp(1))
  yhat <- phi(yhat, mycount, sum(freq), size = size)
  if (!is.null(lambda)) yhat <- yhat + lambda - mycount * log(lambda)

  phat <- myfreq / sum(myfreq)
  ci_width <- qnorm(1 - (1 - conf_level) / 2) *
    sqrt(1 - phat) / sqrt(myfreq - (0.25 * phat + 0.47) * sqrt(myfreq))

  plot_data <- cbind(count, freq, NA, NA, NA, NA, NA)
  plot_data[myindex, 3:7] <- cbind(y, yhat, ci_width, yhat - ci_width, yhat + ci_width)
  plot_data <- as.data.frame(plot_data)
  names(plot_data) <- c(
    "Counts", "Freq", "Metameter", "CI.center",
    "CI.width", "CI.lower", "CI.upper"
  )
  # Fitted "perfect distribution" line evaluated at every count so the column
  # aligns with plot_data even when some counts have zero frequency.
  plot_data$y_line <- predict(fm, newdata = data.frame(mycount = plot_data$Counts))

  if (anyNA(ci_width)) {
    cli::cli_warn(c(
      "Some confidence intervals are undefined and were dropped.",
      "i" = "This can happen for very small frequencies."
    ))
  }

  tidy <- plot_data
  names(tidy) <- c(
    "count", "freq", "metameter", "ci_center",
    "ci_width", "ci_lower", "ci_upper", "y_line"
  )

  estimates <- list(
    type = type,
    par_name = par_name,
    par_estim = par_estim,
    par_ml = unname(par_ml),
    slope = unname(coef(fm)[2]),
    intercept = unname(coef(fm)[1]),
    size = size,
    lambda = lambda,
    conf_level = conf_level
  )

  structure(
    list(
      type = type,
      call = cl,
      data = plot_data,
      tidy = tidy,
      estimates = estimates,
      model = fm,
      size = size,
      conf_level = conf_level,
      n_distinct = n_distinct,
      degenerate = degenerate
    ),
    class = "plotness_fit"
  )
}

#' Validate `plotness()` / `plotness_fit()` input
#'
#' Internal helper; runs every guard once so the plot and data paths share
#' identical validation. All conditions are class-tagged for testability.
#'
#' @param count,freq Numeric vectors of counts and their frequencies.
#' @param type Resolved distribution string.
#' @param size,lambda,conf_level As in [plotness_fit()].
#' @param call Environment to report errors against.
#' @return `invisible(NULL)`, called for its side effects (errors/warnings).
#' @keywords internal
#' @importFrom rlang caller_env
check_plotness_input <- function(count, freq, type, size, lambda, conf_level,
                                 call = caller_env()) {
  if (!is.numeric(conf_level) || length(conf_level) != 1 ||
    is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    cli::cli_abort(
      "{.arg conf_level} must be a single number in (0, 1), not {.val {conf_level}}.",
      class = "plotness_error_conf_level", call = call
    )
  }

  n_na <- sum(is.na(count)) + sum(is.na(freq))
  if (n_na > 0) {
    cli::cli_abort(
      c(
        "{.arg x} contains {n_na} missing value{?s}.",
        "i" = "Remove or impute them before plotting."
      ),
      class = "plotness_error_na", call = call
    )
  }

  if (any(count < 0)) {
    cli::cli_abort(
      "Counts must be non-negative; found {.val {count[count < 0]}}.",
      class = "plotness_error_negative_count", call = call
    )
  }

  noninteger <- abs(count - round(count)) > 1e-8
  if (any(noninteger)) {
    cli::cli_abort(
      "Counts must be whole numbers; found {.val {count[noninteger]}}.",
      class = "plotness_error_noninteger_count", call = call
    )
  }

  if (any(freq < 0)) {
    cli::cli_abort(
      "Frequencies must be non-negative.",
      class = "plotness_error_negative_freq", call = call
    )
  }

  if (!is.null(size)) {
    if (!is.numeric(size) || length(size) != 1 || is.na(size) ||
      size <= 0 || abs(size - round(size)) > 1e-8) {
      cli::cli_abort(
        "{.arg size} must be a single positive whole number.",
        class = "plotness_error_size", call = call
      )
    }
    if (type == "binomial" && any(count > size)) {
      cli::cli_abort(
        c(
          "All counts must be <= {.arg size} ({size}) for a binomialness plot.",
          "x" = "The largest observed count is {max(count)}.",
          "i" = "Increase {.arg size} or check the data."
        ),
        class = "plotness_error_count_gt_size", call = call
      )
    }
  }

  if (!is.null(lambda)) {
    if (!is.numeric(lambda) || length(lambda) != 1 || !is.finite(lambda) ||
      lambda <= 0) {
      cli::cli_abort(
        "{.arg lambda} must be a single positive number.",
        class = "plotness_error_lambda", call = call
      )
    }
    if (type != "poisson") {
      cli::cli_warn(
        "{.arg lambda} is ignored for {.val {type}}.",
        class = "plotness_warning_lambda"
      )
    }
  }

  n_distinct <- length(unique(count[freq > 0]))
  if (n_distinct <= 1) {
    cli::cli_abort(
      c(
        "Cannot fit a line through fewer than two distinct counts.",
        "i" = "A distributionness plot needs at least two distinct observed counts."
      ),
      class = "plotness_error_degenerate", call = call
    )
  }
  if (n_distinct == 2) {
    cli::cli_warn(
      "With only two distinct counts the fit is exact and the plot is not diagnostic.",
      class = "plotness_warning_degenerate"
    )
  }

  invisible(NULL)
}
