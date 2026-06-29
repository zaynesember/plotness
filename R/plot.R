# Okabe-Ito based palette: colourblind-safe and order-independent.
plotness_pal <- c(observed = "#0072B2", line = "#D55E00", ci = "grey40")

# How the regression slope maps to the parameter estimate, per distribution.
rel_label <- function(object) {
  switch(object$type,
    poisson = if (is.null(object$estimates$lambda)) {
      "exp(slope)"
    } else {
      "exp(slope) x lambda"
    },
    binomial = "inv.logit(slope)",
    nbinomial = "1-exp(slope)"
  )
}

# Build the fit-summary subtitle from the (single source of truth) estimates.
fit_subtitle <- function(object) {
  e <- object$estimates
  paste0(
    "slope = ", round(e$slope, 3),
    ",  intercept = ", round(e$intercept, 3), "\n",
    e$par_name, ": ML = ", round(e$par_ml, 3),
    ",  ", rel_label(object), " = ", round(e$par_estim, 3)
  )
}

#' Plot a `plotness_fit` object
#'
#' Draws the Poissonness / binomialness / negative binomialness plot from a
#' fitted [plotness_fit()] object. This is the [ggplot2::autoplot()] method, so
#' `autoplot(fit)` and `plotness(x, ...)` produce the same figure.
#'
#' @param object A [plotness_fit()] object.
#' @param ... Unused; present for `autoplot()` method compatibility.
#' @inheritParams plotness
#' @return A [ggplot2::ggplot] object. Because it is an ordinary ggplot, you can
#'   append further layers or a different theme with `+`.
#' @seealso [plotness_fit()], [plotness()]
#' @examples
#' fit <- plotness_fit(rpois(150, 8), type = "poisson")
#' ggplot2::autoplot(fit)
#' ggplot2::autoplot(fit, conf_int = FALSE, label = FALSE)
#' @importFrom ggplot2 autoplot ggplot aes geom_line geom_point geom_errorbar
#' @importFrom ggplot2 scale_colour_manual labs theme_bw theme coord_cartesian
#' @importFrom rlang .data %||%
#' @exportS3Method ggplot2::autoplot
autoplot.plotness_fit <- function(object, ..., title = NULL, legend = TRUE,
                                  label = TRUE, conf_int = TRUE, xlim = NULL,
                                  xlab = "Number of occurrences",
                                  ylab = "Distribution metameter") {
  df <- object$data
  title <- title %||% switch(object$type,
    poisson = "Poissonness Plot",
    binomial = "Binomialness Plot",
    nbinomial = "Negative Binomialness Plot"
  )

  p <- ggplot(df, aes(x = .data$Counts))

  # Confidence intervals first (drawn underneath), fixed grey, not a series.
  if (isTRUE(conf_int)) {
    p <- p + geom_errorbar(
      aes(ymin = .data$CI.lower, ymax = .data$CI.upper),
      colour = plotness_pal[["ci"]], width = 0.2, na.rm = TRUE
    )
  }

  p <- p +
    geom_line(
      aes(y = .data$y_line, colour = "Perfect distribution"),
      linewidth = 0.75, key_glyph = "point"
    ) +
    geom_point(
      aes(y = .data$Metameter, colour = "Observed distribution"),
      na.rm = TRUE, key_glyph = "point"
    ) +
    scale_colour_manual(values = c(
      "Observed distribution" = unname(plotness_pal[["observed"]]),
      "Perfect distribution" = unname(plotness_pal[["line"]])
    )) +
    labs(x = xlab, y = ylab, title = title, colour = "") +
    theme_bw()

  if (isTRUE(label)) {
    p <- p + labs(subtitle = fit_subtitle(object))
  }
  p <- p + theme(legend.position = if (isTRUE(legend)) "top" else "none")
  if (!is.null(xlim)) {
    p <- p + coord_cartesian(xlim = xlim)
  }
  p
}

#' @rdname plotness_fit
#' @param x A `plotness_fit` object.
#' @param ... Ignored.
#' @export
print.plotness_fit <- function(x, ...) {
  e <- x$estimates # nolint: object_usage_linter. (used in cli/glue strings below)
  cli::cli_text("{.cls plotness_fit}: {.val {x$type}}")
  cli::cli_text("distinct counts: {x$n_distinct}")
  cli::cli_text(
    "slope = {round(e$slope, 3)}, intercept = {round(e$intercept, 3)}"
  )
  cli::cli_text(
    "{e$par_name}: ML = {round(e$par_ml, 3)}, {rel_label(x)} = {round(e$par_estim, 3)}"
  )
  if (isTRUE(x$degenerate)) {
    cli::cli_alert_warning(
      "Only two distinct counts: the fit is exact, not diagnostic."
    )
  }
  invisible(x)
}

#' @rdname plotness_fit
#' @param object A `plotness_fit` object.
#' @return `coef()` returns a named numeric vector `c(intercept, slope)`.
#' @export
coef.plotness_fit <- function(object, ...) {
  c(intercept = object$estimates$intercept, slope = object$estimates$slope)
}

#' @rdname plotness_fit
#' @param legacy Logical. If `TRUE` (default) return the CamelCase data frame
#'   used by the plot; if `FALSE` return a tidy snake_case version.
#' @return `as.data.frame()` returns the per-count plotting data frame.
#' @export
as.data.frame.plotness_fit <- function(x, ..., legacy = TRUE) {
  if (isTRUE(legacy)) x$data else x$tidy
}
