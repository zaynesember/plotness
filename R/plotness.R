#' Poissonness, Binomialness, and Negative Binomialness Plots
#'
#' @description
#' Creates a diagnostic distribution plot or its underlying data frame. Adapted
#' from the \pkg{vcd} package's `distplot()` to produce \pkg{ggplot2} graphics
#' and to allow returning only the data needed for the plot, so the user's
#' preferred graphics package can be used.
#'
#' `plotness()` is a convenience wrapper: it fits the model with
#' [plotness_fit()] and then either draws it (via
#' [autoplot()][autoplot.plotness_fit]) or returns the plotting data.
#'
#' @param x either a vector of counts, a 1-way table of frequencies of counts or
#'  a data frame or matrix with frequencies in the first column and the
#'  corresponding counts in the second column.
#' @param type a character string indicating the distribution.
#' @param size the size argument for the binomial and negative binomial
#'  distribution. If set to NULL and type is "binomial", then size is taken to
#'  be the maximum count. If set to NULL and type is "nbinomial", then size is
#'  estimated from the data.
#' @param lambda parameter of the poisson distribution. If type is "poisson"
#'  and lambda is specified a leveled Poissonness plot is produced.
#' @param plot logical. When `TRUE` a ggplot is returned; otherwise a data frame
#'  containing data necessary for a plot is returned.
#' @param title a character string to be used as the title of the plot. If set
#'  to `NULL` a generic title is produced based on the distribution type.
#' @param legend logical. Should a legend be plotted?
#' @param label logical. Should fit data be printed on the plot?
#' @param conf_int logical. Should confidence intervals be plotted?
#' @param conf_level confidence level for confidence intervals.
#' @param xlim limits for the x-axis.
#' @param xlab a label for the x-axis.
#' @param ylab a label for the y-axis.
#' @return If `plot = TRUE`, a [ggplot2::ggplot] object; otherwise a data frame
#'  containing the data necessary to produce a plot. For programmatic access to
#'  the fitted parameter estimates, use [plotness_fit()].
#' @seealso [plotness_fit()] for the underlying fit object and estimates, and
#'  [autoplot.plotness_fit()] for the plotting method.
#' @examples
#' plotness(rpois(15, 10), type = "poisson")
#' plotness(rbinom(15, 10, prob = 0.5), type = "binomial")
#' plotness(rnbinom(15, 10, prob = 0.5), type = "nbinomial")
#' plotness(rpois(15, 10), type = "poisson", plot = FALSE)
#' @export
plotness <- function(x, type = c("poisson", "binomial", "nbinomial"),
                     size = NULL, lambda = NULL, plot = TRUE, title = NULL,
                     legend = TRUE, label = TRUE, conf_int = TRUE,
                     conf_level = 0.95, xlim = NULL,
                     xlab = "Number of occurrences",
                     ylab = "Distribution metameter") {
  type <- match.arg(type)
  fit <- plotness_fit(
    x,
    type = type, size = size, lambda = lambda, conf_level = conf_level
  )
  if (isTRUE(plot)) {
    autoplot(
      fit,
      title = title, legend = legend, label = label,
      conf_int = conf_int, xlim = xlim, xlab = xlab, ylab = ylab
    )
  } else {
    fit$data
  }
}
