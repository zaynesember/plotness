#' @keywords internal
"_PACKAGE"

# Re-export the ggplot2 autoplot generic so `autoplot(fit)` works for users who
# have only attached plotness.
#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot
