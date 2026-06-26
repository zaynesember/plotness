# Column names referenced inside ggplot2 aes() are evaluated via data masking,
# which R CMD check cannot see. Declare them so the check is clean.
utils::globalVariables(c(
  "Counts", "Metameter", "y_line", "CI.lower", "CI.upper"
))
