#' Non-Centrality Index
#'
#' @param object A lavaan object that was fitted with a MLM estimator/
#'
#'

MNCI <- function(object) {
  ## McDonald, R. P. (1989). An index of goodness-of-fit based on noncentrality.
  ##  Journal of Classification, 6(1), 97-103. doi: 10.1007/bf01908590
  fit <- lavaan::inspect(object, "fit") # lavaan's default output
  chisq <- unlist(fit["chisq.scaled"]) # model Chi-square
  df <- unlist(fit["df.scaled"]) # model df
  n <- object@SampleStats@ntotal
  ncp <- max(chisq - df, 0) # non-centrality parameter
  d <- ncp / (n - 1) # scaled non-centrality parameter
  # McDonald's non-centrality index
  Mc <- exp((d) * -.5)
  # Print the value
  Mc
}
