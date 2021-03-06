#' Gamma Hat from MLM fitted lavaan object
#'
#' @param object A lavaan output object that was fitted with a MLM estimator
#'
#' @return gamma_hat_scaled
#' @export
#'
#' @examples
#' gamma_hat_scaled(lavaan1)
gamma_hat_scaled <- function(object) {
  fit <- lavInspect(object, "fit")
  p <- length(lavaan::lavNames(object, type = "ov", group = 1))
  nParam <- fit["npar"]
  ngroup <- lavInspect(object, "ngroups")
  n <- lavInspect(object, "ntotal")
  n <- n - ngroup
  gammaHatScaled <- p / (p + 2 * ((fit["chisq.scaled"] -
    fit["df.scaled"]) / n))
  adjGammaHatScaled <- 1 - (((ngroup * p * (p + 1)) / 2) / fit["df.scaled"]) *
    (1 - gammaHatScaled)
  adjGammaHatScaled
}
