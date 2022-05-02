#' Gamma Hat from MLM fitted lavaan object
#'
#' @param object A lavaan output object that was fitted with a MLM estimator
#'
#' @return gamma_hat_scaled
#' @export
#'
#' @examples
#' model_voice <- "voice =~ voice1m + voice2m + voice3m"
#' cfa_out <- lavaan::cfa(model_voice, example, group = "country")
#' gamma_hat_scaled(cfa_out)
gamma_hat_scaled <- function(object) {
  fit <- lavaan::lavInspect(object, "fit")
  p <- length(lavaan::lavNames(object, type = "ov", group = 1))
  nParam <- fit["npar"]
  ngroup <- lavaan::lavInspect(object, "ngroups")
  n <- lavaan::lavInspect(object, "ntotal")
  n <- n - ngroup
  gammaHatScaled <- p / (p + 2 * ((fit["chisq.scaled"] -
    fit["df.scaled"]) / n))
  adjGammaHatScaled <- 1 - (((ngroup * p * (p + 1)) / 2) / fit["df.scaled"]) *
    (1 - gammaHatScaled)
  adjGammaHatScaled
}
