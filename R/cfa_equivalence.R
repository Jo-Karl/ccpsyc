
# Helper ------------------------------------------------------------------

#' One-step equivalence testing
#' The function allows for a simple one step test of configural, metric, and
#' scalar equivalence between multiple groups.
#'
#' @param x CFA model identical to models provided to lavaan.
#' @param dat A data frame or tibble containing the raw data for the
#'            specified model.
#' @param dx A character string that indicates the column of dat that contains
#'           the grouping variable. e.g "country"
#' @param standart_lv A boolean that indicates whether the latent variables
#'                    should be standardised.
#' @param orthogonal A boolean that indicates whether the latent variables
#'                    should be orthogonal.
#' @param estim A string indicating the estimator to be used MLM for complete data and MLR for incomplete data. Defaults to MLM
#' @return Returns a data frame with the fit indices for each model and delta
#'         values comparing the different levels of equivalence.
#' @export equival
#' @examples
#' model_voice <- "voice =~ voice1 + voice2 + voice3"
#' equival(x = model_voice, dat = example, dx = "country")
equival <- function(x, dat, group, standart_lv = TRUE, orthog = TRUE, estim = "MLM") {
  PD_1 <- lavaan::cfa(x,
    data = eval(substitute(dat),
      envir = .GlobalEnv,
      enclos = parent.frame()
    ),
    estimator = estim,
    group = paste(group),
    std.lv = standart_lv,
    orthogonal = orthog
  )

  PD_2 <- lavaan::cfa(x,
    data = eval(substitute(dat),
      envir = .GlobalEnv,
      enclos = parent.frame()
    ),
    estimator = estim, group = paste(group),
    group.equal = c("loadings"),
    std.lv = standart_lv,
    orthogonal = orthog
  )

  PD_3 <- lavaan::cfa(x,
    data = eval(substitute(dat),
      envir = .GlobalEnv,
      enclos = parent.frame()
    ),
    estimator = estim, group = paste(group),
    group.equal = c("loadings", "intercepts"),
    std.lv = standart_lv,
    orthogonal = orthog
  )

  pd_conf <- lavaan::fitmeasures(PD_1)


  pd_metr <- lavaan::fitmeasures(PD_2)


  pd_scal <- lavaan::fitmeasures(PD_3)


  PD <- data.frame(
    CFI = c(
      pd_conf["cfi.robust"][[1]],
      pd_metr["cfi.robust"][[1]],
      pd_scal["cfi.robust"][[1]]
    ),
    RMSEA = c(
      pd_conf["rmsea.robust"][[1]],
      pd_metr["rmsea.robust"][[1]],
      pd_scal["rmsea.robust"][[1]]
    ),
    LC = c(
      pd_conf["rmsea.ci.lower.robust"][[1]],
      pd_metr["rmsea.ci.lower.robust"][[1]],
      pd_scal["rmsea.ci.lower.robust"][[1]]
    ),
    UC = c(
      pd_conf["rmsea.ci.upper.robust"][[1]],
      pd_metr["rmsea.ci.upper.robust"][[1]],
      pd_scal["rmsea.ci.upper.robust"][[1]]
    ),
    SRMR = c(
      pd_conf["srmr_bentler"][[1]],
      pd_metr["srmr_bentler"][[1]],
      pd_scal["srmr_bentler"][[1]]
    ),
    GAMMA = c(
      gamma_hat_scaled(PD_1),
      gamma_hat_scaled(PD_2),
      gamma_hat_scaled(PD_3)
    ),
    MNCI = c(
      MNCI(PD_1),
      MNCI(PD_2),
      MNCI(PD_3)
    ),
    DELTACFI = c(
      NA, pd_conf["cfi.robust"][[1]] -
        pd_metr["cfi.robust"][[1]],
      pd_metr["cfi.robust"][[1]] -
        pd_scal["cfi.robust"][[1]]
    ),
    DELTARMSEA = c(
      NA, pd_conf["rmsea.robust"][[1]] -
        pd_metr["rmsea.robust"][[1]],
      pd_metr["rmsea.robust"][[1]] -
        pd_scal["rmsea.robust"][[1]]
    ),
    DELTAGAMMA = c(
      NA, gamma_hat_scaled(PD_1) -
        gamma_hat_scaled(PD_2),
      gamma_hat_scaled(PD_2) -
        gamma_hat_scaled(PD_3)
    ),
    DELTAMNCI = c(
      NA, MNCI(PD_1) -
        MNCI(PD_2),
      MNCI(PD_2) -
        MNCI(PD_3)
    ),
    row.names = c("Configural", "Metric", "Scalar")
  )
  return(PD)
}
