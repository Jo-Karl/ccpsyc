
  # Helper ------------------------------------------------------------------

  #' Title
  #'
  #' @param object
  #'
  #' @return
  #' @export equival
  #' @import lavaan
  #' @examples
  gamma_hat_scaled <- function (object) {
      fit <- lavInspect(object, "fit")
      p <- length(lavaan::lavNames(object, type = "ov", group = 1))
      nParam <- fit["npar"]
      ngroup <- lavInspect(object, "ngroups")
      n <- lavInspect(object, "ntotal")
      n <- n - ngroup
              gammaHatScaled <- p/(p + 2 * ((fit["chisq.scaled"] -
                  fit["df.scaled"])/n))
              adjGammaHatScaled <- 1 - (((ngroup * p * (p + 1))/2)/fit["df.scaled"]) *
                  (1 - gammaHatScaled)
  		adjGammaHatScaled
  }


  MNCI <- function (object){
  ## McDonald, R. P. (1989). An index of goodness-of-fit based on noncentrality.
  ##  Journal of Classification, 6(1), 97-103. doi: 10.1007/bf01908590
  fit <- inspect(object, "fit") # lavaan's default output
  chisq <- unlist(fit["chisq.scaled"]) #model Chi-square
  df <- unlist(fit["df.scaled"]) # model df
  n <- object@SampleStats@ntotal
  ncp <- max(chisq - df,0)  #non-centrality parameter
  d <- ncp/(n-1)  #scaled non-centrality parameter
  #McDonald's non-centrality index
  Mc <- exp((d)*-.5)
  #Print the value
  Mc
  }

  # Main --------------------------------------------------------------------


  equival <- function(x, dat, group){

  PD_1 <- lavaan::cfa(x, data = eval(substitute(dat), envir = .GlobalEnv,
                                     enclos = parent.frame()),
                      estimator = "MLM",
                      group = paste(group))

  PD_2 <- lavaan::cfa(x, data = eval(substitute(dat), envir = .GlobalEnv,
                                     enclos = parent.frame()),
                      estimator = "MLM", group = paste(group),
                      group.equal = c("loadings"))

  PD_3 <- lavaan::cfa(x, data = eval(substitute(dat), envir = .GlobalEnv,
                                     enclos = parent.frame()),
                      estimator = "MLM", group = paste(group),
                      group.equal = c("loadings", "intercepts"))

  pd_conf <- lavaan::fitmeasures(PD_1)


  pd_metr <- lavaan::fitmeasures(PD_2)


  pd_scal <- lavaan::fitmeasures(PD_3)


  PD <- data.frame(CFI = c(pd_conf["cfi.robust"][[1]],
                           pd_metr["cfi.robust"][[1]],
                           pd_scal["cfi.robust"][[1]]
                           ),
                   RMSEA = c(pd_conf["rmsea.robust"][[1]],
                             pd_metr["rmsea.robust"][[1]],
                             pd_scal["rmsea.robust"][[1]]
                             ),
                   LC = c(pd_conf["rmsea.ci.lower.robust"][[1]],
                          pd_metr["rmsea.ci.lower.robust"][[1]],
                          pd_scal["rmsea.ci.lower.robust"][[1]]
                          ),
                   UC = c(pd_conf["rmsea.ci.upper.robust"][[1]],
                          pd_metr["rmsea.ci.upper.robust"][[1]],
                          pd_scal["rmsea.ci.upper.robust"][[1]]
                          ),
                   SRMR = c(pd_conf["srmr_bentler"][[1]],
                            pd_metr["srmr_bentler"][[1]],
                            pd_scal["srmr_bentler"][[1]]
                            ),
                   GAMMA = c(gamma_hat_scaled(PD_1),
                             gamma_hat_scaled(PD_2),
                             gamma_hat_scaled(PD_3)
                             ),
                   MNCI = c(MNCI(PD_1),
                            MNCI(PD_2),
                            MNCI(PD_3)
                            ),
                   DELTACFI = c(NA, pd_conf["cfi.robust"][[1]] -
                                    pd_metr["cfi.robust"][[1]],
                                    pd_metr["cfi.robust"][[1]] -
                                    pd_scal["cfi.robust"][[1]]),
                   DELTARMSEA = c(NA, pd_conf["rmsea.robust"][[1]] -
                                      pd_metr["rmsea.robust"][[1]],
                                      pd_metr["rmsea.robust"][[1]] -
                                      pd_scal["rmsea.robust"][[1]]),
                   DELTAGAMMA = c(NA, gamma_hat_scaled(PD_1) -
                                      gamma_hat_scaled(PD_2),
                                      gamma_hat_scaled(PD_2) -
                                      gamma_hat_scaled(PD_3)),
                   DELTAMNCI = c(NA, MNCI(PD_1) -
                                     MNCI(PD_2),
                                     MNCI(PD_2) -
                                     MNCI(PD_3)),
                   row.names = c("Configural", "Metric", "Scalar")
                   )
  return(PD)
  }





