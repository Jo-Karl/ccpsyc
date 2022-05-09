
####################################################################################
################# Invariance Effect Sizes############################################
################## Johannes A. Karl #################################################
################### Based on code by Heather Gunn####################################
# https://www.tandfonline.com/doi/full/10.1080/10705511.2019.1689507?journalCode=hsem20
####################################################################################
#' Invariance Effect Sizes
#'
#' @param df Multi-group dataframe
#' @param items vector of items for the target construct
#' @param group string defining grouping variable
#' @param intercept_fix Which item should have a fixed intercept defaults to the first item
#' @param lowerLV Lower range of latent variable tested
#' @param upperLV Upper range of latent variable tested
#' @param nodewidth Steps tested
#' @param ... Passes on to lavaan CFA functions
#' @return Returns a dataframe with a row for each item comprising the uni-factorial solution and one column for each invariance effect size.
#'  A detailed interpretation of each effect size is provided in  \href{https://www.tandfonline.com/doi/full/10.1080/10705511.2019.1689507?journalCode=hsem20}{Gunn et al. (2019)}.
#' @export invariance_eff
#'

invariance_eff <- function(df, items, group, nodewidth = .01, intercept_fix = 1, lowerLV = -10, upperLV = 10, ...) {
  ## Declaring global variables to avoid notes by CRAN
  op <- ustart <- label <- NULL

  sum_item <- paste(items, collapse = " + ")
  main <- paste("f1", sum_item, sep = " =~ ")
  fact_intercept <- "f1 ~ c(selc1, selc2) * 1"
  item_int_zero <- paste(items[[intercept_fix]], "~ 0 * 1")
  cfa_model <- paste(main, fact_intercept, sep = "\n")
  cfa_model_f <- paste(cfa_model, item_int_zero, sep = "\n")
  ### Things to fix :
  ### Allow users to select groups
  ### Allow for LV switching
  ### throw an error if model does not converge
  ### Potentially try to fix different factor loadings if convergence fails
  ### Allow users to switch between aproximation of SD and providing variable names for SD
  ### Integrate dMacs signed.
  ### Remove reduntand lines of code, some computations are repeated, some things are named that do not need to be named. (For example the indicators)

  ### Calculating models
  fit.cfa <- lavaan::cfa(main, data = df, group = group, std.lv = T, ...)
  outNoninvariance_a <- lavaan::cfa(main, data = df, std.lv = T, meanstructure = T, ...)
  outNoninvariance <- lavaan::cfa(main, data = df, group = group, std.lv = T, ...)
  outNoninvariance_f <- lavaan::cfa(cfa_model_f,
    data = df, group = group, std.lv = T,
    meanstructure = T, ...
  )

  par_cfa <- lavaan::parametertable(outNoninvariance)
  par_cfa_a <- lavaan::parametertable(outNoninvariance_a)
  par_cfa_f <- lavaan::parametertable(outNoninvariance_f)
  ### Defining paramters
  p <- outNoninvariance@pta$nvar[[1]]
  loading1 <- dplyr::filter(par_cfa, group == "1", op == "=~")$est
  loading2 <- dplyr::filter(par_cfa, group == "2", op == "=~")$est
  loadingS <- dplyr::filter(par_cfa_a, op == "=~")$est
  intercept1 <- dplyr::filter(par_cfa, group == "1", op == "~1", is.na(ustart))$est
  intercept2 <- dplyr::filter(par_cfa, group == "2", op == "~1", is.na(ustart))$est
  stdev1 <- round(dplyr::filter(par_cfa, group == "1", op == "~1", is.na(ustart))$se * sqrt(outNoninvariance@Data@nobs[[1]]), 3)
  stdev2 <- round(dplyr::filter(par_cfa, group == "2", op == "~1", is.na(ustart))$se * sqrt(outNoninvariance@Data@nobs[[2]]), 3)
  fmean1 <- dplyr::filter(par_cfa_f, group == "1", label == "selc1")$est
  fmean2 <- dplyr::filter(par_cfa_f, group == "2", label == "selc2")$est
  fsd1 <- dplyr::filter(par_cfa_f, group == "1", label == "selc1")$se * sqrt(outNoninvariance@Data@nobs[[1]])
  fsd2 <- dplyr::filter(par_cfa_f, group == "2", label == "selc2")$se * sqrt(outNoninvariance@Data@nobs[[2]])
  N1 <- outNoninvariance@Data@nobs[[1]]
  N2 <- outNoninvariance@Data@nobs[[2]]
  interceptS <- dplyr::filter(par_cfa_a, op == "~1", is.na(ustart))$est

  ### Calculation UDI SDI
  # Define evaluation of latent variable
  LV <- seq(lowerLV, upperLV, nodewidth)

  # Create empty matrices for future arrays, matrices, etc.
  DiffExpItem12 <- matrix(NA, length(LV), p)
  pdfLV2 <- matrix(NA, length(LV), 1)
  SDI2numerator <- matrix(NA, length(LV), p)
  UDI2numerator <- matrix(NA, length(LV), p)
  SDI2 <- matrix(NA, p, 1)
  UDI2 <- matrix(NA, p, 1)

  # Calculate SDI2 and UDI2
  for (j in 1:p) {
    for (k in 1:length(LV)) {
      # Calculate difference in expected indicator scores between groups 1 and 2
      DiffExpItem12[k, j] <- (intercept1[j] - intercept2[j]) + (loading1[j] - loading2[j]) * LV[k]
      # probability density function for sample estimate of group 2 latent variable distribution
      pdfLV2[k] <- stats::dnorm(LV[k], mean = fmean2, sd = fsd2)

      # Multiply by latent variable distribution to calculate individual data point in numerator
      SDI2numerator[k, j] <- DiffExpItem12[k, j] * pdfLV2[k] * nodewidth
      UDI2numerator[k, j] <- abs(SDI2numerator[k, j])
    }
    # Sum across range of latent variable using quadrature to calculate numerator & divide by denominator
    SDI2[j] <- sum(SDI2numerator[, j]) / stdev2[j]
    UDI2[j] <- sum(UDI2numerator[, j]) / stdev2[j]
  }





  # Create proportions
  prop1 <- N1 / (N1 + N2)
  prop2 <- N2 / (N1 + N2)

  # Create empty matrices for future arrays, matrices, etc.
  DiffExpItem1S <- matrix(NA, length(LV), p)
  DiffExpItemS2 <- matrix(NA, length(LV), p)
  pdfLV1 <- matrix(NA, length(LV), 1)
  pdfLV2 <- matrix(NA, length(LV), 1)
  WSDInumerator1 <- matrix(NA, length(LV), p)
  WSDInumerator2 <- matrix(NA, length(LV), p)
  WUDInumerator1 <- matrix(NA, length(LV), p)
  WUDInumerator2 <- matrix(NA, length(LV), p)
  WSDI <- matrix(NA, p, 1)
  WUDI <- matrix(NA, p, 1)

  for (j in 1:p) {
    for (k in 1:length(LV)) {
      # Calculate difference in expected indicator scores
      DiffExpItem1S[k, j] <- (intercept1[j] - interceptS[j]) + (loading1[j] - loadingS[j]) * LV[k]
      DiffExpItemS2[k, j] <- (interceptS[j] - intercept2[j]) + (loadingS[j] - loading2[j]) * LV[k]
      # probability density functions for sample estimate of latent variable distribution
      pdfLV1[k] <- stats::dnorm(LV[k], mean = fmean1, sd = fsd1)
      pdfLV2[k] <- stats::dnorm(LV[k], mean = fmean2, sd = fsd2)

      # Multiply by latent variable distribution to calculate individual data point in numerator
      WSDInumerator1[k, j] <- DiffExpItem1S[k, j] * pdfLV1[k] * nodewidth
      WSDInumerator2[k, j] <- DiffExpItemS2[k, j] * pdfLV2[k] * nodewidth
      WUDInumerator1[k, j] <- abs(WSDInumerator1[k, j])
      WUDInumerator2[k, j] <- abs(WSDInumerator2[k, j])
    }

    # Sum across range of latent variable using quadrature to calculate numerator
    # Divide by denominator
    WSDI[j] <- prop1 * sum(WSDInumerator1[, j]) / stdev1[j] +
      prop2 * sum(WSDInumerator2[, j]) / stdev2[j]
    WUDI[j] <- prop1 * sum(WUDInumerator1[, j]) / stdev1[j] +
      prop2 * sum(WUDInumerator2[, j]) / stdev2[j]
  }
  ## This is the section that computes dmacs.

  nitems <- length(names(lavaan::lavInspect(fit.cfa, what = "rsquare")[[1]]))

  cfa_minmax <- function(fit.cfa) {
    dt <- lavaan::inspect(fit.cfa, what = "data")
    latentMin <- min(dt[[1]]) - 1
    latentMax <- max(dt[[1]]) + 1
    out <- cbind(as.numeric(latentMin), as.numeric(latentMax))
    return(out)
  }
  reference_load <- lavaan::inspect(fit.cfa, what = "est")[[1]]$lambda
  focal_load <- lavaan::inspect(fit.cfa, what = "est")[[2]]$lambda
  reference_intrcp <- lavaan::inspect(fit.cfa, what = "est")[[1]]$nu
  focal_intrcp <- lavaan::inspect(fit.cfa, what = "est")[[2]]$nu
  pool.sd <- function(fit.cfa) {
    cfa.se <- dplyr::filter(lavaan::partable(fit.cfa), op == "~1", is.na(ustart))
    cfa.n <- lavaan::lavInspect(fit.cfa, what = "nobs")
    l <- list()
    test <- length(names(lavaan::lavInspect(fit.cfa, what = "rsquare")[[1]]))

    for (i in 1:test) {
      grp1 <- dplyr::filter(cfa.se, group == "1")$se[i] * sqrt(cfa.n[1])
      grp2 <- dplyr::filter(cfa.se, group == "2")$se[i] * sqrt(cfa.n[2])
      numerator <- ((cfa.n[1] - 1) * grp1 + (cfa.n[2] -
        1) * grp2)
      denominator <- (cfa.n[1] - 1) + (cfa.n[2] - 1)
      pooled.sd <- numerator / denominator
      l[[paste("item", i)]] <- pooled.sd
    }
    result <- matrix(unlist(l), nrow = test, byrow = TRUE)
    return(result)
  }
  pld_sd <- pool.sd(fit.cfa)
  fcl_lt_vrnc <- lavaan::inspect(fit.cfa, what = "est")[[2]]$psi
  l <- list()
  rowlab <- c()
  for (i in c(1:nitems)) {
    focal.fn <- function(x) {
      mpr <- focal_intrcp[i] + focal_load[i] * x
      return(mpr)
    }
    reference.fn <- function(x) {
      mpr <- reference_intrcp[i] + reference_load[i] *
        x
      return(mpr)
    }
    diff.fn <- function(x, i = i) {
      d <- ((reference.fn(x) - focal.fn(x))^2) * stats::dnorm(x,
        mean = 0, sd = sqrt(fcl_lt_vrnc)
      )
      return(d)
    }
    dMACS <- round((1 / pld_sd[i]) * sqrt(stats::integrate(diff.fn,
      lower = cfa_minmax(fit.cfa)[, 1], upper = cfa_minmax(fit.cfa)[
        ,
        2
      ]
    )$value), 4)
    l[[length(l) + 1]] <- dMACS
    rowlab[[length(rowlab) + 1]] <- paste("Item", i)
  }



  return(data.frame(
    items = items,
    SDI2 = round(SDI2, 4),
    UDI2 = round(UDI2, 4),
    WSDI = round(WSDI, 4),
    WUDI = round(WUDI, 4),
    dmacs = unlist(l)
  ))
}
