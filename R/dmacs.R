#' Computes dMACS
#'
#' @param fit.cfa Lavaan output object with two groups and a single
#' factor.
#' @param group1 String for first group in the grouping factor
#' @param group2 String for second group in the grouping factor
#'
#' @return Returns dMACS for each item.
#' @importFrom rlang .data
#' @export dMACS
#'
#' @examples
#' dMACS
dMACS <- function(fit.cfa, group1, group2) {

  rsquare_data <- lavaan::lavInspect(fit.cfa, what = "rsquare")
  nitems <- length(names(rsquare_data[[1]]))

  cfa_minmax <- function(fit.cfa) {
    dt <- lavaan::inspect(fit.cfa, what = "data")
    latent_min <- min(dt[[1]]) - 1
    latent_max <- max(dt[[1]]) + 1
    return(cbind(as.numeric(latent_min), as.numeric(latent_max)))
  }

  est_data <- lavaan::inspect(fit.cfa, what = "est")
  reference_load <- est_data[[1]]$lambda
  focal_load <- est_data[[2]]$lambda
  reference_intrcp <- est_data[[1]]$nu
  focal_intrcp <- est_data[[2]]$nu

  pool_sd <- function(fit.cfa) {
    cfa_se <- lavaan::lavInspect(fit.cfa, what = "se")
    cfa_n <- lavaan::lavInspect(fit.cfa, what = "nobs")

    result_list <- vector("list", nitems)

    for (i in seq_len(nitems)) {
      grp1_sd <- cfa_se[[group1]]$nu[i] * sqrt(cfa_n[1])
      grp2_sd <- cfa_se[[group2]]$nu[i] * sqrt(cfa_n[2])
      numerator <- (cfa_n[1] - 1) * grp1_sd + (cfa_n[2] - 1) * grp2_sd
      denominator <- (cfa_n[1] - 1) + (cfa_n[2] - 1)
      pooled_sd <- numerator / denominator
      result_list[[i]] <- pooled_sd
    }

    return(matrix(unlist(result_list), nrow = nitems, byrow = TRUE))
  }

  pooled_sd <- pool_sd(fit.cfa)
  focal_latent_variance <- est_data[[2]]$psi
  minmax_values <- cfa_minmax(fit.cfa)

  dmacs_results <- vector("list", nitems)
  row_labels <- character(nitems)

  for (i in seq_len(nitems)) {
    focal_fn <- function(x) {
      return(focal_intrcp[i] + focal_load[i] * x)
    }

    reference_fn <- function(x) {
      return(reference_intrcp[i] + reference_load[i] * x)
    }

    diff_fn <- function(x) {
      diff_val <- (reference_fn(x) - focal_fn(x))^2
      normal_density <- stats::dnorm(x, mean = 0, sd = sqrt(focal_latent_variance))
      return(diff_val * normal_density)
    }

    integration_result <- stats::integrate(
      diff_fn,
      lower = minmax_values[1],
      upper = minmax_values[2]
    )

    dmacs_value <- round(
      (1 / pooled_sd[i]) * sqrt(integration_result$value),
      3
    )

    dmacs_results[[i]] <- dmacs_value
    row_labels[i] <- paste("Item", i)
  }

  result_matrix <- matrix(
    unlist(dmacs_results),
    nrow = nitems,
    dimnames = list(row_labels, "dMAC")
  )

  return(result_matrix)
}
