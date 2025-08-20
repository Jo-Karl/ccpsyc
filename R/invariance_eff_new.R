
####################################################################################
################# Invariance Effect Sizes############################################
################## Johannes A. Karl #################################################
################### Based on code by Heather Gunn####################################
# https://www.tandfonline.com/doi/full/10.1080/10705511.2019.1689507?journalCode=hsem20
####################################################################################

#' Calculate Differential Item Functioning Indices (This will slowly replace invariance_eff)
#'
#' Computes Signed and Unsigned Differential Item Functioning (SDIF/UDIF) indices
#' along with differential MACS (dMACS) for multi-group confirmatory factor analysis.
#'
#' @param df A data frame containing the item responses and grouping variable
#' @param items Character vector of item column names to include in the analysis
#' @param group Character string specifying the grouping variable column name
#' @param nodewidth Numeric value for numerical integration step size. Default is 0.01.
#'   \strong{Tips for setting nodewidth:}
#'   \itemize{
#'     \item Smaller values (0.001-0.01) = Higher accuracy, slower computation
#'     \item Larger values (0.02-0.1) = Lower accuracy, faster computation
#'     \item For most applications, 0.01 provides good balance
#'     \item Use smaller values (0.005) for high-stakes decisions
#'     \item Use larger values (0.02-0.05) for exploratory analyses
#'   }
#' @param intercept_fix Integer specifying which item intercept to fix to 0 for identification. Default is 1
#' @param lower_lv Numeric value for lower bound of latent variable range. Default is -10.
#'   \strong{Tips for setting lower_lv:}
#'   \itemize{
#'     \item Should cover the practical range of latent variable values
#'     \item Too narrow = Missing important DIF effects in tails
#'     \item Too wide = Unnecessary computation
#'     \item For typical applications: -4 to -6 covers ~99.9% of normal distribution
#'     \item For extreme groups or highly variable data: -8 to -12
#'     \item Consider your data's latent variable distribution from initial CFA
#'   }
#' @param upper_lv Numeric value for upper bound of latent variable range. Default is 10.
#'   \strong{Tips for setting upper_lv:}
#'   \itemize{
#'     \item Should mirror lower_lv for symmetric coverage
#'     \item For typical applications: +4 to +6 covers ~99.9% of normal distribution
#'     \item For extreme groups or highly variable data: +8 to +12
#'     \item Ensure range covers both groups' latent variable distributions
#'     \item Check factor score ranges from preliminary analysis to guide choice
#'   }
#' @param reference_group Character string or numeric value specifying the reference group.
#'   If NULL, uses the first level of the grouping variable
#' @param focal_group Character string or numeric value specifying the focal group.
#'   If NULL, uses the second level of the grouping variable
#' @param convergence_check Logical indicating whether to check model convergence. Default is TRUE
#' @param ... Additional arguments passed to lavaan::cfa()
#'
#' @return A data frame with columns:
#'   \item{items}{Item names}
#'   \item{SDI2}{Signed Differential Item Functioning index for focal group}
#'   \item{UDI2}{Unsigned Differential Item Functioning index for focal group}
#'   \item{WSDI}{Weighted Signed Differential Item Functioning index}
#'   \item{WUDI}{Weighted Unsigned Differential Item Functioning index}
#'   \item{dmacs}{Differential MACS values}
#'
#' @details
#' \strong{Computational Considerations:}
#' The computational time is approximately proportional to:
#' \code{(upper_lv - lower_lv) / nodewidth * number_of_items}
#'
#' For a typical analysis with 10 items:
#' \itemize{
#'   \item Fast setting: range=[-4,4], nodewidth=0.02 → ~4,000 calculations
#'   \item Balanced setting: range=[-6,6], nodewidth=0.01 → ~12,000 calculations
#'   \item Precise setting: range=[-8,8], nodewidth=0.005 → ~32,000 calculations
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with default settings
#' results <- calculate_dif_indices(data, items = c("item1", "item2", "item3"),
#'                                  group = "group_var")
#'
#' # High precision analysis (slower but more accurate)
#' results_precise <- calculate_dif_indices(data, items = c("item1", "item2", "item3"),
#'                                          group = "group_var", nodewidth = 0.005,
#'                                          lower_lv = -6, upper_lv = 6)
#'
#' # Fast exploratory analysis
#' results_fast <- calculate_dif_indices(data, items = c("item1", "item2", "item3"),
#'                                       group = "group_var", nodewidth = 0.02,
#'                                       lower_lv = -4, upper_lv = 4)
#'
#' # Custom group specification
#' results_custom <- calculate_dif_indices(data, items = c("item1", "item2", "item3"),
#'                                         group = "group_var",
#'                                         reference_group = "control",
#'                                         focal_group = "treatment")
#' }
#'
#' @export
calculate_dif_indices <- function(df, items, group, nodewidth = 0.01, intercept_fix = 1,
                                  lower_lv = -10, upper_lv = 10, reference_group = NULL,
                                  focal_group = NULL, convergence_check = TRUE, ...) {

  # Declare global variables to avoid CRAN notes
  op <- ustart <- label <- est <- se <- NULL

  # Input validation
  if (!is.data.frame(df)) {
    stop("df must be a data frame")
  }
  if (!all(items %in% names(df))) {
    stop("All items must be column names in df")
  }
  if (!group %in% names(df)) {
    stop("group must be a column name in df")
  }
  if (intercept_fix < 1 || intercept_fix > length(items)) {
    stop("intercept_fix must be between 1 and the number of items")
  }

  # Set up group levels and handle filtering
  all_group_levels <- unique(df[[group]])

  # If reference_group and focal_group are specified, filter data to only include those groups
  if (!is.null(reference_group) && !is.null(focal_group)) {
    # Validate that specified groups exist in the data
    if (!reference_group %in% all_group_levels) {
      stop("reference_group '", reference_group, "' not found in grouping variable")
    }
    if (!focal_group %in% all_group_levels) {
      stop("focal_group '", focal_group, "' not found in grouping variable")
    }
    if (reference_group == focal_group) {
      stop("reference_group and focal_group must be different")
    }

    # Filter data to only include the two specified groups
    df <- df[df[[group]] %in% c(reference_group, focal_group), ]

    # Remove any factor levels that are no longer present
    if (is.factor(df[[group]])) {
      df[[group]] <- droplevels(df[[group]])
    }

    group_levels <- c(reference_group, focal_group)

  } else {
    # No specific groups provided - use all available groups
    group_levels <- sort(unique(df[[group]]))

    # Check that we have exactly 2 groups
    if (length(group_levels) != 2) {
      stop("Grouping variable must have exactly 2 levels, or specify reference_group and focal_group. ",
           "Found groups: ", paste(all_group_levels, collapse = ", "))
    }

    # Set default reference and focal groups
    reference_group <- group_levels[1]
    focal_group <- group_levels[2]
  }

  # Provide integration parameter guidance
  lv_range <- upper_lv - lower_lv
  n_integration_points <- lv_range / nodewidth
  if (n_integration_points > 10000) {
    warning("Large number of integration points (", round(n_integration_points),
            "). Consider increasing nodewidth or reducing latent variable range for faster computation.")
  }
  if (lv_range < 6) {
    warning("Narrow latent variable range (", lv_range,
            "). Consider expanding range to capture tail effects if groups have extreme scores.")
  }

  # Build CFA model specifications
  sum_item <- paste(items, collapse = " + ")
  main_model <- paste("f1", sum_item, sep = " =~ ")
  fact_intercept <- "f1 ~ c(ref_int, foc_int) * 1"
  item_int_zero <- paste(items[intercept_fix], "~ 0 * 1")
  cfa_model_with_intercept <- paste(main_model, fact_intercept, sep = "\n")
  cfa_model_final <- paste(cfa_model_with_intercept, item_int_zero, sep = "\n")

  # Fit CFA models
  tryCatch({
    fit_cfa <- lavaan::cfa(main_model, data = df, group = group, std.lv = TRUE, ...)
    fit_pooled <- lavaan::cfa(main_model, data = df, std.lv = TRUE, meanstructure = TRUE, ...)
    fit_configural <- lavaan::cfa(main_model, data = df, group = group, std.lv = TRUE, ...)
    fit_intercept <- lavaan::cfa(cfa_model_final, data = df, group = group,
                                 std.lv = TRUE, meanstructure = TRUE, ...)
  }, error = function(e) {
    stop("Model fitting failed: ", e$message)
  })

  # Check convergence if requested
  if (convergence_check) {
    models <- list(fit_cfa, fit_pooled, fit_configural, fit_intercept)
    model_names <- c("fit_cfa", "fit_pooled", "fit_configural", "fit_intercept")

    for (i in seq_along(models)) {
      if (!lavaan::lavInspect(models[[i]], "converged")) {
        stop("Model ", model_names[i], " did not converge")
      }
    }
  }

  # Extract parameter tables
  par_configural <- lavaan::parametertable(fit_configural)
  par_pooled <- lavaan::parametertable(fit_pooled)
  par_intercept <- lavaan::parametertable(fit_intercept)

  # Get number of items
  n_items <- length(items)

  # Map group labels to lavaan's internal group numbers
  group_mapping <- lavaan::lavInspect(fit_configural, "group.label")

  # Find which internal group number corresponds to our reference and focal groups
  ref_group_num <- which(group_mapping == reference_group)
  foc_group_num <- which(group_mapping == focal_group)

  if (length(ref_group_num) == 0) {
    stop("Reference group '", reference_group, "' not found in fitted model")
  }
  if (length(foc_group_num) == 0) {
    stop("Focal group '", focal_group, "' not found in fitted model")
  }

  # Convert to character for parameter table filtering
  ref_group_label <- as.character(ref_group_num)
  foc_group_label <- as.character(foc_group_num)

  # Extract parameters using base R
  loading_ref <- par_configural[par_configural$group == ref_group_label &
                                  par_configural$op == "=~", "est"]
  loading_foc <- par_configural[par_configural$group == foc_group_label &
                                  par_configural$op == "=~", "est"]
  loading_pooled <- par_pooled[par_pooled$op == "=~", "est"]

  intercept_ref <- par_configural[par_configural$group == ref_group_label &
                                    par_configural$op == "~1" &
                                    is.na(par_configural$ustart), "est"]
  intercept_foc <- par_configural[par_configural$group == foc_group_label &
                                    par_configural$op == "~1" &
                                    is.na(par_configural$ustart), "est"]
  intercept_pooled <- par_pooled[par_pooled$op == "~1" &
                                   is.na(par_pooled$ustart), "est"]

  # Extract standard errors and sample sizes
  se_ref <- par_configural[par_configural$group == ref_group_label &
                             par_configural$op == "~1" &
                             is.na(par_configural$ustart), "se"]
  se_foc <- par_configural[par_configural$group == foc_group_label &
                             par_configural$op == "~1" &
                             is.na(par_configural$ustart), "se"]

  n_obs <- lavaan::lavInspect(fit_configural, "nobs")
  n_ref <- n_obs[[ref_group_num]]
  n_foc <- n_obs[[foc_group_num]]

  stdev_ref <- round(se_ref * sqrt(n_ref), 3)
  stdev_foc <- round(se_foc * sqrt(n_foc), 3)

  # Extract factor means and standard deviations
  fmean_ref <- par_intercept[par_intercept$group == ref_group_label &
                               par_intercept$label == "ref_int", "est"]
  fmean_foc <- par_intercept[par_intercept$group == foc_group_label &
                               par_intercept$label == "foc_int", "est"]

  fsd_ref <- par_intercept[par_intercept$group == ref_group_label &
                             par_intercept$label == "ref_int", "se"] * sqrt(n_ref)
  fsd_foc <- par_intercept[par_intercept$group == foc_group_label &
                             par_intercept$label == "foc_int", "se"] * sqrt(n_foc)

  # Calculate proportions
  prop_ref <- n_ref / (n_ref + n_foc)
  prop_foc <- n_foc / (n_ref + n_foc)

  # Define latent variable evaluation points
  lv_points <- seq(lower_lv, upper_lv, nodewidth)
  n_points <- length(lv_points)

  # Calculate DIF indices using vectorized operations
  sdi2 <- numeric(n_items)
  udi2 <- numeric(n_items)
  wsdi <- numeric(n_items)
  wudi <- numeric(n_items)

  for (j in seq_len(n_items)) {
    # Calculate expected score differences
    diff_exp_ref_foc <- (intercept_ref[j] - intercept_foc[j]) +
      (loading_ref[j] - loading_foc[j]) * lv_points
    diff_exp_ref_pool <- (intercept_ref[j] - intercept_pooled[j]) +
      (loading_ref[j] - loading_pooled[j]) * lv_points
    diff_exp_pool_foc <- (intercept_pooled[j] - intercept_foc[j]) +
      (loading_pooled[j] - loading_foc[j]) * lv_points

    # Calculate probability density functions
    pdf_ref <- stats::dnorm(lv_points, mean = fmean_ref, sd = fsd_ref)
    pdf_foc <- stats::dnorm(lv_points, mean = fmean_foc, sd = fsd_foc)

    # Calculate numerators using numerical integration
    sdi2_num <- sum(diff_exp_ref_foc * pdf_foc * nodewidth)
    udi2_num <- sum(abs(diff_exp_ref_foc) * pdf_foc * nodewidth)

    wsdi_num1 <- sum(diff_exp_ref_pool * pdf_ref * nodewidth)
    wsdi_num2 <- sum(diff_exp_pool_foc * pdf_foc * nodewidth)
    wudi_num1 <- sum(abs(diff_exp_ref_pool) * pdf_ref * nodewidth)
    wudi_num2 <- sum(abs(diff_exp_pool_foc) * pdf_foc * nodewidth)

    # Calculate final indices
    sdi2[j] <- sdi2_num / stdev_foc[j]
    udi2[j] <- udi2_num / stdev_foc[j]
    wsdi[j] <- prop_ref * wsdi_num1 / stdev_ref[j] + prop_foc * wsdi_num2 / stdev_foc[j]
    wudi[j] <- prop_ref * wudi_num1 / stdev_ref[j] + prop_foc * wudi_num2 / stdev_foc[j]
  }

  # Calculate dMACS indices
  dmacs_values <- calculate_dmacs_internal(fit_cfa, ref_group_num, foc_group_num)

  # Return results
  data.frame(
    items = items,
    SDI2 = round(sdi2, 4),
    UDI2 = round(udi2, 4),
    WSDI = round(wsdi, 4),
    WUDI = round(wudi, 4),
    dmacs = round(dmacs_values, 4),
    stringsAsFactors = FALSE
  )
}

#' Internal function to calculate dMACS values
#' @param fit_cfa Fitted lavaan CFA object
#' @param ref_group_num Numeric reference group identifier (1 or 2)
#' @param foc_group_num Numeric focal group identifier (1 or 2)
#' @return Vector of dMACS values
#' @keywords internal
calculate_dmacs_internal <- function(fit_cfa, ref_group_num, foc_group_num) {

  # Get number of items
  rsquare_data <- lavaan::lavInspect(fit_cfa, what = "rsquare")
  n_items <- length(names(rsquare_data[[1]]))

  # Calculate min/max values for integration
  dt <- lavaan::inspect(fit_cfa, what = "data")
  latent_min <- min(dt[[1]]) - 1
  latent_max <- max(dt[[1]]) + 1

  # Extract parameter estimates using numeric group indices
  est_data <- lavaan::inspect(fit_cfa, what = "est")
  ref_loadings <- est_data[[ref_group_num]]$lambda
  foc_loadings <- est_data[[foc_group_num]]$lambda
  ref_intercepts <- est_data[[ref_group_num]]$nu
  foc_intercepts <- est_data[[foc_group_num]]$nu

  # Calculate pooled standard deviations
  param_table <- lavaan::parametertable(fit_cfa)
  param_se <- param_table[param_table$op == "~1" & is.na(param_table$ustart), ]

  n_obs <- lavaan::lavInspect(fit_cfa, "nobs")
  n_ref <- n_obs[[ref_group_num]]
  n_foc <- n_obs[[foc_group_num]]

  # Convert group numbers to character for parameter table filtering
  ref_group_char <- as.character(ref_group_num)
  foc_group_char <- as.character(foc_group_num)

  pooled_sd <- numeric(n_items)
  for (i in seq_len(n_items)) {
    se_ref_item <- param_se[param_se$group == ref_group_char, "se"][i] * sqrt(n_ref)
    se_foc_item <- param_se[param_se$group == foc_group_char, "se"][i] * sqrt(n_foc)

    numerator <- (n_ref - 1) * se_ref_item + (n_foc - 1) * se_foc_item
    denominator <- (n_ref - 1) + (n_foc - 1)
    pooled_sd[i] <- numerator / denominator
  }

  # Extract focal group latent variable variance
  foc_latent_var <- est_data[[foc_group_num]]$psi

  # Calculate dMACS for each item
  dmacs_values <- numeric(n_items)

  for (i in seq_len(n_items)) {
    # Define item response functions
    ref_fn <- function(x) ref_intercepts[i] + ref_loadings[i] * x
    foc_fn <- function(x) foc_intercepts[i] + foc_loadings[i] * x

    # Define difference function for integration
    diff_fn <- function(x) {
      diff_val <- (ref_fn(x) - foc_fn(x))^2
      normal_density <- stats::dnorm(x, mean = 0, sd = sqrt(foc_latent_var))
      return(diff_val * normal_density)
    }

    # Perform numerical integration
    integration_result <- stats::integrate(diff_fn, lower = latent_min, upper = latent_max)
    dmacs_values[i] <- (1 / pooled_sd[i]) * sqrt(integration_result$value)
  }

  return(dmacs_values)
}
