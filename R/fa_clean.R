#' Function to quickly organize and clear psych factor loadings
#'
#' @param psych_fa Output from the psych package, can be either fa or principal with at least two dimensions
#' @param cutoff Desired cutoff below which loadings are omitted defaults to .40
#' @param dbl_dist Desired distance between highest and second highest loading for an item to remove double loadings, defaults to .20
#' @param key_file Optional: Either a .csv or .xlsx file with at least two columns:
#'  1 labeled item containing the item labels as in the data frame,
#'  2 a column labeled wording containing the item wording.
#' @param cleaned If true (default), only the cleaned solution with a message for descriptive stats are returned.
#'  If false the function returns a list of data frames one cleaned and one showing all in-between steps
#'
#' @return clean This column contains the assignment after removing NAs and double loadings
#' @return dir This column contains the direction (positive or negative) of the highest loading.
#' @export clearing_fa
#' @importFrom rlang .data
#' @examples
#' library(psych)
#' fa_solution <- fa(example[c(paste0("help", 1:6, "m"), c(paste0("voice", 1:5, "m")))], nfactors = 2)
#' clearing_fa(fa_solution)
clearing_fa <-
  function(psych_fa, cutoff = .40, dbl_dist = .20, key_file = NULL, cleaned = TRUE) {
    ### This section declares variables globally to avoid CRAN nots about non-binded vars
    clean <- X1 <- X2 <- NULL
    loadings <-
      tibble::rownames_to_column(data.frame(unclass(psych_fa$loadings)), "item")
    loadings_s <- loadings


    loadings$auto <- colnames(abs(loadings[2:(ncol(loadings_s))]))[apply(abs(loadings[2:(ncol(loadings_s))]), 1, which.max)]




    loadings$auto[!apply(abs(loadings[2:(ncol(loadings_s))]), 1, max) >= cutoff] <- "Low"


    ordered <- apply(abs(loadings[2:(ncol(loadings_s))]), 1, function(x) {
        sort(t(x), decreasing = T)
    }) %>%
      t() %>%
      data.frame() %>%
      dplyr::mutate(distance = X1 - X2)

    df_full <- cbind(loadings, ordered)

    df_full$clean <- dplyr::if_else(df_full$distance >= dbl_dist, df_full$auto, "Double")
    n_na <- sum((df_full$auto == "Low"))
    n_double <- sum(df_full$clean == "Double")

    df_full$dir <- apply(df_full[2:(ncol(loadings_s))], 1, function(x) dplyr::if_else(x[which.max(abs(x))] < 0, "neg", "pos"))
    df_full$max_load <- apply(df_full[2:(ncol(loadings_s))], 1, function(x) max(abs(x)))
    if (!is.null(key_file) & is.character(key_file)) {
      if (grepl(".xlsx$", key_file)) {
        key_xlsx <- xlsx::read.xlsx(key_file, sheetIndex = 1)
        df_full <- dplyr::left_join(key_xlsx, df_full, by = "item")
      } else if (grepl(".csv$", key_file)) {
        key_csv <- readr::read_csv(key_file)
        df_full <- dplyr::left_join(key_csv, df_full, by = "item")
      } else {
        stop("Please either use a .csv or .xlsx file")
      }
    }
    if (!is.null(key_file) & is.character(key_file)) {
      df_load <- df_full[c("item", "wording", "clean", "dir", "max_load")]
    } else {
      df_load <- df_full[c("item", "clean", "dir", "max_load")]
    }
    if (cleaned) {
      df_load$clean[df_load$clean %in% c("Double", "Low")] <- NA
      cleaned_df <- dplyr::filter(df_load, !is.na(clean))
      message(paste0(
        "Of the ", nrow(df_load), " items ", nrow(cleaned_df), " (", round((nrow(cleaned_df) / nrow(df_load)) * 100, 2), "%) ", "were retained. ", n_na, " (", round((n_na / nrow(df_load)) * 100, 2), "%) ", "showed no substantative loadings and ",
        n_double, " (", round((n_double / nrow(df_load)) * 100, 2), "%) ", "showed double loadings"
      ))
      return(cleaned_df)
    } else {
      l_df <- list(
        filterd = dplyr::filter(df_load, !is.na(clean)),
        all = df_full
      )
      return(l_df)
    }
  }
