#' Multi-group reliability table
#'
#' @param df_s The full dataframe with all groups and items.
#' @param measure_list A named list of vectors containing the item names.
#' The format should be list(measure_name1 = c('Item1', 'Item2', 'Item3'), measure_name2 = c('Item1', 'Item2', 'Item3'))
#' @param group Grouping variable in the dataset as string for example "country"
#' @param digitn Controls the amount of digits shown in the output
#' @param seed Seed for the bootstrapped confidence intervals
#' @export mg_rel_table
#' @return Returns a formatted dataframe with the reliability of all constructs by group


mg_rel_table <- function(df_s, measure_list, group, digitn = 3, seed = 2711) {
  warning("To use this function, please install a version of userfriednlyscience from the CRAN archive")
  if (length(names(measure_list)) == 0) {
    stop("Please provide the measure list as named list.
         The format should be list(measure_name1 = c('Item1', 'Item2', 'Item3'),
         measure_name2 = c('Item1', 'Item2', 'Item3'))")
  }
  df_split <- split.data.frame(df_s, df_s[[group]])
  rel_outer <-
    lapply(measure_list, function(y) {
      rel_list <-
        lapply(df_split, function(x) {
          inter <- ufs::scaleStructure(x,
            items = y,
            digits = digitn,
            interval.type = "normal-theory",
            bootstrapSeed = seed
          )
          out <- dplyr::tibble(
            alpha = paste0(
              format(round(inter$output$cronbach.alpha, digitn), nsmall = digitn),
              "[",
              format(round(inter$output$alpha.ci[1], digitn), nsmall = digitn),
              ", ",
              format(round(inter$output$alpha.ci[2], digitn), nsmall = digitn), "]"
            ),
            omega = paste0(
              format(round(inter$output$omega, digitn), nsmall = digitn),
              "[",
              format(round(inter$output$omega.ci[1], digitn), nsmall = digitn),
              ", ",
              format(round(inter$output$omega.ci[2], digitn), nsmall = digitn),
              "]"
            ),
            glb = paste0(format(round(inter$output$glb, digitn), nsmall = digitn)),
            h = paste0(format(round(inter$output$coefficientH, digitn), nsmall = digitn)),
            percent_pos = inter$intermediate$cor.proPos * 100
          )
          ### Removing leading 0s in the crudest possible way.
          out$alpha <- gsub(x = out$alpha, pattern = "0\\.", replacement = "\\.")
          out$omega <- gsub(x = out$omega, pattern = "0\\.", replacement = "\\.")
          out$glb <- gsub(x = out$glb, pattern = "0\\.", replacement = "\\.")
          out$h <- gsub(x = out$h, pattern = "0\\.", replacement = "\\.")
          out
        })
      do.call(rbind, rel_list) %>%
        cbind(country = names(df_split), .data)
    })
  joined <-
    do.call(rbind, rel_outer) %>%
    cbind(measure = rep(names(measure_list), each = length(names(df_split))), .data)
  rownames(joined) <- NULL
  joined
}
