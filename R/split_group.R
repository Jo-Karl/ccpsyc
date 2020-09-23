#' Split by groups
#'
#' @param df Dataframe
#' @param group Variable from the dataset that defines the groups
#' @param named TRUE/FALSE argument wheter the resulting list should be named
#' @param name.list Supply a list of names same length as number of groups
#' @return Returns a list of dataframes with only the selected groups
#' @export
splitgroup <- function(df, group, named = FALSE, name.list = NA) {
  A <- split(df, df[, paste(group)])
  B <- lapply(seq_along(A), function(x) {
    as.data.frame(A[[x]])[, paste(df[-which(names(df) %in% c(paste(group)))])]
  })
  if (named == TRUE) {
    names(B) <- name.list
  }
  return(B)
}
