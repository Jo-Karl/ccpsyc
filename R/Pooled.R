#' Creating a Pan-Cultural Loading Matrix
#'
#' @param df A data frame contains the variables for the exploratory factor
#'           analysis and the grouping variable.
#' @param group The name of the column tht cointains the grouping supplied as a
#'              string.
#' @param nfactors The number of factors expected.
#'
#' @return returns a Pan-Cultural loading matrix.
#' @export pancultural
#'
#' @examples
#' pancultural(example, "country", 5)
pancultural <- function(df, group, nfactors) {
  if (nfactors > ncol(df) - 1) {
    print("Please specify a lower number of dimensions")
  } else {
    df[[group]] <- as.numeric(factor(df[[group]]))
    PAN <- psych::statsBy(
      data = df, group = paste(group),
      cors = TRUE, cor = "cor", method = "pearson"
    )
    PAN_PCA <- psych::principal(PAN$pooled,
      nfactors = paste(nfactors),
      rotate = "varimax"
    )
    return(PAN_PCA$loadings)
  }
}
