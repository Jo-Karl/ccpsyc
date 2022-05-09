#' Examening chisquare improvement if paths are unconstrained.
#' The function returns the paths to be unconstrained based on chisquare change.
#' Adjusted P-values are calculated based on iterative Bonferroni corrections.
#' @param lavaan.fit Model fitted with lavaan
#' @param ... Arguments passed to lavTestScore
#' @param ndigit Number of digits to round chi and p to
#' @param exp_p Expected p-value
#' @author Maksim Rudnev
#' @export release_bonferroni
#' @return Returns a dataframe representing a Bonferroni corrected version of \code{\link{lavTestScore.clean}}.
#' @importFrom rlang .data

release_bonferroni <- function(lavaan.fit, ndigit = 3, exp_p = .05, ...) {
  lvts <- lavaan::lavTestScore(lavaan.fit, ...)
  plabel <- NULL
  for (lvts.part in names(lvts)[names(lvts) %in% c("uni", "cumulative")]) {
    partab.a <- lavaan::partable(lavaan.fit)[, c(c("lhs", "op", "rhs", "group", "plabel"))] %>%
      dplyr::filter(.data, plabel != "")

    names(partab.a)[1:3] <- c("one", "two", "three")

    out <- merge(as.data.frame(lvts[[lvts.part]]),
      partab.a,
      by.x = c("lhs"), by.y = c("plabel"),
      all.x = T
    )
    out2 <- merge(out,
      partab.a,
      by.x = c("rhs"), by.y = c("plabel"),
      all.x = T, suffixes = c(".lhs", ".rhs")
    )

    out2$group.lhs <- factor(out2$group.lhs, levels = 1:length(lavaan::lavInspect(lavaan.fit, "group.label")), labels = lavaan::lavInspect(lavaan.fit, "group.label"))
    out2$group.rhs <- factor(out2$group.rhs, levels = 1:length(lavaan::lavInspect(lavaan.fit, "group.label")), labels = lavaan::lavInspect(lavaan.fit, "group.label"))

    out3 <- data.frame(
      Term = paste(out2$one.lhs, out2$two.lhs, out2$three.lhs, sep = ""),
      Group1 = out2$group.lhs,
      Group2 = out2$group.rhs,
      Chi.square = round(out2$X2, ndigit), df = out2$df, p.value = round(out2$p.value, ndigit),
      "." = format(as.character(sapply(out2$p.value, function(x) ifelse(x > 0.05, "", ifelse(x > 0.01, "*", ifelse(x > 0.001, "**", "***"))))), justify = "left")
    )

    lvts[[lvts.part]] <- out3
    if (lvts.part == "uni") attr(lvts[[lvts.part]], "header") <- "Chi-square test of releasing single constraints, equivalent to modification indices"
    if (lvts.part == "cumulative") attr(lvts[[lvts.part]], "header") <- "Chi-square test of releasing multiple constraints at the same time"
    class(lvts[[lvts.part]]) <- c("lavaan.data.frame", "data.frame")
  }


  if (any(names(lvts) == c("epc"))) {
    lvts[["epc"]]$group <- factor(lvts[["epc"]]$group,
      levels = 1:length(lavaan::lavInspect(lavaan.fit, "group.label")),
      labels = lavaan::lavInspect(lavaan.fit, "group.label")
    )
  }




  df <- data.frame(
    Term = NULL,
    Group1 = NULL,
    Group2 = NULL,
    Chi.square = NULL,
    df = NULL,
    p.value = NULL
  )
  lvts$uni$p.value_adj <- lvts$uni$p.value
  while (min(lvts$uni$p.value_adj) <= exp_p) {
    out <- lvts$uni[which.max(lvts$uni$Chi.square), ]
    lvts$uni <- lvts$uni[-which.max(lvts$uni$Chi.square), ]
    df <- rbind(df, out)
    lvts$uni$p.value_adj <- lvts$uni$p.value * (nrow(df) + 1)
  }
  return(df)
}
