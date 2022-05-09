#' Get more comprehensible output from lavTestScore
#'
#' @param lavaan.fit Model fitted with lavaan
#' @param ndigit Defines the rounding
#' @param ... Arguments passed to lavTestScore
#' @author Maksim Rudnev
#' @return Returns a dataframe which contains one row for each constrained parameter in the model together with a chi-square test indicating whether the parameter significantly differs between groups.
#' This is a cleaned version identical to \code{\link[lavaan:lavTestScore]{lavTestScore}}.
#' @export lavTestScore.clean

lavTestScore.clean <- function(lavaan.fit, ndigit = 3, ...) {
  plabel <- NULL
  lvts <- lavaan::lavTestScore(lavaan.fit, ...)

  for (lvts.part in names(lvts)[names(lvts) %in% c("uni", "cumulative")]) {
    partab.a <- lavaan::partable(lavaan.fit)[, c(c("lhs", "op", "rhs", "group", "plabel"))] %>%
      dplyr::filter(plabel != "")

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


  return(lvts)
}
