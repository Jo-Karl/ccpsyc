iterative_unconstrain <- function(x, exp_p = .05){
  df <- data.frame(Term = NULL,
                   Group1 = NULL,
                   Group2 = NULL,
                   Chi.square = NULL,
                   df = NULL,
                   p.value = NULL)
  x$uni$p.value_adj <- x$uni$p.value
  while (min(x$uni$p.value_adj) <= exp_p) {
    out <- x$uni[which.max(x$uni$Chi.square),]
    x$uni <- x$uni[-which.max(x$uni$Chi.square),]
    df <- rbind(df, out)
    x$uni$p.value_adj <- x$uni$p.value * (nrow(df) + 1)
  }
  return(df)
}
