#' Procrustes rotation function, returning Tucker's Phi
#'
#' @param loading A correlation matrix to be rotated towards a target
#' @param norm    A correlation matrix that is the goal of the rotation
#' @param rotated A TRUE/FALSE operator indicating if the rotated matrix should
#'                be returned in addition to Tucker's Phi
#' @param digits The number of digits to be displayed in the output, defaults to 2
#'
#' @return Returns Tuckers Phi evaluating the congruence of the loading matrix
#'         to the normative matrix
#' @export prost
prost <- function(loading, norm, rotated = FALSE, digits = 2) {
  if (rotated == TRUE) {
    rotated <- MCMCpack::procrustes(loading, norm)
    nx <- dim(rotated$X.new)[2]
    ny <- dim(norm)[2]
    cross <- t(norm) %*% rotated$X.new
    sumsx <- sqrt(1 / diag(t(rotated$X.new) %*% rotated$X.new))
    sumsy <- sqrt(1 / diag(t(norm) %*% norm))
    result <- matrix(rep(0, nx * ny), ncol = nx)
    result <- round(sumsy * (cross * rep(sumsx, each = ny)), digits)
    congruence <- diag(t(result))
    lin.corr <- stats::cor(rotated$X.new, norm)
    lin <- round(lin.corr[col(lin.corr) == row(lin.corr)], digits)
    congruence.list <- list(
      rotated.matrix = rotated, tuckers.phi = congruence,
      correlation = lin
    )
  } else {
    rotated <- MCMCpack::procrustes(loading, norm)
    nx <- dim(rotated$X.new)[2]
    ny <- dim(norm)[2]
    cross <- t(norm) %*% rotated$X.new
    sumsx <- sqrt(1 / diag(t(rotated$X.new) %*% rotated$X.new))
    sumsy <- sqrt(1 / diag(t(norm) %*% norm))
    result <- matrix(rep(0, nx * ny), ncol = nx)
    result <- round(sumsy * (cross * rep(sumsx, each = ny)), digits)
    congruence <- diag(t(result))
    lin.corr <-  stats::cor(rotated$X.new, norm)
    lin <- round(lin.corr[col(lin.corr) == row(lin.corr)], digits)
    congruence.list <- list(
      tuckers.phi = congruence,
      correlation = lin
    )
  }
  return(congruence.list)
}
