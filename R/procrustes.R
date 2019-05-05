#' Procrustes rotation function, returning Tucker's Phi
#'
#' @param loading A correlation matrix to be rotated towards a target
#' @param norm    A correlation matrix that is the goal of the rotation
#' @param rotated A TRUE/FALSE operator indicating if the rotated matrix should
#'                be returned in addition to Tucker's Phi
#'
#' @return Returns Tuckers Phi evaluating the congruence of the loading matrix
#'         to the normative matrix
#' @import MCMCpack
#' @import psych
#' @export
prost <- function(loading, norm, rotated = FALSE, MCMC = TRUE) {
  if (MCMC == TRUE){
    if (rotated == TRUE){
      rotated <- MCMCpack::procrustes(loading, norm)
      cong.corr <- psych::factor.congruence(rotated$X.new, norm)
      congruence <- cong.corr[col(cong.corr) == row(cong.corr)]
      lin.corr <- cor(rotated$X.new, norm)
      lin <- round(lin.corr[col(lin.corr) == row(lin.corr)],2)
      congruence.list <- list(rotated.matrix = rotated, tuckers.phi = congruence,
                              correlation = lin)
    }else{
      rotated <- MCMCpack::procrustes(loading, norm)
      cong.corr <- psych::factor.congruence(rotated$X.new, norm)
      congruence <- cong.corr[col(cong.corr) == row(cong.corr)]
      lin.corr <- cor(rotated$X.new, norm)
      lin <- round(lin.corr[col(lin.corr) == row(lin.corr)],2)
      congruence.list <- list(tuckers.phi = congruence,
                              correlation = lin)
    }
  }else{
    if (rotated == TRUE){
      rotated <- psych::TargetQ(loading, Target = list(norm))
      cong.corr <- psych::factor.congruence(rotated$loadings, norm)
      congruence <- cong.corr[col(cong.corr) == row(cong.corr)]
      lin.corr <- cor(rotated$loadings, norm)
      lin <- round(lin.corr[col(lin.corr) == row(lin.corr)],2)
      congruence.list <- list(rotated.matrix = rotated, tuckers.phi = congruence,
                              correlation = lin)
    }else{
      rotated <- psych::TargetQ(loading, Target = list(norm))
      cong.corr <- psych::factor.congruence(rotated$loadings, norm)
      congruence <- cong.corr[col(cong.corr) == row(cong.corr)]
      lin.corr <- cor(rotated$loadings, norm)
      lin <- round(lin.corr[col(lin.corr) == row(lin.corr)],2)
      congruence.list <- list(tuckers.phi = congruence,
                              correlation = lin)
    }
  }
  return (congruence.list)
}
