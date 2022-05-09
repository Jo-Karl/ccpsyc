#' Bootstrapped pairwise differences in psychometric function of groups.
#'
#' @param n Number of bootstraps
#' @param n_sample Number of participants to sample
#' @param df Data to resample
#' @param items Items to resample for the model as vector of strings
#' @param group String variable indicating grouping variable
#' @param eff_sizes Effect sizes to be returned
#' @param seed Seed for replicability
#'
#' @return Returns a dataframe with the bootstrapped effect sizes based on the invariance_eff function in this package for two country comparisons.
#' @export boot_inv_eff
#'
#' @examples
#' \donttest{
#' two_country <- dplyr::filter(example, country %in% c("NZ" , "BRA"))
#' boot_inv_eff(n = 10,
#'              n_sample = 200, df = two_country, group = "country",
#'               items = paste0("voice",1:3, "m"))
#'               }
boot_inv_eff <- function(n, n_sample, df, items, group, eff_sizes = c("SDI2", "UDI2", "WSDI", "WUDI", "dmacs"), seed = 2711){
  set.seed(seed)
  t_list <-
    lapply (seq_len(n), function(x){
      bootdat <- df %>% dplyr::group_by(eval(parse(text = group))) %>% dplyr::sample_n(n_sample, replace = T)
      boot_test <- invariance_eff(df = bootdat, items = items, group = group, nodewidth = .01)
      boot_test$boot <-  paste0(x)
      boot_test
    })


  selected_effs <-
    lapply(eff_sizes, function(eff_size){
      eff_m <- do.call(rbind, lapply(t_list, function(x) x[[eff_size]]))

      eff_out <-
        data.frame(
          item = items,
          av = colMeans(eff_m),
          l = t(apply(eff_m, 2, stats::quantile, probs = c(0.025, 0.975)))[,1],
          h = t(apply(eff_m, 2, stats::quantile, probs = c(0.025, 0.975)))[,2],
          eff = eff_size
        )
    })

  do.call(rbind, selected_effs)

}
