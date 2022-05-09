#' Improving boot effectsize output
#'
#' @param x The output of a bootstrapped invariance effect call
#'
#' @return A formatted tibble with all effect sizes reported by boot_inv_eff from this package and significant determined by 95% CIs either crossing 0 or .30
#' @export format_boot_inv_eff
#'
#' @examples
#' \donttest{
#' two_country <- dplyr::filter(example, country %in% c("NZ" , "BRA"))
#' boot_inv_eff_result <- boot_inv_eff(n = 10,
#'                         n_sample = 200, df = two_country, group = "country",
#'                         items = paste0("voice",1:3, "m"))
#' format_boot_inv_eff(boot_inv_eff_result)
#'               }
format_boot_inv_eff <-
  function(x){
    eff <- h <- l <- av <- item <- sig <- res <- NULL
    dplyr::mutate(x, sig = dplyr::if_else(eff %in% c("SDI2", "WSDI"),
                                          dplyr::if_else(h > 0 & l < 0, "", "*"),
                                          dplyr::if_else(h > .30 & l < .30 | av < .30, "", "*") )) %>%
      dplyr::transmute(item = item,
                eff = eff,
                res = paste0(gsub(x = format(round(av, 3), nsmall = 3), pattern = "0\\.", replacement = "\\."),
                             "[",
                             gsub(x = format(round(l, 3), nsmall = 3), pattern = "0\\.", replacement = "\\."),
                             ", ",
                             gsub(x = format(round(h, 3), nsmall = 3), pattern = "0\\.", replacement = "\\."),
                             "]",
                             sig)) %>%
      tidyr::pivot_wider(id_cols = item, names_from = eff, values_from = res)
  }
