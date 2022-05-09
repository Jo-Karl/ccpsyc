
#' Pairwise Effect sizes of similarities and difference in the psychometric structure between multiple groups
#'
#' @param df Multi-Group data frame
#' @param group String variable indicating the grouping variable
#' @param items Vector of strings indicating items for the uni-factorial construct
#' @param eff_sizes Effect sizes to be returned
#'
#' @return The function returns a list of dataframes with the first reporting the averaged results per item
#' and the second reporting the pairwise comparisons.
#' @export multi_group_eff
#'
#' @examples
#' \donttest{
#' example_s <- dplyr::filter(example, country %in% c("NZ", "BRA"))
#' multi_group_eff(df = example, group = "country", items = paste0("voice",1:3, "m"))
#' }
multi_group_eff <-
  function(df, group, items, eff_sizes = c("SDI2", "UDI2", "WSDI", "WUDI", "dmacs")){
    comb_m <- RcppAlgos::comboGeneral(unique(df[[group]]),2, repetition = F)
    list_names <- paste0(comb_m[,1],"-", comb_m[,2])

    df_out <- lapply(lapply(seq_len(nrow(comb_m)), function(x) comb_m[x,]), function(x){
      df[df[[group]] %in% x,]
    })

    eff_out_l <- lapply(df_out, function(x){
      invariance_eff(x, items, group)
    })
    names(eff_out_l) <- list_names

    selected_effs <-
      lapply(eff_sizes, function(eff_size){
        eff_m <-
          do.call(rbind,lapply(eff_out_l, function(x) x[[eff_size]]))

        eff_out <-
          data.frame(
            item = items,
            av = colMeans(eff_m),
            l_ci = t(apply(eff_m, 2, stats::quantile, probs = c(0.025, 0.975)))[,1],
            h_ci = t(apply(eff_m, 2, stats::quantile, probs = c(0.025, 0.975)))[,2],
            sd = apply(eff_m, 2, stats::sd),
            eff = eff_size
          )
      })

    total_df <- do.call(rbind, selected_effs)

    return(list(average = total_df, individual = eff_out_l))
  }
