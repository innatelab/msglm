# convert conditionXeffect matrix to a frame and calculate
# condition2effect weights
#' @export
conditionXeffect_frame <- function(conditionXeffect_mtx, effects_df) {
  effects_df <- maybe_rename(effects_df, c("prior_tau" = "tau"))
  tibble::as_tibble(as.table(conditionXeffect_mtx)) %>%
    dplyr::filter(n != 0) %>% dplyr::rename(mult = n) %>%
    # estimate the weight of a condition in an effect
    dplyr::inner_join(dplyr::select(effects_df, effect, prior_tau)) %>%
    dplyr::mutate(w = prior_tau * abs(mult)) %>%
    dplyr::group_by_at(names(dimnames(conditionXeffect_mtx))[[1]]) %>%
    dplyr::mutate(eff_w = w/sum(w)) %>%
    dplyr::group_by_at(names(dimnames(conditionXeffect_mtx))[[2]]) %>%
    dplyr::mutate(cond_w = eff_w/sum(eff_w)) %>%
    dplyr::ungroup()
}

# extracts the value of given factor from the effect label
#' @export
effect_factor <- function(effects, factor_name, factor_levels, default = factor_levels[1]) {
    # FIXME factor levels are not allowed to start with '_' (to workaround factor with names being prefixes of the other factors)
    factor_vals <- str_match(effects, paste0("(?:^|:)", factor_name, "([^:_][^:]*)(?:$|:)"))[,2]
    factor_vals[is.na(factor_vals)] <- default
    factor(factor_vals, levels=factor_levels)
}

