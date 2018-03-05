require(stringr)

# convert conditionXeffect matrix to a frame and calculate
# condition2effect weights
conditionXeffect_frame <- function(conditionXeffect_mtx, effects_df) {
  as.data.frame(as.table(conditionXeffect_mtx), stringsAsFactors = FALSE) %>%
    dplyr::filter(Freq != 0) %>% dplyr::rename(mult = Freq) %>%
    # estimate the weight of a condition in an effect
    dplyr::inner_join(dplyr::select(effects_df, effect, tau)) %>%
    dplyr::mutate(w = tau * mult) %>%
    dplyr::group_by_(.dots=names(dimnames(conditionXeffect_mtx))[[1]]) %>%
    dplyr::mutate(eff_w = w/sum(w)) %>%
    dplyr::group_by_(.dots=names(dimnames(conditionXeffect_mtx))[[2]]) %>%
    dplyr::mutate(cond_w = eff_w/sum(eff_w)) %>%
    dplyr::ungroup()
}

# compose replicate effects design matrix
replicate_effects_matrix <- function(mschannels_df, msrun_col="msrun", cond_col="condition", replicate_col="msrun") {
  # compose replicate effects
  # ~ (replicate %in% condition) doesn't work because of assymetrical design (diff. N of replicates)
  mschannels_df <- dplyr::select_at(mschannels_df, c("msrun"=msrun_col, "condition"=cond_col, "replicate"=replicate_col))
  msrunXreplEffect_pre.df <- dplyr::group_by(mschannels_df, condition) %>%
    dplyr::do({
      if (n_distinct(.$replicate) <= 1) {
        # return stub to keep msrun
        warning(.$condition[[1]], ": single replicate (", .$msrun[[1]], "), no replicate effects")
        data.frame(msrun = as.character(.$msrun),
                   repl_effect = paste0("__ignore_replicate__", .$condition[[1]]),
                   Freq = 0.0,
                   stringsAsFactors = FALSE)
      } else {
        cond_msruns.df <- dplyr::mutate(., replicate=factor(replicate))
        mtx <- model.matrix(~ replicate,
                          cond_msruns.df,
                          contrasts.arg = list("replicate" = "contr.sum"))
        mtx <- mtx[, colnames(mtx) != "(Intercept)"]
        dimnames(mtx) <- list(msrun = cond_msruns.df$msrun,
                              repl_effect = paste0(cond_col, cond_msruns.df$condition[1], ':', colnames(mtx)))
        as.data.frame(as.table(mtx)) %>% dplyr::filter(Freq != 0.0) %>%
          dplyr::mutate(msrun = as.character(msrun),
                        repl_effect = as.character(repl_effect))
      }
    }) %>% dplyr::ungroup()
  repl_col <- paste0(replicate_col, "_effect")
  res <- frame2matrix(msrunXreplEffect_pre.df %>%
                      dplyr::rename_at("repl_effect", funs(function(x) repl_col)) %>%
                      dplyr::rename_at("msrun", funs(function(x) msrun_col)),
                      msrun_col, repl_col, "Freq")
  res[, !str_detect(colnames(res), "__ignore_replicate__"), drop=FALSE]
}

# extracts the value of given factor from the effect label
effect_factor <- function(effects, factor_name, factor_levels, default = factor_levels[1]) {
    # FIXME factor levels are not allowed to start with '_' (to workaround factor with names being prefixes of the other factors)
    factor_vals <- str_match(effects, paste0("(?:^|:)", factor_name, "([^:_][^:]*)(?:$|:)"))[,2]
    factor_vals[is.na(factor_vals)] <- default
    factor(factor_vals, levels=factor_levels)
}

