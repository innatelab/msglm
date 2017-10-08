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
      if (nrow(.) <= 1) {
        # return stub to keep msrun
        warning(.$condition[[1]], ": single replicate (", .$msrun[[1]], "), no replicate effects")
        data.frame(msrun = as.character(.$msrun[[1]]),
                   repl_effect = paste0("__ignore_replicate__", .$msrun[[1]]),
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
  frame2matrix(dplyr::filter(msrunXreplEffect_pre.df, !str_detect(repl_effect, "__ignore_replicate__")) %>%
                             dplyr::rename_at("repl_effect", funs(function(x) repl_col)) %>%
                             dplyr::rename_at("msrun", funs(function(x) msrun_col)),
                             msrun_col, repl_col, "Freq")
}
