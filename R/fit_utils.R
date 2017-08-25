require(rstan)
require(stringr)

stan_models_path <- file.path(base_scripts_path, "msglm/inst/stan_models") # FIXME
msglm.stan_model <- stan_model(file.path(stan_models_path, "msglm.stan"), "msglm", save_dso = TRUE, auto_write = TRUE)
msglm_local.stan_model <- stan_model(file.path(stan_models_path, "msglm_local.stan"), "msglm_local", save_dso = TRUE, auto_write = TRUE)
msglm_normalize.stan_model <- stan_model(file.path(stan_models_path , "msglm_normalize.stan"), "msglm_normalize", save_dso = TRUE, auto_write = TRUE)

prepare_observations <- function(data.mtx)
{
  data.vec <- as.vector(data.mtx)
  data.vec[ data.vec == 0 ] <- NA
  data.obs.mask <- !is.na( data.vec )
  data.idx <- ifelse( data.obs.mask, cumsum( data.obs.mask ), -cumsum( !data.obs.mask ) )
  data.idx.arr <- matrix( data.idx, ncol = ncol( data.mtx ) )
  list( Nobserved = sum(data.obs.mask),
        oData = array( data.vec[ data.obs.mask ] ),
        dataIndex = as.array(data.idx.arr)
  )
}

obj2effect <- function(expXeff, obj_class, objs, data2obj, data2exp) {
  effXexp <- t(expXeff)
  expXeff_mask <- apply(effXexp, c(2, 1), function(x) x != 0.0)
  objXexp_mask <- matrix(FALSE, nrow = length(objs), ncol = ncol(effXexp))
  if (length(data2obj)>0) { for (i in 1:length(data2obj)) {
    objXexp_mask[data2obj[[i]], data2exp[[i]]] <- TRUE
  } }
  objXeff_mask <- objXexp_mask %*% expXeff_mask
  dnams <- list(objs, dimnames(effXexp)[[1]])
  names(dnams) <- c(obj_class, names(dimnames(effXexp))[[1]])
  dimnames(objXeff_mask) <- dnams
  if (length(objXeff_mask) > 0L) {
    res <- as.data.frame(as.table(t(objXeff_mask))) %>% dplyr::filter(Freq != 0) %>% dplyr::select(-Freq)
  } else {
    res <- data.frame() %>% dplyr::mutate(x = character(), y = character())
    names(res) <- c(obj_class, names(dnams)[[2]])
  }
  res[,paste(colnames(res), collapse='_')] <- seq_len(nrow(res))
  res
}

effects_cumsum <- function(df, group_col) {
  n_effects <- df %>% dplyr::arrange_(group_col) %>% dplyr::group_by_(group_col) %>%
    dplyr::summarize(n_effects = n()) %>% .$n_effects
  if (length(n_effects)>0) {
    res <- cumsum(as.integer(n_effects))
    c(1L, res+1L)
  } else {
    return(0L)
  }
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
  msrunXreplEffect_wide.df <- reshape(msrunXreplEffect_pre.df, direction="wide",
                                      idvar="msrun", timevar="repl_effect", v.names="Freq") %>%
    dplyr::select(-condition, -starts_with("Freq.__ignore_replicate__")) %>%
    dplyr::mutate_at(vars(starts_with("Freq")), funs(if_else(is.na(.), 0.0, .)))
  res <- dplyr::select(msrunXreplEffect_wide.df, -msrun) %>% as.matrix()
  mXre_dims <- list(msrun = msrunXreplEffect_wide.df$msrun,
                    repl_effect = str_replace(colnames(res), "^Freq\\.", ""))
  names(mXre_dims) <- c(msrun_col, paste0(replicate_col, "_effect"))
  dimnames(res) <- mXre_dims
  res
}

calc_contrasts <- function(vars_results, vars_info, dims_info, contrastXmetacondition, conditionXmetacondition.df, contrastXcondition.df,
                           mschannel_col = "mschannel_ix",
                           condition_agg_col = "condition", var_names = c('iaction_labu', 'iaction_labu_replCI', 'obs_labu'),
                           obj_dims = 'glm_object_ix', var_info_cols = 'glm_object_ix', val_trans = NULL,
                           condition.quantiles_lhs = c(0, 1), condition.quantiles_rhs = c(0, 1)) {
  for (vars_category in names(vars_results)) {
    vars_cat_subset_info <- vars_info[[vars_category]]
    vars_cat_subset_info$names <- intersect(vars_cat_subset_info$names, var_names)
    if ( length(vars_cat_subset_info$names) > 0 ) {
      message('Filtering ', vars_category, ' to use for contrasts...')
      cond_stats.df <- dplyr::inner_join(vars_results[[vars_category]]$stats,
                                         dplyr::mutate(contrastXcondition.df, is_lhs = weight > 0)) %>%
        dplyr::group_by_(.dots=c(obj_dims, "contrast", "metacondition", condition_agg_col, "is_lhs")) %>%
        dplyr::group_by_(.dots = condition_agg_col, add=TRUE) %>%
        dplyr::summarize(cond_mean = mean(mean)) %>% dplyr::ungroup()

      cond_agg_stats.df <- dplyr::group_by_(cond_stats.df, .dots=c(obj_dims, "contrast", "metacondition", "is_lhs")) %>%
        dplyr::group_by_(.dots = condition_agg_col, add=TRUE) %>%
        dplyr::summarize(cond_max_mean = max(cond_mean),
                         cond_min_mean = min(cond_mean)) %>%
        dplyr::group_by_(.dots=c(obj_dims, "contrast", "metacondition", "is_lhs")) %>%
        dplyr::mutate(cond_max_qtile = cume_dist(cond_max_mean) - 1/n(),
                      cond_min_qtile = cume_dist(cond_min_mean) - 1/n()) %>% dplyr::ungroup() %>%
        dplyr::select_(.dots = c(obj_dims, "contrast", "metacondition", condition_agg_col,
                                 "is_lhs", "cond_min_qtile", "cond_max_qtile", "cond_max_mean", "cond_min_mean"))

      cond_stats.df <- dplyr::inner_join(cond_stats.df, cond_agg_stats.df) %>%
        dplyr::mutate(is_accepted = (is_lhs & cond_min_qtile >= condition.quantiles_lhs[[1]] & cond_max_qtile <= condition.quantiles_lhs[[2]]) |
                                    (!is_lhs & cond_min_qtile >= condition.quantiles_rhs[[1]] & cond_max_qtile <= condition.quantiles_rhs[[2]]))
      cond_stats.df <- dplyr::filter(cond_stats.df, is_accepted)
      message( 'Calculating contrasts for ', vars_category, ' variables...' )
      samples.df <- dplyr::semi_join(vars_results[[vars_category]]$samples, cond_stats.df) # exclude unused conditions
      experiment_col <- if (mschannel_col %in% colnames(samples.df)) mschannel_col
                        else if ('condition' %in% colnames(samples.df)) 'condition'
                        else 'condition'
      metacondition2experiments.df <- dplyr::inner_join(dplyr::semi_join(conditionXmetacondition.df, cond_stats.df),
                                                        dplyr::distinct(dplyr::select_(samples.df,
                                                                                       .dots=unique(c("condition", experiment_col)))))
      if (nrow(samples.df) == 0) {
        # no samples
        vars_results[[vars_category]]$contrast_stats <- NULL
      } else {
        # adjust samples w.r.t. condition shift
        if ('condition_shift' %in% colnames(samples.df)) {
          for (shift_name in vars_cat_subset_info$names) {
            samples.df[,shift_name] <- samples.df[,shift_name] - samples.df$condition_shift
          }
        }
        contrast_stats.df <- vars_contrast_stats(samples.df, vars_cat_subset_info$names,
                                                 group_cols = obj_dims, condition_col = 'metacondition', experiment_col = experiment_col,
                                                 condition2experiments.df = metacondition2experiments.df,
                                                 contrastXcondition = contrastXmetacondition, val_trans = val_trans)
        var_info.df <- vars_results[[vars_category]]$stats %>%
          dplyr::inner_join(conditionXmetacondition.df) %>%
          dplyr::select_(.dots = setdiff(var_info_cols, c('index_observation', 'iaction_id', 'condition', condition_agg_col, experiment_col))) %>%
          dplyr::distinct()
        # inject contrast statistics into vars_results
        vars_results[[vars_category]]$contrast_stats <- left_join(var_info.df, contrast_stats.df) %>%
          dplyr::mutate(p_value = 2*pmin(prob_nonpos, prob_nonneg)
                        # fake P-values to fit in the plot
                        #prob_nonpos_fake = pmax(if_else(mean < 0, 1E-300, 1E-100), rgamma(n(), shape=1E-2, scale=1E-150), prob_nonpos),
                        #prob_nonneg_fake = pmax(if_else(mean > 0, 1E-300, 1E-100), rgamma(n(), shape=1E-2, scale=1E-150), prob_nonneg),
                        #p_value_fake = 2*pmin(prob_nonpos_fake, prob_nonneg_fake)
          ) %>%
          #dplyr::select(-index_observation, -msrun, -msrun_ix) %>% dplyr::distinct() %>%
          left_join(contrastXmetacondition.df %>% dplyr::filter(contrast %in% rownames(contrastXmetacondition)) %>% dplyr::rename(contrast_weight = weight) %>%
                      dplyr::mutate(metacondition_reported = case_when(.$contrast_type == "filtering" ~ "lhs",
                                                                       .$contrast_type == "comparison" ~ "enriched",
                                                                       TRUE ~ NA_character_))) %>%
          dplyr::filter((metacondition_reported == "lhs" & contrast_weight > 0) |
                        (metacondition_reported == "enriched" & contrast_weight * `50%` > 0)) %>%
          dplyr::select(-metacondition_reported, -contrast_type) %>%
          dplyr::inner_join(dims_info$object)
      }
    }
  }
  return(vars_results)
}

