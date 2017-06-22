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

normalize_experiments <- function(norm_model, def_norm_data, df,
                                  quant_col = "intensity", obj_col = "protgroup_id",
                                  exp_col = "mschannel_ix", cond_col = "bait_orf",
                                  cond_group_col = NULL, exp_shifts = NULL, shift_col="shift",
                                  missing_exp.ratio=0.1, missing_cond.ratio=if (exp_col==cond_col) missing_exp.ratio else 0.0,
                                  method = c("optimizing", "mcmc", "vb"), center_shifts=TRUE,
                                  max_quant.ratio=NA, exp_shift.min = -0.5, Rhat_max = 1.1,
                                  max_objs=1000L,
                                  mcmc.iter=2000L, mcmc.chains=4, mcmc.thin=4, mcmc.adapt_delta=0.9,
                                  vb.iter=100000L,
                                  verbose=FALSE)
{
    method <- match.arg(method)
    if (is.null(cond_group_col)) {
      df$`__all__` <- TRUE
      cond_group_col <- "__all__"
    }
    # prepare "standartized" data frame
    df_std <- df[,c(obj_col, exp_col, cond_col, cond_group_col, quant_col)]
    names(df_std) <-c("obj", "experiment", "condition", "cond_group", "quant")
    if (is.null(exp_shifts)) {
      # no shifts by default
      exp_shifts <- data.frame(experiment = sort(unique(df_std$experiment)),
                               shift = rep.int(0.0, n_distinct(df_std$experiment)),
                               stringsAsFactors = TRUE)
    } else {
      exp_shifts <- exp_shifts[,c(exp_col, shift_col)]
      colnames(exp_shifts) <- c("experiment", "shift")
    }
    if (!is.na(exp_shift.min) && is.finite(exp_shift.min)) {
      exp_shifts <- dplyr::mutate(exp_shifts, is_used = shift >= exp_shift.min)
      if (!all(exp_shifts$is_used)) {
        warning(sum(!exp_shifts$is_used), " experiment(s) with low shifts excluded from ", cond_col, " normalization: ",
                paste0(exp_shifts$experiment[!exp_shifts$is_used], collapse=" "))
        exp_shifts <- dplyr::filter(exp_shifts, is_used)
      }
      df_std <- dplyr::inner_join(df_std, exp_shifts %>% dplyr::select(experiment, shift)) %>%
        dplyr::mutate(norm_quant = quant*exp(-shift)) %>% dplyr::select(-shift)
    } else {
      exp_shifts <- dplyr::mutate(exp_shifts, is_used = TRUE)
      df_std$norm_quant <- df_std$quant
    }
    exp_shifts$experiment <- as.character(exp_shifts$experiment)
    obj_group_stats <- df_std %>% dplyr::filter(!is.na(quant)) %>%
      dplyr::group_by(obj, cond_group) %>%
      dplyr::summarise(n_experiments = n_distinct(experiment),
                       n_conditions = n_distinct(condition)) %>%
      dplyr::ungroup()
    group_stats <- dplyr::group_by(obj_group_stats, cond_group) %>%
      dplyr::summarise(n_max_experiments = max(n_experiments),
                       n_max_conditions = max(n_conditions))
    valid_objs <- dplyr::inner_join(obj_group_stats, group_stats) %>%
      dplyr::filter(n_experiments > 1L & n_experiments >= (1-missing_exp.ratio)*n_max_experiments &
                    n_conditions >= (1-missing_cond.ratio)*n_max_conditions)
    message(n_distinct(group_stats$cond_group), " condition group(s)")
    res <- dplyr::group_by(valid_objs, cond_group) %>% do({
        df_group <- dplyr::inner_join(df_std, .)
        group_objs <- df_group %>%
          dplyr::group_by(obj, n_max_experiments, n_max_conditions, cond_group) %>%
          dplyr::summarise(quant_med = median(norm_quant, na.rm=TRUE),
                           quant_min = min(norm_quant, na.rm=TRUE),
                           quant_max = max(norm_quant, na.rm=TRUE),
                           sd = sd(norm_quant, na.rm=TRUE)) %>% dplyr::ungroup() %>%
          dplyr::mutate(sd_rel = sd/quant_med,
                        sd_quantile = percent_rank(if_else(!is.na(sd_rel), sd_rel, Inf)),
                        is_valid = is.na(max_quant.ratio) |
                                   (quant_min*max_quant.ratio > quant_med &
                                    quant_max < max_quant.ratio*quant_med))
        # max sure at least 25% of group objects or max_objs (if specified) are valid for normalization
        if (sum(group_objs$is_valid) < pmin(0.25*nrow(group_objs),
                                            ifelse(max_objs > 0L, max_objs, Inf))) {
            group_objs$is_valid <- group_objs$sd_quantile <= 0.25
        }
        valid_group_objs <- dplyr::filter(group_objs, is_valid)
        sel_group_objs <- if (max_objs > 0L && nrow(valid_group_objs) > max_objs) {
          sample_n(valid_group_objs, max_objs)
        } else {
          valid_group_objs
        }
        df_group <- dplyr::semi_join(df_group, sel_group_objs)
        norm_data <- def_norm_data
        norm_data$Nexperiments <- sel_group_objs$n_max_experiments[1]
        norm_data$Nconditions <- sel_group_objs$n_max_conditions[1]
        norm_data$Nobjects <- nrow(sel_group_objs)
        message("Normalizing group '", sel_group_objs$cond_group[1], "', ",
                norm_data$Nobjects, " object(s), ",
                norm_data$Nconditions, " condition(s), ",
                norm_data$Nexperiments, " experiment(s)")
        if (nrow(df_group) == 0L) {
            # degenerated case
            warning("No valid observations in group ", .$cond_group[1])
            res <- data.frame(condition = as.character(df_group$condition[1]),
                              shift = 0.0,
                              stringsAsFactors = FALSE)
            colnames(res) <- c(cond_col, "shift")
            return (res)
        } else if (norm_data$Nconditions == 1L) {
            # another degenerated case
            res <- data.frame(condition = as.character(df_group$condition[1]),
                              shift = 0.0,
                              stringsAsFactors = FALSE)
            colnames(res) <- c(cond_col, "shift")
            return (res)
        }
        #if (any(is.na(df_group$quant))) {
        #  stop("Non-quanted observations detected")
        #}
        exp_df <- dplyr::select(df_group, experiment, condition) %>%
          dplyr::distinct() %>% dplyr::arrange(condition, experiment) %>%
          dplyr::mutate(experiment = as.character(experiment),
                        condition= as.character(condition)) %>%
          dplyr::left_join(exp_shifts) %>%
          dplyr::mutate(experiment = factor(experiment, levels=experiment),
                        condition = factor(condition, levels=unique(condition))) %>%
          dplyr::arrange(as.integer(experiment))
        if (any(is.na(exp_df$shift))) {
          warning("No shifts for experiments: ",
                  paste0(exp_df$experiment[is.na(exp_df$shift)], collapse=" "))
        }
        norm_data$experiment2condition <- as.array(as.integer(exp_df$condition))
        norm_data$experiment_shift <- as.array(exp_df$shift)
        df_group <- dplyr::mutate(df_group,
                                  experiment = factor(as.character(experiment),
                                                      levels=levels(exp_df$experiment)))
        df_group <- dplyr::left_join(tidyr::expand(df_group, obj, experiment), df_group) %>%
            dplyr::arrange(as.integer(obj), as.integer(experiment)) %>%
            dplyr::mutate(safe_quant = if_else(is.na(quant), 0.0, quant))

        norm_data$qData <- matrix(df_group$safe_quant, ncol = nrow(exp_df), byrow=TRUE) #nrow=nrow(valid_group_objs)
        message("Running Stan optimization...")
        out_params <- c("data_sigma", "condition_sigma", "condition_shift")
        if (method == 'optimizing') {
          norm_fit <- optimizing(norm_model, norm_data, algorithm="LBFGS",
                                 init=list(condition_sigma=1.0, condition_shift0=as.array(rep.int(0.0, norm_data$Nconditions-1L))),
                                 history_size=10L)
          cond_shift_pars <- norm_fit$par[str_detect(names(norm_fit$par), "^condition_shift\\[\\d+\\]$")]
          cond_shift_ixs <- as.integer(str_match(names(cond_shift_pars), "\\[(\\d+)\\]$")[,2])
          res <- data.frame(condition = levels(exp_df$condition)[cond_shift_ixs],
                            shift = as.numeric(cond_shift_pars),
                            stringsAsFactors=FALSE)
        } else if (method == 'mcmc') {
          norm_fit <- sampling(norm_model, norm_data, chains=mcmc.chains, iter=mcmc.iter, thin=mcmc.thin,
                               pars=out_params, include=TRUE,
                               init=function() list(condition_sigma=1.0, condition_shift0=as.array(rep.int(0.0, norm_data$Nconditions-1L))),
                               control = list(adapt_delta=mcmc.adapt_delta))
          norm_fit_stat <- monitor(norm_fit)
          cond_shift_mask <- str_detect(rownames(norm_fit_stat), "^condition_shift\\[\\d+\\]$")
          nonconv_mask <- norm_fit_stat[cond_shift_mask, 'Rhat'] > Rhat_max
          if (any(nonconv_mask)) {
            warning("Rhat>", Rhat_max, " for ", sum(nonconv_mask), " shift(s)")
          }
          n_eff_min <- 0.1*mcmc.iter
          nonconv_mask <- norm_fit_stat[cond_shift_mask, 'n_eff'] < n_eff_min
          if (any(nonconv_mask)) {
            warning("n_eff<", n_eff_min, " for ", sum(nonconv_mask), " shift(s)")
          }
          cond_shift_pars <- norm_fit_stat[cond_shift_mask, 'mean']
          cond_shift_ixs <- as.integer(str_match(rownames(norm_fit_stat)[cond_shift_mask], "\\[(\\d+)\\]$")[,2])
          res <- data.frame(condition = levels(exp_df$condition)[cond_shift_ixs],
                            shift = as.numeric(cond_shift_pars),
                            Rhat = norm_fit_stat[cond_shift_mask, 'Rhat'],
                            n_eff = norm_fit_stat[cond_shift_mask, 'n_eff'],
                            stringsAsFactors=FALSE)
        } else if (method == 'vb') {
          norm_fit <- vb(norm_model, norm_data, iter=vb.iter,
                         pars=out_params, include=TRUE,
                         init=function() list(condition_sigma=1.0, condition_shift0=as.array(rep.int(0.0, norm_data$Nconditions-1L))) )
          norm_fit_stat <- monitor(norm_fit)
          cond_shift_mask <- str_detect(rownames(norm_fit_stat), "^condition_shift\\[\\d+\\]$")
          cond_shift_pars <- norm_fit_stat[cond_shift_mask, 'mean']
          cond_shift_ixs <- as.integer(str_match(rownames(norm_fit_stat)[cond_shift_mask], "\\[(\\d+)\\]$")[,2])
          res <- data.frame(condition = levels(exp_df$condition)[cond_shift_ixs],
                            shift = as.numeric(cond_shift_pars),
                            stringsAsFactors=FALSE)
        } else {
          stop('Unknown method ', method)
        }
        if (center_shifts) {
          res$shift <- res$shift - median(res$shift)
        }
        col_renames <- "condition"
        names(col_renames) <- cond_col
        res <- dplyr::rename_(res, .dots=col_renames) %>%
          dplyr::mutate(
            method = method,
            n_objects = norm_data$Nobjects,
            n_experiments = norm_data$Nexperiments)
    }) %>% dplyr::ungroup()
    if (cond_group_col == "__all__") {
      res$cond_group <- NULL
    } else {
      new_colnames <- colnames(res)
      new_colnames[new_colnames == "cond_group"] <- cond_group_col
      colnames(res) <- new_colnames
    }
    res
}

calc_contrasts <- function(vars_results, vars_info, dims_info, contrastXmetacondition, conditionXmetacondition.df, contrastXcondition.df,
                           mschannel_col = "mschannel_ix",
                           condition_agg_col = "condition", var_names = c('iaction_labu', 'obs_labu'),
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

process_model_fit_chunk <- function(file_ix, strip_samples=FALSE, strip_stats=FALSE) {
  fit_model.file <- fit_model.files.df[file_ix, 'filename']
  message('Loading ', fit_model.file, '...')
  tmp.env <- new.env(parent = baseenv())
  load(file.path(fit_model_path, fit_model.file), envir = tmp.env)
  if (strip_samples || strip_stats) {
    tmp.env$vars_results <- lapply( tmp.env$vars_results, function( var_results ) {
      if (strip_samples) { var_results$samples <- NULL }
      if (strip_stats) { var_results$stats <- NULL }
      return ( var_results )
    } )
  }
  gc()
  res <- mget(ls(envir=tmp.env), envir=tmp.env)
  return ( res )
}

process_msglm_chunk <- function(file_ix, sel_contrasts = NULL, condition_shifts = NULL,
                                strip_samples=FALSE, strip_stats=FALSE,
                                val_trans = NULL, condition.reported = "lhs",
                                condition.quantiles_lhs = c(0, 1), condition.quantiles_rhs = c(0, 1),
                                condition_agg_col = "condition") {
  fit_model.file <- fit_model.files.df[ file_ix, 'filename' ]
  message( 'Loading ', fit_model.file, '...' )
  tmp.env <- new.env(parent = baseenv())
  load(file.path(fit_model_path, fit_model.file ), envir = tmp.env)
  if (!is.null(sel_contrasts)) {
    # FIXME do once per assembly
    sel_contrastXmetacondition <- contrastXmetacondition[sel_contrasts, , drop=FALSE]
    sel_contrastXmetacondition <- sel_contrastXmetacondition[, colSums(abs(sel_contrastXmetacondition)) != 0, drop=FALSE]
    sel_conditionXmetacondition.df <- dplyr::filter(conditionXmetacondition.df,
                                                    metacondition %in% colnames(sel_contrastXmetacondition))
    if (!is.null(condition_shifts)) {
      sel_conditionXmetacondition.df <- dplyr::inner_join(sel_conditionXmetacondition.df, condition_shifts)
    }
    sel_contrastXcondition.df <- dplyr::semi_join(contrastXcondition.df, sel_conditionXmetacondition.df) %>%
        dplyr::filter(contrast %in% sel_contrasts)
    tmp.env$vars_results <- calc_contrasts(tmp.env$vars_results, tmp.env$vars_info, tmp.env$dims_info,
                                           sel_contrastXmetacondition, sel_conditionXmetacondition.df,
                                           sel_contrastXcondition.df, val_trans = val_trans,
                                           condition_agg_col = condition_agg_col,
                                           condition.reported = condition.reported,
                                           condition.quantiles_lhs = condition.quantiles_lhs,
                                           condition.quantiles_rhs = condition.quantiles_rhs)
  }
  if (strip_samples || strip_stats) {
    tmp.env$vars_results <- lapply( tmp.env$vars_results, function( var_results ) {
      if (strip_samples) { var_results$samples <- NULL }
      if (strip_stats) { var_results$stats <- NULL }
      return ( var_results )
    } )
  }
  gc()
  res <- mget(ls(envir=tmp.env), envir=tmp.env)
  return ( res )
}

join_msglm_reports <- function(section, reports, type) {
  message( 'Assembling joint ', type, ' report for ', section, '...' )
  res <- join_report_frames(reports, frame_extractor = function(report) {
    res <- report$vars_results[[section]][[type]]
    if (!is.null(res) && !('majority_protein_acs' %in% colnames(res))
        && "majority_protein_acs" %in% colnames(report$model_data$protgroups)) {
      res$majority_protein_acs <- rep_len(report$model_data$protgroups$majority_protein_acs[1], nrow(res))
      res$protgroup_id <- rep_len(report$model_data$protgroups$protgroup_id[1], nrow(res))
    }
    if (!is.null(res) && !('pepmod_id' %in% colnames(res))
        && "pepmod_id" %in% colnames(report$model_data$pepmods)) {
      res$pepmod_id <- rep_len(report$model_data$pepmods$pepmod_id[1], nrow(res))
    }
    res
  },
  NULL, global_vars = c())
  return(res)
}
