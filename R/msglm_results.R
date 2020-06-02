#' @export
msglm.prepare_dims_info <- function(model_data, object_cols = NULL)
{
  is_glmm <- "mixeffects" %in% names(model_data)
  xaction_ix_col <- if (is_glmm) "glm_supaction_ix" else "glm_iaction_ix"
  xdition_ix_col <- if (is_glmm) "supcondition_ix" else "condition_ix"
  xdition_col <-  if (is_glmm) "supcondition" else "condition"

  objs_df <- model_data$objects
  if (!is.null(object_cols)) {
    objs_df <- dplyr::select(objs_df, !!!unique(c("glm_object_ix", object_cols)))
  }
  res <- list(iteration = NULL,
    msrun = dplyr::select(model_data$mschannels, msrun_ix, msrun, any_of(c("mschannel", "mschannel_ix", "mstag", "condition", "supcondition"))),
    observation = dplyr::select(model_data$msdata, glm_observation_ix, !!xaction_ix_col, glm_object_ix,
                                !!xdition_ix_col, !!xdition_col, msrun, msrun_ix, any_of(c("mschannel", "mschannel_ix", "mstag"))) %>%
        dplyr::distinct() %>%
        dplyr::inner_join(objs_df),
    object = model_data$objects, # use full object information
    object_effect = model_data$object_effects %>%
        dplyr::mutate(glm_object_ix = as.integer(glm_object_ix)) %>%
        dplyr::inner_join(objs_df) %>%
        dplyr::inner_join(model_data$effects) %>%
        maybe_rename(c("prior_mean" = "mean", "prior_tau" = "tau")),
    object_batch_effect = model_data$object_batch_effects %>%
        dplyr::mutate(glm_object_ix = as.integer(glm_object_ix)) %>%
        dplyr::inner_join(objs_df)
  )
  if (!rlang::has_name(res$object_effect, "prior_mean")) { # set the default min to 0
    message("object_effects$prior_mean missing, setting to 0")
    res$object_effect$prior_mean <- 0.0
  }
  res$object_effect <- mutate(res$object_effect,
                              prior_mean_log2=prior_mean/log(2))
  if ("subobjects" %in% names(model_data)) {
    res$subobject <- dplyr::select(model_data$subobjects, glm_object_ix, glm_subobject_ix,
                                   any_of(c("protregroup_id", "protgroup_id", "pepmod_id", "pepmodstate_id", "charge")))
    if ("suo_subbatch_effects" %in% names(model_data)) {
      res$subobject_subbatch_effect <- model_data$suo_subbatch_effects %>%
        dplyr::mutate(glm_subobject_ix = as.integer(glm_subobject_ix)) %>%
        dplyr::left_join(res$subobject)
    }
  }
  if ("msproto_ix" %in% colnames(model_data$mschannels)) {
    res$msprotocol <- dplyr::select(model_data$mschannels, msproto_ix,
                                    any_of("instrument")) %>%
      dplyr::distinct()
  }
  res$iaction <- dplyr::select(model_data$interactions, glm_iaction_ix,
                               glm_object_ix, iaction_id, condition_ix, condition, is_virtual) %>%
    dplyr::inner_join(objs_df)
  if (is_glmm) {
    res$object_mixeffect <- dplyr::mutate(model_data$mixeffects, tmp="a") %>%
      dplyr::left_join(dplyr::mutate(objs_df, tmp="a")) %>%
                       dplyr::select(-tmp)
    res$object_mixeffect <- maybe_rename(res$object_mixeffect, c("prior_mean" = "mean", "prior_tau" = "tau"))
    if (!rlang::has_name(res$object_mixeffect, "prior_mean")) { # set the default min to 0
      message("mixeffects$prior_mean missing, setting to 0")
      res$object_mixeffect$prior_mean <- 0.0
    }
    res$object_mixeffect <- mutate(res$object_mixeffect,
                                   prior_mean_log2=prior_mean/log(2))
    res$object_mixcoef <- dplyr::mutate(model_data$mixcoefs, tmp="a") %>%
      dplyr::left_join(dplyr::mutate(objs_df, tmp="a")) %>%
      dplyr::select(-tmp)
    res$supaction <- dplyr::select(model_data$superactions, glm_supaction_ix,
                                   glm_object_ix, supaction_id, supcondition_ix, supcondition, is_virtual) %>%
        dplyr::inner_join(objs_df)
  }
  return(res)
}

#' @export
vars_effect_pvalue <- function(samples.df, vars_cat_info, dim_info,
                               tail = c("both", "negative", "positive"))
{
  tail = match.arg(tail)
  group_cols <- paste0('index_', vars_cat_info$dims)
  #print(str(samples.df))
  #print(str(group_cols))
  p_value_all_samples <- function(samples) {
      missing_samples = setdiff(vars_cat_info$names, names(samples))
      if (length(missing_samples) > 0) {
        warning("Missing samples for vars: ", paste0(missing_samples, collapse=" "))
      }
      dplyr::bind_rows(lapply(vars_cat_info$names[vars_cat_info$names %in% names(samples)], function(col) {
          tibble(var = col,
                 p_value = pvalue_not_zero(rlang::as_double(samples[[col]]) - (if (rlang::has_name(samples, "prior_mean")) samples$prior_mean else 0.0), tail = tail))
          }))
  }
  p_value.df <- samples.df %>% dplyr::group_by_at(group_cols) %>% dplyr::do(p_value_all_samples(.))
  if ('fraction' %in% group_cols) {
    p_value.df$fraction <- as.integer(p_value.df$fraction)
  }
  return (p_value.df)
}

#' @export
vars_contrast_stats <- function(samples.df, var_names, group_cols,
                                condition_col, contrastXcondition,
                                experiment_col = condition_col, condition2experiments.df = NULL,
                                nsteps = 100, maxBandwidth = NA,
                                quant.probs = c(0.025, 0.25, 0.50, 0.75, 0.975),
                                val_trans = NULL)
{
  n_max_samples <- n_distinct(samples.df$unpermuted_ix)
  if (n_max_samples == 0L) {
    warning("No samples for ", paste0(var_names, collapse=" "), " ", condition_col)
    return (tibble())
  }
  if (is.na(condition_col)) {
    # no conditions specified, assume the same condition for all experiments
    # FIXME checking user-defined condition2experiments.df and constrastXcondition
    condition2experiments.df <- dplyr::distinct(dplyr::select(samples.df, !!experiment_col)) %>%
      dplyr::mutate(tmp_condition = '__all__')
    samples.df$tmp_condition <- "__all__"
    condition_col <- "tmp_condition"
    # no contrastXcondition matrix, assume contrast===condition,
    # i.e. combine all the experiments and calculate the summary statistics
    if (is.null(contrastXcondition)) {
      contrastXcondition <- diag(1L)
      dimnames(contrastXcondition) <- list(contrast = c("identity"), condition = "__all__")
    }
    else if (ncol(contrastXcondition) != 1) {
      stop("No condition specified, contrastXcondition matrix should have Nx1 size")
    }
  } else if (experiment_col != condition_col) {
    if (condition_col %in% colnames(samples.df)) {
      if (!is.null(condition2experiments.df)) {
        stop("Ambiguous experiment design specification: condition2experiment.df provided, but ", condition_col, " is also present in the samples.df")
      }
      condition2experiments.df <- dplyr::distinct(dplyr::select(samples.df, !!condition_col, !!experiment_col))
    } else {
      if (is.null(condition2experiments.df)) {
        stop("No experiment design: provide condition2experiment.df or ", condition_col, " in the samples.df")
      }
      condition2experiments.df <- dplyr::select(condition2experiments.df, !!condition_col, !!experiment_col)
    }
  } else { # condition_col == experiment_col
    condition2experiments.df <- dplyr::distinct(dplyr::select(samples.df, !!condition_col))
  }

  contrast_stats_all_samples <- function(samples) {
    #print(str(samples))
    # recode experiment indices to match what is in the samples
    samples_grouped <- dplyr::group_by_at(samples, experiment_col)
    samples_stats <- samples_grouped %>%
      dplyr::summarise(n_samples = n_distinct(unpermuted_ix)) %>% dplyr::ungroup()
    n_min_samples <- min(samples_stats$n_samples)
    if (n_distinct(samples_stats$n_samples) > 1L) {
      warning("Different number of samples per experiment: ",
              paste(unique(samples_stats$n_samples), collapse=" "),
              ", sub-sampling")
      samples <- samples_grouped %>% dplyr::sample_n(size = n_min_samples) %>% dplyr::ungroup()
    }
    cur_experiments = unique(as.character(samples[[experiment_col]]))
    cur_cond2expr.df <- condition2experiments.df[as.character(condition2experiments.df[[experiment_col]]) %in% cur_experiments, , drop=FALSE]
    cur_cond2expr.df[[experiment_col]] <- factor(as.character(cur_cond2expr.df[[experiment_col]]),
                                                 levels=cur_experiments)
    # recode condition indices
    cur_contrastXcondition <- contrastXcondition[,colnames(contrastXcondition) %in% unique(as.character(cur_cond2expr.df[[condition_col]])),drop=FALSE]
    cur_cond2expr.df[[condition_col]] <- factor(as.character(cur_cond2expr.df[[condition_col]]),
                                                levels=colnames(cur_contrastXcondition))

    dplyr::bind_rows(lapply(var_names, function(var_col) {
      sample_vals <- samples[[var_col]]
      if (!is.null(val_trans)) {
        if (val_trans == "exp") {
          # FIXME that doesn't work for assembling per-protein reports
          sample_vals <- exp(sample_vals-median(sample_vals))
        }
        else {
          stop("Unknown val_trans=\"", val_trans, "\". Supported transforms: exp")
        }
      }
      samples.arr <- matrix(sample_vals, nrow=n_min_samples, ncol=nrow(samples_stats))

      if (rownames(contrastXcondition) != "identity" && exists("contrasts.df") && rlang::has_name(contrasts.df, "offset")) {
        contrast_offsets <- set_names(contrasts.df$offset,
                                      contrasts.df$contrast)
        contrast_offsets <- contrast_offsets[rownames(cur_contrastXcondition)]
      } else {
        contrast_offsets <- rep.int(0.0, nrow(cur_contrastXcondition))
      }
      #print(str(samples.arr))
      res <- ContrastStatistics(
        samples.arr, as.integer(cur_cond2expr.df[[experiment_col]]),
        as.integer(cur_cond2expr.df[[condition_col]]),
        cur_contrastXcondition,
        contrast_offsets,
        nsteps = nsteps, maxBandwidth = maxBandwidth,
        quant_probs = quant.probs) %>%
        as_tibble() %>%
        dplyr::mutate(var = var_col)
      # TODO se_mean, n_eff, Rhat
      return (res)
    } ) )
  }
  contrast_stats.df <- dplyr::arrange_at(samples.df,
                                         c(group_cols, experiment_col, "iteration", "chain")) %>%
    dplyr::group_by_at(group_cols) %>%
    dplyr::do(contrast_stats_all_samples(.)) %>%
    dplyr::mutate(mean_log2 = mean/log(2),
                  median_log2 = `50%`/log(2),
                  sd_log2 = sd/log(2))
  if ('fraction' %in% group_cols) {
    contrast_stats.df$fraction <- as.integer(contrast_stats.df$fraction)
  }
  if (!is.null(dimnames(contrastXcondition)[[1]])) {
    contrast_col <- names(dimnames(contrastXcondition))[[1]]
    if (is.null(contrast_col)) {
      contrast_col <- "contrast"
    }
    contrast_stats.df[[contrast_col]] <- dimnames(contrastXcondition)[[1]][contrast_stats.df$index_contrast]
  }
  return (contrast_stats.df)
}

.attach_dim_info.by_array_index <- function(samples.df, stan_samples, dim_info)
{
    indexes.mtx <- if (nrow(samples.df) > 0L) {
      do.call(cbind, lapply(seq_along(dim_info), function(dim_ix) {
        last_col <- colnames(samples.df)[[ncol(samples.df)]]
        slice.index(stan_samples[[last_col]], dim_ix+1) %>% as.vector()
      }))
    } else { # degenerated case
      matrix(integer(0L), ncol=length(dim_info))
    }
    colnames(indexes.mtx) <- paste0('index_', names(dim_info))
    indexes.df <- as_tibble(indexes.mtx)
    samples.df <- dplyr::bind_cols(samples.df, indexes.df)
    for (dim_ix in seq_along(dim_info)) {
      dim_name <- names(dim_info)[[dim_ix]]
      if (!is.null(dim_info[[dim_ix]])) {
        message("Adding ", dim_name, " information")
        samples.df <- dplyr::bind_cols(samples.df, dim_info[[dim_ix]][indexes.df[[dim_ix]], , drop=FALSE])
      } else {
        warning("No ", dim_name, " information provided")
      }
    }
    return (samples.df)
}

.attach_dim_info.by_var_name <- function(data.df, dim_info)
{
    res.df <- data.df
    res_index <- extract_index(res.df$var_name, ndim=length(dim_info))
    dims.df <- tibble(ix = seq_along(dim_info),
                      name = names(dim_info)) %>%
      dplyr::group_by(name) %>%
      dplyr::mutate(local_ix = row_number(),
                    suffix = if (n() > 1) paste0('.', local_ix) else '') %>%
      dplyr::ungroup()
    indexes.df <- tibble(!!!purrr::set_names(lapply(seq_along(dim_info), function(dim_ix) {
      res_index[, dim_ix]
    }), paste0('index_', names(dim_info))))

    for (dim_ix in seq_len(nrow(dims.df))) {
        if (dim_ix > ncol(res_index)) {
            stop('Data contain less dimensions than the dim_info')
        }
        dim_name <- dims.df$name[[dim_ix]]
        if (dims.df$local_ix[[dim_ix]] > 1) {
          warning('Dimension ', dim_name, ' is not unique, using suffixes')
        }
        if (!is.null(dim_info[[dim_ix]])) {
            res.df[[dim_name]] <- NULL
            maxdimix = if (nrow(res_index)>0L) max(res_index[,dim_ix]) else 0L
            if (nrow(dim_info[[dim_ix]]) != maxdimix) {
                stop('Dimension #', dim_ix, '(', dim_name,') info contains ', nrow(dim_info[[dim_ix]]),
                     ' elements, the data contain ', maxdimix)
            }
            res_dim_info.df <- dim_info[[dim_ix]][res_index[,dim_ix], , drop=FALSE]
            if (dims.df$suffix[[dim_ix]] != '') {
              colnames(res_dim_info.df) <- paste0(colnames(res_dim_info.df), dims.df$suffix[dim_ix])
            }
            res.df <- dplyr::bind_cols(res.df, res_dim_info.df)
        }
    }
    return (dplyr::bind_cols(indexes.df, res.df))
}

stan_samples_frame <- function(stan_samples, var_names, var_dims) {
    # extract only samples of selected variables
    # convert arrays of samples of selected variables into data frames,
    # add back single copy of dimensions and real iteration information
    if (any(!(var_names %in% names(stan_samples)))) {
      warning("Variable(s) samples not found: ", paste0(setdiff(var_names, names(stan_samples)), collapse=", "))
    }
    avail_var_names <- intersect(var_names, names(stan_samples))
    samples.df <- tibble(!!!purrr::set_names(lapply(avail_var_names, function(var_name) {
      as.vector(stan_samples[[var_name]])
    }), avail_var_names))
    if ('iter_info' %in% names(attributes(stan_samples))) {
      iter_info.df <- attr(stan_samples, 'iter_info')
      samples.df <- dplyr::bind_cols(iter_info.df[rep_len(seq_len(nrow(iter_info.df)), nrow(samples.df)), ],
                                     samples.df)
    } else {
      warning("No iter_info found")
    }

    # add additional dimension information
    if (length(var_dims) > 0) {
        samples.df <- .attach_dim_info.by_array_index(samples.df, stan_samples, var_dims)
    }
    return (samples.df)
}

vars_statistics <- function(vars_category, stan_stats, stan_samples, vars_info, dim_info) {
    message('Extracting ', vars_category, ' variables...')
    vars_cat_info <- vars_info[[vars_category]]
    samples.df <- stan_samples_frame(stan_samples, vars_cat_info$names, dim_info[vars_cat_info$dims])
    # process convergence information
    stats.df <- dplyr::filter(stan_stats, str_detect(var_name,
                                                     paste0('^(', paste0(vars_cat_info$names, collapse='|'), ')(\\[|$)')))
    stats.df$var <- extract_var(stats.df$var_name)
    stats.df <- dplyr::mutate(stats.df,
                              mean_log2 = mean/log(2),
                              median_log2 = `50%`/log(2),
                              sd_log2 = sd/log(2))
    if (length(vars_cat_info$dims) > 0) {
        stats.df <- .attach_dim_info.by_var_name(stats.df, dim_info[vars_cat_info$dims])
    }
    return (list(samples=samples.df, stats=stats.df))
}

vars_opt_convert <- function(vars_category, opt_results, vars_info, dim_info) {
    message('Extracting ', vars_category, ' variables...')
    vars_cat_info <- vars_info[[vars_category]]

    # extract only samples of selected variables
    # convert arrays of samples of selected variables into data frames,
    # add back single copy of dimensions and real iteration information
    res_mask <- str_detect(names(opt_results$par),
                           paste0('^(', paste0(vars_cat_info$names, collapse='|'), ')(\\[|$)'))
    if (!any(res_mask)) {
        # checkf if the variable is degenerated
        if (prod(sapply(dim_info[vars_cat_info$dims], nrow)) > 0) {
          stop('Variables ', paste0("'", vars_cat_info$names, "'", collapse = ", "),
               ' not found in the opt. results')
        } else {
          # skip the degenerated variable
          return (tibble())
        }
    }
    res.df <- tibble(var_name = names(opt_results$par)[res_mask],
                     mean = opt_results$par[res_mask]) %>%
        dplyr::mutate(var = extract_var(var_name))

    # add additional dimension information
    if (length(vars_cat_info$dims) > 0) {
        miss_mask <- !(vars_cat_info$dims %in% names(dim_info))
        if (any(miss_mask)) {
            stop(sum(miss_mask), ' dimension info are missing: ',
                 paste0(vars_cat_info$dims[miss_mask], collapse = ' '))
        }
        res.df <- .attach_dim_info.by_var_name(res.df, dim_info[vars_cat_info$dims])
    }
    return (res.df)
}

#' @export
prepare_contrast_quantiles <- function(contrastXmetacondition) {
  sel_contrasts.df <- dplyr::filter(contrastXmetacondition, !is.na(quantile_min) & !is.na(quantile_max))
  sel_contrasts_lhs.df <- dplyr::filter(sel_contrasts.df, weight > 0)
  sel_contrasts_rhs.df <- dplyr::filter(sel_contrasts.df, weight < 0)

  list(lhs = rlang::set_names(lapply(seq_along(sel_contrasts_lhs.df$contrast), function(i) c(sel_contrasts_lhs.df$quantile_min[[i]], sel_contrasts_lhs.df$quantile_max[[i]])),
                              rlang::as_character(sel_contrasts_lhs.df$contrast)),
       rhs = rlang::set_names(lapply(seq_along(sel_contrasts_rhs.df$contrast), function(i) c(sel_contrasts_rhs.df$quantile_min[[i]], sel_contrasts_rhs.df$quantile_max[[i]])),
                              rlang::as_character(sel_contrasts_rhs.df$contrast)))
}

calc_contrasts <- function(vars_results, vars_info, dims_info,
                           contrastXmetacondition, conditionXmetacondition.df, contrasts.df,
                           contrastXcondition.df = NULL,
                           mschannel_col = "mschannel_ix",
                           condition_agg_col = "condition", var_names = c('iaction_labu', 'iaction_labu_replCI', 'obs_labu'),
                           obj_dim = "object",
                           obj_id_cols = 'glm_object_ix', var_info_cols = obj_id_cols, val_trans = NULL,
                           condition.quantiles_lhs = c(0, 1), condition.quantiles_rhs = c(0, 1)) {
  contrast_col <- names(dimnames(contrastXmetacondition))[[1]]
  condition_col <- condition_agg_col
  metacondition_col <- names(dimnames(contrastXmetacondition))[[2]]
  contrastXmetacondition.df <- as_tibble(as.table(contrastXmetacondition)) %>%
    dplyr::filter(n != 0) %>% dplyr::rename(weight=n) %>%
    dplyr::inner_join(contrasts.df)
  if (is.null(contrastXcondition.df)) { # autogenerate
    contrastXcondition.df <- contrastXmetacondition.df %>%
      dplyr::inner_join(conditionXmetacondition.df) %>%
      dplyr::arrange_at(c(contrast_col, "contrast_type", metacondition_col, condition_col))
  } else {
    message("Using existing contrastXcondition")
  }
  if (!rlang::has_name(contrastXcondition.df, "is_preserved_condition")) {
    contrastXcondition.df$is_preserved_condition <- FALSE
  }
  for (vars_category in names(vars_results)) {
    vars_cat_subset_info <- vars_info[[vars_category]]
    vars_cat_subset_info$names <- intersect(vars_cat_subset_info$names, var_names)
    if (length(vars_cat_subset_info$names) > 0) {
      message('Filtering ', vars_category, ' to use for contrasts...')
      if (!rlang::has_name(vars_results[[vars_category]]$stats, "condition")) {
        next # FIXME skipping supcondition-related vars
      }
      cond_stats.df <- dplyr::inner_join(vars_results[[vars_category]]$stats,
                                         dplyr::mutate(contrastXcondition.df, is_lhs = weight > 0)) %>%
        dplyr::group_by_at(c(obj_id_cols, contrast_col, metacondition_col,
                             condition_agg_col, "is_lhs", "is_preserved_condition")) %>%
        dplyr::summarize(cond_mean = mean(mean)) %>% dplyr::ungroup()
      cond_agg_stats.df <- dplyr::group_by_at(cond_stats.df, c(obj_id_cols, contrast_col,
                                              metacondition_col, condition_agg_col, "is_lhs")) %>%
        dplyr::summarize(cond_max_mean = max(cond_mean),
                         cond_min_mean = min(cond_mean)) %>%
        dplyr::group_by_at(c(obj_id_cols, contrast_col, metacondition_col, "is_lhs")) %>%
        dplyr::mutate(cond_max_qtile = cume_dist(cond_max_mean) - 1/n(),
                      cond_min_qtile = cume_dist(cond_min_mean) - 1/n()) %>% dplyr::ungroup() %>%
        dplyr::select(!!!c(obj_id_cols, contrast_col, metacondition_col, condition_agg_col),
                      is_lhs, cond_min_qtile, cond_max_qtile, cond_max_mean, cond_min_mean)
      # compose threshold dataframe
      contrast_quantile_thresholds <- function(cond_qtls, is_lhs) {
        res <- tibble(contrast = unique(cond_agg_stats.df[[contrast_col]]),
                      cond_qtile.min_thresh = 0, cond_qtile.max_thresh = 1,
                      is_lhs = is_lhs)
        if (nrow(res) > 0) {
          if (is.list(cond_qtls)) {
            rownames(res) <- res$contrast
            res[names(cond_qtls), 'cond_qtile.min_thresh'] <- sapply(cond_qtls, function(qtl) qtl[[1]])
            rownames(res) <- res$contrast
            res[names(cond_qtls), 'cond_qtile.max_thresh'] <- sapply(cond_qtls, function(qtl) qtl[[2]])
          } else {
            res$cond_qtile.min_thresh <- cond_qtls[[1]]
            res$cond_qtile.max_thresh <- cond_qtls[[2]]
          }
        }
        return(res)
      }
      contr_qtl_thresh.df <- dplyr::bind_rows(contrast_quantile_thresholds(condition.quantiles_lhs, is_lhs = TRUE),
                                              contrast_quantile_thresholds(condition.quantiles_rhs, is_lhs = FALSE))
      cond_stats.df <- dplyr::inner_join(cond_stats.df, cond_agg_stats.df) %>%
        dplyr::left_join(contr_qtl_thresh.df) %>%
        dplyr::mutate(is_accepted = (cond_min_qtile >= cond_qtile.min_thresh) &
                                    ((cond_max_qtile <= cond_qtile.max_thresh) |
                                     (!is_lhs & is_preserved_condition)))
      cond_stats.df <- dplyr::filter(cond_stats.df, is_accepted)
      message('Calculating contrasts for ', vars_category,
              ' variables (', paste0(vars_cat_subset_info$names, collapse=' '), ')...')
      samples.df <- dplyr::semi_join(vars_results[[vars_category]]$samples, cond_stats.df) # exclude unused conditions
      experiment_col <- if (mschannel_col %in% colnames(samples.df)) mschannel_col
                        else if (condition_col %in% colnames(samples.df)) condition_col
                        else condition_col
      metacondition2experiments.df <- dplyr::inner_join(dplyr::semi_join(dplyr::select_at(conditionXmetacondition.df, c(condition_col, metacondition_col)),
                                                                         cond_stats.df),
                                                        dplyr::distinct(dplyr::select(samples.df,
                                                                                      !!!unique(c(condition_col, experiment_col)))))
      if (nrow(samples.df) == 0) {
        # no samples
        vars_results[[vars_category]]$contrast_stats <- NULL
      } else {
        invalid_contrast_types <- !(contrastXmetacondition.df$contrast_type %in% c("filter", "filtering", "comparison"))
        if (any(invalid_contrast_types)) {
          stop("Unsupported contrast types: ",
               paste0(unique(contrastXmetacondition.df$contrast_type[invalid_contrast_types]),
                      collapse=" "))
        }
        # include only the object columns that are included in the processed vars_category
        objs_df <- dplyr::select(dims_info[[obj_dim]], !!!intersect(colnames(vars_results[[vars_category]]$stats),
                                                                    colnames(dims_info[[obj_dim]])))
        # adjust samples w.r.t. condition shift
        if ('condition_shift' %in% colnames(samples.df)) {
          for (shift_name in vars_cat_subset_info$names) {
            samples.df[, shift_name] <- samples.df[, shift_name] - samples.df$condition_shift
          }
        }
        contrast_stats.df <- vars_contrast_stats(samples.df, vars_cat_subset_info$names,
                                                 group_cols = obj_id_cols, condition_col = metacondition_col,
                                                 experiment_col = experiment_col,
                                                 condition2experiments.df = metacondition2experiments.df,
                                                 contrastXcondition = contrastXmetacondition, val_trans = val_trans)
        var_info.df <- vars_results[[vars_category]]$stats %>%
          dplyr::inner_join(conditionXmetacondition.df) %>%
          dplyr::select(!!!setdiff(var_info_cols,
                                   c("index_observation", "iaction_id",
                                     condition_col, condition_agg_col, experiment_col))) %>%
          dplyr::distinct()
        # inject contrast statistics into vars_results
        vars_results[[vars_category]]$contrast_stats <- dplyr::left_join(var_info.df, contrast_stats.df) %>%
          dplyr::mutate(p_value = 2*pmin(prob_nonpos, prob_nonneg),
                        # fake P-values to fit in the plot
                        #prob_nonpos_fake = pmax(if_else(mean < 0, 1E-300, 1E-100), rgamma(n(), shape=1E-2, scale=1E-150), prob_nonpos),
                        #prob_nonneg_fake = pmax(if_else(mean > 0, 1E-300, 1E-100), rgamma(n(), shape=1E-2, scale=1E-150), prob_nonneg),
                        #p_value_fake = 2*pmin(prob_nonpos_fake, prob_nonneg_fake)
          ) %>%
          #dplyr::select(-index_observation, -msrun, -msrun_ix) %>% dplyr::distinct() %>%
          dplyr::left_join(group_by(cond_stats.df, contrast, is_lhs) %>%
                           dplyr::summarise(conditions = str_c(condition, collapse=" ")) %>%
                           dplyr::group_by(contrast) %>%
                           dplyr::summarise(conditions_lhs = conditions[is_lhs],
                                            conditions_rhs = conditions[!is_lhs]) %>%
                           dplyr::ungroup()) %>%
          dplyr::left_join(contrastXmetacondition.df %>% dplyr::filter(contrast %in% rownames(contrastXmetacondition)) %>%
                           dplyr::rename(contrast_weight = weight) %>%
                           dplyr::mutate(metacondition_reported =
                                           case_when(contrast_type %in% c("filter", "filtering") ~ "lhs",
                                                     contrast_type == "comparison" ~ "enriched",
                                                     TRUE ~ NA_character_))) %>%
          dplyr::filter((metacondition_reported == "lhs" & contrast_weight > 0) |
                        (metacondition_reported == "enriched" & contrast_weight * `50%` > 0)) %>%
          dplyr::select(-metacondition_reported, -contrast_type) %>%
          dplyr::inner_join(objs_df)
      }
    }
  }
  return(vars_results)
}

calc_contrasts_subset <- function(vars_results, vars_info, dims_info,
                                  contrasts = NULL, condition_shifts = NULL,
                                  val_trans = NULL,
                                  condition.reported = "lhs",
                                  condition.quantiles_lhs = c(0, 1), condition.quantiles_rhs = c(0, 1),
                                  condition_agg_col = "condition")
{
    # FIXME do once per assembly
    sel_contrastXmetacondition <- contrastXmetacondition[contrasts, , drop=FALSE]
    sel_contrastXmetacondition <- sel_contrastXmetacondition[, colSums(abs(sel_contrastXmetacondition)) != 0, drop=FALSE]
    sel_conditionXmetacondition.df <- dplyr::filter(conditionXmetacondition.df,
                                                    metacondition %in% colnames(sel_contrastXmetacondition))
    if (!is.null(condition_shifts)) {
      sel_conditionXmetacondition.df <- dplyr::inner_join(sel_conditionXmetacondition.df, condition_shifts)
    }
    sel_contrasts.df <- dplyr::filter(contrasts.df, contrast %in% contrasts)

    calc_contrasts(vars_results, vars_info, dims_info,
                   sel_contrastXmetacondition, sel_conditionXmetacondition.df,
                   sel_contrasts.df, val_trans = val_trans,
                   condition_agg_col = condition_agg_col,
                   condition.reported = condition.reported,
                   condition.quantiles_lhs = condition.quantiles_lhs,
                   condition.quantiles_rhs = condition.quantiles_rhs)
}

#' @export
default_contrast_vars <- function(vars_info) {
  unlist(lapply(vars_info, function(vi) {
    rel_dims <- names(dims_info)[sapply(names(dims_info), function(dname) any(str_detect(colnames(dims_info[[dname]]), "^(msrun|mschannel|condition|action)")))]
    if (any(vi$dims %in% rel_dims)) str_subset(vi$names, "(?:_labu)(?:_replCI)?$") else c()
  }))
}

#' @export
process.stan_fit <- function(msglm.stan_fit, dims_info,
                             vars_info = attr(msglm.stan_fit, "msglm_vars_info"),
                             mschannel_col = "msrun_ix",
                             effect_vars = unlist(lapply(vars_info, function(vi) str_subset(vi$names, "_(?:mix)?effect(?:_replCI)?$"))),
                             contrast_vars = default_contrast_vars(vars_info),
                             condition.quantiles_lhs = c(0, 1), condition.quantiles_rhs = c(0, 1),
                             keep.samples=FALSE, min.iteration=NA, chains=NA, verbose=FALSE)
{
  message('Extracting MCMC samples...')
  all_pars <- unlist(sapply(vars_info, function(vi) vi$names))
  msglm.stan_samples <- stan.extract_samples(msglm.stan_fit,
                                             pars=all_pars, min.iteration=min.iteration, chains=chains,
                                             permuted=TRUE)

  message('Computing parameters statistics...')
  msglm.stan_stats <- msglm.stan_fit %>%
    stan.extract_samples(pars = all_pars, min.iteration = min.iteration, chains=chains) %>%
    monitor(print = FALSE) %>% as.data.frame
  msglm.stan_stats$var_name <- rownames(msglm.stan_stats)

  message('Composing results...')
  res <- lapply(names(vars_info), vars_statistics, msglm.stan_stats, msglm.stan_samples,
                vars_info, dims_info)
  names(res) <- names(vars_info)

  # add interaction CI with respect to observations variability
  is_glmm <- "supactions" %in% names(res)
  if (is_glmm) {
    # FIXME what to do for glmm? there's no interaction observations
  } else {
    res$iactions_obsCI <- list(stats = vars_contrast_stats(res$observations$samples,
                                                           c('obs_labu'),
                                                           c('glm_iaction_ix', 'glm_iaction_ix', 'condition'),
                                                           condition_col = NA, contrastXcondition = NULL, mschannel_col) %>%
                                      dplyr::ungroup() %>%
                                      dplyr::select(-index_contrast, -contrast))
  }

  message('Calculating P-values..')
  local({
    for (ctg in names(res)) {
      # unset prior mean for non-effect vars (e.g. obj_effect_sigma)
      if (rlang::has_name(res[[ctg]]$stats, "prior_mean")) {
        res[[ctg]]$stats <<- dplyr::mutate(res[[ctg]]$stats,
                                           prior_mean=if_else(var %in% effect_vars, prior_mean, NA_real_),
                                           prior_mean_log2=if_else(var %in% effect_vars, prior_mean_log2, NA_real_))
      }
      ctg_subset_info <- vars_info[[ctg]]
      ctg_subset_info$names <- intersect(ctg_subset_info$names, effect_vars)

      if (length(ctg_subset_info$names) > 0 & !is.null(res[[ctg]]$samples)) {
        message('Calculating P-values for ', ctg, ' variables...')
        p_value.df <- vars_effect_pvalue(res[[ctg]]$samples, ctg_subset_info, dims_info)
        #print(str(p_value.df))
        #print(str(stats.df))
        res[[ctg]]$stats <<- dplyr::left_join(res[[ctg]]$stats %>%
                             dplyr::mutate(p_value = NULL), p_value.df)
      }
  }})

  message("Calculating contrasts...")
  res <- calc_contrasts(res, vars_info, dims_info,
                        contrastXmetacondition.mtx, conditionXmetacondition.df,
                        dplyr::distinct(dplyr::select(contrastXmetacondition.df, contrast, contrast_type)),
                        contrastXcondition.df = if(exists('contrastXcondition.df')){contrastXcondition.df}else{NULL},
                        var_names = contrast_vars,
                        mschannel_col = mschannel_col,
                        condition.quantiles_lhs = condition.quantiles_lhs,
                        condition.quantiles_rhs = condition.quantiles_rhs)

  message("Removing MCMC samples...")
  if (is.character(keep.samples)) {
    samples_mask <- !(names(res) %in% keep.samples)
  } else {
    samples_mask <- rep.int(!keep.samples, length(res))
  }
  for (ctg_ix in which(samples_mask)) {
    if (verbose) message("Removing ", names(res)[[ctg_ix]], " samples")
    res[[ctg_ix]]$samples <- NULL
  }

  return (res)
}
