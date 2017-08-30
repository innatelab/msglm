require(rstan)
require(insilicoMop)

msglm.prepare_dims_info <- function(model_data, object_cols = NULL)
{
  objs_df <- model_data$objects
  if (!is.null(object_cols)) {
    objs_df <- dplyr::select_(objs_df, .dots=unique(c("glm_object_ix", object_cols)))
  }
  list(iteration = NULL,
    effect = data.frame(effect = all_effects, stringsAsFactors = FALSE),
    repl_effect = data.frame(repl_effect = all_repl_effects, stringsAsFactors = FALSE),
    batch_effect = data.frame(batch_effect = all_batch_effects, stringsAsFactors = FALSE),
    condition = data.frame(condition = all_conditions, index_condition = seq_along(all_conditions),
                           stringsAsFactors = FALSE),
    msrun = dplyr::select(model_data$mschannels, msrun_ix, msrun, condition),
    iaction = dplyr::select(model_data$interactions, glm_iaction_ix, glm_object_ix, iaction_id, condition_ix, condition) %>%
        dplyr::inner_join(objs_df),
    observation = dplyr::select(model_data$ms_data, glm_iaction_ix, glm_object_ix, condition_ix, condition, msrun, msrun_ix) %>%
        dplyr::inner_join(objs_df),
    object = model_data$objects, # use full object information
    object_effect = model_data$object_effects %>%
        dplyr::mutate(glm_object_ix = as.integer(glm_object_ix)) %>%
        dplyr::inner_join(objs_df),
    object_repl_effect = model_data$object_repl_effects %>%
        dplyr::mutate(glm_object_ix = as.integer(glm_object_ix)) %>%
        dplyr::inner_join(objs_df),
    object_batch_effect = model_data$object_batch_effects %>%
        dplyr::mutate(glm_object_ix = as.integer(glm_object_ix)) %>%
        dplyr::inner_join(objs_df)
  )
}

vars_effect_pvalue <- function(samples.df, vars_cat_info, dim_info, tail = c("both", "negative", "positive"))
{
   tail = match.arg(tail)
   group_cols <- paste0( 'index_', vars_cat_info$dims )
   #print( str(samples.df) )
   #print( str(group_cols ) )
   p_value_all_samples <- function( samples ) {
       do.call( rbind, lapply( vars_cat_info$names, function(col) {
           data.frame( var = col,
                       p_value = pvalue_not_zero( samples[[col]], tail = tail ),
                       stringsAsFactors = FALSE )
           } ) )
   }
   p_value.df <- samples.df %>%
        group_by_( .dots = group_cols ) %>%
        do( p_value_all_samples(.) )
   if ('fraction' %in% group_cols) {
     p_value.df$fraction <- as.integer( p_value.df$fraction )
   }
   return ( p_value.df )
}

vars_contrast_stats <- function(samples.df, var_names, group_cols,
                                condition_col, contrastXcondition,
                                experiment_col = condition_col, condition2experiments.df = NULL,
                                nsteps = 100, maxBandwidth = NA,
                                quant.probs = c(0.025, 0.25, 0.50, 0.75, 0.975),
                                val_trans = NULL)
{
  n_samples <- n_distinct(samples.df$unpermuted_ix)
  if (n_samples == 0L) {
    warning("No samples for ", paste0(var_names, collapse=" "), " ", condition_col)
    return(data.frame())
  }
  if (is.na(condition_col)) {
    # no conditions specified, assume the same condition for all experiments
    # FIXME checking user-defined condition2experiments.df and constrastXcondition
    condition2experiments.df <- dplyr::distinct(dplyr::select_(samples.df, .dots=c("experiment" = experiment_col))) %>%
      dplyr::mutate(condition = '__all__')
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
      condition2experiments.df <- dplyr::distinct(dplyr::select_(samples.df, .dots=c("condition" = condition_col, "experiment" = experiment_col)))
    } else {
      if (is.null(condition2experiments.df)) {
        stop("No experiment design: provide condition2experiment.df or ", condition_col, " in the samples.df")
      }
      condition2experiments.df <- dplyr::select_(condition2experiments.df, .dots=c("condition" = condition_col, "experiment" = experiment_col))
    }
  } else {
    condition2experiments.df <- dplyr::distinct(dplyr::select_(samples.df, .dots=c("condition" = condition_col))) %>%
      dplyr::mutate(experiment = condition)
  }

  contrast_stats_all_samples <- function( samples ) {
    #print(str(samples))
    # recode experiment indices to match what is in the samples
    cur_experiments = unique(samples[[experiment_col]])
    cur_cond2expr.df <- dplyr::filter(condition2experiments.df, experiment %in% cur_experiments) %>%
      dplyr::mutate(experiment = factor(as.character(experiment), levels=cur_experiments))
    # recode condition indices
    cur_contrastXcondition <- contrastXcondition[,colnames(contrastXcondition) %in% unique(cur_cond2expr.df$condition),drop=FALSE]
    cur_cond2expr.df <- dplyr::mutate(cur_cond2expr.df,
      condition = factor(as.character(condition), levels=colnames(cur_contrastXcondition)))

    bind_rows(lapply(var_names, function(var_col) {
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
      samples.arr <- matrix(sample_vals, nrow=n_samples)

      #print(str(samples.arr))
      res <- insilicoMop:::ContrastStatistics( 
        samples.arr, as.integer(cur_cond2expr.df$experiment), as.integer(cur_cond2expr.df$condition),
        cur_contrastXcondition,
        nsteps = nsteps, maxBandwidth = maxBandwidth,
        quant_probs = quant.probs ) %>%
        as.data.frame(stringsAsFactors = FALSE, optional = TRUE) %>%
        mutate(var = var_col)
      # TODO se_mean, n_eff, Rhat
      return(res)
    } ) )
  }
  contrast_stats.df <- dplyr::arrange_(samples.df,
                                       .dots = c(group_cols, experiment_col, 'iteration', 'chain') ) %>%
    group_by_( .dots = group_cols ) %>%
    do( contrast_stats_all_samples(.) )
  if ( 'fraction' %in% group_cols ) {
    contrast_stats.df$fraction <- as.integer( contrast_stats.df$fraction )
  }
  if ( !is.null(dimnames(contrastXcondition)[[1]]) ) {
    contrast_stats.df$contrast <- dimnames(contrastXcondition)[[1]][contrast_stats.df$index_contrast]
  }
  return ( contrast_stats.df )
}

.attach_dim_info.by_array_index <- function( samples.df, stan_samples, dim_info )
{
    indexes.df <- do.call( cbind, lapply( seq_along(dim_info), function(dim_ix) {
      slice.index( stan_samples[[ colnames(samples.df)[[ncol(samples.df)]] ]], dim_ix+1 ) %>% as.vector()
    } ) ) %>% as.data.frame()
    colnames(indexes.df) <- paste0( 'index_', names(dim_info) )
    samples.df <- cbind( samples.df, indexes.df )
    for ( dim_ix in seq_along(dim_info) ) {
      dim_name <- names(dim_info)[[dim_ix]]
      if ( !is.null(dim_info[[dim_ix]]) ) {
        message( "Adding ", dim_name, " information" )
        samples.df <- cbind( samples.df, dim_info[[dim_ix]][ indexes.df[[dim_ix]], , drop=FALSE ] )
      } else {
        warning( "No ", dim_name, " information provided" )
      }
    }
    return ( samples.df )
}

.attach_dim_info.by_var_name <- function( data.df, dim_info )
{
    res.df <- data.df
    res_index <- extract_index( res.df$var_name )
    dims.df <- data.frame( ix = seq_along(dim_info),
                           name = names(dim_info),
                           stringsAsFactors = FALSE ) %>%
      group_by( name ) %>%
      mutate( local_ix = row_number(), suffix = if ( n() > 1 ) paste0( '.', local_ix ) else '' ) %>%
      ungroup()
    indexes.df <- do.call( cbind, lapply( seq_along(dim_info), function(dim_ix) {
      res_index[,dim_ix]
    } ) )
    colnames(indexes.df) <- paste0( 'index_', names(dim_info) )

    for ( dim_ix in seq_len(nrow(dims.df)) ) {
        if ( dim_ix > ncol(res_index) ) {
            stop( 'Data contain less dimensions than the dim_info' )
        }
        dim_name <- dims.df$name[[dim_ix]]
        if ( dims.df$local_ix[[dim_ix]] > 1 ) {
          warning( 'Dimension ', dim_name, ' is not unique, using suffixes' )
        }
        if ( !is.null(dim_info[[dim_ix]]) ) {
            res.df[[dim_name]] <- NULL
            if ( nrow(dim_info[[dim_ix]]) != max(res_index[,dim_ix]) ) {
                stop( 'Dimension ', dim_ix, ' info contains ', nrow(dim_info[[dim_ix]]),
                      ' elements, the data contain ', max(res_index[,dim_ix]) )
            }
            res_dim_info.df <- dim_info[[dim_ix]][ res_index[,dim_ix], , drop=FALSE ]
            if ( dims.df$suffix[[dim_ix]] != '' ) {
              colnames(res_dim_info.df) <- paste0( colnames(res_dim_info.df), dims.df$suffix[dim_ix] )
            }
            res.df <- cbind( res.df, res_dim_info.df )
        }
    }
    return ( cbind( indexes.df, res.df ) )
}

stan_samples_frame <- function(stan_samples, var_names, var_dims) {
    # extract only samples of selected variables
    # convert arrays of samples of selected variables into data frames,
    # add back single copy of dimensions and real iteration information
    if (any(!(var_names %in% names(stan_samples)))) {
      warning("Variable(s) samples not found: ", paste0(setdiff(var_names, names(stan_samples)), collapse=", "))
    }
    avail_var_names <- intersect(var_names, names(stan_samples))
    samples.df <- do.call(data.frame, lapply(avail_var_names, function(var_name) {
      as.vector(stan_samples[[var_name]])
    }))
    if (ncol(samples.df) > 0) {
        colnames(samples.df) <- avail_var_names
        if ('iter_info' %in% names(attributes(stan_samples))) {
          samples.df <- cbind(attr(stan_samples, 'iter_info'), samples.df)
        } else {
          warning("No iter_info found")
        }

        # add additional dimension information
        if (length(var_dims) > 0) {
            samples.df <- .attach_dim_info.by_array_index(samples.df, stan_samples, var_dims)
        }
    }
}

vars_statistics <- function(vars_category, stan_stats, stan_samples, vars_info, dim_info) {
    message( 'Extracting ', vars_category, ' variables...' )
    vars_cat_info <- vars_info[[vars_category]]
    samples.df <- stan_samples_frame(stan_samples, vars_cat_info$names, dim_info[vars_cat_info$dims])
    # process convergence information
    stats.df <- subset( stan_stats, grepl( paste0( '^(', paste0( vars_cat_info$names, collapse='|' ), ')(\\[|$)' ), var_name ) )
    stats.df$var <- extract_var( stats.df$var_name )
    if ( nrow(stats.df) > 0 && length(vars_cat_info$dims) > 0 ) {
        stats.df <- .attach_dim_info.by_var_name( stats.df, dim_info[vars_cat_info$dims] )
    }
    return ( list( samples=samples.df, stats=stats.df ) )
}

vars_opt_convert <- function(vars_category, opt_results, vars_info, dim_info) {
    message( 'Extracting ', vars_category, ' variables...' )
    vars_cat_info <- vars_info[[vars_category]]

    # extract only samples of selected variables
    # convert arrays of samples of selected variables into data frames,
    # add back single copy of dimensions and real iteration information
    res_mask <- grepl( paste0( '^(', paste0( vars_cat_info$names, collapse='|' ), ')(\\[|$)' ), names(opt_results$par) )
    if ( !any(res_mask) ) {
        # checkf if the variable is degenerated
        if ( prod( sapply( dim_info[vars_cat_info$dims], nrow ) ) > 0 ) {
          stop( 'Variables ', paste0( "'", vars_cat_info$names, "'", collapse = ", " ), ' not found in the opt. results' ) 
        } else {
          # skip the degenerated variable
          return ( data.frame() )
        }
    }
    res.df <- data.frame(
            var_name = names(opt_results$par)[res_mask],
            mean = opt_results$par[ res_mask ],
            stringsAsFactors = FALSE ) %>%
      mutate( var = extract_var( var_name ) )

    # add additional dimension information
    if ( length(vars_cat_info$dims) > 0 ) {
        miss_mask <- !(vars_cat_info$dims %in% names(dim_info))
        if ( any(miss_mask) ) {
            stop( sum(miss_mask), ' dimension info are missing: ',
                  paste0( vars_cat_info$dims[miss_mask], collapse = ' ' ) )
        }
        res.df <- .attach_dim_info.by_var_name( res.df, dim_info[vars_cat_info$dims] )
    }
    return ( res.df )
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
        # include only the object columns that are included in the processed vars_category
        objs_df <- dplyr::select_(dims_info$object, .dots=intersect(colnames(vars_results[[vars_category]]$stats),
                                                                    colnames(dims_info$object)))
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
    sel_contrastXcondition.df <- dplyr::semi_join(contrastXcondition.df, sel_conditionXmetacondition.df) %>%
        dplyr::filter(contrast %in% contrasts)
    calc_contrasts(vars_results, vars_info, dims_info,
                   sel_contrastXmetacondition, sel_conditionXmetacondition.df,
                   sel_contrastXcondition.df, val_trans = val_trans,
                   condition_agg_col = condition_agg_col,
                   condition.reported = condition.reported,
                   condition.quantiles_lhs = condition.quantiles_lhs,
                   condition.quantiles_rhs = condition.quantiles_rhs)
}

process.stan_fit <- function(msglm.stan_fit, dims_info, keep.samples=FALSE, verbose=FALSE)
{
  message( 'Extracting MCMC samples...' )
  msglm.stan_samples <- stan.extract_samples(msglm.stan_fit,
                                             pars = unlist(sapply(msglm.vars_info, function(vi) vi$names)),
                                             min.iteration = min.iteration,
                                             permuted = TRUE)

  message('Computing parameters statistics...')
  msglm.stan_stats <- msglm.stan_fit %>%
    stan.extract_samples(pars = unlist(sapply(msglm.vars_info, function(vi) vi$names)),
                         min.iteration = min.iteration) %>%
    monitor(print = FALSE) %>% as.data.frame
  msglm.stan_stats$var_name <- rownames(msglm.stan_stats)

  message( 'Composing results...' )
  res <- lapply(names(msglm.vars_info), vars_statistics, msglm.stan_stats, msglm.stan_samples,
                msglm.vars_info, dims_info)
  names(res) <- names(msglm.vars_info)

  # add interaction CI with respect to observations variability
  res$iactions_obsCI <- list(stats = vars_contrast_stats(res$observations$samples,
                                                         c('obs_labu'),
                                                         c('glm_iaction_ix', 'protgroup_id', 'condition'),
                                                         condition_col = NA, contrastXcondition = NULL, 'msrun_ix') %>%
                                    dplyr::ungroup() %>%
                                    dplyr::select(-index_contrast, -contrast))

  message( 'Calculating P-values..' )
  local({
    for (ctg in names(res)) {
      ctg_subset_info <- msglm.vars_info[[ctg]]
      ctg_subset_info$names <- intersect(ctg_subset_info$names, c('obj_effect', 'obj_effect_replCI', 'obj_batch_effect_unscaled'))

      if (length(ctg_subset_info$names) > 0 & !is.null(res[[ctg]]$samples)) {
        message( 'Calculating P-values for ', ctg, ' variables...' )
        p_value.df <- vars_effect_pvalue(res[[ctg]]$samples, ctg_subset_info, dims_info)
        #print(str(p_value.df))
        #print(str(stats.df))
        res[[ctg]]$stats <<- left_join(res[[ctg]]$stats %>% mutate(p_value = NULL), p_value.df)
      }
  }})

  message("Calculating contrasts...")
  res <- calc_contrasts(res, msglm.vars_info[c("iactions", "observations", "object_batch_shifts")], dims_info,
                        contrastXmetacondition.mtx, conditionXmetacondition.df, contrastXcondition.df,
                        mschannel_col = 'msrun')

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
