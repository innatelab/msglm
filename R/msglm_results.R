#' Generate the dimensions information for the MSGLM report.
#'
#' @param model_data MSGLM model data that is used to get the dimensions info.
#' @param object_cols Non-standard object information columns to
#'                    include in the reports.
#'
#' @returns The named list of data frames with the dimension information.
#'          The element names are the names of the dimensions,
#'          Each data frame has ``index_<dimname>`` key,
#'          which matches the index along the corresponding dimension in
#'          the MSGLM model.
#'
#' @seealso [process.stan_fit()]
#'
#' @export
msglm_dims <- function(model_data)
{
  model_def <- model_data$model_def
  # FIXME use model_data class
  is_glmm <- "mixeffects" %in% names(model_data)
  xaction_ix_col <- if (is_glmm) "glm_supaction_ix" else "index_interaction"
  xdition_ix_col <- if (is_glmm) "supcondition_ix" else "index_condition"
  xdition_col <-  if (is_glmm) "supcondition" else "condition"

  objs_df <- model_data$objects
  res <- list(iteration = NULL,
    mschannel = model_data$mschannels,
    observation = model_data$observations %>%
        dplyr::distinct() %>%
        dplyr::inner_join(objs_df, by=c("index_object", "object_id")),
    object = model_data$objects, # use full object information
    object_effect = model_data$object_effects %>%
        dplyr::inner_join(objs_df, by=c("index_object", "object_id")) %>%
        dplyr::inner_join(model_def$effects, by=c("index_effect", "effect")),
    object_batch_effect = model_data$object_batch_effects %>%
        dplyr::inner_join(objs_df, by=c("index_object", "object_id"))
  )
  if ("subobjects" %in% names(model_data)) {
    res$subobject <- model_data$subobjects
    res$subobject_batch_effect <- model_data$subobject_batch_effects %>%
        dplyr::left_join(res$subobject, by = c('index_subobject', 'subobject_id'))
  }
  res$mscalib <- dplyr::select(model_data$mschannels, index_mscalib,
                               any_of("instrument")) %>%
    dplyr::distinct()
  res$interaction <- dplyr::select(model_data$interactions, index_interaction,
                                   index_object, iaction_id, index_condition, condition, is_virtual) %>%
    dplyr::inner_join(objs_df, by="index_object")
  if (is_glmm) {
    res$object_mixeffect <- dplyr::mutate(model_data$mixeffects, tmp="a") %>%
      dplyr::left_join(dplyr::mutate(objs_df, tmp="a")) %>%
                       dplyr::select(-tmp)
    res$object_mixeffect <- maybe_rename(res$object_mixeffect, c("prior_mean" = "mean", "prior_tau" = "tau"))
    if (!rlang::has_name(res$object_mixeffect, "prior_mean")) { # set the default min to 0
      message("mixeffects$prior_mean missing, setting to 0")
      res$object_mixeffect$prior_mean <- 0.0
    }
    res$object_mixcoef <- dplyr::mutate(model_data$mixcoefs, tmp="a") %>%
      dplyr::left_join(dplyr::mutate(objs_df, tmp="a")) %>%
      dplyr::select(-tmp)
    res$supaction <- dplyr::select(model_data$superactions, glm_supaction_ix,
                                   index_object, supaction_id, supcondition_ix, supcondition, is_virtual) %>%
        dplyr::inner_join(objs_df)
  }
  res <- lapply(setNames(names(res), names(res)), function(dimname) {
    dim_df <- res[[dimname]]
    if (!is.null(dim_df)) {
      # check index correctness
      index_col <- paste0('index_', dimname)
      checkmate::assert_set_equal(dim_df[[index_col]], seq_len(nrow(dim_df)), ordered=TRUE)
      # remove unused dimensions (everything except index_object and index_interaction)
      dim_df <- dplyr::select(dim_df, -(starts_with("index_") & !any_of(c(index_col, "index_interaction", "index_object", "index_subobject"))))
    }
    return(dim_df)
  })
  return(res)
}

#' Calculate contrasts (linear combination of model effects) using the posterior
#' MCMC draws.
#'
#' The procedure allows advanced calculation of contrasts based on contrast-specific
#' grouping of the variables (provided by `vargroups`). In that case, the
#' contrasts for all possible combinations of variables in the groups are being
#' calculated and combined together to form the posterior distribution of the
#' contrast.
#' For example, consider the group \eqn{A} contains variables \eqn{a_1} and \eqn{a_2},
#' and the group \eqn{B} is \eqn{\{b_1, b_2\}}. For the \eqn{A - B} contrast,
#' the procedure will calculate 4 linear combinations: \eqn{a_1 - b_1}, \eqn{a_1 - b_2},
#' \eqn{a_2 - b_1} and \eqn{a_2 - b_2}, then it will combine the resulting 4 sets of
#' MCMC draws and proceed with calculating the posterior summary.
#'
#' @param vars_draws MCMC draws for the model variables to use for contrast
#'     calculation in `posterior::draws_array` 3D array format.
#' @param vars_stats pre-calculated summary statistics for MCMC draws
#' @param vargroups grouped data frame that defines how variables are grouped.
#'     The variables are identified by `index_varspec` column, which corresponds
#'     to the index of the variable in `vars_draws` array.
#'     The existing grouping defines how the variables would be grouped.
#'     The grouping should contain the column having the same name as
#'     the columns dimension name of the `vargroupXcontrast` matrix.
#'     For contrast-specific grouping of the variables the grouping
#'     can optionally include the `contrast` column, and its values should
#'     match the `contrast` column of `contrasts` data frame.
#' @param vargroupXcontrast the matrix specifying contrasts. The rows are
#'       contrasts (the row names should match the `contrast` column of `contrasts`),
#'       the columns are variable groups (the column names should match the
#'       corresponding column of `vargroups`)
#' @param contrasts data frame with all contrasts. The `contrast` is the
#'       obligatory contrast identifier, the `offset_col` column specifies
#'       what the contrast should be compared with for significance testing.
#' @param offset_col (defaults to `offset`) column of the `contrasts` frame with the offsets that
#'       should be used for significance testing (i.e. the resuling
#'       contrast is adjusted by `-offset`)
#' @param tail which comparison to do, one of `both` (the default),
#'             `negative` (\eqn{P(X \leq t)}) or `positive` (\eqn{P(X \geq t)}),
#'             where \eqn{t} is the `prior_mean` of the variable \eqn{X}.
#'             For `both` the double of the minimal of the two p-values is given.
#' @param summary also calculate summary statistic (mean, quartiles etc)
#' @param method which method for P-value calculation to use:
#'               * `draws` (approximate contrast distribution by MCMC samples; slow)
#'               * `normal` (approximate contrast distribution by normal distribution; fast)
#' @param nsteps (defaults to 100) how many bins to use for the calculation of p-values
#' @param maxBandwidth constrain the rule-of-thumb bandwidth for the posterior distribution
#'       if it is above the specified limit
#'
#' @export
vars_contrast_stats <- function(vars_draws, vars_stats, vargroups,
                                vargroupXcontrast, contrasts,
                                offset_col="offset",
                                tail = c("both", "negative", "positive"),
                                summary=TRUE,
                                method=c("draws", "normal"),
                                nsteps = 100L, maxBandwidth = NA_real_,
                                mlog10pvalue_threshold = 10.0,
                                mlog10pvalue_hard_threshold_factor = 3.0)
{
  method <- match.arg(method)
  tail <- match.arg(tail)

  if (dplyr::n_groups(vargroups) < nrow(vargroupXcontrast)) {
    stop("Number of vargroups (", dplyr::n_groups(vargroups),
         ") less than the rows of vargroupXcontrast (", nrow(vargroupXcontrast), ")")
  }
  contrasts <- dplyr::mutate(contrasts,
                             `__contrast_ix__` = match(as.character(contrast), colnames(vargroupXcontrast)))
  if (any(is.na(contrasts$`__contrast_ix__`))) {
    stop(sum(is.na(contrasts$`__contrast_ix__`)), " contrasts(s) not defined: ",
         paste0(dplyr::filter(contrasts, is.na(`__contrast_ix__`)), collapse=", "))
  }

  vargroupid_col = names(dimnames(vargroupXcontrast))[[1]]
  if (is.na(vargroupid_col)) {
    stop("No name for columns dimension of vargroupXcontrast found")
  }
  if (!rlang::has_name(vargroups, vargroupid_col)) {
    stop("No vargroup identifier column (", vargroupid_col, ") found in vargroups data frame")
  }
  # add contrast_ix to vargroup (each vargroup is allowed to have contrast-specific composition)
  if (!rlang::has_name(vargroups, 'contrast')) {
    # vargroups are not contrast-specific, add (redundant) contrast information
    vargroups <- dplyr::full_join(vargroups, dplyr::select(contrasts, `__contrast_ix__`), by = character())
  } else {
    vargroups <- dplyr::inner_join(vargroups, dplyr::select(contrasts, contrast, `__contrast_ix__`), by = "contrast")
  }
  vargroups <- dplyr::mutate(vargroups, `__vargroup_ix__` = match(as.character(!!sym(vargroupid_col)),
                                                                  rownames(vargroupXcontrast))) %>%
    dplyr::group_by(`__vargroup_ix__`, `__contrast_ix__`, .add=TRUE)
  vargroup_info.df <- dplyr::group_keys(vargroups)
  if (!rlang::has_name(vargroup_info.df, vargroupid_col)) {
    stop("No vargroup identifier column (", vargroupid_col, ") found in vargroups grouping variables")
  }
  if (any(is.na(vargroup_info.df$`__vargroup_ix__`))) {
    stop(sum(is.na(vargroup_info.df$`__vargroup_ix__`)), " vargroup(s) not defined: ",
         paste0(dplyr::filter(vargroup_info.df, is.na(`__vargroup_ix__`)), collapse=", "))
  }
  contrast_offsets <- rlang::set_names(dplyr::pull(contrasts, !!offset_col), contrasts$contrast)
  contrast_stats <- if (method == "draws") {
        ContrastStatistics_draws(vars_draws,
            vargroups$index_varspec, vargroups$`__vargroup_ix__`, vargroups$`__contrast_ix__`,
            vargroupXcontrast, contrast_offsets[colnames(vargroupXcontrast)],
            nsteps = nsteps, maxBandwidth = maxBandwidth,
            mlog10pvalue_threshold = mlog10pvalue_threshold,
            mlog10pvalue_hard_threshold_factor = mlog10pvalue_hard_threshold_factor,
            summaryfun = if (summary) function(draws) {
                posterior::as_draws_array(draws) %>%
                posterior::summarise_draws(posterior_summary_metrics) %>%
                dplyr::select(-variable)
            } else NULL)
  } else if (method == "normal") {
        summaryfun = if (summary) function(perm_means, perm_variances) {
          mean_means <- mean(perm_means)
          mean_sds <- mean(sqrt(perm_variances))
          qtls = qnorm(c(0.025, 0.25, 0.75, 0.975), mean_means, mean_sds)
          tibble(mean = mean_means,
                 median = mean_means,
                 sd = mean_sds,
                 mad = mean_sds*qnorm(0.75),
                 sd_of_means = sd(perm_means),
                 q2.5 = qtls[[1]],
                 q25 = qtls[[2]],
                 q75 = qtls[[3]],
                 q97.5 = qtls[[4]])
        } else NULL
        ContrastStatistics_normal(vars_stats$mean, vars_stats$sd,
            vargroups$index_varspec, vargroups$`__vargroup_ix__`, vargroups$`__contrast_ix__`,
            vargroupXcontrast, contrast_offsets[colnames(vargroupXcontrast)],
            mlog10pvalue_threshold = mlog10pvalue_threshold,
            mlog10pvalue_hard_threshold_factor = mlog10pvalue_hard_threshold_factor,
            summaryfun = summaryfun)
  }
  res <- dplyr::left_join(dplyr::select(contrasts, -dplyr::matches("^[lr]hs_quantile")),
                          contrast_stats, by = "__contrast_ix__") %>%
    dplyr::select(-`__contrast_ix__`) %>%
    dplyr::mutate(p_value = dplyr::case_when(
                                  tail == "negative" ~ prob_nonpos,
                                  tail == "positive" ~ prob_nonneg,
                                  # 2x correction as both tails are tested
                                  tail == "both" ~ 2*pmin(0.5, prob_nonneg, prob_nonpos),
                                  TRUE ~ NA_real_))
  return (res)
}

vars_identity_contrast_stats <- function(vars_draws, vars_stats,
                                         varspecs, groups_info,
                                         group_idcol, group_cols=character(),
                                         offset_col='offset',
                                         ...) {
    groups_df <- dplyr::select(groups_info, !!group_idcol, any_of(offset_col)) %>% dplyr::distinct()
    if (!rlang::has_name(groups_df, offset_col)) {
      groups_df <- dplyr::mutate(groups_df, !!sym(offset_col):=0)
    }
    group_ids <- dplyr::pull(groups_df, !!group_idcol)
    groups_diag <- diag(nrow = length(group_ids), ncol = length(group_ids))
    dimnames(groups_diag) <- list(group_ids, group_ids) %>%
        rlang::set_names(c(group_idcol, "contrast"))
    res <- vars_contrast_stats(vars_draws, vars_stats,
                        dplyr::inner_join(varspecs, groups_info,
                                          by=c("var_index",
                                               intersect(colnames(varspecs), colnames(groups_info))) %>%
                                              unique()) %>%
                        dplyr::group_by_at(c(group_idcol, group_cols)),
                        vargroupXcontrast = groups_diag,
                        contrasts = dplyr::rename(groups_df, contrast=!!group_idcol),
                        offset_col=offset_col,
                        ...) %>%
                        dplyr::rename(!!group_idcol := contrast)
    if (!rlang::has_name(groups_df, offset_col)) {
      res <- dplyr::select(res, -!!sym(offset_col))
    }
    return(res)
}

#' Compare the posterior of each variable in varspecs with its `prior_mean`
#' and return the p-value for the significance of the difference (w.r.t `tail`)
#'
#' @param vars_draws MCMC draws from the posterior in `posterior::draws_array` format
#' @param varspecs Stan variables (with indices specified) to calculate P-values for
#' @param ... p-value calculate options to pass to [vars_contrast_stats()]
#'
#' @returns data frame with the p-value per each varspec
#' @seealso [vars_contrast_stats()]
#' @export
vars_pvalues <- function(vars_draws, vars_stats, varspecs, ...){
  varspecs <- dplyr::mutate(varspecs,
                            has_mcmc_draws = varspec %in% dimnames(vars_draws)$variable)
  if (!all(varspecs$has_mcmc_draws)) {
    warning("Missing MCMC samples for vars: ",
            paste0(varspecs$varspec[!varspecs$has_mcmc_draws], collapse=" "))
  }

  res_df <- vars_identity_contrast_stats(vars_draws, vars_stats,
              dplyr::filter(varspecs, has_mcmc_draws), varspecs,
              group_idcol = "varspec", offset_col = "prior_mean",
              summary=FALSE, ...)
  return(res_df)
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
    res.df <- tibble(varspec = names(opt_results$par)[res_mask],
                     mean = opt_results$par[res_mask]) %>%
        dplyr::mutate(var = extract_var(varspec))

    # add additional dimension information
    if (length(vars_cat_info$dims) > 0) {
        miss_mask <- !(vars_cat_info$dims %in% names(dim_info))
        if (any(miss_mask)) {
            stop(sum(miss_mask), ' dimension info are missing: ',
                 paste0(vars_cat_info$dims[miss_mask], collapse = ' '))
        }
        res.df <- .attach_dim_info_by_varspec(res.df, dim_info[vars_cat_info$dims])
    }
    return (res.df)
}

#
append_contrasts_stats <- function(vars_results, standraws, stanstats, varspecs,
        metaconditionXcontrast, contrasts.df, conditionXcontrast.df,
        condition_agg_col = "condition", # filtering is based on pregrouping quantiles in metacondition using this column (per contrast)
        object_cols = 'index_object', metacondition_cols = c(),
        group_cols = c(), method = c("draws", "normal")
){
  method <- match.arg(method)
  metacondition_col <- names(dimnames(metaconditionXcontrast))[[1]]
  contrast_col <- names(dimnames(metaconditionXcontrast))[[2]]
  condition_col <- 'condition'
  cutoff_quantiles.df <- dplyr::bind_rows(
    dplyr::transmute(contrasts.df, contrast, is_lhs = TRUE,
                     qtile_min = lhs_quantile_min,
                     qtile_max = lhs_quantile_max),
    dplyr::transmute(contrasts.df, contrast, is_lhs = FALSE,
                     qtile_min = rhs_quantile_min,
                     qtile_max = rhs_quantile_max))

  # generate contrast reports
  contrast_results <- dplyr::group_by(varspecs$spec_info, category) %>%
  dplyr::group_map(.keep = TRUE, function(cat_varspecs, cat) {
    message('  * ', cat$category,
            ' variables: ', paste0(unique(cat_varspecs$var), collapse=', '), '...')
    cat_info <- varspecs$cats_info[[cat$category]]
    cat_vars_stats <- dplyr::semi_join(vars_results[[cat$category]]$stats,
                                       dplyr::select(cat_varspecs, varspec), by = "varspec")
    cat_cols <- intersect(colnames(cat_vars_stats), colnames(cat_info))

    conditionXcontrast_pregroup_stats.df <- dplyr::inner_join(cat_vars_stats,
                                                              dplyr::mutate(conditionXcontrast.df, is_lhs = weight > 0),
                                                              by = condition_col) %>%
      dplyr::group_by_at(c("var", group_cols, contrast_col, metacondition_col, condition_agg_col,
                           "is_lhs", "is_preserved_condition")) %>%
      dplyr::summarize(var_pregroup_max = max(mean), var_pregroup_min = min(mean)) %>%
      dplyr::group_by_at(c("var", group_cols, contrast_col, metacondition_col, "is_lhs")) %>%
      dplyr::mutate(var_pregroup_max_qtile = cume_dist(var_pregroup_max) - 1/n(),
                    var_pregroup_min_qtile = cume_dist(var_pregroup_min) - 1/n()) %>%
      dplyr::ungroup() %>% dplyr::inner_join(cutoff_quantiles.df, by = c('contrast', 'is_lhs')) %>%
      dplyr::mutate(is_accepted = (var_pregroup_min_qtile >= qtile_min) &
                                  ((var_pregroup_max_qtile <= qtile_max) |
                                  # FIXME use contrast_type for is_lhs check?
                                   (!is_lhs & is_preserved_condition)))

    vargroups.df <- dplyr::inner_join(cat_varspecs, cat_info, by = "var_index") %>%
      dplyr::inner_join(dplyr::filter(conditionXcontrast_pregroup_stats.df, is_accepted), by = c("var", group_cols, condition_col)) %>%
      dplyr::select_at(c("var", "index_varspec", "category", group_cols, contrast_col, metacondition_col, condition_col)) %>%
      dplyr::group_by_at(c("var", "category", group_cols)) %>%
      dplyr::group_modify(~ vars_contrast_stats(standraws, stanstats,
                              vargroups = dplyr::group_by_at(.x, c(metacondition_col, contrast_col)),
                              vargroupXcontrast = metaconditionXcontrast,
                              contrasts = contrasts.df, method = method))
    return (vargroups.df)
  })
  # add contrast reports to the var_results
  for (cat_contrasts in contrast_results) {
    cat <- cat_contrasts$category[[1]]
    cat_results <- vars_results[[cat]]
    cat_results$contrast_stats <- cat_contrasts
    vars_results[[cat]] <- cat_results
  }
  return(vars_results)
}

#' @export
default_contrast_vars <- function(vars_info) {
  intersect(c('iaction_labu', 'iaction_labu_replCI'),
            unlist(lapply(vars_info, function(vi) vi$names )))
  # obs_labu is skipped because it's much more expensive to compute and Rhat statistics doesn't make sense
}

#' @export
process.stan_fit <- function(msglm.stan_fit, model_data, dims_info = msglm_dims(model_data),
                             vars_info = attr(msglm.stan_fit, "msglm_vars_info"),
                             effect_vars = unlist(lapply(vars_info, function(vi) str_subset(vi$names, "_(?:mix)?effect(?:_replCI)?$"))),
                             contrast_vars = default_contrast_vars(vars_info), contrast_group_cols = c(),
                             condition_agg_col = "condition", object_cols = setdiff(colnames(dims_info$object), "index_object"),
                             min.iteration=NA, chains=NA, verbose=model_data$model_def$verbose,
                             contrast_method = c("draws", "normal"))
{
  contrast_method <- match.arg(contrast_method)
  model_def <- model_data$model_def
  message('Extracting MCMC samples...')
  all_vars <- unlist(sapply(vars_info, function(vi) vi$names))
  avail_vars = intersect(all_vars, msglm.stan_fit$metadata()$stan_variables) # some variables, e.g. suo_batch_shifts might be empty
  # TODO: use draws_rvars because it's faster?
  msglm.stan_draws <- msglm.stan_fit$draws(variables = avail_vars, format="draws_array"
                       #inc_warmup = !is.na(min.iteration)
  )
  if (!is.na(chains)) {
    msglm.stan_draws <- posterior::subset_draws(msglm.stan_draws, chain = chains)
  }
  if (!is.na(min.iteration)) {
    msglm.stan_draws <- posterior::subset_draws(msglm.stan_draws,
                                                iteration = min.iteration:(stanfit$metadata()$iter_sampling))
  }
  message('Collecting variables information...')
  varspecs <- process_varspecs(dimnames(msglm.stan_draws)$variable, vars_info, dims_info)

  message('Computing variables summary statistics...')
  msglm.stan_stats <- msglm.stan_fit$summary(variables = avail_vars, posterior_summary_metrics) %>%
    # FIXME, min.iteration = min.iteration, chains=chains) %>%
    dplyr::rename(varspec = variable)

  message('Composing results...')
  avail_cats = unique(varspecs$spec_info$category)
  res <- lapply(avail_cats, function(cat) {
    message('  * ', cat, ' variables...')
    # process convergence information
    cat_info.df <- varspecs$cats_info[[cat]]
    cat_varspecs.df <- dplyr::filter(varspecs$spec_info, category == cat)
    res.df <- msglm.stan_stats %>%
      dplyr::inner_join(cat_varspecs.df, by = "varspec") %>%
      dplyr::left_join(cat_info.df, by = "var_index") %>%
      dplyr::select(-category)
    cat_eff_varspecs.df <- dplyr::filter(cat_varspecs.df, var %in% effect_vars) %>%
      dplyr::left_join(dplyr::select(cat_info.df, var_index, any_of("prior_mean")), by = "var_index")
    if (nrow(cat_eff_varspecs.df) > 0) {
      message('    - calculating P-values for: ',
              paste0(unique(cat_eff_varspecs.df$var), collapse=', '), '...')
      p_values.df <- vars_pvalues(msglm.stan_draws, msglm.stan_stats,
                                  cat_eff_varspecs.df, method=contrast_method)
      res.df <- dplyr::left_join(res.df, dplyr::select(p_values.df, -any_of("prior_mean")), by="varspec") %>%
        dplyr::select(-index_varspec, -var_index)
    }
    list(stats = res.df)
  })
  names(res) <- avail_cats

  # add interaction CI with respect to observations variability
  is_glmm <- "supactions" %in% names(res)
  if (is_glmm) {
    # FIXME what to do for glmm? there's no interaction observations
  } else if ('obs_labu' %in% all_vars) {
    # we combine posteriors from all observations of a given interaction to get yet another
    # estimate of interaction abundance
    # we reuse(abuse) the contrast calculation for that -- just to group the appropriate draws
    # and get the summary statistics, but we don't need the contrasts
    message("  * obs_labu aggregate statistics...")
    res$iactions_obsCI <- list(stats = vars_identity_contrast_stats(
            msglm.stan_draws, msglm.stan_stats,
            dplyr::filter(varspecs$spec_info, var == 'obs_labu'),
            varspecs$cats_info$observations,
            group_idcol = "index_interaction",
            method=contrast_method))
  }

  message("Calculating contrasts...")
  contrasts.df <- model_def$contrasts
  # set default quantile thresholds, if not defined
  if (!rlang::has_name(contrasts.df, "lhs_quantile_min")) {
    contrasts.df$lhs_quantile_min <- 0.0
  }
  if (!rlang::has_name(contrasts.df, "lhs_quantile_max")) {
    contrasts.df$lhs_quantile_max <- 1.0
  }
  if (!rlang::has_name(contrasts.df, "rhs_quantile_min")) {
    contrasts.df$rhs_quantile_min <- 0.0
  }
  if (!rlang::has_name(contrasts.df, "rhs_quantile_max")) {
    contrasts.df$rhs_quantile_max <- 1.0
  }

  metaconditionXcontrast.df <- as.data.frame.table(model_def$metaconditionXcontrast, responseName="weight") %>%
    dplyr::filter(weight != 0) %>%
    dplyr::inner_join(contrasts.df, by="contrast")
  conditionXcontrast.df <- as.data.frame.table(model_def$conditionXmetacondition, responseName="part_of") %>%
    dplyr::filter(part_of) %>% dplyr::select(-part_of) %>%
    dplyr::inner_join(dplyr::select(metaconditionXcontrast.df, metacondition, contrast, weight, contrast_type),
                      by="metacondition")
  if (rlang::has_name(model_def, "conditionXcontrast")) {
    old_nrow <- nrow(conditionXcontrast.df)
    join_cols <- intersect(intersect(colnames(conditionXcontrast.df), colnames(model_def$conditionXcontrast)),
                           c("contrast", "contrast_type", "condition", "metacondition", "weight"))
    conditionXcontrast.df <- dplyr::left_join(conditionXcontrast.df, model_def$conditionXcontrast, by=join_cols)
    if (old_nrow != nrow(conditionXcontrast.df)) {
      stop("Incompatible model_def$conditionXcontrast, check that it matches metaconditionXcondition and metaconditionXcontrast")
    }
  }
  if (!rlang::has_name(conditionXcontrast.df, "is_preserved_condition")) {
    conditionXcontrast.df$is_preserved_condition <- FALSE
  }
  conditionXcontrast.df <- dplyr::mutate(conditionXcontrast.df,
                                         is_preserved_condition = coalesce(is_preserved_condition, FALSE)) %>%
      dplyr::arrange(contrast, contrast_type, metacondition, condition)

  # subset the varspec for those that could be used for contrasts calculation
  contrast_varspecs = list(spec_info = dplyr::filter(varspecs$spec_info, var %in% contrast_vars))
  contrast_varspecs$cats_info = varspecs$cats_info[unique(contrast_varspecs$spec_info$category)]
  # filter for vars that have condition associated with them
  #contrast_varspecs$cats_info = varspecs$cats_info[sapply(varspecs$cats_info, function(df) rlang::has_name(df, condition_agg_col)]
  #contrast_varspecs$spec_info = dplyr::filter(contrast_varspecs$spec_info, category %in% names(contrast_varspecs$cats_info))
  res <- append_contrasts_stats(res, msglm.stan_draws, msglm.stan_stats, contrast_varspecs,
            model_def$metaconditionXcontrast, contrasts.df, conditionXcontrast.df,
            condition_agg_col = condition_agg_col,
            object_cols = object_cols, metacondition_cols = c(),
            group_cols = contrast_group_cols,
            method = contrast_method)
  return (res)
}
