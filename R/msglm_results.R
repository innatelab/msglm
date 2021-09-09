#' Generate the dimensions information for the MSGLM report.
#'
#' @param model_data MSGLM model data that is used to get the dimensions info.
#' @param object_cols Non-standard object information columns to
#'                    include in the reports.
#' @returns The named list of data frames with the dimension information.
#'          The element names are the names of the dimensions,
#'          Each data frame has `index_<dimname>` key,
#'          which matches the index along the corresponding dimension in
#'          the MSGLM model.
#'
#' @seealso process.stan_fit
#'
#' @export
msglm_dims <- function(model_data)
{
  is_glmm <- "mixeffects" %in% names(model_data)
  xaction_ix_col <- if (is_glmm) "glm_supaction_ix" else "index_interaction"
  xdition_ix_col <- if (is_glmm) "supcondition_ix" else "index_condition"
  xdition_col <-  if (is_glmm) "supcondition" else "condition"

  objs_df <- model_data$objects
  if (!is.null(object_cols)) {
    objs_df <- dplyr::select(objs_df, !!!unique(c("index_object", object_cols)))
  }
  res <- list(iteration = NULL,
    msrun = dplyr::select(model_data$mschannels, index_msrun, msrun, any_of(c("mschannel", "index_ mschannel", "mstag", "condition", "supcondition"))),
    observation = dplyr::select(model_data$msdata, index_observation, !!xaction_ix_col, index_object,
                                !!xdition_ix_col, !!xdition_col, msrun, index_msrun, any_of(c("mschannel", "index_ mschannel", "mstag"))) %>%
        dplyr::distinct() %>%
        dplyr::inner_join(objs_df),
    object = model_data$objects, # use full object information
    object_effect = model_data$object_effects %>%
        dplyr::mutate(index_object = as.integer(index_object)) %>%
        dplyr::inner_join(objs_df) %>%
        dplyr::inner_join(model_data$effects) %>%
        maybe_rename(c("prior_mean" = "mean", "prior_tau" = "tau")),
    object_batch_effect = model_data$object_batch_effects %>%
        dplyr::mutate(index_object = as.integer(index_object)) %>%
        dplyr::inner_join(objs_df)
  )
  if (!rlang::has_name(res$object_effect, "prior_mean")) { # set the default min to 0
    message("object_effects$prior_mean missing, setting to 0")
    res$object_effect$prior_mean <- 0.0
  }
  res$object_effect <- mutate(res$object_effect,
                              prior_mean_log2=prior_mean/log(2))
  if ("subobjects" %in% names(model_data)) {
    res$subobject <- dplyr::select(model_data$subobjects, index_object, index_subobject,
                                   any_of(c("protregroup_id", "protgroup_id", "pepmod_id", "pepmodstate_id", "charge")))
    if ("suo_subbatch_effects" %in% names(model_data)) {
      res$subobject_subbatch_effect <- model_data$suo_subbatch_effects %>%
        dplyr::mutate(index_subobject = as.integer(index_subobject)) %>%
        dplyr::left_join(res$subobject)
    }
  }
  if ("index_mscalib" %in% colnames(model_data$mschannels)) {
    res$mscalib <- dplyr::select(model_data$mschannels, index_mscalib,
                                 any_of("instrument")) %>%
      dplyr::distinct()
  }
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
    res$object_mixeffect <- mutate(res$object_mixeffect,
                                   prior_mean_log2=prior_mean/log(2))
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

#' Compare the posterior of each variable in varspecs with its `prior_mean`
#' and return the p-value for the significance of the difference (w.r.t `tail`)
#'
#' @param vars_draws MCMC draws from the posterior in `posterior::draws_array` format
#' @param varspecs variables to calculate P-values for
#' @param tail which comparison to do, one of `both` (the default),
#'             `negative` (``P(x \leq t)``) or `positive` (``P(x \geq t)``),
#'             where ``t`` is the `prior_mean` of the variable.
#'             For `both` the double of the minimal of the two p-values is given.
#' @param nsteps (defaults to 100) how many bins to use for the calculation of p-values
#' @param maxBandwidth constrain the rule-of-thumb bandwidth for the posterior distribution
#'       if it is above the specified limit
#' @returns data frame with the p-value per each varspec
#' @export
vars_pvalues <- function(vars_draws, varspecs, tail = c("both", "negative", "positive"),
                         nsteps = 100L, maxBandwidth = NA_real_,
                         mlog10pvalue_threshold = 10.0,
                         mlog10pvalue_hard_threshold_factor = 3.0
){
  tail = match.arg(tail)
  if (!rlang::has_name(varspecs, "prior_mean")) {
    varspecs$prior_mean <- 0.0
  }
  res.df <- dplyr::mutate(varspecs,
                          has_mcmc_draws = varspec %in% dimnames(vars_draws)$variable) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(p_value = if (has_mcmc_draws) {
                        pvalue_compare(rlang::as_double(posterior::subset_draws(vars_draws, variable = varspec)),
                                       prior_mean, tail = tail, nsteps = nsteps, bandwidth = maxBandwidth,
                                       mlog10_threshold = mlog10pvalue_threshold,
                                       mlog10_hard_threshold_factor = mlog10pvalue_hard_threshold_factor)
                      } else { NA_real_ }) %>%
    dplyr::ungroup()
  if (!all(res.df$has_mcmc_draws)) {
    warning("Missing MCMC samples for vars: ", paste0(res.df$varspec[!res.df$has_mcmc_draws], collapse=" "))
  }
  res.df <- dplyr::select(res.df, -has_mcmc_draws)
  return(res.df)
}

#' Calculate contrasts (linear combination of model effects) using the posterior
#' MCMC draws.
#'
#' The procedure allows advanced calculation of contrasts based on contrast-specific
#' grouping of the variables (provided by `vargroups`). In that case, the
#' contrasts for all possible combinations of variables in the groups are being
#' calculated and combined together to form the posterior distribution of the
#' contrast.
#' For example, consider the group ``A`` contains variables ``a_1`` and ``a_2``,
#' and the group ``B`` is ``{b_1, b_2}``. For the ``A - B`` contrast,
#' the procedure will calculate 4 linear combinations: ``a_1 - b_1``, ``a_1 - b_2``,
#' ``a_2 - b_1`` and ``a_2 - b_2``, then it will combine the resulting 4 sets of
#' MCMC draws and proceed with calculating the posterior summary.
#'
#' @param vars_draws MCMC draws for the model variables to use for contrast
#'     calculation in `posterior::draws_array` 3D array format.
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
#'       obligatory contrast identifier, the obligatory `offset` specifies
#'       what the contrast should be compared with for significance testing.
#' @param nsteps (defaults to 100) how many bins to use for the calculation of p-values
#' @param maxBandwidth constrain the rule-of-thumb bandwidth for the posterior distribution
#'       if it is above the specified limit
#'
#' @export
vars_contrast_stats <- function(vars_draws, vargroups,
                                vargroupXcontrast, contrasts,
                                nsteps = 100L, maxBandwidth = NA_real_,
                                mlog10pvalue_threshold = 10.0,
                                mlog10pvalue_hard_threshold_factor = 3.0)
{
  if (n_groups(vargroups) < nrow(vargroupXcontrast)) {
    stop("Number of vargroups (", n_groups(vargroups),
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
  contrast_offsets <- rlang::set_names(contrasts$offset, contrasts$contrast)
  res <- dplyr::left_join(dplyr::select(contrasts, -dplyr::matches("^[lr]hs_quantile")),
        ContrastStatistics(vars_draws,
                            vargroups$index_varspec, vargroups$`__vargroup_ix__`, vargroups$`__contrast_ix__`,
                            vargroupXcontrast, contrast_offsets[colnames(vargroupXcontrast)],
                            nsteps = nsteps, maxBandwidth = maxBandwidth,
                            mlog10pvalue_threshold = mlog10pvalue_threshold,
                            mlog10pvalue_hard_threshold_factor = mlog10pvalue_hard_threshold_factor,
                            summaryfun = function(draws) {
                              posterior::as_draws_array(draws) %>%
                              posterior::summarise_draws(
                                  c("mean", "median", "sd", "mad", "quantiles_ci",
                                    "rhat", "ess_bulk", "ess_tail"))
                            }), by = "__contrast_ix__") %>%
    dplyr::select(-variable, -`__contrast_ix__`) %>%
    dplyr::mutate(mean_log2 = mean/log(2),
                  median_log2 = median/log(2),
                  sd_log2 = sd/log(2),
                  p_value = 2*pmin(prob_nonpos, prob_nonneg, 0.5))
  return (res)
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
append_contrasts_stats <- function(vars_results, standraws, varspecs,
        metaconditionXcontrast, contrasts.df, conditionXcontrast.df,
        condition_agg_col = "condition", # filtering is based on pregrouping quantiles in metacondition using this column (per contrast)
        object_cols = 'index_object', metacondition_cols = c(),
        group_cols = c()
){
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
    cat_vars_stats <- dplyr::semi_join(vars_results[[cat$category]]$stats, dplyr::select(cat_varspecs, varspec), by = "varspec")
    cat_cols <- intersect(colnames(cat_vars_stats), colnames(cat_info))

    conditionXcontrast_pregroup_stats.df <- dplyr::inner_join(cat_vars_stats,
                                                              dplyr::mutate(conditionXcontrast.df, is_lhs = weight > 0)) %>%
      dplyr::group_by_at(c("var", group_cols, contrast_col, metacondition_col, condition_agg_col,
                           "is_lhs", "is_preserved_condition")) %>%
      dplyr::summarize(var_pregroup_max = max(mean), var_pregroup_min = min(mean)) %>%
      dplyr::group_by_at(c("var", group_cols, contrast_col, metacondition_col, "is_lhs")) %>%
      dplyr::mutate(var_pregroup_max_qtile = cume_dist(var_pregroup_max) - 1/n(),
                    var_pregroup_min_qtile = cume_dist(var_pregroup_min) - 1/n()) %>%
      dplyr::ungroup() %>% dplyr::inner_join(cutoff_quantiles.df) %>%
      dplyr::mutate(is_accepted = (var_pregroup_min_qtile >= qtile_min) &
                                  ((var_pregroup_max_qtile <= qtile_max) |
                                  # FIXME use contrast_type for is_lhs check?
                                   (!is_lhs & is_preserved_condition)))

    vargroups.df <- dplyr::inner_join(cat_varspecs, cat_info, by = "var_index") %>%
      dplyr::inner_join(dplyr::filter(conditionXcontrast_pregroup_stats.df, is_accepted)) %>%
      dplyr::select_at(c("var", "index_varspec", "category", group_cols, contrast_col, metacondition_col, condition_col)) %>%
      dplyr::group_by_at(c("var", "category", group_cols)) %>%
      dplyr::group_modify(~ vars_contrast_stats(standraws,
                              vargroups = dplyr::group_by_at(.x, c(metacondition_col, contrast_col)),
                              vargroupXcontrast = metaconditionXcontrast,
                              contrasts = contrasts.df))
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

# quantiles for symmetric 50% and 95% credible intervals
#' @export
quantiles_ci <- function(x) { posterior::quantile2(x, probs=c(0.025, 0.25, 0.75, 0.975)) }

#' @export
process.stan_fit <- function(msglm.stan_fit, dims_info,
                             vars_info = attr(msglm.stan_fit, "msglm_vars_info"),
                             effect_vars = unlist(lapply(vars_info, function(vi) str_subset(vi$names, "_(?:mix)?effect(?:_replCI)?$"))),
                             contrast_vars = default_contrast_vars(vars_info), contrast_group_cols = c(),
                             condition_agg_col = "condition", object_cols = setdiff(colnames(dims_info$object), "index_object"),
                             min.iteration=NA, chains=NA, verbose=FALSE)
{
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
  msglm.stan_stats <- msglm.stan_fit$summary(variables = avail_vars,
                                             c("mean", "median", "sd", "mad", "quantiles_ci",
                                               "rhat", "ess_bulk", "ess_tail")) %>%
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
      dplyr::select(-category) %>%
      dplyr::mutate(mean_log2 = mean/log(2),
                    median_log2 = median/log(2),
                    sd_log2 = sd/log(2))
    cat_eff_varspecs.df <- dplyr::filter(cat_varspecs.df, var %in% effect_vars) %>%
      dplyr::left_join(dplyr::select(cat_info.df, var_index, any_of("prior_mean")), by = "var_index")
    if (nrow(cat_eff_varspecs.df) > 0) {
      message('    - calculating P-values for: ',
              paste0(unique(cat_eff_varspecs.df$var), collapse=', '), '...')
      p_values.df <- vars_pvalues(msglm.stan_draws, dplyr::select(cat_eff_varspecs.df, varspec, any_of("prior_mean")))
      res.df <- dplyr::left_join(res.df, p_values.df) %>%
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
    iactions.df <- dplyr::select(varspecs$cats_info$observations, index_interaction) %>% dplyr::distinct()
    iactions_diag <- diag(nrow = nrow(iactions.df), ncol = nrow(iactions.df))
    dimnames(iactions_diag) <- list(index_interaction = iactions.df$index_interaction,
                                    contrast = iactions.df$index_interaction)
    res$iactions_obsCI <- list(stats = vars_contrast_stats(msglm.stan_draws,
                                                           dplyr::filter(varspecs$spec_info, var == 'obs_labu') %>%
                                                           dplyr::inner_join(varspecs$cats_info$observations, by="var_index") %>%
                                                           dplyr::group_by(condition, index_object, index_interaction),
                                                           vargroupXcontrast = iactions_diag,
                                                           contrasts = dplyr::mutate(iactions.df, contrast=index_interaction, offset=0)) %>%
                                      # actually, we don't need contrasts
                                      dplyr::select(-contrast, -offset))
  }

  message("Calculating contrasts...")
  contrasts.df <- rlang::env_get(nm="contrasts.df", default=dplyr::distinct(dplyr::transmute(contrastXmetacondition.df, contrast, contrast_type, offset=0)), inherit=TRUE)
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

  metaconditionXcontrast.df <- as.data.frame.table(metaconditionXcontrast.mtx, responseName="weight") %>%
    dplyr::filter(weight != 0) %>%
    dplyr::inner_join(contrasts.df)
  metacondition_col <- names(dimnames(metaconditionXcontrast.mtx))[[1]]
  contrast_col <- names(dimnames(metaconditionXcontrast.mtx))[[2]]
  # get or create conditionXcontrast
  conditionXcontrast.df <- rlang::env_get(nm="conditionXcontrast.df",
      default = metaconditionXcontrast.df %>%
            dplyr::inner_join(conditionXmetacondition.df) %>%
            dplyr::arrange_at(c(contrast_col, "contrast_type", metacondition_col, 'condition')))
  if (!rlang::has_name(conditionXcontrast.df, "is_preserved_condition")) {
    conditionXcontrast.df$is_preserved_condition <- FALSE
  }

  # subset the varspec for those that could be used for contrasts calculation
  contrast_varspecs = list(spec_info = dplyr::filter(varspecs$spec_info, var %in% contrast_vars))
  contrast_varspecs$cats_info = varspecs$cats_info[unique(contrast_varspecs$spec_info$category)]
  # filter for vars that have condition associated with them
  #contrast_varspecs$cats_info = varspecs$cats_info[sapply(varspecs$cats_info, function(df) rlang::has_name(df, condition_agg_col)]
  #contrast_varspecs$spec_info = dplyr::filter(contrast_varspecs$spec_info, category %in% names(contrast_varspecs$cats_info))
  res <- append_contrasts_stats(res, msglm.stan_draws, contrast_varspecs,
            metaconditionXcontrast.mtx, contrasts.df, conditionXcontrast.df,
            condition_agg_col = condition_agg_col,
            object_cols = object_cols, metacondition_cols = c(),
            group_cols = contrast_group_cols)
  return (res)
}
