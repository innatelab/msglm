# use experimental design matrices and add
# object-specific object/replicate/batch effects information to model_data
prepare_expanded_effects <- function(model_data, verbose=model_data$model_def$verbose)
{
  model_def <- model_data$model_def

  # FIXME use model_def class information
  is_glmm <- rlang::has_name(model_def, "supXcond.mtxs") &&
             rlang::has_name(model_def, "mixeffects.df")
  if (is_glmm) {
    if (verbose) message("Detected mixeffects data for GLMM model")
    supXcond.mtxs <- model_def$supXcond.mtxs
    mixcoefXeff.mtx <- model_def$mixcoefXeff.mtx
    # "GLMM" mixing model
    if (any(names(supXcond.mtxs) != rownames(mixcoefXeff.mtx))) {
      stop("Mismatch between supXcond.mtxs matrix names and mix coefficient names")
    }
    mixeffects_df <- model_def$mixeffects
    if (any(mixeffects_df$mixeffect != colnames(mixcoefXeff.mtx))) {
      stop("Mismatch between mixcoefXeff.mtx matrix colnames and mix effect names")
    }
    for (i in seq_along(supXcond.mtxs)) {
      if (any(colnames(supXcond.mtxs[[i]]) != model_data$conditions$condition)) {
        stop("Mismatch between supXcond.mtxs[[", names(supXcond.mtxs)[[i]],
             "]] matrix colnames and condition names")
      }
      if (any(rownames(supXcond.mtxs[[i]]) != model_data$supconditions$supcondition)) {
        stop("Mismatch between supXcond.mtxs[[", names(supXcond.mtxs)[[i]],
             "]] matrix rownames and supcondition names")
      }
    }
    model_data$mixeffects <- mixeffects_df
    model_data$mixcoefs <- tibble(mixcoef = rownames(mixcoefXeff.mtx))
    model_data$mixcoefXeff <- mixcoefXeff.mtx
    model_data$supXcond <- supXcond.mtxs
    conditions <- levels(model_data$conditions$condition)
    sactXiact.mtxs <- lapply(model_data$supXcond, function(supXcond.mtx){
      res <- extrude_matrix(supXcond.mtx,
                            model_data$superactions$index_object,
                            model_data$superactions$index_supcondition,
                            blockdim = "index_object", rowixdim = "index_supaction",
                            coldim = "interaction_id", corecoldim = "condition_id")
      names(dimnames(res$mtx)) <- c("sact", "iact")
      res$mtx <- res$mtx[, colSums(abs(res$mtx)) != 0.0, drop=FALSE]
      res$df <- dplyr::rename(res$df, condition=eff, iaction_id=objeff)
      res$iaction_df <- dplyr::rename(res$objeff_df, condition=eff,
                                      iaction_id=objeff,
                                      index_object=obj)
      res$objeff_df <- NULL
      return(res)
    }) %>% pivot_wider(id_cols = )
    # collect all used interaction ids to use as factor levels
    iaction_ids = unique(unlist(lapply(sactXiact.mtxs, function(x) x$iaction_df$iaction_id)))
    # convert condition and iaction_id into factor
    sactXiact.mtxs <- lapply(sactXiact.mtxs, function(x) {
      x$df <- mutate(x$df,
                     condition = factor(condition, levels=conditions),
                     iaction_id = factor(iaction_id, levels=iaction_ids))
      x$iaction_df <- mutate(x$iaction_df,
                             condition = factor(condition, levels=conditions),
                             iaction_id = factor(iaction_id, levels=iaction_ids))
      return(x)
    })
    model_data$interactions <- dplyr::distinct(dplyr::bind_rows(lapply(sactXiact.mtxs, function(x) x$iaction_df))) %>%
      dplyr::mutate(index_condition = as.integer(condition)) %>%
      dplyr::mutate(index_interaction = row_number(),
             is_virtual = FALSE)
    model_data$mixtions <- dplyr::bind_rows(lapply(names(sactXiact.mtxs), function(mixcoef){
      tibble::tibble(mixcoef = mixcoef,
                     iaction_id = factor(colnames(sactXiact.mtxs[[mixcoef]]$mtx), levels=iaction_ids))
    })) %>%
      dplyr::mutate(mixtion_ix = row_number(),
                    mixtion = paste0(mixcoef, " ", iaction_id),
                    mixcoef = factor(mixcoef, levels=rownames(mixcoefXeff.mtx)),
                    mixcoef_ix = replace_na(as.integer(mixcoef), 0L)) %>%
    dplyr::left_join(dplyr::select(model_data$interactions, iaction_id, index_interaction)) %>%
    dplyr::arrange(mixtion_ix)
    model_data$supactXmixt <- bind_cols(lapply(sactXiact.mtxs, function(x) x$mtx))
    dimnames(model_data$supactXmixt) <- list(superaction = model_data$superactions$supaction_id,
                                             mixtion = model_data$mixtions$mixtion)
  } else {
    iaction_ids <- unique(model_data$interactions$iaction_id)
    obs_ids <- unique(model_data$observations$observation_id)
  }

  iactXobjeff_df <- dplyr::full_join(matrix2frame(model_def$conditionXeffect, row_col="condition", col_col = "effect"),
                                     dplyr::select(model_data$objects, index_object, object_id),
                                     by = character()) %>%
    dplyr::inner_join(dplyr::select(model_def$effects, effect, index_effect), by="effect") %>%
    dplyr::inner_join(dplyr::select(model_data$interactions, index_interaction, index_object, condition),
                      by=c("index_object", "condition")) %>%
    dplyr::mutate(object_effect = paste0(effect, '@', object_id))
  model_data$object_effects <- dplyr::select(iactXobjeff_df, index_object, object_id,
                                             index_effect, effect, object_effect) %>%
    dplyr::distinct() %>%
    dplyr::arrange(index_object, index_effect) %>%
    dplyr::mutate(index_object_effect = row_number())
  model_data$iactXobjeff <- frame2matrix(iactXobjeff_df, row_col="index_interaction", col_col="object_effect",
                                         rows=model_data$interactions$index_interaction,
                                         cols=model_data$object_effects$object_effect)

  if (rlang::has_name(model_def, "msexperimentXeffect")) {
    if (verbose) message("Using model_def$msexperimentXeffect for per-MS experiment design matrix")
    msexp_dim <- names(dimnames(model_def$msexperimentXeffect))[[1]]
    obsXobjeff_df <- dplyr::full_join(matrix2frame(model_def$msexperimentXeffect, row_col="condition", col_col = msexp_dim),
                                      dplyr::select(model_data$objects, index_object, object_id),
                                      by = character()) %>%
      dplyr::inner_join(dplyr::select(model_data$observations, index_observation, index_object, index_mschannel, !!sym(msexp_dim)),
                        by=c("index_object", msexp_dim)) %>%
      dplyr::mutate(object_effect = paste0(effect, '@', object_id))
  } else {
    obsXobjeff_df <- dplyr::inner_join(iactXobjeff_df,
                                       dplyr::select(model_data$observations, index_interaction, index_observation),
                                       by="index_interaction")
  }
  model_data$obsXobjeff <- frame2matrix(obsXobjeff_df, row_col="index_observation", col_col="object_effect",
                                        rows=model_data$observations$index_observation,
                                        cols=model_data$object_effects$object_effect)

  if (rlang::has_name(model_def, "msexperimentXbatchEffect")) {
    msexpXbatchEffect <- model_def$msexperimentXbatchEffect
    msexp_idcol <- names(dimnames(msexpXbatchEffect))[[1]]
    batch_effects_df <- model_def$batch_effects
    if (verbose) message(nrow(batch_effects_df), " batch effects on ", msexp_idcol, " level defined")
  } else {
    if (verbose) message("No batch effects defined")
    msexp_idcol <- "mschannel"
    msexpXbatchEffect <- constant_matrix(0, list(mschannel = model_data$mschannels$mschannel,
                                                 batch_effect = character(0)))
    batch_effects_df <- tibble::tibble(batch_effect = character(0),
                                       index_batch_effect = integer(0),
                                       is_positive = logical(0))
  }
  obsXobjbatcheff_df <- dplyr::full_join(matrix2frame(msexpXbatchEffect, row_col=msexp_idcol, col_col = "batch_effect"),
                                         dplyr::select(model_data$objects, index_object, object_id, object_label),
                                         by = character()) %>%
    dplyr::inner_join(dplyr::select(model_data$observations, !!sym(msexp_idcol), index_object, index_observation), by=c("index_object", msexp_idcol)) %>%
    dplyr::inner_join(dplyr::select(batch_effects_df, batch_effect, index_batch_effect), by="batch_effect") %>%
    dplyr::mutate(object_batch_effect = paste0(batch_effect, '@', object_id))
  model_data$object_batch_effects <- dplyr::select(obsXobjbatcheff_df,
                                                   object_batch_effect,
                                                   index_batch_effect, batch_effect,
                                                   index_object, object_id) %>%
    dplyr::distinct() %>%
    dplyr::arrange(index_object, index_batch_effect) %>%
    dplyr::mutate(index_object_batch_effect = row_number())
  model_data$obsXobjbatcheff <- frame2matrix(obsXobjbatcheff_df, row_col="index_observation", col_col="object_batch_effect",
                                             rows=model_data$observations$index_observation,
                                             cols=model_data$object_batch_effects$object_batch_effect)
  if ("subobjects" %in% names(model_data)) {
    if (rlang::has_name(model_def, "msexperimentXquantBatchEffect")) {
      msexpXquantBatchEffect <- model_def$msexperimentXquantBatchEffect
      msexp_idcol <- names(dimnames(msexpXquantBatchEffect))[[1]]
      quant_batch_effects_df <- model_def$quant_batch_effects
      if (verbose) message(nrow(quant_batch_effects_df), " quantification batch effects on ", msexp_idcol, " level defined")
    } else {
      if (verbose) message("No quantification-level batch effects defined")
      msexp_idcol <- "mschannel"
      msexpXquantBatchEffect <- constant_matrix(0, list(mschannel = model_data$mschannels$mschannel,
                                                        quant_batch_effect = character(0)))
      quant_batch_effects_df <- tibble::tibble(quant_batch_effect = character(0),
                                               index_quant_batch_effect = integer(0),
                                               is_positive = logical(0))
    }
    # remove reference subobject (one with the smallest index of EACH object in the model) from quant batch effects
    subobjs_df <- dplyr::group_by(model_data$subobjects, index_object) %>%
      dplyr::filter(index_subobject > min(index_subobject)) %>%
      dplyr::ungroup()
    subobsXsubobjbatcheff_df <- dplyr::full_join(matrix2frame(msexpXquantBatchEffect, row_col=msexp_idcol, col_col="quant_batch_effect"),
                                                 dplyr::select(subobjs_df, index_object, index_subobject, subobject_id),
                                                 by = character()) %>%
      dplyr::inner_join(dplyr::select(model_data$msdata, !!sym(msexp_idcol), index_object, index_subobject, index_subobservation), by=c("index_object", "index_subobject", msexp_idcol)) %>%
      dplyr::inner_join(dplyr::select(quant_batch_effects_df, quant_batch_effect, index_quant_batch_effect), by="quant_batch_effect") %>%
      dplyr::mutate(subobject_batch_effect = paste0(quant_batch_effect, '@', subobject_id)) %>%
      dplyr::arrange(index_quant_batch_effect, index_subobservation)
    model_data$subobject_batch_effects <- dplyr::select(subobsXsubobjbatcheff_df,
                                                        subobject_batch_effect,
                                                        index_quant_batch_effect, quant_batch_effect,
                                                        index_subobject, subobject_id) %>%
      dplyr::distinct() %>%
      dplyr::arrange(index_subobject, index_quant_batch_effect) %>%
      dplyr::mutate(index_subobject_batch_effect = row_number())
    model_data$subobsXsubobjbatcheff <- frame2matrix(subobsXsubobjbatcheff_df,
                                                     row_col="index_subobservation",
                                                     col_col="subobject_batch_effect",
                                                     cols = model_data$subobject_batch_effects$subobject_batch_effect)
  }

  return(model_data)
}

#' @export
# FIXME more checks/control over the columns of intensities_df/stats_df
impute_intensities <- function(intensities_df, stats_df, log2_mean_offset=-1.8, log2_sd_scale=0.3){
  res <- dplyr::inner_join(intensities_df, stats_df) %>%
    dplyr::mutate(intensity_imputed = if_else(is.na(intensity),
                                              2^(rnorm(n(), mean=log2_intensity.mean + log2_intensity.sd * log2_mean_offset,
                                                       sd=log2_intensity.sd * log2_sd_scale)),
                                              intensity)) %>%
    dplyr::ungroup()
  # don't take stats_df columns along
  dplyr::select(res, one_of(colnames(intensities_df)), intensity_imputed)
}

#' @export
cluster_msprofiles <- function(msdata, msrun_stats, obj_col="pepmodstate_id", msrun_col="msrun", nclu=4) {
  # create matrix of intensities
  intensities.df <- tidyr::expand(msdata, !!!rlang::syms(c(obj_col, msrun_col))) %>%
    dplyr::left_join(dplyr::select(msdata, any_of(c(obj_col, msrun_col, "intensity"))),
                     by=c(obj_col, msrun_col)) %>%
    impute_intensities(msrun_stats) %>%
    dplyr::arrange_at(c(obj_col, msrun_col))
  # handle trivial cases
  if (n_distinct(intensities.df[[obj_col]]) == 1L ||
      n_distinct(intensities.df[[msrun_col]]) == 1L) {
    return(tibble(!!obj_col := unique(intensities.df[[obj_col]]),
                  profile_cluster = 1L,
                  nsimilar_profiles = 1L))
  }
  obj_stats.df <- group_by(intensities.df, !!!rlang::syms(obj_col)) %>%
    summarise(n_quants = sum(!is.na(intensity))) %>%
    dplyr::ungroup()
  # add a bit of noise to avoid zero variance
  intensities.mtx <- matrix(log2(pmax(intensities.df$intensity_imputed + rnorm(nrow(intensities.df)), 0)),
                            ncol = n_distinct(intensities.df[[obj_col]]),
                            dimnames = list(msrun = unique(intensities.df[[msrun_col]]),
                                            object = unique(intensities.df[[obj_col]])))
  obj.pca <- stats::prcomp(intensities.mtx, scale.=TRUE)
  # create object feature matrix
  obj.pca_featmtx <- obj.pca$rotation * crossprod(t(rep.int(1, nrow(obj.pca$rotation))),
                                                  summary(obj.pca)$importance[2,])

  res <- tibble(!!obj_col := parse_integer(rownames(obj.pca_featmtx)),
                tmp_profile_cluster = stats::cutree(hclust(dist(obj.pca_featmtx), method="single"),
                                                    min(c(nclu, nrow(obj.pca_featmtx), ncol(obj.pca_featmtx)))))
  # assign profile_cluster indices from largest to smallest clusters
  res_clustats <- dplyr::inner_join(res, obj_stats.df, by=obj_col) %>%
    dplyr::group_by(tmp_profile_cluster) %>%
    dplyr::summarise(nsimilar_profiles = n(),
                     n_quants_cluster_total = sum(n_quants, na.rm=TRUE),
                     n_quants_cluster_median = median(n_quants, na.rm=TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(n_quants_cluster_total), desc(nsimilar_profiles)) %>%
    dplyr::mutate(profile_cluster = row_number())

  return(inner_join(res, res_clustats, by="tmp_profile_cluster") %>%
           dplyr::select(-tmp_profile_cluster))
}

# annotate msdata (full subobject X msrun Cartesian product)
# as reliable (column name) according to
# the specificity of observing the quant object in the MS experiment groups
# defined by `model_def$msdata_specificity_msexp_group_cols`
# or by co-occurrence of distinct quant subobjects in the given MS experiment group
# defined by `model_def$msdata_cooccurrence_msexp_group_cols`
annotate_msdata <- function(msdata_df, model_def, verbose = model_def$verbose,
                            specificity_pvalue = 1E-3,
                            specificity_msexp_group_cols = 'condition',
                            specificity_quantobject_group_cols = NULL,
                            observation_fdr = 0.01,
                            cooccurrence_msexp_group_cols = 'msrun',
                            cooccurrence_pvalue = 1E-3) {
  if (verbose) message("Identifying reliable quantifications")
  msdata_df <- dplyr::mutate(msdata_df, is_observed = !is.na(intensity))

  spec_quantobj_group_cols <- c('object_id', specificity_quantobject_group_cols %||%
                                       intersect(colnames(msdata_df),
                                                 c('subobject_id', 'object_id'))) %>% unique()
  spec_id_cols <- c(specificity_msexp_group_cols, spec_quantobj_group_cols)

  # calculate probabilities that observations of an object are specific to the given interaction
  spec_stats_df <- dplyr::group_by_at(msdata_df, spec_id_cols) %>%
    dplyr::summarise(nms_observed = sum(is_observed), nms_missed = sum(!is_observed),
                     .groups = "drop") %>%
    dplyr::group_by_at(spec_quantobj_group_cols) %>%
    dplyr::mutate(nms_observed_total = sum(nms_observed), nms_missed_total = sum(nms_missed)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(is_reliable = phyper(nms_observed-1, nms_observed_total, nms_missed_total,
                                       nms_observed + nms_missed, lower.tail=FALSE) <= specificity_pvalue)

  msdata_df <- dplyr::left_join(msdata_df,
      dplyr::select_at(spec_stats_df, c(spec_id_cols, 'is_reliable')),
      by = spec_id_cols)
  if (rlang::has_name(msdata_df, 'subobject_id')) {
    cooccur_id_cols <- unique(c('object_id', cooccurrence_msexp_group_cols))

    # calculate probabilities that all quantitations of subobjects in a given observation are false discoveries
    cooccur_stats_df <- dplyr::group_by_at(msdata_df, cooccur_id_cols) %>%
      dplyr::summarise(nsubobj_observed = n_distinct(subobject_id[is_observed]),
                       nsubobj_missed = n_distinct(subobject_id[!is_observed]),
                       .groups = "drop") %>%
      dplyr::mutate(is_cooccurring = pbinom(nsubobj_observed - 1L, nsubobj_observed + nsubobj_missed,
                                            observation_fdr, lower.tail=FALSE) <= cooccurrence_pvalue)

    msdata_df <- dplyr::left_join(msdata_df,
        dplyr::select_at(cooccur_stats_df, c(cooccur_id_cols, 'is_cooccurring')),
        by = cooccur_id_cols) %>%
        dplyr::mutate(is_reliable = is_reliable | is_cooccurring,
                      is_cooccurring = NULL)
  }
  msdata_df <- dplyr::mutate(msdata_df, is_reliable = !is.na(intensity) & is_reliable)

  return(msdata_df)
}

# select relevant msdata and subobjects from msdata_df
# and define model_data$subobjects and model_data$observations
prepare_msdata <- function(model_data, msdata, verbose = model_data$model_def$verbose, ...) {
  model_def <- model_data$model_def
  intensities_dfname <- paste0(model_def$quantobject, "_",
                               if (any(coalesce(model_data$mschannels$mstag, "F") == "F")) "" else "tag",
                               "intensities")
  if (!rlang::has_name(msdata, intensities_dfname)) {
    stop("No intensities (", intensities_dfname, " data frame) found in msdata")
  }
  intensities_df <- msdata[[intensities_dfname]]
  msexp_idcol <- intersect(c("mschannel", "msrun"), colnames(intensities_df))[[1]]
  modelobj_idcol <- paste0(model_def$modelobject, "_id")
  quantobj_idcol <- paste0(model_def$quantobject, "_id")
  intensities_df <- dplyr::select(intensities_df, !!sym(quantobj_idcol), !!sym(msexp_idcol), intensity)


  if (model_def$modelobject == model_def$quantobject) {
    # modelobj is quanted directly
    msdata_df <- dplyr::left_join(model_data$observations,
                                  intensities_df, by = c(quantobj_idcol, msexp_idcol)) %>%
        annotate_msdata(model_def) %>%
        dplyr::arrange(index_observation)
  } else if (model_def$quantobject == "pepmodstate") {
    # quant specific pepmodstates of modelobj
    modelobj2pepmod_df <- msdata[[paste0(model_def$modelobj, "2pepmod")]]
    subobjs_df <- dplyr::inner_join(modelobj2pepmod_df,
                                    dplyr::select(model_data$objects, !!sym(modelobj_idcol), object_id, index_object),
                                    by=modelobj_idcol) %>%
      dplyr::filter(is_specific) %>%
      dplyr::inner_join(dplyr::select(msdata$pepmodstates, pepmod_id, charge, pepmodstate_id), by="pepmod_id") %>%
      dplyr::mutate(subobject_id = pepmodstate_id)
    if (nrow(subobjs_df) == 0L) stop("No specific ", model_def$quantobject, "s found for ", modelobj_idcol, "=", model_data$modelobj_id)
    if (verbose) message(nrow(subobjs_df), " specific ", model_def$quantobject, "(s) found")
    msdata_df <- dplyr::inner_join(dplyr::select(model_data$observations, index_mschannel, index_observation, index_object, object_id),
                                   dplyr::select(subobjs_df, index_object, subobject_id), by="index_object") %>%
        dplyr::inner_join(model_data$mschannels, by = c("index_mschannel")) %>%
        dplyr::left_join(intensities_df, by=c(subobject_id="pepmodstate_id", msexp_idcol))
    if (all(is.na(msdata_df$intensity))) stop("No quantifications for ", nrow(subobjs_df), " specific ",
                                              model_def$quantobject, "(s) of ", modelobj_idcol, "=",
                                              model_data$object_id)
    msdata_df <- annotate_msdata(msdata_df, model_def, verbose=verbose, ...)

    # arrange pepmodstates by object, by profile cluster and by the number of quantitations
    max_subobjects <- model_def$max_subobjects %||% 20L
    subobject_group_size <- model_def$subobject_group_size %||% (max_subobjects %/% 2)
    subobj_stats_df <- msdata_df %>%
      dplyr::group_by(index_object, subobject_id) %>%
      dplyr::summarise(n_quants = sum(!is.na(intensity)),
                       intensity_med = median(intensity, na.rm=TRUE),
                       .groups = "drop") %>%
      dplyr::inner_join(
        dplyr::group_by(msdata_df, index_object) %>%
        dplyr::group_modify(~ cluster_msprofiles(.x, msdata[[paste0(msexp_idcol, "_pepmodstate_stats")]],
                                                 obj_col='subobject_id', msrun_col=msexp_idcol)) %>%
          dplyr::ungroup(),
        by = c("subobject_id", "index_object")) %>%
      dplyr::left_join(dplyr::select(subobjs_df, subobject_id, pepmod_id, is_specific, charge), by = 'subobject_id') %>%
      dplyr::arrange(index_object, profile_cluster, desc(is_specific), desc(n_quants), desc(intensity_med),
                     pepmod_id, charge) %>%
      dplyr::group_by(index_object, profile_cluster) %>%
      dplyr::mutate(index_subobject_group = row_number() %/% subobject_group_size, # put objects within cluster into groups of 20
                    index_subobject_local = row_number() %% subobject_group_size) %>%
      dplyr::ungroup()
    # take the first group of 10 objects from each cluster, then continue with the second group etc
    model_data$subobjects <- subobj_stats_df %>%
      dplyr::arrange(index_object, index_subobject_group, profile_cluster, index_subobject_local) %>%
      dplyr::mutate(index_subobject = row_number()) %>%
      dplyr::filter(index_subobject <= max_subobjects) # remove less abundant subobjects of rich objects
    if (verbose) message(nrow(model_data$subobjects), " ", model_def$quantobject, "(s) from ",
                         n_distinct(model_data$subobjects$profile_cluster), " cluster(s) selected")
    msdata_df <- dplyr::inner_join(msdata_df,
                                   dplyr::select(model_data$subobjects, index_object, index_subobject, subobject_id),
                                   by=c("index_object", "subobject_id")) %>%
      dplyr::arrange(index_object, index_observation, index_subobject) %>%
      dplyr::mutate(index_subobservation = row_number())
  } else {
     stop("Unsupported combination of modelobject=", model_def$modelobject,
         " and quantobject=", model_def$quantobject)
  }
  # separately index quantifications and missing data
  model_data$msdata <- mutate(msdata_df,
                              index_qdata = if_else(is_observed, cumsum(is_observed), NA_integer_),
                              index_mdata = if_else(!is_observed, cumsum(!is_observed), NA_integer_))
  message(nrow(model_data$msdata), " ", if_else(rlang::has_name(model_data$msdata, 'index_subobservation'),
                                                'subobservation', 'observation'), '(s) of ',
          n_distinct(model_data$msdata$index_object), ' ', model_def$modelobject, '(s)',
          if_else(rlang::has_name(model_data$msdata, 'index_subobject'),
                  paste0(' with ', n_distinct(model_data$msdata$index_subobject), ' ',
                         model_def$quantobject, '(s)'), ''), ': ',
          sum(model_data$msdata$is_observed), " quantitation(s) (",
          sum(model_data$msdata$is_reliable), " reliable), ",
          sum(!model_data$msdata$is_observed), " missed")
  return(model_data)
}

#' Title
#'
#' @param model_def
#' @param msdata
#' @param object_ids
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
msglm_data <- function(model_def, msdata, object_ids, verbose = model_def$verbose,
                       msexperiment_extra_cols = character(0),
                       msexperiment_shift_col = "total_shift", ...) {
  model_data <- list(model_def = model_def, object_id = object_ids)
  modelobj <- model_def$modelobject
  modelobj_idcol <- str_c(modelobj, "_id")

  if (!rlang::has_name(msdata, paste0(modelobj, "s"))) {
    stop("No model object (", modelobj, ") information found in MS data")
  }
  modelobjs_df <- dplyr::mutate(msdata[[paste0(modelobj, "s")]],
                                object_id = !!sym(modelobj_idcol),
                                object_label = !!sym(str_c(modelobj, "_label")))
  modelobj_cols <- model_def$object_cols %||%
      (c(intersect(c("majority_protein_acs", "protein_acs",
                     "gene_names", "protein_names"), colnames(modelobjs_df)),
         str_subset(colnames(modelobjs_df), "^is_")) %>% unique())
  if (verbose) {
    message("Model object (", modelobj, ") columns to use: ", paste0(modelobj_cols, collapse=", "))
  }
  model_data$objects <- dplyr::filter(modelobjs_df, object_id %in% object_ids) %>%
      dplyr::select_at(c("object_id", "object_label", modelobj_idcol,
                         modelobj_cols)) %>%
      dplyr::arrange(object_id) %>%
      dplyr::mutate(index_object = row_number())
  missing_obj_ids <- setdiff(object_ids, unique(modelobjs_df$object_id))
  if (length(missing_obj_ids) > 0) {
    stop("Objects not found: ", paste0(missing_obj_ids, ", "))
  }

  # FIXME support experiments with mschannels
  if (rlang::has_name(msdata, "mschannels")) {
    if (verbose) message("Using msdata$mschannels for experimental design")
    checkmate::assert_tibble(msdata$mschannels)
    checkmate::assert_names(names(msdata$mschannels), must.include = c("mschannel", "msrun", "mstag", "condition"))
    mschannels_df <- msdata$mschannels
  } else if (rlang::has_name(msdata, "msruns")) {
    if (verbose) message("Using msdata$msruns for experimental design")
    checkmate::assert_tibble(msdata$msruns)
    checkmate::assert_names(names(msdata$msruns), must.include = c("msrun", "condition"))
    mschannels_df <- dplyr::mutate(msdata$msruns, mstag = NA_character_, mschannel = msrun)
  } else {
    stop("Cannot find MS experiments information (msdata$mschannels or msdata$msruns)")
  }
  checkmate::assert_subset(mschannels_df$condition, model_def$conditions$condition)
  # add normalization shifts to mschannels
  if (rlang::has_name(msdata, "mschannel_shifts")) {
    checkmate::assert_tibble(msdata$mschannel_shifts)
    checkmate::assert_names(names(msdata$mschannel_shifts), must.include = c("mschannel", msexperiment_shift_col))
    if (verbose) message("Using msdata$mschannel_shifts$", msexperiment_shift_col, " for MS experiments normalization")
    checkmate::assert_set_equal(msdata$mschannel_shifts$mschannel, mschannels_df$mschannel)
    mschannels_df <- dplyr::inner_join(mschannels_df, dplyr::select_at(msdata$mschannel_shifts,
                                                                       c("msrun", mschannel_shift=msexperiment_shift_col)), by="mschannel")
    checkmate::assert_numeric(msdata$mschannel_shifts[[msexperiment_shift_col]], any.missing = FALSE, finite=TRUE, lower=-20, upper=20)
  } else if (rlang::has_name(msdata, "msrun_shifts")) {
    checkmate::assert_tibble(msdata$msrun_shifts)
    checkmate::assert_names(names(msdata$msrun_shifts), must.include = c("msrun", msexperiment_shift_col))
    if (verbose) message("Using msdata$msrun_shifts$", msexperiment_shift_col, " for MS experiments normalization")
    checkmate::assert_set_equal(msdata$msrun_shifts$msrun, mschannels_df$msrun)
    checkmate::assert_numeric(msdata$msrun_shifts[[msexperiment_shift_col]], any.missing = FALSE, finite=TRUE, lower=-20, upper=20)
    mschannels_df <- dplyr::inner_join(mschannels_df,
                                       dplyr::select_at(msdata$msrun_shifts,
                                                        c("msrun", mschannel_shift=msexperiment_shift_col)),
                                       by="msrun")
  } else {
    stop("Cannot find MS experiments normalization information (msdata$mschannel_shifts or msdata$msrun_shifts)")
  }
  checkmate::assert_set_equal(unique(mschannels_df$condition),
                              dplyr::filter(model_def$conditions, !is_virtual)$condition)
  model_data$mschannels <- mschannels_df %>%
    dplyr::inner_join(dplyr::select(model_def$conditions, condition, index_condition, any_of("mstag")),
                      by=intersect(c('condition', 'mstag'), colnames(model_def$conditions))) %>%
    dplyr::arrange(index_condition, msrun, mstag) %>%
    dplyr::mutate(index_mschannel = row_number(),
                  index_msrun = as.integer(factor(msrun, levels=unique(msrun))),
                  index_mscalib = 1L) # FIXME there could be multiple protocols
  if (rlang::has_name(msdata, "msexperimentXeffect")) {
    msexp_dim <- names(dimnames(model_def$msexperimentXeffect))
    checkmate::assert_character(mschannels_df[[msexp_dim]], unique=TRUE)
    checkmate::assert_set_equal(mschannels_df[[msexp_dim]], rownames(model_def$msexperimentXeffect))
  }

  # all objects X all conditions (including virtual)
  # FIXME rename objectXcondition
  model_data$interactions <- dplyr::full_join(dplyr::select(model_data$objects, object_id, index_object),
                                              dplyr::select(model_def$conditions, index_condition, is_virtual, condition),
                                              by=character()) %>%
    dplyr::mutate(iaction_id = paste(condition, object_id)) %>%
    dplyr::arrange(index_condition, index_object) %>%
    dplyr::mutate(index_interaction = row_number())

  # all objects X all MS channels (only those with actual experiments)
  # FIXME rename objectXmschannel
  model_data$observations <- dplyr::inner_join(model_data$interactions, model_data$mschannels,
                                               by = c("index_condition", "condition")) %>%
    dplyr::arrange(index_object, index_interaction, index_mschannel) %>%
    dplyr::mutate(index_observation = row_number(),
                  observation_id = str_c(object_id, '_', mschannel))

  model_data$quantobj_mscalib <- msdata[[paste0(quantobj, "_mscalib")]]
  model_data$quantobj_labu_shift <- msdata[[paste0(quantobj, '_labu_shift')]]
  model_data$quantobj_labu_min <- msdata[[paste0(quantobj, '_labu_min')]]

  model_data <- prepare_msdata(model_data, msdata, verbose=verbose, ...)
  model_data <- prepare_expanded_effects(model_data, verbose=verbose)

  return (model_data)
}
