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

  msprobeXeffect_name <- "msprobeXeffect"
  if (rlang::has_name(model_def, "msprobeXeffect")) {
    if (verbose) message("Using model_def$", msprobeXeffect_name, " for per-", msprobe, " design matrix")
    msprb_dim <- names(dimnames(model_def[[msprobeXeffect_name]]))[[1]]
    if (msprb_dim != msprb_idcol) {
      stop("Name of model_def$", msprobeXeffect_name, " rows dimension (", msprb_dim, ") inconsistent with ", msprb_idcol)
    }
    obsXobjeff_df <- dplyr::full_join(matrix2frame(model_def[[msprobeXeffect_name]], row_col="condition", col_col = msprb_idcol),
                                      dplyr::select(model_data$objects, index_object, object_id),
                                      by = character()) %>%
      dplyr::inner_join(dplyr::select(model_data$observations, index_observation, index_object, index_mschannel, !!sym(msprb_idcol)),
                        by=c("index_object", msprb_idcol)) %>%
      dplyr::mutate(object_effect = paste0(effect, '@', object_id))
  } else {
    obsXobjeff_df <- dplyr::inner_join(iactXobjeff_df,
                                       dplyr::select(model_data$observations, index_interaction, index_observation),
                                       by="index_interaction")
  }
  model_data$obsXobjeff <- frame2matrix(obsXobjeff_df, row_col="index_observation", col_col="object_effect",
                                        rows=model_data$observations$index_observation,
                                        cols=model_data$object_effects$object_effect)

  if (rlang::has_name(model_def, "msprobeXbatchEffect")) {
    msprobeXbatchEffect <- model_def$msprobeXbatchEffect
    batch_effects_df <- model_def$batch_effects
    if (verbose) message(nrow(batch_effects_df), " batch effects on ",
                         msprb, " level defined")
  } else {
    if (verbose) message("No batch effects defined")
    msprobeXbatchEffect <- constant_matrix(0, list(msprobe = character(0),
                                                   batch_effect = character(0)))
    batch_effects_df <- tibble::tibble(batch_effect = character(0),
                                       index_batch_effect = integer(0),
                                       is_positive = logical(0))
  }
  obsXobjbatcheff_df <- dplyr::full_join(matrix2frame(msprobeXbatchEffect, row_col="msprobe", col_col="batch_effect"),
                                         dplyr::select(model_data$objects, index_object, object_id, object_label),
                                         by = character()) %>%
    dplyr::inner_join(dplyr::select(batch_effects_df, batch_effect, index_batch_effect), by="batch_effect") %>%
    dplyr::left_join(dplyr::select(model_data$observations, msprobe, index_object, index_observation),
                     by=c("index_object", "msprobe")) %>%
    dplyr::mutate(object_batch_effect = paste0(batch_effect, '@', object_id))
  model_data$object_batch_effects <- dplyr::select(obsXobjbatcheff_df,
                                                   object_batch_effect,
                                                   index_batch_effect, batch_effect,
                                                   index_object, object_id) %>%
    dplyr::distinct() %>%
    dplyr::arrange(index_object, index_batch_effect) %>%
    dplyr::mutate(index_object_batch_effect = row_number())
  model_data$obsXobjbatcheff <- frame2matrix(obsXobjbatcheff_df, row_col="index_observation", col_col="object_batch_effect",
                                             rows=if (nrow(obsXobjbatcheff_df) > 0) {
                                                    model_data$observations$index_observation
                                                  } else integer(0),
                                             cols=model_data$object_batch_effects$object_batch_effect)
  if (rlang::has_name(model_data, "subobjects")) {
    if (rlang::has_name(model_def, "mschannelXquantBatchEffect")) {
      mschannelXquantBatchEffect <- model_def$mschannelXquantBatchEffect
      mschan <- model_data$msentities[['mschannel']]
      mschan_idcol <- mschan
      quant_batch_effects_df <- model_def$quant_batch_effects
      if (verbose) message(nrow(quant_batch_effects_df),
                           " quantification batch effects on ", mschan_idcol, " level defined")
    } else {
      if (verbose) message("No quantification-level batch effects defined")
      mschannelXquantBatchEffect <- constant_matrix(0, list(mschannel = character(0),
                                                            quant_batch_effect = character(0)))
      quant_batch_effects_df <- tibble::tibble(quant_batch_effect = character(0),
                                               index_quant_batch_effect = integer(0),
                                               is_positive = logical(0))
    }
    # remove reference subobject (one with the smallest index of EACH object in the model) from quant batch effects
    subobjs_df <- dplyr::group_by(model_data$subobjects, index_object) %>%
      dplyr::filter(index_subobject > min(index_subobject)) %>%
      dplyr::ungroup()
    subobsXsubobjbatcheff_df <- dplyr::full_join(matrix2frame(mschannelXquantBatchEffect,
                                                              row_col="mschannel", col_col="quant_batch_effect"),
                                                 dplyr::select(subobjs_df, index_object, index_subobject, subobject_id),
                                                 by = character()) %>%
      dplyr::inner_join(dplyr::select(quant_batch_effects_df, quant_batch_effect, index_quant_batch_effect), by="quant_batch_effect") %>%
      dplyr::left_join(dplyr::select(model_data$msdata, mschannel, index_object, index_subobject, index_subobservation),
                       by=c("index_object", "index_subobject", "mschannel")) %>%
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
                                                     rows = if (nrow(subobsXsubobjbatcheff_df)>0) {
                                                              model_data$msdata$index_subobservation
                                                            } else integer(0),
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
cluster_msprofiles <- function(msdata, mschannel_stats,
                               obj_col=msdata$msentities[['quantobject']],
                               mschannel_col=msdata$msentities[['mschannel']], nclu=4) {
  # create matrix of intensities
  objs.df <- dplyr::select_at(msdata, obj_col) %>%
    dplyr::distinct() %>% dplyr::arrange(!!sym(obj_col)) %>%
    dplyr::mutate(`__index_msobject__` = row_number())
  intensities.df <- tidyr::expand(msdata, !!!rlang::syms(c(obj_col, mschannel_col))) %>%
    dplyr::left_join(dplyr::select(msdata, any_of(c(obj_col, mschannel_col, "intensity"))),
                     by=c(obj_col, mschannel_col)) %>%
    impute_intensities(mschannel_stats) %>%
    dplyr::inner_join(objs.df, by=obj_col) %>%
    dplyr::arrange(`__index_msobject__`, !!sym(mschannel_col))
  # handle trivial cases
  if (n_distinct(intensities.df[[obj_col]]) == 1L ||
      n_distinct(intensities.df[[mschannel_col]]) == 1L) {
    return(tibble(!!obj_col := unique(intensities.df[[obj_col]]),
                  profile_cluster = 1L,
                  nsimilar_profiles = 1L))
  }
  obj_stats.df <- group_by(intensities.df, `__index_msobject__`) %>%
    summarise(n_quants = sum(!is.na(intensity))) %>%
    dplyr::ungroup()
  # add a bit of noise to avoid zero variance
  intensities.mtx <- matrix(log2(pmax(intensities.df$intensity_imputed + rnorm(nrow(intensities.df)), 0)),
                            ncol = nrow(objs.df),
                            dimnames = list(mschannel = unique(intensities.df[[mschannel_col]]),
                                            `__index_msobject__` = NULL))
  obj.pca <- stats::prcomp(intensities.mtx, scale.=TRUE)
  # create object feature matrix
  obj.pca_featmtx <- obj.pca$rotation * crossprod(t(rep.int(1, nrow(obj.pca$rotation))),
                                                  summary(obj.pca)$importance[2,])

  res <- tibble(`__index_msobject__` = seq_len(nrow(obj.pca_featmtx)),
                tmp_profile_cluster = stats::cutree(hclust(dist(obj.pca_featmtx), method="single"),
                                                    min(c(nclu, nrow(obj.pca_featmtx), ncol(obj.pca_featmtx))))) %>%
    dplyr::inner_join(objs.df, by="__index_msobject__")
  # assign profile_cluster indices from largest to smallest clusters
  res_clustats <- dplyr::inner_join(res, obj_stats.df, by="__index_msobject__") %>%
    dplyr::group_by(tmp_profile_cluster) %>%
    dplyr::summarise(nsimilar_profiles = n(),
                     n_quants_cluster_total = sum(n_quants, na.rm=TRUE),
                     n_quants_cluster_median = median(n_quants, na.rm=TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(n_quants_cluster_total), desc(nsimilar_profiles)) %>%
    dplyr::mutate(profile_cluster = row_number())

  return(inner_join(res, res_clustats, by="tmp_profile_cluster") %>%
           dplyr::select(-tmp_profile_cluster, -`__index_msobject__`))
}

# annotate msdata (full subobject X msrun Cartesian product)
# as reliable (column name) according to
# the specificity of observing the quant object in the MS experiment groups
# defined by `msprobe$spec_msexp_group`
# or by co-occurrence of distinct quant subobjects in the given MS experiment group
# defined by `msprobe$cooccur_msexp_group`
annotate_msdata <- function(msdata_df, model_def, verbose = model_def$verbose,
                            specificity_pvalue = 1E-3,
                            specificity_quantobject_group_cols = NULL,
                            observation_fdr = 0.01,
                            cooccurrence_pvalue = 1E-3) {
  if (verbose) message("Identifying reliable quantifications")
  msdata_df <- dplyr::mutate(msdata_df, is_observed = !is.na(intensity))

  spec_quantobj_group_cols <- c('object_id', specificity_quantobject_group_cols %||%
                                       intersect(colnames(msdata_df),
                                                 c('subobject_id', 'object_id'))) %>% unique()
  spec_id_cols <- c("spec_msexp_group", spec_quantobj_group_cols)

  # calculate probabilities that observations of an object are specific to the given interaction
  spec_stats_df <- dplyr::group_by(msdata_df, spec_msexp_group, !!!syms(spec_quantobj_group_cols)) %>%
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
    # calculate probabilities that all quantitations of subobjects in a given observation are false discoveries
    cooccur_stats_df <- dplyr::group_by(msdata_df, object_id, cooccur_msexp_group) %>%
      dplyr::summarise(nsubobj_observed = n_distinct(subobject_id[is_observed]),
                       nsubobj_missed = n_distinct(subobject_id[!is_observed]),
                       .groups = "drop") %>%
      dplyr::mutate(is_cooccurring = pbinom(nsubobj_observed - 1L, nsubobj_observed + nsubobj_missed,
                                            observation_fdr, lower.tail=FALSE) <= cooccurrence_pvalue)

    msdata_df <- dplyr::left_join(msdata_df,
        dplyr::select(cooccur_stats_df, object_id, cooccur_msexp_group, is_cooccurring),
        by = c("object_id", "cooccur_msexp_group")) %>%
        dplyr::mutate(is_reliable = is_reliable | is_cooccurring,
                      is_cooccurring = NULL)
  }
  msdata_df <- dplyr::mutate(msdata_df, is_reliable = !is.na(intensity) & is_reliable)

  return(msdata_df)
}

# select relevant msdata and subobjects from msdata_df
# and define model_data$subobjects and model_data$observations
prepare_msdata <- function(model_data, msdata, verbose = model_data$model_def$verbose,
                           max_subobjects = 20L,
                           specificity_quantobject_group_cols = NULL,
                           ...) {
  model_def <- model_data$model_def
  msprb <- msdata$msentities[['msprobe']]
  msprb_idcol <- msprb
  mschan <- msdata$msentities[['mschannel']]
  mschan_idcol <- mschan
  modelobj <- msdata$msentities[['modelobject']]
  modelobj_idcol <- paste0(modelobj, "_id")
  quantobj <- msdata$msentities[['quantobject']]
  quantobj_idcol <- paste0(quantobj, "_id")
  intensities_dfname <- paste0(quantobj, "_intensities")
  if (!rlang::has_name(msdata, intensities_dfname)) {
    stop("No intensities (", intensities_dfname, " data frame) found in msdata")
  }

  if (modelobj == quantobj) {
    # modelobj is quanted directly
    intensities_df <- dplyr::select_at(msdata[[intensities_dfname]],
                                       c(object_id=modelobj_idcol, msprobe=msprb_idcol, "intensity"))
    msdata_df <- dplyr::left_join(model_data$observations,
                                  intensities_df, by = c("object_id", "msprobe")) %>%
        annotate_msdata(model_def) %>%
        dplyr::arrange(index_observation)
  } else {
    # quant specific quantobjects of modelobj
    modelobj2quantobj_df <- msdata[[paste0(modelobj, "2", quantobj)]]
    subobjs_df <- dplyr::inner_join(modelobj2quantobj_df,
                                    dplyr::select(model_data$objects, !!sym(modelobj_idcol), object_id, index_object),
                                    by=modelobj_idcol) %>%
      dplyr::filter(is_specific) %>%
      dplyr::inner_join(dplyr::select(msdata[[paste0(quantobj, "s")]], !!sym(quantobj_idcol),
                                      any_of(c("msfraction", "charge"))),
                        by=quantobj_idcol) %>%
      dplyr::mutate(subobject_id = !!sym(quantobj_idcol))
    if (nrow(subobjs_df) == 0L) stop("No specific ", model_def$quantobject, "s found for ", modelobj_idcol, "=", model_data$modelobj_id)
    if (verbose) message(nrow(subobjs_df), " specific ", model_def$quantobject, "(s) found")
    intensities_df <- dplyr::select_at(msdata[[intensities_dfname]],
                                       c(subobject_id = quantobj_idcol, mschannel=mschan_idcol, "intensity"))
    msdata_df <- dplyr::inner_join(dplyr::select(model_data$observations, index_msprobe, index_observation, index_object, object_id),
                                   dplyr::select(subobjs_df, index_object, subobject_id, any_of("msfraction")), by="index_object") %>%
        dplyr::inner_join(dplyr::select(model_data$mschannels, index_mschannel, mschannel, index_msprobe, any_of("msfraction")),
                          by = c("index_msprobe", intersect("msfraction", colnames(subobjs_df)))) %>%
        dplyr::inner_join(dplyr::select(model_data$msprobe, index_msprobe, msprobe,
                                        spec_msexp_group, cooccur_msexp_group,
                                        any_of(specificity_quantobject_group_cols)),
                          by=c("index_msprobe")) %>%
        dplyr::left_join(intensities_df, by=c("subobject_id", "mschannel"))
    if (all(is.na(msdata_df$intensity))) stop("No quantifications for ", nrow(subobjs_df), " specific ",
                                              quantobj, "(s) of ", modelobj_idcol, "=", model_data$object_id)
    msdata_df <- annotate_msdata(msdata_df, model_def, verbose=verbose,
                                 specificity_quantobject_group_cols = specificity_quantobject_group_cols,
                                 ...)

    # arrange pepmodstates by object, by profile cluster and by the number of quantitations
    subobject_group_size <- model_def$subobject_group_size %||% (max_subobjects %/% 2)
    subobj_info_cols <- intersect(colnames(subobjs_df), c("pepmod_id", "msfraction", "charge"))

    subobj_stats_df <- msdata_df %>%
      dplyr::group_by(index_object, subobject_id) %>%
      dplyr::summarise(n_quants = sum(!is.na(intensity)),
                       intensity_med = median(intensity, na.rm=TRUE),
                       .groups = "drop") %>%
      dplyr::inner_join(
        dplyr::group_by(msdata_df, index_object) %>%
        dplyr::group_modify(~ cluster_msprofiles(.x, dplyr::rename(msdata[[paste0(mschan, "_", quantobj, "_stats")]],
                                                                   mschannel = !!sym(mschan_idcol)),
                                                 obj_col='subobject_id', mschannel_col="mschannel")) %>%
        dplyr::ungroup(),
        by = c("subobject_id", "index_object")) %>%
      dplyr::left_join(dplyr::select(subobjs_df, object_id, subobject_id, is_specific,
                                     any_of(quantobj_idcol), !!!syms(subobj_info_cols)),
                       by = 'subobject_id') %>%
      dplyr::arrange(index_object, profile_cluster,
                     desc(is_specific), desc(n_quants), desc(intensity_med),
                     !!!syms(subobj_info_cols)) %>%
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
  #} else {
  #   stop("Unsupported combination of modelobject=", modelobj,
  #       " and quantobject=", quantobj)
  }
  # separately index quantifications and missing data
  model_data$msdata <- mutate(msdata_df,
                              index_qdata = if_else(is_observed, cumsum(is_observed), NA_integer_),
                              index_mdata = if_else(!is_observed, cumsum(!is_observed), NA_integer_))
  message(nrow(model_data$msdata), " ",
          if (rlang::has_name(model_data$msdata, 'index_subobservation'))
            'subobservation' else 'observation', '(s) of ',
          n_distinct(model_data$msdata$index_object), ' ', modelobj, '(s)',
          if (rlang::has_name(model_data$msdata, 'index_subobject'))
            paste0(' with ', n_distinct(model_data$msdata$index_subobject), ' ',
                   quantobj, '(s)') else '', ': ',
          sum(model_data$msdata$is_observed), " quantitation(s) (",
          sum(model_data$msdata$is_reliable), " reliable), ",
          sum(!model_data$msdata$is_observed), " missed")
  return(model_data)
}

#' Prepare *MSGLM* model input data for specified objects.
#'
#' @param model_def *msglm_model* object with MSGLM model definition
#' @param msdata *msglm_data_collection* object with all MS data
#' @param object_ids vector of *model objects* IDs to analyze
#' @param verbose
#'
#' @return object of *msglm_model_data* class
#' @export
#'
#' @examples
#' @seealso [msglm_model()], [import_msglm_data()], [to_standata()]
msglm_data <- function(model_def, msdata, object_ids, verbose = model_def$verbose,
                       mschannel_extra_cols = character(0),
                       mschannel_shift_col = paste0("total_", msdata$msentities[['mschannel']], "_shift"),
                       max_subobjects = 20L,
                       specificity_msexp_group_cols = msdata$msentities[['condition']],
                       cooccurrence_msexp_group_cols = msdata$msentities[['msprobe']],
                       ...) {
  checkmate::assert_class(model_def, "msglm_model")
  checkmate::assert_class(msdata, "msglm_data_collection")
  model_data <- list(model_def = model_def, object_id = object_ids,
                     msentities = msdata$msentities)
  modelobj <- msdata$msentities[['modelobject']]
  modelobj_idcol <- paste0(modelobj, "_id")

  quantobj <- msdata$msentities[['quantobject']]

  if (!rlang::has_name(msdata, paste0(modelobj, "s"))) {
    stop("No model object (", modelobj, ") information found in MS data")
  }
  modelobjs_df <- msdata$modelobjects
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

  msprb <- msdata$msentities[['msprobe']]
  if (is.na(msprb)) stop("No msprobe entity specified in msdata")
  mschan <- msdata$msentities[['mschannel']]
  if (is.na(mschan)) stop("No mschannel entity specified in msdata")
  condition <- msdata$msentities[['condition']]
  if (is.na(condition)) stop("No condition entity specified in msdata")
  msprb_idcol <- msprb
  mschan_idcol <- mschan
  mschan_parent = if_else(msprb == mschan, condition, msprb)

  msrun <- msdata$msentities[['msrun']]
  msexp <- msdata$msentities[['msexperiment']]

  msprbs_df <- msdata[[paste0(msprb,'s')]]
  if (verbose) message("Using msdata$", msprb, "s for MS probes information")
  checkmate::assert_data_frame(msprbs_df)
  msprb_cols <- c(msprobe = msprb_idcol,
                  condition = condition,
                  mstag = msdata$msentities[['mstag']])
  msprbs_orig_df <- msdata[[paste0(msprb,'s')]]
  checkmate::assert_names(colnames(msprbs_orig_df),
                          must.include = Filter(Negate(is.na), msprb_cols))
  msprbs_df <- dplyr::select_at(msprbs_orig_df, Filter(Negate(is.na), msprb_cols))
  if (!is.null(specificity_msexp_group_cols)) {
    msprbs_df <- dplyr::bind_cols(msprbs_df, dplyr::transmute(msprbs_orig_df,
                                  spec_msexp_group = factor(paste(!!!syms(specificity_msexp_group_cols)))))
  } else {
    msprbs_df <- dplyr::mutate(msprbs_df, spec_msexp_group = NA_integer_)
  }
  if (!is.null(cooccurrence_msexp_group_cols)) {
    msprbs_df <- dplyr::bind_cols(msprbs_df, dplyr::transmute(msprbs_orig_df,
                                  cooccur_msexp_group = factor(paste(!!!syms(cooccurrence_msexp_group_cols)))))
  } else {
    msprbs_df <- dplyr::mutate(msprbs_df, cooccur_msexp_group = NA_integer_)
  }
  checkmate::assert_character(as.character(msprbs_df$msprobe), unique=TRUE, any.missing=FALSE)
  checkmate::assert_subset(as.character(msprbs_df$condition), model_def$conditions$condition)
  checkmate::assert_set_equal(unique(as.character(msprbs_df$condition)),
                              dplyr::filter(model_def$conditions, !is_virtual)$condition)
  # fill missing columns with NA
  for (col in setdiff(names(msprb_cols), colnames(msprbs_df))) {
    if (verbose) message("Adding empty ", col, " column to MS probe information")
    msprbs_df[[col]] <- NA_character_
  }

  if (mschan != msprb) {
    if (!rlang::has_name(msdata, paste0(mschan, "s"))) {
      stop("MS channels data msdata$", mschan, " not found")
    }
    if (verbose) message("Using msdata$", mschan, "s for MS channels information")
    mschans_df <- msdata[[paste0(mschan,"s")]]
    checkmate::assert_data_frame(mschans_df)
    mschan_cols <- c(mschannel = mschan_idcol, msprobe = msprb_idcol,
                     msrun = msdata$msentities[['msrun']],
                     msfraction = msdata$msentities[['msfraction']],
                     mstag = msdata$msentities[['mstag']])
    checkmate::assert_names(colnames(mschans_df),
                            must.include = Filter(Negate(is.na), mschan_cols))
    mschans_df <- dplyr::select_at(mschans_df, Filter(Negate(is.na), mschan_cols))
    checkmate::assert_character(as.character(mschans_df$mschannel), unique=TRUE, any.missing=FALSE)
    checkmate::assert_set_equal(unique(as.character(mschans_df$msprobe)),
                                as.character(msprbs_df$msprobe))
    # fill missing columns with NA
    for (col in setdiff(names(mschan_cols), colnames(mschans_df))) {
      if (verbose) message("Adding empty ", col, " column to MS channel information")
      mschans_df[[col]] <- NA_character_
    }
  } else {
    mschans_df <- dplyr::mutate(msprbs_df, mschannel = msprobe)
  }
  # add normalization shifts to msprobes
  mschan_shifts_dfname <- NA_character_
  for (kv in list(c("mschannel", mschan), c("msprobe", msprb), c("msrun", msrun), c("msexperiment", msexp))) {
    msentity <- kv[[1]]
    msentity_name <- kv[[2]]
    if (!is.na(msentity_name) && rlang::has_name(msdata, paste0(msentity_name, "_shifts"))) {
      mschan_shifts_dfname <- paste0(msentity_name, "_shifts")
      mschan_shifts_idcol <- msentity_name
      mschan_shift_refcol <- msentity
      break
    }
  }
  if (is.na(mschan_shifts_dfname)) {
    stop("Cannot find MS experiments normalization information ",
         "(msdata$", mschan, "_shifts or msdata$", msdata$msentities[['msrun']], "_shifts etc)")
  }
  if (verbose) message("Using msdata$", mschan_shifts_dfname, "$", mschannel_shift_col, " for MS experiments normalization")
  mschan_shifts_df <- msdata[[mschan_shifts_dfname]]
  checkmate::assert_data_frame(mschan_shifts_df, .var.name = paste0("msdata$", mschan_shifts_dfname))
  checkmate::assert_names(names(mschan_shifts_df), must.include = c(mschan_shifts_idcol, mschannel_shift_col))
  checkmate::assert_set_equal(as.character(mschan_shifts_df[[mschan_shifts_idcol]]),
                              as.character(mschans_df[[mschan_shift_refcol]]))
  mschan_shifts_df = dplyr::select_at(mschan_shifts_df,
                                      c(mschan_shifts_idcol, mschannel_shift=mschannel_shift_col))
  checkmate::assert_numeric(mschan_shifts_df$mschannel_shift, any.missing = FALSE, finite=TRUE, lower=-20, upper=20)
  mschans_df <- dplyr::inner_join(mschans_df, mschan_shifts_df,
                                  by=rlang::set_names(mschan_shifts_idcol, mschan_shift_refcol))
  msprbs_df <- msprbs_df %>%
    dplyr::inner_join(dplyr::select(model_def$conditions, condition, index_condition, any_of("mstag")),
                      by=intersect(c('condition', 'mstag'), colnames(model_def$conditions))) %>%
    dplyr::arrange(dplyr::across(c(index_condition, any_of(c("mstag", "msexperiment")), msprobe)))
  if (rlang::has_name(model_def, "msprobeXeffect")) {
    msprb_dim <- names(dimnames(model_def$msprobeXeffect))
    if (msprb_dim != msprb_idcol) {
      stop("MS experiments object (", msprb_idcol, ") does not match",
           " the name of model_def$msprobeXeffect rows dimension (",
           msprb_dim, ")")
    }
    msprobes_order <- rownames(model_def$msprobeXeffect)
  } else {
    msprobes_order <- NULL
  }
  # fix the order of msprobes
  msprbs_df <- ensure_primary_index_column(
        msprbs_df, 'index_msprobe',
        id_col="msprobe", ids_ordered=msprobes_order,
        create=TRUE, paste0("msdata$", msprb, "s"))
  # fix the order of mschannels
  mschan_import_cols <- c("msprobe", "index_msprobe", "index_condition")
  if (msprb != mschan) {
    mschan_import_cols <- append(mschan_import_cols, c("condition"))
  }
  mschans_df <- dplyr::inner_join(mschans_df, dplyr::select_at(msprbs_df, mschan_import_cols),
                                  by='msprobe') %>%
      dplyr::arrange(dplyr::across(c(index_msprobe, any_of(c("mstag", "msfraction", "msrun")))))
  model_data$msprobes <- msprbs_df
  model_data$mschannels <- dplyr::mutate(mschans_df,
      index_msrun = match(msrun, unique(msrun)),
      index_mschannel = row_number(),
      index_mscalib = 1L) # FIXME there could be multiple mscalib per dataset

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
  model_data$observations <- dplyr::inner_join(model_data$interactions, model_data$msprobes,
                                               by = c("index_condition", "condition")) %>%
    dplyr::arrange(index_object, index_interaction, index_msprobe) %>%
    dplyr::mutate(index_observation = row_number(),
                  observation_id = paste0(object_id, '_', msprobe))

  model_data$quantobj_mscalib <- msdata[[paste0(quantobj, "_mscalib")]]
  model_data$quantobj_labu_shift <- msdata[[paste0(quantobj, '_labu_shift')]]
  model_data$quantobj_labu_min <- msdata[[paste0(quantobj, '_labu_min')]]
  logbase <- logintensityBase(model_data$quantobj_mscalib, silent=TRUE)
  if (logbase != 2) {
    warning("msdata$", quantobj, "_mscalib logintensityBase=", logbase,
            " converting mscalib model and log-intensities to log2-based ones")
    model_data$quantobj_mscalib <- convert_logintensityBase(model_data$quantobj_mscalib, new_base=2)
    #k <- log(logbase) / log(2)
    #model_data$quantobj_labu_shift <- model_data$quantobj_labu_shift * k
    #model_data$quantobj_labu_min <- model_data$quantobj_labu_min * k
  }

  model_data <- prepare_msdata(model_data, msdata, max_subobjects=max_subobjects,
                               verbose=verbose, ...)
  model_data <- prepare_expanded_effects(model_data, verbose=verbose)

  return (structure(model_data, class="msglm_model_data"))
}
