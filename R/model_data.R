# use experimental design matrices and add
# object-specific object/replicate/batch effects information to model_data
prepare_expanded_effects <- function(model_data, verbose=model_data$model_def$verbose)
{
  model_def <- model_data$model_def

  # FIXME use model_def class information
  is_glmm <- rlang::has_name(model_def, "mixtureXcondition.mtxs") &&
             rlang::has_name(model_def, "mixeffects.df")
  if (is_glmm) {
    if (verbose) message("Detected mixeffects data for GLMM model")
    mixtureXcondition.mtxs <- model_def$mixtureXcondition.mtxs
    mixcoefXeff.mtx <- model_def$mixcoefXeff.mtx
    # "GLMM" mixing model
    if (any(names(mixtureXcondition.mtxs) != rownames(mixcoefXeff.mtx))) {
      stop("Mismatch between mixtureXcondition.mtxs matrix names and mix coefficient names")
    }
    mixeffects_df <- model_def$mixeffects
    if (any(mixeffects_df$mixeffect != colnames(mixcoefXeff.mtx))) {
      stop("Mismatch between mixcoefXeff.mtx matrix colnames and mix effect names")
    }
    for (i in seq_along(mixtureXcondition.mtxs)) {
      if (any(colnames(mixtureXcondition.mtxs[[i]]) != model_data$conditions$condition)) {
        stop("Mismatch between mixtureXcondition.mtxs[[", names(mixtureXcondition.mtxs)[[i]],
             "]] matrix colnames and condition names")
      }
      if (any(rownames(mixtureXcondition.mtxs[[i]]) != model_data$mixtures$mixture)) {
        stop("Mismatch between mixtureXcondition.mtxs[[", names(mixtureXcondition.mtxs)[[i]],
             "]] matrix rownames and mixture names")
      }
    }
    model_data$mixeffects <- mixeffects_df
    model_data$mixcoefs <- tibble(mixcoef = rownames(mixcoefXeff.mtx))
    model_data$mixcoefXeff <- mixcoefXeff.mtx
    model_data$mixtureXcondition <- mixtureXcondition.mtxs
    conditions <- levels(model_data$conditions$condition)
    obj_mixXcond.mtxs <- lapply(model_data$mixtureXcondition, function(mixtureXcondition.mtx){
      res <- extrude_matrix(mixtureXcondition.mtx,
                            model_data$object_mixtures$index_object,
                            model_data$object_mixtures$index_mixture,
                            blockdim = "index_object", rowixdim = "index_object_mixture",
                            coldim = "object_condition_id", corecoldim = "condition_id")
      names(dimnames(res$mtx)) <- c("sact", "iact")
      res$mtx <- res$mtx[, colSums(abs(res$mtx)) != 0.0, drop=FALSE]
      res$df <- dplyr::rename(res$df, condition=eff, objcond_id=objeff)
      res$objconds_df <- dplyr::rename(res$objeff_df, condition=eff,
                                      objcond_id=objeff,
                                      index_object=obj)
      res$objeff_df <- NULL
      return(res)
    }) %>% pivot_wider(id_cols = )
    # collect all used object_condition ids to use as factor levels
    objcond_ids = unique(unlist(lapply(obj_mixXcond.mtxs, function(x) x$objconds_df$objcond_id)))
    # convert condition and objcond_id into factor
    obj_mixXcond.mtxs <- lapply(obj_mixXcond.mtxs, function(x) {
      x$df <- mutate(x$df,
                     condition = factor(condition, levels=conditions),
                     objcond_id = factor(objcond_id, levels=objcond_ids))
      x$objconds_df <- mutate(x$objconds_df,
                             condition = factor(condition, levels=conditions),
                             objcond_id = factor(objcond_id, levels=objcond_ids))
      return(x)
    })
    model_data$object_conditions <- dplyr::distinct(dplyr::bind_rows(lapply(obj_mixXcond.mtxs, function(x) x$objconds_df))) %>%
      dplyr::mutate(index_condition = as.integer(condition)) %>%
      dplyr::mutate(index_object_condition = row_number(),
             is_virtual = FALSE)
    model_data$mixtions <- dplyr::bind_rows(lapply(names(obj_mixXcond.mtxs), function(mixcoef){
      tibble::tibble(mixcoef = mixcoef,
                     objcond_id = factor(colnames(obj_mixXcond.mtxs[[mixcoef]]$mtx), levels=objcond_ids))
    })) %>%
      dplyr::mutate(mixtion_ix = row_number(),
                    object_premix = paste0(mixcoef, " ", objcond_id),
                    mixcoef = factor(mixcoef, levels=rownames(mixcoefXeff.mtx)),
                    mixcoef_ix = replace_na(as.integer(mixcoef), 0L)) %>%
    dplyr::left_join(dplyr::select(model_data$object_conditions, objcond_id, index_object_condition)) %>%
    dplyr::arrange(mixtion_ix)
    model_data$object_mixtureXpremix <- bind_cols(lapply(obj_mixXcond.mtxs, function(x) x$mtx))
    dimnames(model_data$object_mixtureXpremix) <- list(object_mixture = model_data$object_mixtures$object_mixture_id,
                                                       object_premix = model_data$object_premixes$object_premix)
  } else {
    objcond_ids <- unique(model_data$object_conditions$objcond_id)
    objprobe_ids <- unique(model_data$object_msprobes$objprobe_id)
  }

  obj_condXeff_df <- dplyr::full_join(matrix2frame(model_def$conditionXeffect, row_col="condition", col_col = "effect"),
                                     dplyr::select(model_data$objects, index_object, object_id),
                                     by = character()) %>%
    dplyr::inner_join(dplyr::select(model_def$effects, effect, index_effect), by="effect") %>%
    dplyr::inner_join(dplyr::select(model_data$object_conditions, index_object_condition, index_object, condition),
                      by=c("index_object", "condition")) %>%
    dplyr::mutate(object_effect = paste0(effect, '@', object_id))
  model_data$object_effects <- dplyr::select(obj_condXeff_df, index_object, object_id,
                                             index_effect, effect, object_effect) %>%
    dplyr::distinct() %>%
    dplyr::arrange(index_object, index_effect) %>%
    dplyr::mutate(index_object_effect = row_number())
  model_data$object_conditionXeffect <- frame2matrix(obj_condXeff_df, row_col="index_object_condition", col_col="object_effect",
                                          rows=model_data$object_conditions$index_object_condition,
                                          cols=model_data$object_effects$object_effect)

  msprobeXeffect_name <- "msprobeXeffect"
  msprb <- model_data$msentities[['msprobe']]
  msprb_idcol <- msprb
  if (rlang::has_name(model_def, "msprobeXeffect")) {
    if (verbose) message("Using model_def$", msprobeXeffect_name, " for per-", msprobe, " design matrix")
    msprb_dim <- names(dimnames(model_def[[msprobeXeffect_name]]))[[1]]
    if (msprb_dim != msprb_idcol) {
      stop("Name of model_def$", msprobeXeffect_name, " rows dimension (", msprb_dim, ") inconsistent with ", msprb_idcol)
    }
    obj_probeXeff_df <- dplyr::full_join(matrix2frame(model_def[[msprobeXeffect_name]], row_col="condition", col_col = msprb_idcol),
                                         dplyr::select(model_data$objects, index_object, object_id),
                                         by = character()) %>%
      dplyr::inner_join(dplyr::select(model_data$object_msprobes, index_object_msprobe, index_object, index_mschannel, !!sym(msprb_idcol)),
                        by=c("index_object", msprb_idcol)) %>%
      dplyr::mutate(object_effect = paste0(effect, '@', object_id))
  } else {
    obj_probeXeff_df <- dplyr::inner_join(obj_condXeff_df,
                                          dplyr::select(model_data$object_msprobes, index_object_condition, index_object_msprobe),
                                          by="index_object_condition")
  }
  model_data$object_msprobeXeffect <- frame2matrix(obj_probeXeff_df, row_col="index_object_msprobe", col_col="object_effect",
                                           rows=model_data$object_msprobes$index_object_msprobe,
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
  obj_probeXbatcheff_df <- dplyr::full_join(matrix2frame(msprobeXbatchEffect, row_col="msprobe", col_col="batch_effect"),
                                            dplyr::select(model_data$objects, index_object, object_id, object_label),
                                            by = character()) %>%
    dplyr::inner_join(dplyr::select(batch_effects_df, batch_effect, index_batch_effect), by="batch_effect") %>%
    dplyr::left_join(dplyr::select(model_data$object_msprobes, msprobe, index_object, index_object_msprobe),
                     by=c("index_object", "msprobe")) %>%
    dplyr::mutate(object_batch_effect = paste0(batch_effect, '@', object_id))
  model_data$object_batch_effects <- dplyr::select(obj_probeXbatcheff_df,
                                                   object_batch_effect,
                                                   index_batch_effect, batch_effect,
                                                   index_object, object_id) %>%
    dplyr::distinct() %>%
    dplyr::arrange(index_object, index_batch_effect) %>%
    dplyr::mutate(index_object_batch_effect = row_number())
  model_data$object_msprobeXbatch_effect <- frame2matrix(obj_probeXbatcheff_df, row_col="index_object_msprobe", col_col="object_batch_effect",
                                                rows=if (nrow(obj_probeXbatcheff_df) > 0) {
                                                    model_data$object_msprobes$index_object_msprobe
                                                  } else integer(0),
                                                cols=model_data$object_batch_effects$object_batch_effect)
  if (rlang::has_name(model_data, "quantobjects")) {
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
    qobj_probeXbatcheff_df <- dplyr::full_join(matrix2frame(mschannelXquantBatchEffect,
                                                            row_col="mschannel", col_col="quant_batch_effect"),
                                               dplyr::select(model_data$quantobjects, index_object, index_quantobject, quantobject_id),
                                               by = character()) %>%
      dplyr::inner_join(dplyr::select(quant_batch_effects_df, quant_batch_effect, index_quant_batch_effect), by="quant_batch_effect") %>%
      dplyr::right_join(dplyr::select(model_data$msdata, mschannel, index_msprobe, index_object, index_quantobject, index_quantobject_msprobe),
                        by=c("index_object", "index_quantobject", "mschannel")) %>%
      # index msdata without batch effect as 0 (needed for ref quant object removal later)
      dplyr::mutate(index_quant_batch_effect = replace_na(index_quant_batch_effect, 0L)) %>%
      dplyr::group_by(index_quantobject) %>%
      # remove batch effects that are present in all msdata of a given quantobject
      # (otherwise it creates redundancy)
      dplyr::group_modify(~{
        nmschans <- n_distinct(.x$mschannel)
        dplyr::group_by(.x, index_quant_batch_effect) %>%
        dplyr::filter(n() < nmschans) %>%
        dplyr::ungroup()
      }) %>%
      # remove reference quantobject for each object (one with the smallest index
      # among those that are detected in all batch_effect groups)
      dplyr::group_by(index_object) %>%
      dplyr::group_modify(~{
        # TODO make igraph dependency optional
        # indentify quant_batch_effect groups that don't share probes between each other
        # (which means that quant_batch_effects should have 1 degree-of-freedom less,
        #  because intensities in these groups are, strictly speaking, not comparable)
        # quant batch effects can share probes if the same biological sample (i.e. probe)
        # was analysed more than once
        probe2batcheff_df <- dplyr::select(.x, index_msprobe, index_quant_batch_effect) %>%
          dplyr::transmute(index_msprobe, index_quant_batch_effect,
                           probe_id = paste0("p", index_msprobe),
                           batcheff_id = paste0("e", index_quant_batch_effect)) %>%
          dplyr::distinct()
        probe2batcheff_graph <- igraph::graph_from_data_frame(
          dplyr::select(probe2batcheff_df, batcheff_id, probe_id), directed = FALSE)
        batcheffs <- as.integer(igraph::components(probe2batcheff_graph)$membership)
        batcheffs_df <- dplyr::select(probe2batcheff_df, batcheff_id, index_quant_batch_effect) %>%
          dplyr::distinct() %>% dplyr::inner_join(
            tibble::as_tibble(batcheffs, rownames = "batcheff_id"),
            by = "batcheff_id") %>%
          dplyr::select(index_quant_batch_effect,
                        index_quant_batch_effect_group = value)
        ngroups <- dplyr::n_distinct(batcheffs_df$index_quant_batch_effect_group)

        if (ngroups > 1L) {
          .x <- dplyr::inner_join(.x, batcheffs_df, by="index_quant_batch_effect")
          # select potential reference objects that are shared by effect groups
          overlap_qobjs_df <- dplyr::group_by(.x, index_quantobject) %>%
            dplyr::summarise(ngroups_qobj = n_distinct(index_quant_batch_effect_group), .groups = "drop") %>%
            dplyr::filter(ngroups_qobj == ngroups)
          if (nrow(overlap_qobjs_df) > 0L) {
            # remove one reference quantobject batch effects from each group to
            # eliminate redundancy
            del_qobj_batch_effs_df <- dplyr::semi_join(.x, overlap_qobjs_df, by ="index_quantobject") %>%
              dplyr::filter(index_quant_batch_effect > 0L) %>%
              dplyr::arrange(index_quant_batch_effect_group, index_quantobject, index_quant_batch_effect) %>%
              dplyr::group_by(index_quant_batch_effect) %>%
              dplyr::filter(row_number() == 1L) %>% dplyr::ungroup()
            .x <- dplyr::anti_join(.x, del_qobj_batch_effs_df,
                                   by =c("index_quant_batch_effect",
                                         "index_quantobject"))
          }
        }
        .x
      }) %>% dplyr::ungroup() %>%
      dplyr::filter(index_quant_batch_effect > 0L) %>%
      dplyr::mutate(quantobject_batch_effect = paste0(quant_batch_effect, '@', quantobject_id)) %>%
      dplyr::arrange(index_quant_batch_effect, index_quantobject_msprobe)
    model_data$quantobject_batch_effects <- dplyr::select(qobj_probeXbatcheff_df,
                                                        quantobject_batch_effect,
                                                        index_quant_batch_effect, quant_batch_effect,
                                                        index_quantobject, quantobject_id) %>%
      dplyr::distinct() %>%
      dplyr::arrange(index_quantobject, index_quant_batch_effect) %>%
      dplyr::mutate(index_quantobject_batch_effect = row_number())
    model_data$quantobject_msprobeXquant_batch_effect <- frame2matrix(
            qobj_probeXbatcheff_df,
            row_col="index_quantobject_msprobe",
            col_col="quantobject_batch_effect",
            rows = if (nrow(qobj_probeXbatcheff_df)>0) {
                    model_data$msdata$index_quantobject_msprobe
                  } else integer(0),
            cols = model_data$quantobject_batch_effects$quantobject_batch_effect)
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
cluster_msprofiles <- function(msdata, mschannel_stats, obj_col, mschannel_col, nclu=4) {
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
  mschannels <- unique(intensities.df[[mschannel_col]])
  if (nrow(objs.df) == 1L || length(mschannels) == 1L) {
    return(tibble(!!obj_col := objs.df[[obj_col]],
                  profile_cluster = 1L,
                  nsimilar_profiles = 1L))
  }
  if (nrow(intensities.df) != length(mschannels)*nrow(objs.df)) {
    stop("Duplicate intensities detected, check the input data")
  }
  obj_stats.df <- group_by(intensities.df, `__index_msobject__`) %>%
    summarise(n_quants = sum(!is.na(intensity))) %>%
    dplyr::ungroup()
  # add a bit of noise to avoid zero variance
  intensities.mtx <- matrix(log2(pmax(if_else(is.finite(intensities.df$intensity_imputed),
                                              intensities.df$intensity_imputed, 0.0), 1)) +
                            rnorm(nrow(intensities.df), sd=0.01),
                            ncol = nrow(objs.df),
                            dimnames = list(mschannel = mschannels, `__index_msobject__` = NULL))
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

# annotate msdata (full quantobject X msrun Cartesian product)
# as reliable (column name) according to
# the specificity of observing the quant object in the MS experiment groups
# defined by `msprobe$spec_msexp_group`
# or by co-occurrence of distinct quant quantobjects in the given MS experiment group
# defined by `msprobe$cooccur_msexp_group`
annotate_msdata <- function(msdata_df, model_def, verbose = model_def$verbose,
                            specificity_pvalue = 1E-3,
                            specificity_quantobject_group_cols = NULL,
                            quantobject_fdr = 0.01,
                            cooccurrence_pvalue = 1E-3) {
  if (verbose) message("Identifying reliable quantifications")
  msdata_df <- dplyr::mutate(msdata_df, is_observed = !is.na(intensity))

  spec_quantobj_group_cols <- c('object_id', specificity_quantobject_group_cols %||%
                                       intersect(colnames(msdata_df),
                                                 c('quantobject_id', 'object_id'))) %>% unique()
  spec_id_cols <- c("spec_msexp_group", spec_quantobj_group_cols)

  # calculate probabilities that object_msprobes of an object are specific to the given object_condition
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
  if (rlang::has_name(msdata_df, 'quantobject_id')) {
    # calculate probabilities that all quantitations of quantobjects in a given object_msprobe are false discoveries
    cooccur_stats_df <- dplyr::group_by(msdata_df, object_id, cooccur_msexp_group) %>%
      dplyr::summarise(nqobj_observed = n_distinct(quantobject_id[is_observed]),
                       .groups = "drop") %>%
      dplyr::inner_join(dplyr::group_by(msdata_df, object_id) %>%
                        dplyr::summarise(nqobj_total = n_distinct(quantobject_id),
                                         .groups = "drop"), by="object_id") %>%
      dplyr::mutate(is_cooccurring = pbinom(nqobj_observed - 1L, nqobj_total,
                                            quantobject_fdr, lower.tail=FALSE) <= cooccurrence_pvalue)

    msdata_df <- dplyr::left_join(msdata_df,
        dplyr::select(cooccur_stats_df, object_id, cooccur_msexp_group, is_cooccurring),
        by = c("object_id", "cooccur_msexp_group")) %>%
        dplyr::mutate(is_reliable = is_reliable | is_cooccurring,
                      is_cooccurring = NULL)
    if (rlang::has_name(msdata_df, "ident_type")) {
      msdata_df <- dplyr::group_by(msdata_df, quantobject_id, spec_msexp_group) %>%
        dplyr::mutate(is_reliable = is_reliable | any(ident_type %in% c("MULTI-MSMS"))) %>%
        dplyr::ungroup()
    }
  }
  msdata_df <- dplyr::mutate(msdata_df, is_reliable = !is.na(intensity) & is_reliable)

  return(msdata_df)
}

# select relevant msdata and quantobjects from msdata_df
# and define model_data$quantobjects and model_data$object_msprobes
prepare_msdata <- function(model_data, msdata, verbose = model_data$model_def$verbose,
                           max_quantobjects = 20L,
                           specificity_quantobject_group_cols = NULL,
                           ...) {
  model_def <- model_data$model_def
  msprb <- msdata$msentities[['msprobe']]
  msprb_idcol <- msprb
  mschan <- msdata$msentities[['mschannel']]
  mschan_idcol <- mschan
  obj <- msdata$msentities[['object']]
  obj_idcol <- paste0(obj, "_id")
  quantobj <- msdata$msentities[['quantobject']]
  quantobj_idcol <- paste0(quantobj, "_id")
  intensities_dfname <- paste0(quantobj, "_intensities")
  if (!rlang::has_name(msdata, intensities_dfname)) {
    stop("No intensities (", intensities_dfname, " data frame) found in msdata")
  }

  if (obj == quantobj) {
    # obj is quanted directly
    intensities_df <- dplyr::select(msdata[[intensities_dfname]],
                                    object_id=!!sym(obj_idcol), msprobe=!!sym(msprb_idcol),
                                    intensity, any_of("ident_type"))
    msdata_df <- dplyr::left_join(model_data$object_msprobes,
                                  intensities_df, by = c("object_id", "msprobe")) %>%
        annotate_msdata(model_def) %>%
        dplyr::arrange(index_object_msprobe)
  } else {
    # quant specific quantobjects of object
    obj2quantobj_df <- msdata[[paste0(obj, "2", quantobj)]]
    qobjs_df <- dplyr::inner_join(obj2quantobj_df,
                                  dplyr::select(model_data$objects, !!sym(obj_idcol), object_id, index_object),
                                  by=obj_idcol) %>%
      dplyr::filter(is_specific) %>%
      dplyr::inner_join(dplyr::select(msdata[[paste0(quantobj, "s")]], !!sym(quantobj_idcol),
                                      any_of(c("msfraction", "charge"))),
                        by=quantobj_idcol) %>%
      dplyr::mutate(quantobject_id = !!sym(quantobj_idcol))
    if (nrow(qobjs_df) == 0L) stop("No specific ", quantobj, "s found for ",
                                   obj_idcol, "=", model_data$object_id)
    if (verbose) message(nrow(qobjs_df), " specific ", quantobj, "(s) found")
    intensities_df <- dplyr::select(msdata[[intensities_dfname]],
                                    quantobject_id=!!sym(quantobj_idcol), mschannel=!!sym(mschan_idcol),
                                    intensity, any_of("ident_type"))
    msdata_df <- dplyr::inner_join(dplyr::select(model_data$object_msprobes, index_msprobe, index_object_msprobe, index_object, object_id),
                                   dplyr::select(qobjs_df, index_object, quantobject_id, any_of("msfraction")), by="index_object") %>%
        dplyr::inner_join(dplyr::select(model_data$mschannels, index_mschannel, mschannel, index_msprobe, index_msprotocol, any_of(c("msfraction", "msprotocol"))),
                          by = c("index_msprobe", intersect("msfraction", colnames(qobjs_df)))) %>%
        dplyr::inner_join(dplyr::select(model_data$msprobe, index_msprobe, msprobe,
                                        spec_msexp_group, cooccur_msexp_group,
                                        any_of(specificity_quantobject_group_cols)),
                          by=c("index_msprobe")) %>%
        dplyr::left_join(intensities_df, by=c("quantobject_id", "mschannel"))
    if (all(is.na(msdata_df$intensity))) stop("No quantifications for ", nrow(qobjs_df), " specific ",
                                              quantobj, "(s) of ", obj_idcol, "=", model_data$object_id)
    if (n_distinct(msdata_df$index_msprotocol) > 1L) {
      # removing missing quantitations if a given quantobject is completely missing in a given msprotocol
      msdata_df <- dplyr::group_by(msdata_df, quantobject_id, index_msprotocol) %>%
        dplyr::filter(any(!is.na(intensity))) %>%
        dplyr::ungroup()
    }
    nqobj_probes <- nrow(dplyr::distinct(msdata_df, quantobject_id, mschannel))
    if (nrow(msdata_df) != nqobj_probes) {
      if (verbose) warning(nrow(msdata_df) - nqobj_probes, " of ", nrow(msdata_df),
                           " ", quantobj, " MS intensities are duplicate, summing duplicate entries")
      msdata_df <- dplyr::group_by(msdata_df, dplyr::across(!any_of(c("intensity", "ident_type")))) %>%
        dplyr::summarise(dplyr::across(intensity, ~sum(.x, na.rm = TRUE)),
                         dplyr::across(any_of("ident_type"), ~min(.x, na.rm = TRUE)),
                         .groups="drop") %>%
        dplyr::mutate(intensity = if_else(intensity != 0, intensity, NA_real_))
    }
    msdata_df <- annotate_msdata(msdata_df, model_def, verbose=verbose,
                                 specificity_quantobject_group_cols = specificity_quantobject_group_cols,
                                 ...)

    # arrange pepmodstates by object, by profile cluster and by the number of quantitations
    quantobject_group_size <- model_def$quantobject_group_size %||% (max_quantobjects %/% 2)
    qobj_info_cols <- intersect(colnames(qobjs_df), c("pepmod_id", "msfraction", "charge"))

    qobj_stats_df <- msdata_df %>%
      dplyr::group_by(index_object, quantobject_id) %>%
      dplyr::summarise(n_quants = sum(!is.na(intensity)),
                       intensity_med = median(intensity, na.rm=TRUE),
                       .groups = "drop") %>%
      dplyr::inner_join(
        dplyr::group_by(msdata_df, index_object) %>%
        dplyr::group_modify(~ cluster_msprofiles(.x, dplyr::rename(msdata[[paste0(mschan, "_", quantobj, "_stats")]],
                                                                   mschannel = !!sym(mschan_idcol)),
                                                 obj_col='quantobject_id', mschannel_col="mschannel")) %>%
        dplyr::ungroup(),
        by = c("quantobject_id", "index_object")) %>%
      dplyr::left_join(dplyr::select(qobjs_df, object_id, quantobject_id, is_specific,
                                     any_of(quantobj_idcol), !!!syms(qobj_info_cols)),
                       by = 'quantobject_id') %>%
      dplyr::arrange(index_object, profile_cluster,
                     desc(is_specific), desc(n_quants), desc(intensity_med),
                     !!!syms(qobj_info_cols)) %>%
      dplyr::group_by(index_object, profile_cluster) %>%
      dplyr::mutate(index_quantobject_group = row_number() %/% quantobject_group_size, # put objects within cluster into groups of 20
                    index_quantobject_local = row_number() %% quantobject_group_size) %>%
      dplyr::ungroup()
    # take the first group of 10 objects from each cluster, then continue with the second group etc
    model_data$quantobjects <- qobj_stats_df %>%
      dplyr::arrange(index_object, index_quantobject_group, profile_cluster, index_quantobject_local) %>%
      dplyr::mutate(index_quantobject = row_number()) %>%
      dplyr::filter(index_quantobject <= max_quantobjects) # remove less abundant quantobjects of rich objects
    if (verbose) message(nrow(model_data$quantobjects), " ", quantobj, "(s) from ",
                         n_distinct(model_data$quantobjects$profile_cluster), " cluster(s) selected")
    msdata_df <- dplyr::inner_join(msdata_df,
                                   dplyr::select(model_data$quantobjects, index_object, index_quantobject, quantobject_id),
                                   by=c("index_object", "quantobject_id")) %>%
      dplyr::arrange(index_object, index_object_msprobe, index_quantobject) %>%
      dplyr::mutate(index_quantobject_msprobe = row_number())
  #} else {
  #   stop("Unsupported combination of object=", obj,
  #       " and quantobject=", quantobj)
  }
  # separately index quantifications and missing data
  model_data$msdata <- mutate(msdata_df,
                              index_qdata = if_else(is_observed, cumsum(is_observed), NA_integer_),
                              index_mdata = if_else(!is_observed, cumsum(!is_observed), NA_integer_))
  message(nrow(model_data$msdata), " ",
          if (rlang::has_name(model_data$msdata, 'index_quantobject_msprobe'))
            quantobj else obj, '-in-msprobe(s) of ',
          n_distinct(model_data$msdata$index_object), ' ', obj, '(s)',
          if (rlang::has_name(model_data$msdata, 'index_quantobject'))
            paste0(' with ', n_distinct(model_data$msdata$index_quantobject), ' ',
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
                       max_quantobjects = 20L,
                       specificity_msexp_group_cols = msdata$msentities[['condition']],
                       cooccurrence_msexp_group_cols = msdata$msentities[['msprobe']],
                       ...) {
  checkmate::assert_class(model_def, "msglm_model")
  checkmate::assert_class(msdata, "msglm_data_collection")
  model_data <- list(model_def = model_def, object_id = object_ids,
                     msentities = msdata$msentities)
  obj <- msdata$msentities[['object']]
  obj_idcol <- paste0(obj, "_id")

  quantobj <- msdata$msentities[['quantobject']]

  if (!rlang::has_name(msdata, paste0(obj, "s"))) {
    stop("No model object (", obj, ") information found in MS data")
  }
  objs_df <- msdata$objects
  obj_cols <- msdata$msentities_extra_columns$object
  # FIXME remove this if
  if (is.null(obj_cols)) {
    if (verbose) warning("msdata$msentities_extra_columns$object not found, resorting to temorary guess")
    obj_cols <- (c(intersect(c("majority_protein_acs", "protein_acs",
                     "gene_names", "protein_names"), colnames(objs_df)),
         str_subset(colnames(objs_df), "^is_")) %>% unique())
    if (verbose) {
      message("Model object (", obj, ") columns to use: ", paste0(obj_cols, collapse=", "))
    }
  }
  model_data$objects <- dplyr::filter(objs_df, object_id %in% object_ids) %>%
      dplyr::select_at(c("object_id", "object_label", obj_idcol,
                         obj_cols) %>% unique()) %>%
      dplyr::arrange(object_id) %>%
      dplyr::mutate(index_object = row_number())
  missing_obj_ids <- setdiff(object_ids, unique(objs_df$object_id))
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
  checkmate::assert_character(as.character(msprbs_df$msprobe), unique=TRUE, any.missing=FALSE)
  checkmate::assert_subset(as.character(msprbs_df$condition), model_def$conditions$condition)
  checkmate::assert_set_equal(unique(as.character(msprbs_df$condition)),
                              dplyr::filter(model_def$conditions, !is_virtual)$condition)
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
                     mstag = msdata$msentities[['mstag']],
                     msprotocol = msdata$msentities[['msprotocol']])
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
      index_msprotocol = if (rlang::has_name(mschans_df, "msprotocol")) {
          match(msprotocol, unique(msprotocol))
        } else {
          1L
        }
    )

  # all objects X all conditions (including virtual)
  # FIXME rename objectXcondition
  model_data$object_conditions <- dplyr::full_join(dplyr::select(model_data$objects, object_id, index_object),
                                              dplyr::select(model_def$conditions, index_condition, is_virtual, condition),
                                              by=character()) %>%
    dplyr::mutate(objcond_id = paste(condition, object_id)) %>%
    dplyr::arrange(index_condition, index_object) %>%
    dplyr::mutate(index_object_condition = row_number())

  # all objects X all MS channels (only those with actual experiments)
  # FIXME rename objectXmschannel
  model_data$object_msprobes <- dplyr::inner_join(model_data$object_conditions, model_data$msprobes,
                                               by = c("index_condition", "condition")) %>%
    dplyr::arrange(index_object, index_object_condition, index_msprobe) %>%
    dplyr::mutate(index_object_msprobe = row_number(),
                  objprobe_id = paste0(object_id, '_', msprobe))

  model_data$quantobject_mscalib <- msdata[[paste0(quantobj, "_mscalib")]]
  model_data$quantobject_labu_min <- msdata[[paste0(quantobj, '_labu_min')]]
  logbase <- logintensityBase(model_data$quantobject_mscalib, silent=TRUE)
  if (logbase != 2) {
    warning("msdata$", quantobj, "_mscalib logintensityBase=", logbase,
            " converting mscalib model and log-intensities to log2-based ones")
    model_data$quantobject_mscalib <- convert_logintensityBase(model_data$quantobject_mscalib, new_base=2)
    #k <- log(logbase) / log(2)
    #model_data$quantobj_labu_min <- model_data$quantobj_labu_min * k
  }

  model_data <- prepare_msdata(model_data, msdata, max_quantobjects=max_quantobjects,
                               verbose=verbose, ...)
  model_data <- prepare_expanded_effects(model_data, verbose=verbose)

  return (structure(model_data, class="msglm_model_data"))
}
