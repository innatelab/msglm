# given a experimentXeffect matrix,
# and mapping of the data to objects and experiments, construct
# a data.frame of all related objectXeffect combinations
#' @export
iactXeffect <- function(expXeff, iact2obj, iact2exp) {
  iact_infos <- lapply(seq_along(iact2obj), function(iact_ix){
    iactXeff <- expXeff[iact2exp[[iact_ix]], , drop=FALSE]
    eff_mask <- colSums(abs(iactXeff), na.rm=TRUE) != 0.0
    if (!any(eff_mask)) return (NULL)
    list(mtx = iactXeff[, eff_mask, drop=FALSE],
         df = tibble(eff = colnames(expXeff)[eff_mask],
                     iact_ix = iact_ix))
  })
  objeffs <- unique(unlist(lapply(seq_along(iact_infos),
                           function(iact_ix) if (!is.null(iact_infos[[iact_ix]])) paste0(colnames(iact_infos[[iact_ix]]$mtx), "@", iact2obj[[iact_ix]]) else c())))
  iact_mask <- !sapply(iact_infos, is.null)
  iact_infos_masked <- iact_infos[iact_mask]
  df <- if (any(iact_mask)) {
      dplyr::bind_rows(lapply(iact_infos_masked, function(iact_info) iact_info$df))
  } else {
      tibble(eff = character(),
             iact_ix = integer())
  }
  df <- df %>%
    dplyr::mutate(obj = iact2obj[iact_ix],
                  objeff = paste0(eff, '@', obj))
  objeff_df <- dplyr::select(df, obj, eff, objeff) %>%
    dplyr::distinct() %>%
    dplyr::mutate(eff = factor(eff, levels=intersect(colnames(expXeff), unique(df$eff)))) %>%
    dplyr::arrange(obj, eff) %>%
    dplyr::mutate(objeff = factor(objeff, levels=objeff))
  df <- dplyr::mutate(df, objeff = factor(objeff, levels = levels(objeff_df$objeff)))
  mtx <- matrix(0.0, ncol = nrow(objeff_df), nrow = length(iact2obj),
                dimnames = list(iact = NULL, objeff = levels(objeff_df$objeff)))
  for (iact_ix in which(iact_mask)) {
    iact_mtx <- iact_infos[[iact_ix]]$mtx
    mtx[iact_ix, paste0(colnames(iact_mtx), "@", iact2obj[iact_ix])] <- iact_mtx
  }
  return(list(mtx = mtx, df = df, objeff_df = objeff_df))
}

# use experimental design matrices and add
# object/replicate/batch effects information to model_data
#' @export
prepare_effects <- function(model_data, underdefined_iactions=FALSE)
{
  is_glmm <- exists("supXcond.mtxs") && exists("mixeffects.df")
  if (is_glmm) {
    message("Detected mixeffects data for GLMM model")
    # "GLMM" mixing model
    if (any(names(supXcond.mtxs) != rownames(mixcoefXeff.mtx))) {
      stop("Mismatch between supXcond.mtxs matrix names and mix coefficient names")
    }
    if (any(mixeffects.df$mixeffect != colnames(mixcoefXeff.mtx))) {
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
    model_data$mixeffects <- mixeffects.df
    model_data$mixcoefs <- tibble(mixcoef = rownames(mixcoefXeff.mtx))
    model_data$mixcoefXeff <- mixcoefXeff.mtx
    model_data$supXcond <- supXcond.mtxs
    conditions <- levels(model_data$conditions$condition)
    sactXiact.mtxs <- lapply(model_data$supXcond, function(supXcond.mtx){
      res <- iactXeffect(supXcond.mtx,
                         model_data$superactions$glm_object_ix,
                         model_data$superactions$supcondition_ix)
      names(dimnames(res$mtx)) <- c("sact", "iact")
      res$mtx <- res$mtx[, colSums(abs(res$mtx)) != 0.0, drop=FALSE]
      res$df <- dplyr::rename(res$df, condition=eff, iaction_id=objeff)
      res$iaction_df <- dplyr::rename(res$objeff_df, condition=eff,
                                      iaction_id=objeff,
                                      glm_object_ix=obj)
      res$objeff_df <- NULL
      return(res)
    })
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
      dplyr::mutate(condition_ix = as.integer(condition)) %>%
      dplyr::mutate(glm_iaction_ix = row_number(),
             is_virtual = FALSE)
    model_data$mixtions <- dplyr::bind_rows(lapply(names(sactXiact.mtxs), function(mixcoef){
      tibble::tibble(mixcoef = mixcoef,
                     iaction_id = factor(colnames(sactXiact.mtxs[[mixcoef]]$mtx), levels=iaction_ids))
    })) %>%
      dplyr::mutate(mixtion_ix = row_number(),
                    mixtion = paste0(mixcoef, " ", iaction_id),
                    mixcoef = factor(mixcoef, levels=rownames(mixcoefXeff.mtx)),
                    mixcoef_ix = replace_na(as.integer(mixcoef), 0L)) %>%
    dplyr::left_join(dplyr::select(model_data$interactions, iaction_id, glm_iaction_ix)) %>%
    dplyr::arrange(mixtion_ix)
    model_data$supactXmixt <- do.call(cbind, lapply(sactXiact.mtxs, function(x) x$mtx))
    dimnames(model_data$supactXmixt) <- list(superaction = model_data$superactions$supaction_id,
                                             mixtion = model_data$mixtions$mixtion)
  } else {
    iaction_ids <- unique(model_data$interaction$iaction_id)
  }

  iactXobjeff <- iactXeffect(conditionXeffect.mtx, model_data$interactions$glm_object_ix,
                             model_data$interactions$condition_ix)
  model_data$object_effects <- iactXobjeff$objeff_df %>%
    dplyr::rename(glm_object_ix = obj,
                  effect = eff, object_effect = objeff) %>%
    dplyr::arrange(object_effect)
  model_data$iactXobjeff <- iactXobjeff$mtx
  dimnames(model_data$iactXobjeff) <- list(interaction = iaction_ids,
                                           objeff = dimnames(model_data$iactXobjeff)[[2]])
  model_data$effects <- effects.df %>%
      dplyr::mutate(effect = factor(effect, levels=levels(model_data$object_effects$effect))) %>%
      dplyr::arrange(effect)

  obsXobjbatcheff <- iactXeffect(msrunXbatchEffect.mtx[unique(model_data$mschannels$msrun), , drop=FALSE],
                                 model_data$observations$glm_object_ix,
                                 model_data$observations$msrun_ix)
  model_data$object_batch_effects <- obsXobjbatcheff$objeff_df %>%
    dplyr::rename(glm_object_ix = obj,
                  batch_effect = eff,
                  object_batch_effect = objeff) %>%
    dplyr::arrange(object_batch_effect)
  model_data$batch_effects <- dplyr::select(batch_effects.df, batch_effect, is_positive) %>%
      dplyr::mutate(batch_effect = factor(batch_effect, levels=levels(model_data$object_batch_effects$batch_effect))) %>%
      dplyr::arrange(batch_effect)
  model_data$obsXobjbatcheff <- obsXobjbatcheff$mtx

  if (underdefined_iactions) {
    # detect proteins that have no quantifications for estimating object_shift
    oe2iact.df <- model_data$object_effects %>%
      dplyr::inner_join(dplyr::mutate(conditionXeffect.df, effect = factor(effect, levels=levels(model_data$object_effects$effect)))) %>%
      dplyr::inner_join(model_data$interactions)
    effectless_iactions.df <- model_data$interactions %>% dplyr::anti_join(oe2iact.df)
    underdefined_objects.df <- dplyr::anti_join(model_data$objects,
                                                effectless_iactions.df %>% dplyr::filter(!is_virtual)) %>%
      dplyr::select(glm_object_ix)

    model_data$objects <- dplyr::mutate(model_data$objects,
                                        is_underdefined = glm_object_ix %in% underdefined_objects.df$glm_object_ix)
  } else {
    model_data$objects <- dplyr::mutate(model_data$objects,
                                        is_underdefined = FALSE)
  }

  return(model_data)
}
