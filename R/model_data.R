# prepare observations array
prepare_observations <- function(data.mtx)
{
  data.vec <- as.vector(data.mtx)
  data.vec[data.vec == 0] <- NA
  data.obs.mask <- !is.na(data.vec)
  data.idx <- ifelse(data.obs.mask, cumsum(data.obs.mask), -cumsum(!data.obs.mask))
  data.idx.arr <- matrix(data.idx, ncol = ncol(data.mtx))
  list(Nobserved = sum(data.obs.mask), # total non-NA  datapoints
       oData = array(data.vec[data.obs.mask]), # only non-NA values
       dataIndex = as.array(data.idx.arr) # either index in oData if observed, or -index of non-observed
  )
}

# given a experimentXeffect matrix,
# and mapping of the data to objects and experiments, construct
# a data.frame of all related objectXeffect combinations
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
      bind_rows(lapply(iact_infos_masked, function(iact_info) iact_info$df))
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
prepare_effects <- function(model_data, underdefined_iactions=FALSE)
{
  iactXobjeff <- iactXeffect(conditionXeffect.mtx, model_data$interactions$glm_object_ix,
                             model_data$interactions$condition_ix)
  model_data$object_effects <- iactXobjeff$objeff_df %>%
    dplyr::rename(glm_object_ix = obj,
                  effect = eff, object_effect = objeff) %>%
    dplyr::arrange(object_effect)
  model_data$iactXobjeff <- iactXobjeff$mtx
  dimnames(model_data$iactXobjeff) <- list(interaction = model_data$interactions$iaction_id,
                                           objeff = dimnames(model_data$iactXobjeff)[[2]])
  model_data$effects <- dplyr::select(effects.df, effect, is_positive) %>%
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
