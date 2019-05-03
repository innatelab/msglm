# prepare observations array
#' @export
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
  is_glmm <- exists("mix_condXact.mtxs") && exists("mix_effects.df")
  if (is_glmm) {
    message("Detected mix_effects data for GLMM model")
    # "GLMM" mixing model
    if (any(names(mix_condXact.mtxs) != c("(Intercept)", as.character(mix_effects.df$mix_effect)))) {
      stop("Mismatch between mix_condXact.mtxs matrix names and mix effect names")
    }
    model_data$mix_effects <- mix_effects.df
    model_data$mix_condXact <- mix_condXact.mtxs
    sactXiact.mtxs <- lapply(model_data$mix_condXact, function(condXact.mtx){
      res <- iactXeffect(condXact.mtx,
                         model_data$superactions$glm_object_ix,
                         model_data$superactions$condition_ix)
      dimnames(res$mtx) <- list(sact=dimnames(res$mtx)[[1]],
                                iact=dimnames(res$mtx)[[2]])
      res$mtx <- res$mtx[, colSums(abs(res$mtx)) != 0.0, drop=FALSE]
      res$df <- dplyr::rename(res$df, act=eff, iact=objeff)
      res$iaction_df <- dplyr::rename(res$objeff_df, action=eff,
                                      iaction_id=objeff,
                                      glm_object_ix=obj) %>%
        dplyr::mutate(condition = action)
      res$objeff_df <- NULL
      return(res)
    })
    model_data$interactions <- dplyr::distinct(dplyr::bind_rows(lapply(sactXiact.mtxs, function(x) x$iaction_df))) %>%
      dplyr::mutate(action_ix = factor(action, levels = colnames(mix_condXact.mtxs[[1]]))) %>%
      dplyr::left_join(dplyr::select(model_data$conditions, condition, condition_ix)) %>%
      dplyr::mutate(glm_iaction_ix = row_number(),
             is_virtual = FALSE)
    model_data$mixtions <- dplyr::bind_rows(lapply(names(sactXiact.mtxs), function(mix_eff){
      tibble::tibble(mix_effect = mix_eff,
                     iaction_id = colnames(sactXiact.mtxs[[mix_eff]]$mtx))
    })) %>%
      dplyr::mutate(mixtion_ix = row_number(),
                    mixtion = paste0(mix_effect, " ", iaction_id),
                    mix_effect = factor(mix_effect, levels=levels(mix_effects.df$mix_effect)),
                    mix_effect_ix = coalesce(as.integer(mix_effect), 0L)) %>%
    dplyr::left_join(dplyr::select(model_data$interactions, iaction_id, glm_iaction_ix)) %>%
    dplyr::arrange(mixtion_ix)
    model_data$supactXmixt <- do.call(cbind, lapply(sactXiact.mtxs, function(x) x$mtx))
    dimnames(model_data$supactXmixt) <- list(superaction = model_data$superactions$supaction_id,
                                             mixtion = model_data$mixtions$mixtion)
  }

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
