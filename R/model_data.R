# prepare observations array
prepare_observations <- function(data.mtx)
{
  data.vec <- as.vector(data.mtx)
  data.vec[ data.vec == 0 ] <- NA
  data.obs.mask <- !is.na( data.vec )
  data.idx <- ifelse( data.obs.mask, cumsum( data.obs.mask ), -cumsum( !data.obs.mask ) )
  data.idx.arr <- matrix( data.idx, ncol = ncol( data.mtx ) )
  list( Nobserved = sum(data.obs.mask), # total non-NA  datapoints
        oData = array( data.vec[ data.obs.mask ] ), # only non-NA values
        dataIndex = as.array(data.idx.arr) # either index in oData if observed, or -index of non-observed
  )
}

# given a experimentXeffect matrix,
# and mapping of the data to objects and experiments, construct
# a data.frame of all related objectXeffect combinations
obj2effect <- function(expXeff, obj_class, objs, data2obj, data2exp) {
  effXexp <- t(expXeff)
  expXeff_mask <- apply(effXexp, c(2, 1), function(x) x != 0.0)
  objXexp_mask <- matrix(FALSE, nrow = length(objs), ncol = ncol(effXexp))
  if (length(data2obj)>0) { for (i in 1:length(data2obj)) {
    objXexp_mask[data2obj[[i]], data2exp[[i]]] <- TRUE
  } }
  objXeff_mask <- objXexp_mask %*% expXeff_mask
  dnams <- list(objs, dimnames(effXexp)[[1]])
  names(dnams) <- c(obj_class, names(dimnames(effXexp))[[1]])
  dimnames(objXeff_mask) <- dnams
  if (length(objXeff_mask) > 0L) {
    res <- as.data.frame(as.table(t(objXeff_mask))) %>% dplyr::filter(Freq != 0) %>% dplyr::select(-Freq)
  } else {
    res <- data.frame() %>% dplyr::mutate(x = character(), y = character())
    names(res) <- c(obj_class, names(dnams)[[2]])
  }
  res[,paste(colnames(res), collapse='_')] <- seq_len(nrow(res))
  res
}

# use experimental design matrices and add
# object/replicate/batch effects information to model_data
prepare_effects <- function(model_data, underdefined_iactions=FALSE)
{
  model_data$object_effects <- obj2effect(conditionXeffect.mtx, "glm_object_ix", model_data$objects$glm_object_ix,
                                          model_data$interactions$glm_object_ix,
                                          model_data$interactions$condition_ix) %>%
    dplyr::mutate(glm_object_ix = as.integer(glm_object_ix))
  model_data$effects <- dplyr::select(effects.df, effect, is_positive) %>%
      dplyr::mutate(effect = factor(effect, levels=levels(model_data$object_effects$effect))) %>%
      dplyr::arrange(effect)

  model_data$object_repl_effects <- obj2effect(msrunXreplEffect.mtx[model_data$mschannels$msrun, , drop=FALSE], "glm_object_ix",
                                               model_data$objects$glm_object_ix,
                                               model_data$ms_data$glm_object_ix,
                                               model_data$ms_data$msrun_ix) %>%
    dplyr::mutate(glm_object_ix = as.integer(glm_object_ix))
  
  model_data$object_batch_effects <- obj2effect(msrunXbatchEffect.mtx[model_data$mschannels$msrun, , drop=FALSE], "glm_object_ix",
                                                model_data$objects$glm_object_ix,
                                                model_data$ms_data$glm_object_ix,
                                                model_data$ms_data$msrun_ix) %>%
    dplyr::mutate(glm_object_ix = as.integer(glm_object_ix))
  
  if (underdefined_iactions) {
    # detect proteins that have no quantifications for estimating object_shift
    oe2iact.df <- model_data$object_effects %>%
      dplyr::inner_join(conditionXeffect.df) %>%
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
