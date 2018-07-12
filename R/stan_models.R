require(rstan)
require(Matrix)

stan_models_path <- file.path(base_scripts_path, "R/msglm/inst/stan_models") # FIXME
msglm.stan_model <- stan_model(file.path(stan_models_path, "msglm.stan"), "msglm", save_dso = TRUE, auto_write = TRUE)
msglm_local.stan_model <- stan_model(file.path(stan_models_path, "msglm_local.stan"), "msglm_local", save_dso = TRUE, auto_write = TRUE)
msglm_local_subobjects.stan_model <- stan_model(file.path(stan_models_path, "msglm_local_subobjects.stan"), "msglm_local_subobjects", save_dso = TRUE, auto_write = TRUE)
msglm_normalize.stan_model <- stan_model(file.path(stan_models_path , "msglm_normalize.stan"), "msglm_normalize", save_dso = TRUE, auto_write = TRUE)

# variables description for msglm_local model
msglm.vars_info <- list(
  global = list( names = c('suo_shift_sigma', 'suo_msproto_shift_sigma'), dims=c() ),#obj_shift_sigma', 'obj_effect_tau'), dims = c() ),
  #batch_effects = list(names = c('batch_effect_sigma'), dims = c('batch_effect')),
  iactions = list(names = c('iact_repl_shift_sigma', 'iaction_labu', 'iaction_labu_replCI'), dims = c('iaction')),
  observations = list(names = c('obs_labu', "obs_repl_shift"), dims = c('observation')),
  subobjects = list(names = c('suo_shift_unscaled'), dims = c('subobject')),
  subobjectXmsprotocol = list(names = c('suo_msproto_shift_unscaled'), dims=c('subobjectXmsprotocol1')),
  objects = list(names = c('obj_base_labu', 'obj_base_labu_replCI', "obj_base_repl_shift_sigma"), dims = c('object')),
  object_effects = list(names = c('obj_effect_sigma', 'obj_effect_repl_shift_sigma', 'obj_effect', 'obj_effect_replCI'),
                        dims = c('object_effect')),
  object_batch_effects = list(names = c('obj_batch_effect'), #'obj_batch_effect_unscaled', 'obj_batch_effect_sigma',
                              dims = c('object_batch_effect'))
)

# get the indices of the first rows in a group
# also get the index+1 of the last row
nrows_cumsum <- function(df, group_col) {
  n_rows <- df %>% dplyr::arrange_(group_col) %>% dplyr::group_by_(group_col) %>%
    dplyr::summarize(n_rows = n()) %>% .$n_rows
  if (length(n_rows)>0) {
    res <- cumsum(as.integer(n_rows))
    c(1L, res+1L)
  } else {
    return(0L)
  }
}

stan.prepare_data <- function(base_input_data, model_data,
                              global_labu_shift = global_protgroup_labu_shift,
                              base_repl_shift_tau=0.1, effect_repl_shift_tau=0.25, batch_tau=0.3)
{
  message('Converting MSGLM model data to Stan-readable format...')
  if (any(as.integer(model_data$effects$effect) != seq_len(nrow(model_data$effects)))) {
    stop("model_data$effects are not ordered")
  }
  if (any(as.integer(model_data$batch_effects$batch_effect) != seq_len(nrow(model_data$batch_effects)))) {
    stop("model_data$batch_effects are not ordered")
  }
  obs_df <- dplyr::select(model_data$ms_data, glm_observation_ix, msrun_ix, glm_object_ix, glm_iaction_ix) %>%
    dplyr::distinct()
  if (any(obs_df$glm_observation_ix != seq_len(nrow(obs_df)))) {
    stop("model_data$ms_data not ordered by observations / have missing observations")
  }
  res <- base_input_data
  res <- c(res, list(
    Niactions = nrow(model_data$interactions),
    Nobservations = nrow(obs_df),
    Nexperiments = n_distinct(model_data$mschannels$msrun_ix),
    Nconditions = nrow(conditionXeffect.mtx),
    Nobjects = n_distinct(model_data$ms_data$glm_object_ix),
    experiment_shift = as.array(model_data$mschannels$model_mschannel_shift),
    observation2experiment = as.array(obs_df$msrun_ix),
    observation2iaction = as.array(obs_df$glm_iaction_ix),
    iaction2obj = as.array(model_data$interactions$glm_object_ix),
    iaction2condition = as.array(model_data$interactions$condition_ix),
    Neffects = ncol(conditionXeffect.mtx),
    effect_is_positive = as.array(as.integer(model_data$effects$is_positive)),
    NobjEffects = nrow(model_data$object_effects),
    obj_effect2effect = as.array(as.integer(model_data$object_effects$effect)),
    NbatchEffects = ncol(msrunXbatchEffect.mtx),
    batch_effect_is_positive = as.array(as.integer(model_data$batch_effects$is_positive)),
    NobjBatchEffects = nrow(model_data$object_batch_effects),
    obj_batch_effect2batch_effect = as.array(as.integer(model_data$object_batch_effects$batch_effect)),
    NunderdefObjs = sum(model_data$objects$is_underdefined),
    underdef_objs = as.array(dplyr::filter(model_data$objects, is_underdefined) %>% .$glm_object_ix),
    Nquanted = sum(!is.na(model_data$ms_data$intensity)),
    Nmissed = sum(is.na(model_data$ms_data$intensity)),
    quant2observation = as.array(model_data$ms_data$glm_observation_ix[!is.na(model_data$ms_data$qdata_ix)]),
    miss2observation = as.array(model_data$ms_data$glm_observation_ix[!is.na(model_data$ms_data$mdata_ix)]),
    qData = as.array(model_data$ms_data$intensity[!is.na(model_data$ms_data$intensity)]),
    global_labu_shift = global_labu_shift,
    effect_tau = effects.df$tau,
    obj_base_repl_shift_tau = base_repl_shift_tau,
    obj_effect_repl_shift_tau = effect_repl_shift_tau,
    obj_batch_effect_tau = batch_tau,
    obj_base_labu_sigma = 4.0,
    #obj_batch_effect_sigma = 0.25,
    underdef_obj_shift = -8.0#,
    #zShift = mean(log(ms_data$protgroup_intensities$intensity), na.rm = TRUE),
    #zScale = 1.0/sd(log(ms_data$protgroup_intensities$intensity), na.rm = TRUE)
  )) %>%
    append_sparse("iactXobjeff", model_data$iactXobjeff) %>%
    append_sparse("obsXobjbatcheff", model_data$obsXobjbatcheff)

  if ("subobjects" %in% names(model_data)) {
    # data have subobjects
    res$Nsubobjects <- nrow(model_data$subobjects)
    res$suo2obj <- as.array(as.integer(model_data$subobjects$glm_object_ix))
    res$quant2suo <- as.array(as.integer(model_data$ms_data$glm_subobject_ix[!is.na(model_data$ms_data$qdata_ix)]))
    res$miss2suo <- as.array(as.integer(model_data$ms_data$glm_subobject_ix[!is.na(model_data$ms_data$mdata_ix)]))
    res$Nmsprotocols <- 0L
    res$experiment2msproto <- integer(0)
  }
  if ('msproto_ix' %in% names(model_data$mschannels)) {
    res$Nmsprotocols <- n_distinct(model_data$mschannels$msproto_ix)
    res$experiment2msproto <- as.array(model_data$mschannels$msproto_ix)
  }
  if ("Nsubobjects" %in% names(res)) {
    message(res$Niactions, " interaction(s) of ", res$Nobjects, " object(s) with ",
            res$Nsubobjects, " subobject(s), ",
            res$Nquanted, " quantitation(s), ",
            res$Nmissed, " missed")
  } else {
    message(res$Niactions, " interaction(s) of ", res$Nobjects, " object(s), ",
            res$Nquanted, " quantitation(s), ",
            res$Nmissed, " missed")
  }
  return(res)
}

stan.sampling <- function(stan_input_data, iter=4000, chains=8, thin=4,
                          adapt_delta=0.9, max_treedepth=11)
{
    message("Running Stan MCMC...")
    vars_info <- msglm.vars_info
    if ("Nsubobjects" %in% names(stan_input_data)) {
      stanmodel <- msglm_local_subobjects.stan_model
    } else {
      stanmodel <- msglm_local.stan_model
      # exclude subobject-related
      vars_info$subobjects <- NULL
      vars_info$subobjectXmsprotocol <- NULL
      vars_info$global$names <- setdiff(vars_info$global$names,
                                        c('suo_shift_sigma', 'suo_msproto_shift_sigma'))
    }
    res <- sampling(stanmodel,
             pars=unlist(lapply(vars_info, function(vi) vi$names)), include=TRUE,
             data = stan_input_data,
             #init = function() { pcp_peaks_glm.generate_init_params(pcp_peaks_glm.model_data) },
             iter = iter, chains = chains, thin = thin,
             control = list(adapt_delta=adapt_delta, max_treedepth=max_treedepth))
    attr(res, "msglm_vars_info") <- vars_info
    return(res)
}
