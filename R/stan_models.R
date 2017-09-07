require(rstan)

stan_models_path <- file.path(base_scripts_path, "msglm/inst/stan_models") # FIXME
msglm.stan_model <- stan_model(file.path(stan_models_path, "msglm.stan"), "msglm", save_dso = TRUE, auto_write = TRUE)
msglm_local.stan_model <- stan_model(file.path(stan_models_path, "msglm_local.stan"), "msglm_local", save_dso = TRUE, auto_write = TRUE)
msglm_normalize.stan_model <- stan_model(file.path(stan_models_path , "msglm_normalize.stan"), "msglm_normalize", save_dso = TRUE, auto_write = TRUE)

# variables description for msglm_local model
msglm.vars_info <- list(
  global = list( names = c(), dims=c() ),#obj_shift_sigma', 'obj_effect_tau'), dims = c() ),
  effects = list( names = c('effect_shift_replCI_sigma'), dims = c('effect') ),
  #batch_effects = list(names = c('obj_batch_effect_sigma'), dims = c('batch_effect')),
  conditions = list(names = c('obj_repl_shift_sigma'), dims = c("condition")),
  iactions = list(names = c('iaction_labu', 'iaction_labu_replCI'), dims = c('iaction')),
  observations = list(names = c('obs_labu'), dims = c('observation')),
  objects = list(names = c('obj_base_labu'), dims = c('object')),
  object_effects = list(names = c('obj_effect_lambda', 'obj_effect', 'obj_effect_replCI'),
                        dims = c('object_effect')),
  object_repl_shifts = list(names = c('objXexp_repl_shift_unscaled'), dims = c('object', 'msrun')),
  object_batch_effects = list(names = c('obj_batch_effect_unscaled'), #'protein_batch_effect_lambda',
                              dims = c('object_batch_effect')),
  object_batch_shifts = list(names = c('objXexp_batch_shift'), dims = c('object', 'msrun'))
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
                              repl_tau=0.25, batch_tau=0.3)
{
  message('Converting MSGLM model data to Stan-readable format...')
  res <- list(
    Niactions = nrow(model_data$interactions),
    Nobservations = nrow(model_data$ms_data),
    Nexperiments = n_distinct(model_data$mschannels$msrun_ix),
    Nconditions = nrow(conditionXeffect.mtx),
    Nobjects = n_distinct(model_data$ms_data$glm_object_ix),
    experiment2condition = as.array(model_data$mschannels$condition_ix),
    effectXcondition = t(conditionXeffect.mtx),
    inv_effectXcondition = t(inv_conditionXeffect.mtx),
    experiment_shift = as.array(model_data$mschannels$model_mschannel_shift),
    batchEffectXexperiment = t(msrunXbatchEffect.mtx),
    replEffectXexperiment = t(msrunXreplEffect.mtx),
    observation2experiment = as.array(model_data$ms_data$msrun_ix),
    observation2iaction = as.array(model_data$ms_data$glm_iaction_ix),
    iaction2obj = as.array(model_data$interactions$glm_object_ix),
    iaction2condition = as.array(model_data$interactions$condition_ix),
    Neffects = ncol(conditionXeffect.mtx),
    NobjEffects = nrow(model_data$object_effects),
    obj_effect_is_positive = model_data$object_effects$is_positive,
    NeffectsPerObjCumsum = as.array(nrows_cumsum(model_data$object_effects, 'glm_object_ix')),
    obj_effect2obj = as.array(as.integer(model_data$object_effects$glm_object_ix)),
    obj_effect2effect = as.array(as.integer(model_data$object_effects$effect)),
    NreplEffects = ncol(msrunXreplEffect.mtx),
    NobjReplEffects = nrow(model_data$object_repl_effects),
    NreplEffectsPerObjCumsum = as.array(nrows_cumsum(model_data$object_repl_effects, 'glm_object_ix')),
    obj_repl_effect2repl_effect = as.array(as.integer(model_data$object_repl_effects$replicate_effect)),
    NbatchEffects = ncol(msrunXbatchEffect.mtx),
    NobjBatchEffects = nrow(model_data$object_batch_effects),
    NbatchEffectsPerObjCumsum = as.array(nrows_cumsum(model_data$object_batch_effects, 'glm_object_ix')),
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
    obj_repl_shift_tau = repl_tau,
    obj_batch_effect_tau = batch_tau,
    obj_base_labu_sigma = 4.0,
    #obj_batch_effect_sigma = 0.25,
    underdef_obj_shift = -8.0#,
    #zShift = mean(log(ms_data$protgroup_intensities$intensity), na.rm = TRUE),
    #zScale = 1.0/sd(log(ms_data$protgroup_intensities$intensity), na.rm = TRUE)
  ) %>% c(base_input_data)
  message(res$Niactions, " interaction(s), ",
          res$Nquanted, " quantitation(s), ",
          res$Nmissed, " missed")
  return(res)
}

stan.sampling <- function(stan_input_data, iter=4000, chains=8, thin=4, adapt_delta=0.9)
{
    message("Running Stan MCMC...")
    sampling(msglm_local.stan_model,
             pars=unlist(lapply(msglm.vars_info, function(vi) vi$names)), include=TRUE,
             data = stan_input_data,
             #init = function() { pcp_peaks_glm.generate_init_params(pcp_peaks_glm.model_data) },
             iter = iter, chains = chains, thin = thin,
             control = list(adapt_delta=adapt_delta, max_treedepth=11))
}
