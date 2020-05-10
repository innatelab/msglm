# variables description for msglm_local model
msglm.vars_info <- list(
  global = list(names=c('suo_shift_sigma', 'effect_slab_c', 'batch_effect_slab_c'),
                dims=c()),#obj_shift_sigma', 'obj_effect_tau'), dims = c() ),
  #batch_effects = list(names = c('batch_effect_sigma'),
  #                     dims = c('batch_effect')),
  iactions = list(names=c('iact_repl_shift_sigma', 'iaction_labu', 'iaction_labu_replCI'),
                  dims=c('iaction')),
  observations = list(names=c('obs_labu', "obs_repl_shift", "obs_batch_shift"),
                      dims=c('observation')),
  subobjects = list(names=c('suo_shift', 'suo_llh'),
                    dims=c('subobject')),
  #subobject_subbatch_shifts = list(names=c('suoxobs_subbatch_shift'),
  #                                 dims=c('subobjectXobservation')),
  subobject_subbatch_effects = list(names=c('suo_subbatch_effect', 'suo_subbatch_effect_sigma'),
                                    dims=c('subobject_subbatch_effect')),
  objects = list(names=c('obj_base_labu', 'obj_base_labu_replCI', "obj_base_repl_shift_sigma"),
                 dims=c('object')),
  object_effects = list(names=c('obj_effect_sigma', 'obj_effect_repl_shift_sigma', 'obj_effect', 'obj_effect_replCI'),
                        dims=c('object_effect')),
  object_batch_effects = list(names=c('obj_batch_effect', 'obj_batch_effect_sigma'), #'obj_batch_effect_unscaled',
                              dims=c('object_batch_effect'))
)

msglmm.vars_info <- msglm.vars_info
msglmm.vars_info$object_mixeffects <- list(names = c("obj_mixeffect", "obj_mixeffect_sigma"),
                                            dims = c("object_mixeffect"))
msglmm.vars_info$object_mixcoefs <- list(names = c("obj_mixcoef"),
                                         dims = c("object_mixcoef"))
msglmm.vars_info$supactions <- list(names = c("supaction_labu", "supact_repl_shift_sigma"),
                                    dims = c("supaction"))

# get the indices of the first rows in a group
# also get the index+1 of the last row
nrows_cumsum <- function(df, group_col) {
  n_rows <- df %>% dplyr::arrange_at(group_col) %>% dplyr::group_by_at(group_col) %>%
    dplyr::summarize(n_rows = n()) %>% .$n_rows
  if (length(n_rows)>0) {
    res <- cumsum(as.integer(n_rows))
    c(1L, res+1L)
  } else {
    return(0L)
  }
}

#' @export
stan.prepare_data <- function(base_input_data, model_data,
                              global_labu_shift = global_protgroup_labu_shift,
                              effect_slab_df = 4, effect_slab_scale = 2.5,
                              obj_labu_min = -10, obj_labu_min_scale = 1,
                              base_repl_shift_tau=0.1, effect_repl_shift_tau=0.25,
                              batch_tau=0.3, batch_slab_df = 4, batch_slab_scale = 2.5, batch_df = 1.0, batch_df2 = 4.0,
                              subbatch_tau=batch_tau)
{
  message('Converting MSGLM model data to Stan-readable format...')
  is_glmm <- "mixeffects" %in% names(model_data)
  if (any(as.integer(model_data$effects$effect) != seq_len(nrow(model_data$effects)))) {
    stop("model_data$effects are not ordered")
  }
  if (any(as.integer(model_data$batch_effects$batch_effect) != seq_len(nrow(model_data$batch_effects)))) {
    stop("model_data$batch_effects are not ordered")
  }
  if (!("mschannel_ix" %in% names(model_data$mschannels))) {
    warn("No mschannel_ix column found, using msrun_ix")
    model_data$mschannels$mschannel_ix <- model_data$mschannels$msrun_ix
    model_data$msdata$mschannel_ix <- model_data$msdata$msrun_ix
  }
  xaction_ix_col <- if (is_glmm) "glm_supaction_ix" else "glm_iaction_ix"
  obs_df <- dplyr::select(model_data$observations, glm_observation_ix, mschannel_ix, msrun_ix, glm_object_ix, !!xaction_ix_col) %>%
    dplyr::distinct()
  if (any(obs_df$glm_observation_ix != seq_len(nrow(obs_df)))) {
    stop("model_data$msdata not ordered by observations / have missing observations")
  }
  if (is_glmm && any(as.integer(model_data$mixeffects$mixeffect) !=
                     seq_len(nrow(model_data$mixeffects)))) {
    stop("model_data$mixeffects are not ordered")
  }
  model_data$effects <- maybe_rename(model_data$effects, c("prior_mean" = "mean", "prior_tau" = "tau"))
  if (is_glmm) {
    model_data$mixeffects <- maybe_rename(model_data$mixeffects, c("prior_mean" = "mean", "prior_tau" = "tau"))
  }
  res <- base_input_data
  res <- c(res, list(
    Nobservations = nrow(obs_df),
    Nexperiments = n_distinct(model_data$mschannels$mschannel_ix),
    Nconditions = nrow(conditionXeffect.mtx),
    Nobjects = n_distinct(model_data$interactions$glm_object_ix),
    experiment_shift = as.array(model_data$mschannels$model_mschannel_shift),
    observation2experiment = as.array(obs_df$mschannel_ix),
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
    Nquanted = sum(!is.na(model_data$msdata$intensity)),
    Nmissed = sum(is.na(model_data$msdata$intensity)),
    quant2observation = as.array(model_data$msdata$glm_observation_ix[!is.na(model_data$msdata$qdata_ix)]),
    miss2observation = as.array(model_data$msdata$glm_observation_ix[!is.na(model_data$msdata$mdata_ix)]),
    qData = as.array(model_data$msdata$intensity[!is.na(model_data$msdata$intensity)]),
    global_labu_shift = global_labu_shift,
    effect_tau = effects.df$prior_tau,
    effect_mean = effects.df$prior_mean,
    effect_df = if (rlang::has_name(effects.df, "prior_df")) {effects.df$prior_df} else {rep.int(1.0, nrow(effects.df))},
    effect_df2 = if (rlang::has_name(effects.df, "prior_df2")) {effects.df$prior_df2} else {rep.int(1.0, nrow(effects.df))},
    effect_slab_df = effect_slab_df,
    effect_slab_scale = effect_slab_scale,
    obj_labu_min = obj_labu_min, obj_labu_min_scale = obj_labu_min_scale,
    obj_base_labu_sigma = 15.0,
    obj_base_repl_shift_tau = base_repl_shift_tau,
    obj_effect_repl_shift_tau = effect_repl_shift_tau,
    obj_batch_effect_tau = batch_tau,
    batch_effect_df = batch_df, batch_effect_df2 = batch_df2,
    batch_effect_slab_df =  batch_slab_df, batch_effect_slab_scale = batch_slab_scale,
    #obj_batch_effect_sigma = 0.25,
    underdef_obj_shift = -8.0#,
    #zShift = mean(log(msdata$protgroup_intensities$intensity), na.rm = TRUE),
    #zScale = 1.0/sd(log(msdata$protgroup_intensities$intensity), na.rm = TRUE)
  )) %>%
    modifyList(matrix2csr("obsXobjbatcheff", model_data$obsXobjbatcheff))

  if (is_glmm) {
    message("Setting GLMM interaction data...")
    res$Nmix <- nrow(model_data$mixcoefXeff)
    res$NmixEffects <- nrow(model_data$mixeffects)
    res$mixeffect_mean <- as.array(model_data$mixeffects$prior_mean)
    res$mixeffect_tau <- as.array(model_data$mixeffects$prior_tau)
    res$mixcoefXeff <- as.matrix(model_data$mixcoefXeff)

    res <- modifyList(res, matrix2csr("iactXobjeff", model_data$iactXobjeff))
    iact_data <- list(Nsupactions = nrow(model_data$superactions),
                      observation2supaction = as.array(obs_df$glm_supaction_ix),
                      supaction2obj = as.array(model_data$superactions$glm_object_ix),
                      Niactions = nrow(model_data$interactions),
                      iaction2obj = as.array(model_data$interactions$glm_object_ix),
                      Nmixtions = nrow(model_data$mixtions),
                      mixt2iact = as.array(model_data$mixtions$glm_iaction_ix),
                      mixt2mix = as.array(model_data$mixtions$mixcoef_ix))
    res <- modifyList(res, iact_data)
    res <- modifyList(res, matrix2csr("supactXmixt", model_data$supactXmixt))
  } else {
    res <- modifyList(res, matrix2csr("iactXobjeff", model_data$iactXobjeff))
    iact_data <- list(Niactions = nrow(model_data$interactions),
                      observation2iaction = as.array(obs_df$glm_iaction_ix),
                      iaction2obj = as.array(model_data$interactions$glm_object_ix))
    res <- modifyList(res, iact_data)
  }
  if ("subobjects" %in% names(model_data)) {
    # data have subobjects
    res$Nsubobjects <- nrow(model_data$subobjects)
    res$suo2obj <- as.array(as.integer(model_data$subobjects$glm_object_ix))
    res$quant2suo <- as.array(as.integer(model_data$msdata$glm_subobject_ix[!is.na(model_data$msdata$qdata_ix)]))
    res$miss2suo <- as.array(as.integer(model_data$msdata$glm_subobject_ix[!is.na(model_data$msdata$mdata_ix)]))
    res$Nmsprotocols <- 0L
    res$experiment2msproto <- integer(0)
    # subobject-specific batch effects
    if ("suo_subbatch_effects" %in% names(model_data)) {
      res$NsubBatchEffects <- ncol(msrunXsubbatchEffect.mtx)
      res$NsuoBatchEffects <- nrow(model_data$suo_subbatch_effects)
      res$suo_subbatch_effect_tau <- subbatch_tau
      res$suo_subbatch_effect2subbatch_effect <- as.array(as.integer(model_data$suo_subbatch_effects$subbatch_effect))
      res$subbatch_effect_is_positive = as.array(as.integer(model_data$subbatch_effects$is_positive))
      res <- modifyList(res, matrix2csr("suoxobsXsuobatcheff", model_data$suoxobsXsuobatcheff))
    }
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

# all models
msglm_model_names <- function() {
  system.file('stan_models', package="msglm", mustWork=FALSE) %>%
  list.files(pattern="\\.stan$") %>% sapply(basename) %>% str_remove("\\.stan$")
}

msglm_stan_model <- function(model_name) {
  stan_models_path <- system.file('stan_models', package="msglm", mustWork=FALSE)
  message("Loading ", model_name, " Stan model")
  return (rstan::stan_model(file.path(stan_models_path, paste0(model_name, ".stan")),
                            model_name, save_dso = TRUE, auto_write = TRUE))
}

#' @export
stan.sampling <- function(stan_input_data, iter=4000, chains=8, thin=4,
                          max_treedepth=12L, stepsize_jitter=0.1, ...)
{
    message("Running Stan MCMC...")
    if ("Nsubobjects" %in% names(stan_input_data)) {
      if ("Nsupactions" %in% names(stan_input_data)) {
        vars_info <- msglmm.vars_info
        stanmodel <- msglm_stan_model("msglmm_local_subobjects")
      } else {
        vars_info <- msglm.vars_info
        stanmodel <- msglm_stan_model("msglm_local_subobjects")
      }
    } else {
      if ("Nsupactions" %in% names(stan_input_data)) {
        vars_info <- msglmm.vars_info
        stanmodel <- msglm_stan_model("msglmm_local")
      } else {
        vars_info <- msglm.vars_info
        stanmodel <- msglm_stan_model("msglm_local")
      }
      # exclude subobject-related
      vars_info$subobjects <- NULL
      vars_info$subobject_subbatch_effects <- NULL
      vars_info$global$names <- setdiff(vars_info$global$names,
                                        c('suo_shift_sigma'))
    }
    control_params <- list(..., max_treedepth=max_treedepth, stepsize_jitter=stepsize_jitter)
    res <- rstan::sampling(stanmodel,
             pars=unlist(lapply(vars_info, function(vi) vi$names)), include=TRUE,
             data = stan_input_data,
             #init = function() { pcp_peaks_glm.generate_init_params(pcp_peaks_glm.model_data) },
             iter = iter, chains = chains, thin = thin,
             control = control_params)
    attr(res, "msglm_vars_info") <- vars_info
    return(res)
}
