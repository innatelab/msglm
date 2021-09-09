# variables description for msglm_local model
msglm.vars_info <- list(
  global = list(names=c('subobj_shift_sigma', 'effect_slab_c'),# 'batch_effect_slab_c'),
                dims=c()),#obj_shift_sigma', 'obj_effect_tau'), dims = c() ),
  #batch_effects = list(names = c('batch_effect_sigma'),
  #                     dims = c('batch_effect')),
  interactions = list(names=c('iact_repl_shift_sigma', 'iaction_labu', 'iaction_labu_replCI'),
                      dims=c('interaction')),
  observations = list(names=c('obs_labu', "obs_repl_shift", "obs_batch_shift"),
                      dims=c('observation')),
  subobjects = list(names=c('subobj_shift', 'subobj_llh'),
                    dims=c('subobject')),
  #subobject_batch_shifts = list(names=c('subobs_batch_shift'),
  #                              dims=c('subobservations')),
  subobject_batch_effects = list(names=c('subobj_batch_effect', 'subobj_batch_effect_sigma'),
                                 dims=c('subobject_batch_effect')),
  objects = list(names=c('obj_base_labu', 'obj_base_labu_replCI'), #"obj_base_repl_shift_sigma"),
                 dims=c('object')),
  object_effects = list(names=c('obj_effect_sigma', 'obj_effect', 'obj_effect_replCI'),# 'obj_effect_repl_shift_sigma'),
                        dims=c('object_effect')),
  object_batch_effects = list(names=c('obj_batch_effect'),#, 'obj_batch_effect_sigma'), #'obj_batch_effect_unscaled',
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

#' Coverts the data and experimental design into Stan format.
#'
#' @param base_input_data input data already in Stan format (e.g. MS noise model parameters)
#' @param model_data the list with MS data an experimental design
#'
#' @param obj_labu_shift the average log-abundance of model objects
#' @param effect_slab_df the *degrees of freedom* for the prior of object effect *slab* regularization parameter
#' @param effect_slab_scale the *scale* parameter for the prior of object effect *slab* regularization parameter
#' @param obj_labu_min minimal object log-abundance
#' @param obj_labu_min_scale how strongly the estimates below *obj_labu_min* would be penalized (smaller = more stringent)
#' @param iact_repl_shift_tau the *tau* parameter for the prior of biochemical log-abundance variation between the replicates
#' @param iact_repl_shift_df the *degrees of freedom* parameter for the prior of biochemical log-abundance variation between the replicates
#' @param hsprior_lambda_a_offset offset from zero for the "a" parameter of horseshoe priors. Helps NUTS integration step
#' @param hsprior_lambda_t_offset offset from zero for the "a" parameter of horseshoe priors. Helps NUTS integration step
#' @param batch_effect_sigma sigma parameter for the normal distribution prior of batch effects
#' @param subbatch_tau
#' @param subbatch_df
#' @param subbatch_c
#' @param empty_observation_sigmoid_scale
#'
#' @export
stan.prepare_data <- function(model_data,
                              effect_slab_df = 4, effect_slab_scale = 2.5,
                              obj_labu_min_scale = 1,
                              iact_repl_shift_tau=0.03, iact_repl_shift_df=4.0,
                              hsprior_lambda_a_offset = 0.05,
                              hsprior_lambda_t_offset = 0.01,
                              batch_effect_sigma=0.5,
                              subbatch_tau=0.3, subbatch_df=4, subbatch_c=10,
                              empty_observation_sigmoid_scale = 1.0,
                              verbose = FALSE)
{
  model_def <- model_data$model_def
  if (verbose) message('Converting MSGLM model data to Stan-readable format...')
  is_glmm <- rlang::has_name(model_def, "mixeffects")
  ensure_primary_index_column(model_def$effects, 'index_effect')
  ensure_primary_index_column(model_data$mschannels, 'index_mschannel')

  xaction_ix_col <- if (is_glmm) "glm_supaction_ix" else "index_interaction"
  obs_df <- dplyr::select(model_data$observations, index_observation, index_mschannel, index_msrun, index_object, !!xaction_ix_col) %>%
    dplyr::distinct()
  if (any(obs_df$index_observation != seq_len(nrow(obs_df)))) {
    stop("model_data$msdata not ordered by observations / have missing observations")
  }
  if (is_glmm) {
    ensure_primary_index_column(model_def$mixeffects, 'index_mixeffect')
    # FIXME remove model_def$mixeffects <- maybe_rename(model_def$mixeffects, c("prior_mean" = "mean", "prior_tau" = "tau"))
  }
  msdata_obs_flags.df <- dplyr::group_by(model_data$msdata, index_observation) %>%
    dplyr::transmute(is_empty_observation = all(is.na(intensity))) %>%
    dplyr::ungroup()
  checkmate::assert_set_equal(msdata_obs_flags.df$index_observation,
                              model_data$msdata$index_observation, ordered=TRUE)

  missing_mask <- is.na(model_data$msdata$intensity)
  res <- list(
    Nobservations = nrow(obs_df),
    Nexperiments = n_distinct(model_data$mschannels$index_mschannel),
    Nconditions = nrow(model_def$conditions),
    Nobjects = nrow(model_data$objects),
    experiment_shift = as.array(model_data$mschannels$mschannel_shift),
    observation2experiment = as.array(obs_df$index_mschannel),

    Neffects = ncol(model_def$conditionXeffect),
    effect_is_positive = as.array(as.integer(model_def$effects$is_positive)),
    effect_tau = model_def$effects$prior_tau,
    effect_mean = model_def$effects$prior_mean,
    effect_df = model_def$effects$prior_df1,
    effect_df2 = model_def$effects$prior_df2,
    effect_slab_df = effect_slab_df,
    effect_slab_scale = effect_slab_scale,

    NobjEffects = nrow(model_data$object_effects),
    obj_effect2effect = as.array(model_data$object_effects$index_effect),

    NbatchEffects = n_distinct(model_data$object_batch_effects$index_batch_effect),
    batch_effect_is_positive = as.array(if (nrow(model_data$object_batch_effects) > 0)
                                        as.integer(model_def$batch_effects$is_positive)
                                        else integer(0)),
    batch_effect_sigma = batch_effect_sigma,

    NobjBatchEffects = nrow(model_data$object_batch_effects),
    obj_batch_effect2batch_effect = as.array(model_data$object_batch_effects$index_batch_effect),

    obj_labu_shift = model_data$quantobj_labu_shift,
    obj_base_labu_sigma = 15.0,
    obj_labu_min = rlang::set_names(model_data$quantobj_labu_min, NULL),
    obj_labu_min_scale = obj_labu_min_scale,

    Nquanted = sum(!missing_mask),
    Nmissed = sum(missing_mask),
    missing_sigmoid_scale = as.array(if_else(msdata_obs_flags.df$is_empty_observation[missing_mask],
                                             empty_observation_sigmoid_scale, 1.0)),
    qData = as.array(model_data$msdata$intensity[!missing_mask]),
    quant_isreliable = as.array(model_data$msdata$is_reliable[!missing_mask]),

    hsprior_lambda_a_offset = hsprior_lambda_a_offset,
    hsprior_lambda_t_offset = hsprior_lambda_t_offset,
    iact_repl_shift_tau = iact_repl_shift_tau, iact_repl_shift_df = iact_repl_shift_df
  ) %>%
    modifyList(model_data$quantobj_mscalib) %>%
    modifyList(matrix2csr("obsXobjbatcheff", model_data$obsXobjbatcheff))

  iact_data <- list(Niactions = nrow(model_data$interactions),
                    iaction2obj = as.array(model_data$interactions$index_object))
  iact_data <- modifyList(iact_data, matrix2csr("iactXobjeff", model_data$iactXobjeff))

  if (is_glmm) {
    message("Setting GLMM interaction data...")
    res$Nmix <- nrow(model_def$mixcoefXeff)
    res$NmixEffects <- nrow(model_def$mixeffects)
    res$mixeffect_mean <- as.array(model_def$mixeffects$prior_mean)
    res$mixeffect_tau <- as.array(model_def$mixeffects$prior_tau)
    res$mixcoefXeff <- as.matrix(model_def$mixcoefXeff)

    iact_data <- modifyList(iact_data,
                 list(Nsupactions = nrow(model_data$superactions),
                      observation2supaction = as.array(obs_df$glm_supaction_ix),
                      supaction2obj = as.array(model_data$superactions$index_object),
                      Nmixtions = nrow(model_data$mixtions),
                      mixt2iact = as.array(model_data$mixtions$index_interaction),
                      mixt2mix = as.array(model_data$mixtions$mixcoef_ix)))
    iact_data <- modifyList(iact_data, matrix2csr("supactXmixt", model_data$supactXmixt))
  } else {
    iact_data <- modifyList(iact_data, matrix2csr("obsXobjeff", model_data$obsXobjeff))
    iact_data$observation2iaction <- as.array(obs_df$index_interaction)
  }
  res <- modifyList(res, iact_data)

  if ("subobjects" %in% names(model_data)) {
    # data have subobjects
    subobj_data <- list(
      Nsubobjects = nrow(model_data$subobjects),
      subobj2obj = as.array(model_data$subobjects$index_object),
      Nsubobservations = nrow(model_data$msdata),
      quant2subobs = as.array(model_data$msdata$index_subobservation[!is.na(model_data$msdata$index_qdata)]),
      miss2subobs = as.array(model_data$msdata$index_subobservation[!is.na(model_data$msdata$index_mdata)]),
      subobs2subobj = as.array(model_data$msdata$index_subobject),
      subobs2obs = as.array(model_data$msdata$index_observation),
      # TODO support different noise models
      Nmsprotocols = 0L,
      experiment2msproto = integer(0),

      # subobject-specific batch effects
      NquantBatchEffects = n_distinct(model_data$subobject_batch_effects$index_quant_batch_effect),
      quant_batch_effect_tau = subbatch_tau,
      quant_batch_effect_df = subbatch_df,
      quant_batch_effect_c = subbatch_c,
      quant_batch_effect_is_positive = as.array(if (nrow(model_data$subobject_batch_effects)>0)
                                                as.integer(model_def$quant_batch_effects$is_positive)
                                                else integer(0)),
      NsubobjBatchEffects = nrow(model_data$subobject_batch_effects),
      subobj_batch_effect2quant_batch_effect = as.array(model_data$subobject_batch_effects$index_quant_batch_effect)
    ) %>%
    modifyList(matrix2csr("subobsXsubobjbatcheff", model_data$subobsXsubobjbatcheff))
    res <- modifyList(res, subobj_data)
  } else {
    res$quant2obs <- as.array(model_data$msdata$index_observation[!is.na(model_data$msdata$index_qdata)])
    res$miss2obs <- as.array(model_data$msdata$index_observation[!is.na(model_data$msdata$index_mdata)])
  }
  if ('index_mscalib' %in% names(model_data$mschannels)) {
    res$Nmsprotocols <- n_distinct(model_data$mschannels$index_mscalib)
    res$experiment2msproto <- as.array(model_data$mschannels$index_mscalib)
  }
  if ("Nsubobjects" %in% names(res)) {
    message(res$Niactions, " interaction(s) of ", res$Nobjects, " object(s) with ",
            res$Nsubobjects, " subobject(s), ",
            res$Nquanted, " quantitation(s) (", sum(res$quant_isreliable), " reliable), ",
            res$Nmissed, " missed (", sum(msdata_obs_flags.df$is_empty_observation), " in empty observations)")
  } else {
    message(res$Niactions, " interaction(s) of ", res$Nobjects, " object(s), ",
            res$Nquanted, " quantitation(s) (", sum(res$quant_isreliable), " reliable), ",
            res$Nmissed, " missed (", sum(msdata_obs_flags.df$is_empty_observation), " in empty observations)")
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
  return (cmdstanr::cmdstan_model(file.path(stan_models_path, paste0(model_name, ".stan"))))
}

#' @export
stan.sampling <- function(stan_input_data, iter=4000L, refresh=100L, chains=8L,
                          max_treedepth=12L, ...)
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
                                        c('subobj_shift_sigma'))
    }
    res <- stanmodel$sample(
              data = stan_input_data,
              #pars=unlist(lapply(vars_info, function(vi) vi$names)), include=TRUE,
              #init = function() { pcp_peaks_glm.generate_init_params(pcp_peaks_glm.model_data) },
              iter_warmup=0.5*iter, iter_sampling=0.5*iter,
              refresh=refresh, chains=chains, parallel_chains=chains,
              max_treedepth=max_treedepth, ...)
    attr(res, "msglm_vars_info") <- vars_info
    return(res)
}
