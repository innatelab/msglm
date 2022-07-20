# variables description for msglm_local model
msglm.vars_info <- list(
  global = list(names=c('effect_slab_c'),# 'batch_effect_slab_c'),
                dims=c()),#'qobj_shift_sigma', 'obj_effect_tau'), dims = c() ),
  #batch_effects = list(names = c('batch_effect_sigma'),
  #                     dims = c('batch_effect')),
  object_conditions = list(names=c('obj_probe_shift_sigma', 'obj_cond_labu', 'obj_cond_labu_replCI'),
                           dims=c('object_condition')),
  object_msprobes = list(names=c('obj_probe_labu', "obj_probe_shift", "obj_probe_batch_shift"),
                         dims=c('object_msprobe')),
  quantobjects = list(names=c('qobj_shift', 'qobj_llh'),
                      dims=c('quantobject')),
  #quantobject_batch_shifts = list(names=c('qobj_batch_shift'),
  #                                dims=c('quantobject_mschannel')),
  quantobject_batch_effects = list(names=c('qobj_batch_effect', 'qobj_batch_effect_sigma'),
                                   dims=c('quantobject_batch_effect')),
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

# convert object into a list of Stan variables for the "data" block
#' @export
to_standata <- function(obj, ...) UseMethod("to_standata")

#' Coverts the data and experimental design into Stan format.
#'
#' @param model_data the list with MS data an experimental design
#'
#' @param effect_slab_df the *degrees of freedom* for the prior of object effect *slab* regularization parameter
#' @param effect_slab_scale the *scale* parameter for the prior of object effect *slab* regularization parameter
#' @param object_labu_min minimal object log-abundance
#' @param object_labu_min_scale how strongly the estimates below *obj_labu_min* would be penalized (smaller = more stringent)
#' @param object_msprobe_shift_tau the *tau* parameter for the prior of biochemical log-abundance variation between the *probes* (replicates) of the same *condition*
#' @param object_msprobe_shift_df the *degrees of freedom* parameter for the prior of biochemical log-abundance variation between the *probes* (replicates) of the same *condition*
#' @param hsprior_lambda_a_offset offset from zero for the "a" parameter of horseshoe priors. Helps NUTS integration step
#' @param hsprior_lambda_t_offset offset from zero for the "a" parameter of horseshoe priors. Helps NUTS integration step
#' @param batch_effect_sigma sigma parameter for the normal distribution prior of batch effects
#' @param quant_batch_tau
#' @param quant_batch_df
#' @param quant_batch_c
#' @param empty_msprobe_sigmoid_scale how match the lower bound prior of object-in-msprobe abundance
#'                                    is relaxed if there are no MS observations of the object in a given MS probe
#'
#' @export
to_standata.msglm_model_data <- function(model_data,
                              effect_slab_df = 4, effect_slab_scale = 2.5,
                              object_labu_min_scale = 1,
                              object_msprobe_shift_tau=0.05, object_msprobe_shift_df=4.0,
                              hsprior_lambda_a_offset = 0.05, hsprior_lambda_t_offset = 0.01,
                              batch_effect_sigma=0.5,
                              quant_batch_tau=1.0, quant_batch_df=2, quant_batch_c=10,
                              quantobject_shift_df=2.0, quantobject_shift_sigma=1.0,
                              empty_msprobe_sigmoid_scale = 1.0,
                              verbose = FALSE)
{
  model_def <- model_data$model_def
  if (verbose) message('Converting MSGLM model data to Stan-readable format...')
  ensure_primary_index_column(model_def$effects, 'index_effect')
  ensure_primary_index_column(model_data$msprobes, 'index_msprobe')
  ensure_primary_index_column(model_data$mschannels, 'index_mschannel')

  is_glmm <- rlang::inherits_all(model_def, "msglmm_model")
  xaction_ix_col <- if (is_glmm) "index_object_mixture" else "index_object_condition"
  objprobes_df <- dplyr::select(model_data$object_msprobes, index_object_msprobe, index_msprobe,
                                index_object, !!xaction_ix_col)
  if (any(objprobes_df$index_object_msprobe != seq_len(nrow(objprobes_df)))) {
    stop("model_data$msdata not ordered by object probe / have missing object probe indices")
  }
  msdata_objprobe_flags.df <- dplyr::group_by(model_data$msdata, index_object_msprobe) %>%
    # FIXME set as non-empty if there are idents
    dplyr::transmute(is_empty_msprobe = all(is.na(intensity))) %>%
    dplyr::ungroup()
  checkmate::assert_set_equal(msdata_objprobe_flags.df$index_object_msprobe,
                              model_data$msdata$index_object_msprobe, ordered=TRUE)

  missing_mask <- is.na(model_data$msdata$intensity)
  res <- list(
    Nconditions = nrow(model_def$conditions),
    Nobjects = nrow(model_data$objects),

    Nprobes = nrow(model_data$msprobes),

    NobjProbes = nrow(objprobes_df),
    obj_probe2probe = as.array(objprobes_df$index_msprobe),

    Nmschannels = n_distinct(model_data$mschannels$index_mschannel),
    mschannel_shift = as.array(model_data$mschannels$mschannel_shift),

    Neffects = ncol(model_def$conditionXeffect),
    effect_is_positive = as.array(as.integer(model_def$effects$is_positive)),
    effect_tau = as.array(model_def$effects$prior_tau),
    effect_mean = as.array(model_def$effects$prior_mean),
    effect_df = as.array(model_def$effects$prior_df1),
    effect_df2 = as.array(model_def$effects$prior_df2),
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

    obj_base_labu_sigma = 15.0,
    obj_labu_min = rlang::set_names(model_data$quantobject_labu_min, NULL),
    obj_labu_min_scale = object_labu_min_scale,

    NobjConditions = nrow(model_data$object_conditions),
    obj_cond2obj = as.array(model_data$object_conditions$index_object),

    Nquanted = sum(!missing_mask),
    Nmissed = sum(missing_mask),
    missing_sigmoid_scale = as.array(if_else(msdata_objprobe_flags.df$is_empty_msprobe[missing_mask],
                                             empty_msprobe_sigmoid_scale, 1.0)),
    qData = as.array(model_data$msdata$intensity[!missing_mask]),
    quant_isreliable = as.array(model_data$msdata$is_reliable[!missing_mask]),

    hsprior_lambda_a_offset = hsprior_lambda_a_offset,
    hsprior_lambda_t_offset = hsprior_lambda_t_offset,

    obj_probe_shift_tau = object_msprobe_shift_tau,
    obj_probe_shift_df = object_msprobe_shift_df
  ) %>%
    modifyList(to_standata(model_data$quantobject_mscalib, silent=!verbose)) %>%
    modifyList(matrix2stancsr(model_data$object_conditionXeffect, "obj_condXeff")) %>%
    modifyList(matrix2stancsr(model_data$object_msprobeXbatch_effect, "obj_probeXbatcheff"))

  res <- modifyList(res, object_msprobeXcondition_to_standata(model_data))

  if (rlang::has_name(model_data, "quantobjects")) {
    # data have quantobjects
    qobj_data <- list(
      Nquantobjects = nrow(model_data$quantobjects),
      quantobj2obj = as.array(model_data$quantobjects$index_object),
      NqobjChannels = nrow(model_data$msdata),
      quant2qobj_channel = as.array(model_data$msdata$index_quantobject_mschannel[!is.na(model_data$msdata$index_qdata)]),
      miss2qobj_channel = as.array(model_data$msdata$index_quantobject_mschannel[!is.na(model_data$msdata$index_mdata)]),
      qobj_channel2quantobj = as.array(model_data$msdata$index_quantobject),
      qobj_channel2mschannel = as.array(model_data$msdata$index_mschannel),
      qobj_channel2obj_probe = as.array(model_data$msdata$index_object_msprobe),
      # TODO support different noise models
      Nmsprotocols = 0L,
      mschannel2msproto = integer(0),

      qobj_shift_df = quantobject_shift_df,
      qobj_shift_sigma = quantobject_shift_sigma,

      # quantobject-specific batch effects
      NquantBatchEffects = n_distinct(model_data$quantobject_batch_effects$index_quant_batch_effect),
      quant_batch_effect_tau = quant_batch_tau,
      quant_batch_effect_df = quant_batch_df,
      quant_batch_effect_c = quant_batch_c,
      quant_batch_effect_is_positive = as.array(if (nrow(model_data$quantobject_batch_effects)>0)
                                                as.integer(model_def$quant_batch_effects$is_positive)
                                                else integer(0)),
      NqobjBatchEffects = nrow(model_data$quantobject_batch_effects),
      qobj_batch_effect2quant_batch_effect = as.array(model_data$quantobject_batch_effects$index_quant_batch_effect)
    ) %>%
    modifyList(matrix2stancsr(model_data$quantobject_mschannelXquant_batch_effect, "qobj_channelXqbatcheff"))
    res <- modifyList(res, qobj_data)
  } else {
    if (rlang::has_name(model_def, "quant_batch_effects")) {
      stop("Object-level model data provides quant_batch_effects, please rerun msglm_data().")
    }
    res$quant2obj_probe <- as.array(model_data$msdata$index_object_msprobe[!is.na(model_data$msdata$index_qdata)])
    res$quant2mschannel <- as.array(model_data$msdata$index_mschannel[!is.na(model_data$msdata$index_qdata)])
    res$miss2obj_probe <- as.array(model_data$msdata$index_object_msprobe[!is.na(model_data$msdata$index_mdata)])
    res$miss2mschannel <- as.array(model_data$msdata$index_mschannel[!is.na(model_data$msdata$index_mdata)])
  }
  if (rlang::has_name(model_data, 'index_msprotocol')) {
    res$Nmsprotocols <- n_distinct(model_data$mschannels$index_msprotocol)
    res$mschannel2msproto <- as.array(model_data$mschannels$index_msprotocol)
  }
  message(res$NobjConditions, " object-in-condition(s) of ", res$Nobjects, " object(s)",
          if (rlang::has_name(res, "Nquantobjects")) paste0(" with ", res$Nquantobjects, " quantobject(s)") else "", ",",
          res$Nquanted, " quantitation(s) (", sum(res$quant_isreliable), " reliable), ",
          res$Nmissed, " missed (", sum(msdata_objprobe_flags.df$is_empty_msprobe), " in empty object probe(s))")
  return(structure(res, class="msglm_standata", msglm_model_data=model_data))
}

object_msprobeXcondition_to_standata <- function(model_data) {
  UseMethod("object_msprobeXcondition_to_standata")
}

object_msprobeXcondition_to_standata.msglm_model_data <- function(model_data) {
  res <- matrix2stancsr(model_data$object_msprobeXeffect, "obj_probeXeff")
  res$obj_probe2obj_cond <- as.array(model_data$object_msprobes$index_object_condition)
  return(res)
}

object_msprobeXcondition_to_standata.msglmm_model_data <- function(model_data) {
  model_def <- model_data$model_def
  ensure_primary_index_column(model_def$mixeffects, 'index_mixeffect')

  # FIXME remove model_def$mixeffects <- maybe_rename(model_def$mixeffects, c("prior_mean" = "mean", "prior_tau" = "tau"))
  message("Setting GLMM object mixtures data...")
  res <- list(
    NobjMixtures = nrow(model_data$object_mixtures),
    obj_probe2obj_mix = as.array(model_data$object_msprobes$index_object_mixture),
    obj_mix2obj = as.array(model_data$object_mixtures$index_object),

    Nmix = nrow(model_def$mixcoefXeff),
    NmixEffects = nrow(model_def$mixeffects),
    mixeffect_mean = as.array(model_def$mixeffects$prior_mean),
    mixeffect_tau = as.array(model_def$mixeffects$prior_tau),
    mixcoefXeff = as.matrix(model_def$mixcoefXeff),

    Nmixtions = nrow(model_data$mixtions),
    mixt2obj_cond = as.array(model_data$mixtions$index_object_condition),
    mixt2mix = as.array(model_data$mixtions$mixcoef_ix)
  )
  return(res)
}

# all models
msglm_model_names <- function() {
  system.file('stan_models', package="msglm", mustWork=FALSE) %>%
  list.files(pattern="\\.stan$") %>% sapply(basename) %>% str_remove("\\.stan$")
}

msglm_stan_model <- function(model_name, ...) {
  stan_models_path <- system.file('stan_models', package="msglm", mustWork=FALSE)
  message("Loading ", model_name, " Stan model")
  return (cmdstanr::cmdstan_model(file.path(stan_models_path, paste0(model_name, ".stan")),
                                  include_paths = stan_models_path, ...))
}

stan_model_name <- function(standata) {
  model_data <- attr(standata, "msglm_model_data")
  model_def <- model_data$model_def
  res <- c(msglm_model = "msglm_local",
           msglmm_model = "msglmm_local")[class(model_def)[[1]]]
  obj <- model_data$msentities[['object']]
  quantobj <- model_data$msentities[['quantobject']]
  if (obj == quantobj) {
    # do nothing
  } else if (quantobj == "pepmodstate") {
    res <- paste0(res, '_quantobjects')
  } else {
     stop("Unsupported combination of object=", obj,
          " and quantobject=", quantobj)
  }
  return(res)
}

#' @param model_data either *msglm_stan_data* or *msglm_model_data* object
#' @seealso [prepare_stan_data()], [msglm_data()]
#' @export
fit_model <- function(model_data, ...) UseMethod("fit_model")

#' @export
fit_model.msglm_model_data <- function(model_data, stanmodel_options=list(), ...) {
  standata <- rlang::exec(to_standata, model_data, !!!stanmodel_options)
  fit_model(standata, ...)
}

#' @export
fit_model.msglm_standata <- function(standata, method = c("mcmc", "variational"),
                                     iter=4000L, refresh=100L, chains=8L,
                                     max_treedepth=12L,
                                     advi_iter=NULL, ...)
{
    method <- match.arg(method)
    # TODO convert "if" into virtual method(s)
    vars_info <- if (rlang::has_name(standata, "Nsupactions")) {
      msglmm.vars_info
    } else {
      msglm.vars_info
    }
    if (!rlang::has_name(standata, "Nquantobjects")) {
      # exclude quantobject-related
      vars_info$quantobjects <- NULL
      vars_info$quantobject_batch_effects <- NULL
      vars_info$global$names <- setdiff(vars_info$global$names,
                                        c('qobj_shift_sigma'))
    }
    message("Running Stan (", method, " mode)...")
    stanmodel <- msglm_stan_model(stan_model_name(standata))
    if (method == "mcmc") {
      res <- stanmodel$sample(
                data = structure(standata, class="list"),
                #pars=unlist(lapply(vars_info, function(vi) vi$names)), include=TRUE,
                iter_warmup=0.5*iter, iter_sampling=0.5*iter,
                refresh=refresh, chains=chains, parallel_chains=chains,
                max_treedepth=max_treedepth, ...)
    } else if (method == "variational") {
      # FIXME take thin into account
      res <- stanmodel$variational(
                data = structure(standata, class="list"),
                iter=advi_iter, output_samples=0.5*iter*chains, ...)
    }
    attr(res, "msglm_vars_info") <- vars_info
    return(res)
}
