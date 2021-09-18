# extracts the value of given factor from the effect label
#' @export
effect_factor <- function(effects, factor_name, factor_levels, default = factor_levels[1]) {
    # FIXME factor levels are not allowed to start with '_' (to workaround factor with names being prefixes of the other factors)
    factor_vals <- str_match(effects, paste0("(?:^|:)", factor_name, "([^:_][^:]*)(?:$|:)"))[,2]
    factor_vals[is.na(factor_vals)] <- default
    factor(factor_vals, levels=factor_levels)
}

#' Define the experimental design of the MSGLM model
#'
#' @param conditionXeffect experimental design matrix,
#'        the rows are conditions, the columns are the effects
#' @param conditions the data frame of conditions
#' @param effects the data frame of effect with optional prior specification (FIXME)
#' @param msexperimentXeffect optional experimental design matrix for
#'        per-MS experiment effect specification. Normally it is not required as `conditionXeffect`
#'        matrix should be enough. It's reserved for the rare cases where the there is a high variability
#'        in biological responses between the biological replicates.
#' @param modelobject name of the MS object, which abundance would be modeled by MSGLM
#' @param quantobject name of the MS object to take the intensity information from
#' @param verbose
#'
#' @return msglm_model object
#' @export
#'
#' @examples
msglm_model <- function(conditionXeffect,
                        conditions, effects,
                        msexperimentXeffect = NULL,
                        modelobject = c("protgroup", "protregroup", "ptmgroup"),
                        quantobject = c("protgroup", "protregroup", "ptmgroup", "pepmodstate"),
                        verbose = FALSE)
{
  modelobject <- match.arg(modelobject)
  quantobject <- match.arg(quantobject)
  if (verbose) message("Initializing MSGLM model for ", modelobject,
                       "s using ", quantobject, " intensities")

  checkmate::assert_matrix(conditionXeffect, mode="numeric",
                           any.missing = FALSE, nrows = nrow(conditions),
                           min.cols = nrow(effects), max.cols=nrow(effects)+1L)
  # remove intercept if it's in the design matrix
  if ("(Intercept)" %in% rownames(conditionXeffect)) {
    if (verbose) warning('Removing (Intercept) effect from the conditionXeffect matrix (always present in the model)')
    conditionXeffect <- conditionXeffect[rownames(conditionXeffect) != "(Intercept)",, .drop=FALSE]
  }
  # FIXME allow row/col permutations, but then fix the row/col order
  checkmate::assert_matrix(conditionXeffect, mode="numeric",
                           any.missing = FALSE,
                           row.names = "unique", col.names = "unique")
  if (verbose) message("  * experimental design with ", nrow(conditionXeffect), " condition(s) and ",
                       ncol(conditionXeffect), " effect(s)")

  # process effects
  checkmate::assert_tibble(effects)
  # TODO remove automatic sounds-like column renames?
  effects <- maybe_rename(effects, c("prior_mean" = "mean", "prior_tau" = "tau"))
  checkmate::assert_names(names(effects), must.include = 'effect')
  effects <- dplyr::mutate(effects, effect = as.character(effect))
  checkmate::assert_character(effects$effect, unique=TRUE)
  # intercept should not be defined in effects
  checkmate::assert_disjunct(effects$effect, "(Intercept)")
  effects <- ensure_primary_index_column(effects, 'index_effect', create=TRUE,
                                         id_col = 'effect', ids_ordered = colnames(conditionXeffect))
  # define default effects prior
  if (rlang::has_name(effects, "prior_mean")) {
    checkmate::assert_numeric(effects$prior_mean, any.missing = FALSE, finite = TRUE)
  } else {
    if (verbose) message('Setting prior_mean=0 for all effects')
    effects <- dplyr::mutate(effects, prior_mean = 0.0)
  }
  if (rlang::has_name(effects, "prior_tau")) {
    checkmate::assert_numeric(effects$prior_tau, lower=0, any.missing = FALSE, finite=TRUE)
  } else {
    if (verbose) message('Setting prior_tau=1 for all effects')
    effects <- dplyr::mutate(effects, prior_tau = 1.0)
  }
  if (rlang::has_name(effects, "prior_df1")) {
    checkmate::assert_numeric(effects$prior_df1, lower=0, any.missing = FALSE, finite=TRUE)
  } else {
    if (verbose) message('Setting prior_df1=1 for all effects')
    effects <- dplyr::mutate(effects, prior_df1 = 1.0)
  }
  if (rlang::has_name(effects, "prior_df2")) {
    checkmate::assert_numeric(effects$prior_df2, lower=0, any.missing = FALSE, finite=TRUE)
  } else {
    if (verbose) message('Setting prior_df2=1 for all effects')
    effects <- dplyr::mutate(effects, prior_df2 = 1.0)
  }
  if (rlang::has_name(effects, "is_positive")) {
    checkmate::assert_logical(effects$is_positive, any.missing = FALSE)
  } else {
    if (verbose) message('Setting is_positive=FALSE for all effects')
    effects <- dplyr::mutate(effects, is_positive=FALSE)
  }

  # process conditions
  checkmate::assert_tibble(conditions)
  checkmate::assert_names(names(conditions), must.include = 'condition')
  conditions <- dplyr::mutate(conditions, condition = as.character(condition))
  checkmate::assert_character(conditions$condition, unique=TRUE)
  conditions <- ensure_primary_index_column(conditions, 'index_condition', create=TRUE,
                                            id_col="condition", ids_ordered = rownames(conditionXeffect))
  if (rlang::has_name(conditions, 'is_virtual')) {
    checkmate::assert_logical(conditions$is_virtual, any.missing = FALSE)
  } else {
    if (verbose) message('Setting is_virtual=FALSE for all conditions')
    conditions <- dplyr::mutate(conditions, is_virtual = FALSE)
  }

  model_def <- list(
    modelobject = modelobject, quantobject = quantobject,
    effects = effects,
    conditions = conditions,
    conditionXeffect = conditionXeffect.mtx,
    verbose = verbose
  )
  if (!is.null(msexperimentXeffect)) {
    if (verbose) message("MS experiment-specific experimental design specified")
    checkmate::assert_matrix(msexperimentXeffect, mode="numeric", any.missing = FALSE)
    msexp_dim <- names(dimnames(msexperimentXeffect))[[1]]
    checkmate::assert_choice(msexp_dim, c("msrun", "mschannel"))
    checkmate::assert_subset(colnames(msexp_dim), choices=effects$effect)
    model_def$msexperimentXeffect <- msexperimentXeffect
  }

  return(model_def)
}

#' Title
#'
#' @param model_def
#' @param mschannelXbatchEffect
#' @param batch_effects
#' @param applies_to
#' @param verbose
#'
#' @return updated msglm_model object
#' @export
#'
#' @examples
set_batch_effects <- function(model_def,
                              msexperimentXbatchEffect,
                              batch_effects = NULL,
                              applies_to = c('modelobject', 'quantobject'),
                              verbose = model_def$verbose
){
  applies_to <- match.arg(applies_to)
  id_col <- c(modelobject="batch_effect", quantobject="quant_batch_effect")[applies_to]
  df_name <- paste0(id_col, 's')
  checkmate::assert_matrix(msexperimentXbatchEffect, mode="numeric", any.missing=FALSE,
                           min.rows = 1, min.cols = 1, row.names = "unique", col.names = "unique")
  msexp_dim <- names(dimnames(msexperimentXbatchEffect))[[1]]
  checkmate::assert_choice(msexp_dim, c("msrun", "mschannel"))
  msexperiments <- rownames(msexperimentXbatchEffect)
  if (is.null(msexperiments)) stop("No names for MS experiments in the rows of the matrix")

  mtx_name <- paste0("msexperimentX",
                     c(modelobject="batchEffect", quantobject="quantBatchEffect")[[applies_to]])
  mtx_batch_effects <- colnames(msexperimentXbatchEffect)

  if (verbose) message("  * ", ncol(msexperimentXbatchEffect), " ", applies_to, "-level batch effect(s) specified")

  # process batch_effects
  if (!is.null(batch_effects)) {
    checkmate::assert_tibble(batch_effects)
    # TODO remove automatic sounds-like column renames?
    batch_effects <- maybe_rename(batch_effects, c("prior_mean" = "mean", "prior_tau" = "tau"))
    checkmate::assert_names(names(batch_effects), must.include = id_col)
    batch_effects <- dplyr::mutate(batch_effects, !!sym(id_col) := as.character(!!sym(id_col)))
    checkmate::assert_character(batch_effects[[id_col]], unique=TRUE)
    batch_effects <- ensure_primary_index_column(batch_effects, paste0('index_', id_col), create=TRUE,
                                                 id_col = id_col, ids_ordered = mtx_batch_effects)
    # define default batch_effects prior
    if (!rlang::has_name(batch_effects, "is_positive")) {
      if (verbose) message('Setting is_positive=FALSE for all ', df_name)
      batch_effects <- dplyr::mutate(batch_effects, is_positive=FALSE)
    } else {
      checkmate::assert_logical(batch_effects$is_positive, any.missing = FALSE)
    }
  } else {
    if (verbose) message('No batch_effect frame provided, using default batch effects settings')
    batch_effects <- tibble::tibble(!!sym("index_", id_col) := seq_along(mtx_batch_effects),
                                    !!sym(id_col) := mtx_batch_effects,
                                    is_positive = FALSE)
  }
  # update model_def
  model_def[[mtx_name]] <- msexperimentXbatchEffect
  model_def[[df_name]] <- batch_effects
  return(model_def)
}

#' Title
#'
#' @param model_def
#' @param metaconditionXcontrast
#' @param conditionXcontrast
#' @param verbose
#'
#' @return updated msglm_model object
#' @export
#'
#' @examples
set_contrasts <- function(model_def,
                          metaconditionXcontrast = NULL,
                          conditionXmetacondition = NULL,
                          contrasts = NULL,
                          conditionXcontrast = NULL,
                          verbose = model_def$verbose
){
  checkmate::assert_matrix(metaconditionXcontrast, mode="numeric", any.missing = FALSE,
                           row.names = "unique", col.names = "unique")
  metacond_dim <- names(dimnames(metaconditionXcontrast))[[1]]
  checkmate::assert_choice(metacond_dim, c("condition", "metacondition"))
  checkmate::assert_true(names(dimnames(metaconditionXcontrast))[[2]] == "contrast")
  if (verbose) message("  * ", ncol(metaconditionXcontrast), " contrast(s) using ",
                       nrow(metaconditionXcontrast), " ", metacond_dim, "s specified")

  if (!is.null(conditionXmetacondition)) {
    checkmate::assert_true(metacond_dim == "metacondition")
    checkmate::assert_matrix(conditionXmetacondition, mode="logical", any.missing = FALSE,
                             ncol = nrow(metaconditionXcontrast),
                             row.names="unique", col.names="unique")
    checkmate::assert_set_equal(names(dimnames(conditionXmetacondition)), c("condition", "metacondition"), ordered = TRUE)
    checkmate::assert_set_equal(colnames(conditionXmetacondition), rownames(metaconditionXcontrast))
    checkmate::assert_subset(rownames(conditionXmetacondition), model_def$conditions$condition)
  } else {
    if (verbose) message("No metaconditions specified, using conditions directly")
    conditionXmetacondition <- constant_matrix(FALSE, dimnames = list(condition = model_def$conditions$condition, metacondition = model_def$conditions$condition))
    # one-to-one correspondence between conditions and metaconditions
    for (cond in model_def$conditions$condition) {
      conditionXmetacondition[cond, cond] <- TRUE
    }
  }
  if (!is.null(contrasts)) {
    checkmate::assert_tibble(contrasts)
    contrasts <- ensure_primary_index_column(contrasts, "index_contrast", id_col = "contrast",
                                             ids_ordered = colnames(metaconditionXcontrast),
                                             create=TRUE)
  } else {
    # default contrast properties
    if (verbose) message("  * No contrast details specified, using the defaults")
    contrasts <- tibble(contrast = colnames(metaconditionXcontrast),
                        index_contrast = seq_len(ncol(metaconditionXcontrast)),
                        offset = 0,
                        contrast_type = "comparison")
  }

  # update model_def
  model_def$conditionXmetacondition <- conditionXmetacondition
  model_def$metaconditionXcontrast <- metaconditionXcontrast
  model_def$contrasts <- contrasts
  model_def$metaconditionXcontrast <- metaconditionXcontrast
  if (!is.null(conditionXcontrast)) {
    if (verbose) message("  * Using condition-level annotations for contrasts")
    checkmate::assert_tibble(conditionXcontrast)
    checkmate::assert_names(colnames(conditionXcontrast),
                            must.include = c("condition", "metacondition", "contrast",
                                             "is_preserved_condition"))
    checkmate::assert_subset(conditionXcontrast$condition, model_def$conditions$condition)
    checkmate::assert_subset(conditionXcontrast$contrast, contrasts$contrast)
    model_def$conditionXcontrast <- conditionXcontrast
  }
  return(model_def)
}
