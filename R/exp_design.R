# extracts the value of given factor from the effect label
#' @export
effect_factor <- function(effects, factor_name, factor_levels, default = factor_levels[1]) {
    # FIXME factor levels are not allowed to start with '_' (to workaround factor with names being prefixes of the other factors)
    factor_vals <- str_match(effects, paste0("(?:^|:)", factor_name, "([^:_][^:]*)(?:$|:)"))[,2]
    factor_vals[is.na(factor_vals)] <- default
    factor(factor_vals, levels=factor_levels)
}

#' Define the experimental design of the MSGLM model.
#'
#' @param conditionXeffect experimental design matrix,
#'        the rows are conditions, the columns are the effects
#' @param conditions the data frame of conditions
#' @param effects the data frame of effect with optional prior specification (FIXME)
#' @param msprobeXeffect optional experimental design matrix for
#'        per-MS probe effect specification. Normally it is not required as `conditionXeffect`
#'        matrix already specifies the study design. It's reserved for the rare cases where
#'        the there is a high variability in the scale of biological responses
#'        between the biological replicates, and the coefficients of `msprobeXeffect`
#'        matrix could be used to "normalize" the scale of effects in each MS probe.
#' @param verbose show extended diagnostic output
#'
#' @return *msglm_model* object
#' @export
#'
#' @examples
#' @seealso [set_batch_effects()], [set_contrasts()]
msglm_model <- function(conditionXeffect,
                        conditions, effects,
                        msprobeXeffect = NULL,
                        verbose = FALSE)
{
  if (verbose) message("Initializing MSGLM model")

  checkmate::assert_matrix(conditionXeffect, mode="numeric",
                           any.missing = FALSE, nrows = nrow(conditions),
                           min.cols = nrow(effects), max.cols=nrow(effects)+1L, min.rows = 1)

  # Check rank of the conditionXeffect matrix
  rank_conditionXeffect <- if (length(conditionXeffect) > 0) Matrix::rankMatrix(conditionXeffect) else 0
  if (rank_conditionXeffect < nrow(effects)) {
    warning('The rank of conditionXeffect matrix (', rank_conditionXeffect,
            ') is lower than the number of effects (', nrow(effects), '), i.e. there are redundant effects')
  }
  if (nrow(effects) < 1) {
    warning('No effects in the experimental design')
  }

  # remove intercept if it's in the design matrix
  if ("(Intercept)" %in% colnames(conditionXeffect)) {
    if (verbose) warning('Removing (Intercept) effect from the conditionXeffect matrix (always present in the model)')
    conditionXeffect <- conditionXeffect[, colnames(conditionXeffect) != "(Intercept)", drop=FALSE]
  }
  # FIXME allow row/col permutations, but then fix the row/col order
  checkmate::assert_matrix(conditionXeffect, mode="numeric",
                           any.missing = FALSE,
                           row.names = "unique", col.names = "unique")
  if (verbose) message("  * experimental design with ", nrow(conditionXeffect), " condition(s) and ",
                       ncol(conditionXeffect), " effect(s)")

  # process effects
  checkmate::assert_data_frame(effects)
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
  checkmate::assert_data_frame(conditions)
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

  model_def <- structure(list(
    effects = effects,
    conditions = conditions,
    conditionXeffect = conditionXeffect,
    verbose = verbose
  ), class="msglm_model")
  if (!is.null(msprobeXeffect)) {
    if (verbose) message("MS probe-specific experimental design specified")
    checkmate::assert_matrix(msprobeXeffect, mode="numeric", any.missing = FALSE)
    msprobe_dim <- names(dimnames(msprobeXeffect))[[1]]
    checkmate::assert_choice(msprobe_dim, c("msexperiment", "msprobe", "msrun", "mschannel"))
    checkmate::assert_subset(colnames(msprobe_dim), choices=effects$effect)
    model_def$msprobeXeffect <- msprobeXeffect
  }

  return(model_def)
}

#' Set *batch effects* for the MSGLM model.
#'
#' The batch effects are those that influence the measured data
#' (so they should be considered when fitting the model),
#' but don't represent the primary interest of the study (so
#' they should be defined separately from the ones in `conditionXeffect`
#' matrix). Batch effects are taken into account when fitting the
#' intensity of *quantobject* to the measurements, but are ignored when
#' inferring the adundance of *object-in-condition* and calculating
#' contrasts.
#'
#' *MSGLM* discriminates the two types of batch effects:
#'   * *object-level* batch effects that happen during sample preparation
#'     and affect the actual content of the MS *probe*. This includes the
#'     use of different sample preparation protocols, reagent batches etc.
#'   * *quantobject-level* batch effects that happen when the sample is being
#'     measured by MS. These effects cannot alter the contents of MS *probe*,
#'     but rather how it is measured by MS: which MS1 peaks (i.e. *quantobjects*)
#'     are identified and their relative intensities. Any alterations to
#'     MS protocol can introduce such batch effects: use of different MS
#'     instrument, LC system, gradient, fractionation, batch of MS-specific
#'     reagents etc.
#'
#' **Note:** While MS data sets with intensities specified at the *object-level*
#' (LFQ, iBAQ etc) may also contain batch effects related to MS protocol,
#' for the MSGLM model inference these batch effects are no different than
#' the ones related to sample preparation. In such cases, these batch effects
#' have to be specified via *object-level* batch effects matrix.
#'
#' @param model_def *msglm_model* object
#' @param batchEffectMatrix experimental design matrix for batch effects:
#'        the columns are the batch effects, and the rows are either *MS probes*
#'        (for *object-level* batch effects) or *MS channels* (for
#'        *quantobject-level* or MS protocol-related batch effects).
#'        MSGLM will check the row names of the name rows dimension to
#'        make sure that the matrix is specified correctly
#' @param batch_effects optional data frame with parameters of individual
#'        batch effects identified by `batch_effect`/`quant_batch_effect`
#'        column. Right now it supports specifying whether the batch
#'        effect is positive with `is_positive` column.
#' @param applies_to which type of batch effect it is: `object` for *object-level*
#'        batch effects, `quantobject` for *quantobject-level* ones.
#' @param verbose provide extended diagnostic output
#'
#' @return updated *msglm_model* object
#' @export
#'
#' @examples
#' @seealso [msglm_model()]
set_batch_effects <- function(model_def,
                              batchEffectMatrix,
                              batch_effects = NULL,
                              applies_to = c('object', 'quantobject'),
                              verbose = model_def$verbose
){
  checkmate::assert_class(model_def, "msglm_model")
  applies_to <- match.arg(applies_to)

  id_col <- c(object="batch_effect", quantobject="quant_batch_effect")[applies_to]
  df_name <- paste0(id_col, 's')
  checkmate::assert_matrix(batchEffectMatrix, mode="numeric", any.missing=FALSE,
                           min.rows = 1, min.cols = 1, row.names = "unique", col.names = "unique")
  msprobe_dim <- names(dimnames(batchEffectMatrix))[[1]]
  msprobe_dimnames_allowed <- c("msrun", "mschannel")
  if (applies_to == "object") msprobe_dimnames_allowed <- append(msprobe_dimnames_allowed, c("msexperiment", "msprobe"))
  checkmate::assert_choice(msprobe_dim, msprobe_dimnames_allowed)
  msprobes <- rownames(batchEffectMatrix)
  if (is.null(msprobes)) stop("No names for MS probes in the rows of the matrix")

  mtx_name <- c(object="msprobeXbatchEffect",
                quantobject="mschannelXquantBatchEffect")[[applies_to]]
  mtx_batch_effects <- colnames(batchEffectMatrix)

  if (verbose) message("  * ", ncol(batchEffectMatrix), " ", applies_to, "-level batch effect(s) specified")

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
    if (verbose) message('No batch_effects frame provided, using default batch effects settings')
    batch_effects <- tibble::tibble(!!sym(paste0("index_", id_col)) := seq_along(mtx_batch_effects),
                                    !!sym(id_col) := mtx_batch_effects,
                                    is_positive = FALSE)
  }
  # update model_def
  if (rlang::has_name(model_def, mtx_name)) warning("Redefining ", mtx_name)
  model_def[[mtx_name]] <- batchEffectMatrix
  if (rlang::has_name(model_def, df_name)) warning("Redefining ", df_name)
  model_def[[df_name]] <- batch_effects
  return(model_def)
}

#' Set *contrasts* for the MSGLM model.
#'
#' TODO description
#'
#' @param model_def *msglm_model* object
#' @param metaconditionXcontrast matrix that defines *contrasts* (matrix columns)
#'        as linear combinations of *metaconditions* (matrix rows)
#' @param conditionXmetacondition optional logical matrix defining *metaconditions*
#'        (matrix columns) as subsets of model's *conditions* (matrix rows). If
#'        this matrix is not specified, it is assumed that *metaconditions* are
#'        identical to *conditions*
#' @param contrasts optional data frame defining the properties of *contrasts*
#' @param conditionXcontrast optional data frame defining how specific *conditions*
#'        are used in individual *contrast*
#' @param verbose provide extended diagnostic output
#'
#' @return updated *msglm_model* object
#' @export
#'
#' @examples
#' @seealso [msglm_model()]
set_contrasts <- function(model_def,
                          metaconditionXcontrast,
                          conditionXmetacondition = NULL,
                          contrasts = NULL,
                          conditionXcontrast = NULL,
                          verbose = model_def$verbose
){
  checkmate::assert_class(model_def, "msglm_model")
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
  if (rlang::has_name(model_def, "conditionXmetacondition")) warning("Redefining conditionXmetacondition")
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
