# FIXME
#' @export
mschannel_quantobj_statistics <- function(msdata) {
  checkmate::assert_class(msdata, "msglm_data_collection")
  quantobj <- msdata$msentities[['quantobject']]
  mschan <- msdata$msentities[['mschannel']]
  msrun <- msdata$msentities[['msrun']]
  mstag <- msdata$msentities[['mstag']]
  quantobjs_dfname <- str_c(quantobj, "s")
  quantobj_idcol = str_c(quantobj, "_id")
  if (!rlang::has_name(msdata, quantobjs_dfname)) {
    stop("msdata contains no `", quantobjs_dfname, " data")
  }
  quantobj_intensities_dfname <- str_c(quantobj, "_intensities")
  if (!rlang::has_name(msdata, quantobj_intensities_dfname)) {
    stop("msdata contains no `", quantobj_intensities_dfname, "` data")
  }
  quantobj_intensities_df <- msdata[[quantobj_intensities_dfname]]
  groupat <- c(msrun=msrun)
  if (!is.na(mstag)) {
    groupat[['mstag']] <- mstag
  }
  quantobj_intensities_df <- dplyr::select_at(quantobj_intensities_df,
      append(groupat, c(quantobj_id=quantobj_idcol, "intensity"))) %>%
    dplyr::filter(!is.na(quantobj_id))
  if (!is.na(mstag)) {
    groupat[['mschannel']] <- mschan
    mschans_df <- dplyr::select_at(msdata[[paste0(mschan,"s")]], groupat) %>% dplyr::distinct()
    quantobj_intensities_df <- dplyr::inner_join(quantobj_intensities_df, mschans_df, by = c("msrun", "mstag"))
  }
  expandat <- append(names(groupat), "quantobj_id")
  intensity_stats_df <- tidyr::expand(quantobj_intensities_df, !!!syms(expandat)) %>%
    dplyr::left_join(quantobj_intensities_df, by=expandat) %>%
    #dplyr::group_by(quantobj_id, condition) %>% dplyr::filter(any(!is.na(intensity))) %>%
    dplyr::group_by_at(names(groupat)) %>%
    summarize(log2_intensity.mean = mean(log2(intensity[!is.na(intensity)])),
              log2_intensity.median = median(log2(intensity[!is.na(intensity)])),
              log2_intensity.sd = sd(log2(intensity[!is.na(intensity)])),
              n = n(),
              n_missing = sum(is.na(intensity)), .groups="drop")

  quantobj_idents_dfname <- str_c(quantobj, "_idents")
  if (!rlang::has_name(msdata, quantobj_idents_dfname)) {
    warning("msdata contains no `", quantobj_idents_dfname, "` data, using intensities instead")
    quantobj_idents_dfname <- quantobj_intensities_dfname
  }
  quantobj_idents_df <- dplyr::select_at(msdata[[quantobj_idents_dfname]],
    c(quantobj_id=quantobj_idcol, msrun=msrun, "ident_type"))
  ident_stats_df <- dplyr::mutate(quantobj_idents_df,
      is_matching = ident_type %in% c("By matching", "MULTI-MATCH", "MULTI-MATCH-MSMS"),
      is_msms = ident_type %in% c("By MS/MS", "MULTI-MSMS", "ISO-MSMS", "MSMS", "MULTI-SECPEP")) %>%
    dplyr::group_by(msrun) %>%
    summarize(n_matching = n_distinct(quantobj_id[is_matching], na.rm = TRUE),
              n_msms = n_distinct(quantobj_id[is_msms], na.rm = TRUE),
              .group = "drop")
  return (dplyr::left_join(intensity_stats_df, ident_stats_df, by="msrun"))
}

#' Collect MS data that would be used by MSGLM models.
#'
#' @param msdata list of various MS-related data frames
#' @param model_def *msglm_model* object
#' @param mscalib *mscalib* MS intensities noise model (on "quantobject" leve)
#' @param verbose if TRUE, produces debugging output
#'
#' @return *msglm_data_collection* object with the data for all potential model objects
#' @export
#'
#' @examples
import_msglm_data <- function(msdata, model_def,
                              mscalib = get(paste0(model_def$quantobject, "_mscalib")),
                              condition = "condition",
                              msexperiment = NA_character_,
                              mschannel = NA_character_,
                              msrun = NA_character_,
                              msfraction = NA_character_,
                              mstag = NA_character_,
                              min_intensity_quantile = 0.001,
                              min_intensity_offset = -5,
                              verbose=model_def$verbose)
{
  checkmate::assert_class(model_def, "msglm_model")
  checkmate::assert_list(msdata)
  res <- structure(list(), class="msglm_data_collection")
  modelobj <- model_def$modelobject
  modelobj_idcol <- paste0(modelobj, "_id")
  quantobj <- model_def$quantobject
  quantobj_idcol <- paste0(quantobj, "_id")

  # FIXME support experiments with mschannels
  if (!is.na(msexperiment)) {
    msexp <- msexperiment
  } else if (rlang::has_name(msdata, "msexperiments")) {
    msexp <- "msexperiment"
  } else {
    msexp <- NA_character_
  }
  if (!is.na(msexp)) {
    msexps_dfname <- paste0(msexp, "s")
    if (verbose) message("Using msdata$", msexps_dfname, " as MS experiments source")
    msexps_df <- msdata[[msexps_dfname]]
    checkmate::assert_data_frame(msexps_df, .var.name = paste0("msdata$", msexps_dfname))
    checkmate::assert_names(names(msexps_df),
                            must.include = c(msexp, condition))
    checkmate::assert_character(as.character(msexps_df[[msexp]]), unique=TRUE, any.missing = FALSE,
                                .var.name = paste0("msdata$", msexps_dfname, "$", msexp))
    mschan_parent <- msexp
  } else {
    if (verbose) message("No explicit MS experiments frame specified")
    msexps_df <- NULL
    msexp <- NA_character_
    mschan_parent <- condition
  }

  mschan_required_cols <- NULL
  if (!is.na(mschannel)) {
    mschan <- mschannel
    if (!is.na(mstag)) {
      if (is.na(msrun)) stop("Specifying mstag column (", mstag, ") also requires specifying msrun column")
      mschan_required_cols <- c(mschan, msrun, mstag, mschan_parent)
    } else {
      if (!is.na(msrun)) stop("Specifying mschannel (", mschan, ") and msrun (", msrun, ") column also requires specifying mstag column")
    }
  } else if (!is.na(msrun)) {
    if (!is.na(mstag)) stop("For data with mstag column (", mstag, "), mschannel column has to be specified")
    mschan <- msrun
  } else if (rlang::has_name(msdata, "mschannels")) {
    mschan <- "mschannel"
    mstag <- mstag %||% "mstag"
    msrun <- msrun %||% "msrun"
    mschan_required_cols <- c(mschan, mstag, msrun, mschan_parent)
  } else if (rlang::has_name(msdata, "msruns")) {
    if (!is.na(mstag)) stop("MS tag column (", mstag, ") specified, but no tagged MS channels found")
    mschan <- "msrun"
    msrun <- "msrun"
  } else {
    stop("Cannot find MS channels information (msdata$mschannels or msdata$msruns)")
  }
  mschan_required_cols <- mschan_required_cols %||% c(mschan, mschan_parent)
  if (verbose) message("Using msdata$", mschan, "s as MS channels source")
  mschans_dfname <- paste0(mschan, "s")
  mschans_df <- msdata[[mschans_dfname]]
  checkmate::assert_data_frame(mschans_df, .var.name = paste0("msdata$", mschans_dfname))
  if (!is.na(msfraction)) {
    msfrac <- mfraction
  } else if (rlang::has_name(mschans_df, "msfraction")) {
    msfrac <- "msfraction"
    if (verbose) message("MS fraction column (", msfrac, ") detected")
  }
  if (!is.na(msfrac)) {
    mschan_required_cols <- c(mschan_required_cols, msfrac)
  }
  checkmate::assert_names(names(mschans_df), must.include = mschan_required_cols)
  if (is.na(msexp)) {
    msexp <- mschan
    if (verbose) message("Using msdata$", msexp, "s as MS experiments source")
  }

  res$msentities <- c(
    quantobject = quantobj,
    modelobject = modelobj,
    condition = condition,
    msexperiment = msexp,
    mschannel = mschan,
    msrun = msrun,
    mstag = mstag,
    msfraction = msfrac)

  msexp_idcol <- msexp
  mschan_idcol <- mschan
  msrun_idcol <- msrun
  if (rlang::has_name(msexps_df, "is_used")) {
    checkmate::assert_logical(msexps_df$is_used, any.missing=FALSE,
                              .var.name=paste0("msdata$", msexp, "s$is_used"))
    if (verbose) message("Importing ", sum(msexps_df$is_used), " used ", msexp,
                         "(s) of ", nrow(msexps_df))
    msexps_df <- dplyr::filter(msexps_df, is_used)
  } else {
    if (verbose) message("Importing all ", nrow(msexps_df), " MS experiment(s)")
  }
  virtual_conditions <- setdiff(unique(msexps_df[[condition]]), model_def$conditions$condition)
  if (length(virtual_conditions) > 0) {
    warning(length(virtual_conditions), ' condition(s) not covered by MS experiments')
  }
  checkmate::assert_subset(as.character(msexps_df[[condition]]), as.character(model_def$conditions$condition),
                           .var.name = paste0("msdata$", msexps_dfname, "$", condition))
  res[[msexps_dfname]] <- msexps_df
  if (mschan != msexp) {
    mschans_df <- dplyr::semi_join(msdata[[mschans_dfname]],
                                  dplyr::select(msexps_df, !!sym(msexp_idcol)),
                                  by=msexp_idcol)
    if (rlang::has_name(mschans_df, "is_used")) {
      checkmate::assert_logical(mschans_df$is_used, any.missing=FALSE,
                                .var.name=paste0("msdata$", mschan, "s$is_used"))
      if (verbose) message("Importing ", sum(mschans_df$is_used), " used ", mschan,
                          "(s) of ", nrow(mschans_df))
      mschans_df <- dplyr::filter(mschans_df, is_used)
    } else {
      if (verbose) message("Importing all ", nrow(mschans_df), " MS channel(s)")
    }
    res[[mschans_dfname]] <- mschans_df
    checkmate::assert_set_equal(mschans_df[[msexp_idcol]], msexps_df[[msexp_idcol]],
                                .var.name = paste0("msdata$", mschans_dfname, "$", msexp_idcol))
  } else {
    mschans_df <- msexps_df
  }

  quantobj_intensities_dfname <- paste0(quantobj, "_intensities")
  if (rlang::has_name(msdata, quantobj_intensities_dfname)) {
    if (verbose) message("Importing ", quantobj_intensities_dfname, "...")
    quantobj_intensities_df <- msdata[[quantobj_intensities_dfname]]
    checkmate::assert_data_frame(quantobj_intensities_df)
    checkmate::assert_names(colnames(quantobj_intensities_df),
                            must.include = c(mschan_idcol, quantobj_idcol, "intensity"),
                            .var.name = paste0("msdata$", quantobj_intensities_dfname))
    quantobj_intensities_df <- dplyr::semi_join(quantobj_intensities_df,
                                                dplyr::select(mschans_df, !!sym(mschan_idcol)),
                                                by=mschan_idcol)
  } else {
    stop("msdata$", quantobj_intensities_dfname, " not found")
  }
  modelobjs_dfname <- paste0(modelobj, "s")
  if (rlang::has_name(msdata, modelobjs_dfname)) {
    modelobjs_df <- msdata[[modelobjs_dfname]]
    if (verbose) message("Found ", nrow(modelobjs_df), " ", modelobj, "(s)")
  } else {
    stop("msdata$", modelobjs_dfname, " not found")
  }
  if (quantobj != modelobj) {
    quantobjs_dfname <- paste0(quantobj, "s")
    if (rlang::has_name(msdata, quantobjs_dfname)) {
      quantobjs_df <- msdata[[quantobjs_dfname]]
      if (verbose) message("Found ", nrow(quantobjs_df), " ", quantobj, "(s)")
    } else {
      stop("msdata$", quantobjs_dfname, " not found")
    }
  } else {
    quantobjs_dfname <- modelobjs_dfname
    quantobjs_df <- modelobjs_df
  }
  modelobj_idents_dfname <- paste0(modelobj, "_idents")
  if (rlang::has_name(msdata, modelobj_idents_dfname)) {
    if (verbose) message("Importing ", modelobj_idents_dfname, "...")
    modelobj_idents_df <- msdata[[modelobj_idents_dfname]]
    checkmate::assert_data_frame(modelobj_idents_df)
    checkmate::assert_names(colnames(modelobj_idents_df),
                            must.include = c(msexp_idcol, modelobj_idcol, "ident_type"),
                            .var.name = paste0("msdata$", modelobj_idents_dfname))
    res[[modelobj_idents_dfname]] <- dplyr::semi_join(modelobj_idents_df,
                                                      dplyr::select(msexps_df, !!sym(msexp_idcol)),
                                                      by=msexp_idcol)
  } else {
    warning("msdata$", modelobj_idents_dfname, " not found")
  }
  if (modelobj == quantobj) {
    quantobj_intensities_df <- dplyr::filter(quantobj_intensities_df, !is.na(intensity))
    if (verbose) message(nrow(quantobj_intensities_df), " of ", nrow(quantobj_intensities_orig_df),
                         " ", quantobj, " intensity measurements selected")
    modelobjs_df <- dplyr::semi_join(modelobjs_df, quantobj_intensities_dfname, by=modelobj_idcol)
    if (verbose) message(nrow(modelobjs_df), " of ", nrow(modelobjs_orig_df),
                         " ", modelobj, "s selected")
  } else if (quantobj == "pepmodstate") {
    modelobj2pepmod_dfname = paste0(modelobj, "2pepmod")
    modelobj2quantobj_dfname = paste0(modelobj, "2", quantobj)
    modelobj2pepmod_df <- msdata[[modelobj2pepmod_dfname]]
    if (!rlang::has_name(msdata, modelobj2quantobj_dfname)) {
      warning("msdata$", modelobj2quantobj_dfname, " not found, trying msdata$", modelobj2pepmod_dfname)
      if (!rlang::has_name(msdata, modelobj2pepmod_dfname)) {
        stop("msdata$", modelobj2pepmod_dfname, " not found")
      }
      modelobj2quantobj_df <- dplyr::inner_join(modelobj2pepmod_df, quantobjs_df, by="pepmod_id") %>%
        dplyr::select(!!syms(modelobj_idcol, quantobj_idcol, "pepmod_id", "is_specific")) %>% dplyr::distinct()
    } else {
      modelobj2quantobj_df <- msdata[[modelobj2quantobj_dfname]]
      modelobj2pepmod_df <- NULL
    }
    checkmate::assert_data_frame(modelobj2quantobj_df, .var.name=paste0("msdata$", modelobj2quantobj_dfname))
    checkmate::assert_names(colnames(modelobj2quantobj_df),
                            must.include=c(paste0(model_def$modelobject, "_id"),
                                           paste0(model_def$quantobject, "_id"),
                                           "is_specific"),
                            .var.name = paste0("msdata$", modelobj2quantobj_dfname))
    if (verbose) message("Filtering for specific and quantified ", quantobj,
                         "s and corresponding ", modelobj, "s...")
    quantobjs_orig_df <- quantobjs_df
    quantobj_intensities_orig_df <- quantobj_intensities_df
    modelobjs_orig_df <- modelobjs_df
    modelobj2quantobj_df <- dplyr::filter(modelobj2quantobj_df, is_specific)
    # TODO require identification in a selected set of conditions
    # (e.g. for AP-MS data should be identified in AP-MS conditions,
    #  whereas identification in the background is not required)
    quantobj_intensities_df <- dplyr::filter(quantobj_intensities_df, !is.na(intensity)) %>%
        dplyr::semi_join(modelobj2quantobj_df, by=quantobj_idcol)
    if (verbose) message("  * ", nrow(quantobj_intensities_df), " of ", nrow(quantobj_intensities_orig_df),
                         " specific ", quantobj, " intensity measurements selected")
    quantobjs_df <- dplyr::semi_join(quantobjs_df, quantobj_intensities_df, by=quantobj_idcol)
    if (verbose) message("  * ", nrow(quantobjs_df), " of ", nrow(quantobjs_orig_df),
                         " specific ", quantobj, "s selected")
    modelobj2quantobj_df <- dplyr::semi_join(modelobj2quantobj_df, quantobjs_df, by=quantobj_idcol)
    modelobjs_df <- dplyr::semi_join(modelobjs_df, modelobj2quantobj_df, by=modelobj_idcol)
    if (verbose) message("  * ", nrow(modelobjs_df), " of ", nrow(modelobjs_orig_df),
                         " ", modelobj, "s have specific ", quantobj, " quantitations")
    res[[modelobj2quantobj_dfname]] <- modelobj2quantobj_df
    res[[quantobjs_dfname]] <- quantobjs_df

    quantobj_idents_dfname <- paste0(quantobj, "_idents")
    if (rlang::has_name(msdata, quantobj_idents_dfname)) {
      if (verbose) message("Importing ", quantobj_idents_dfname, "...")
      quantobj_idents_df <- msdata[[quantobj_idents_dfname]]
      checkmate::assert_data_frame(quantobj_idents_df)
      checkmate::assert_names(colnames(quantobj_idents_df),
                              must.include = c(msrun_idcol, quantobj_idcol, "ident_type"),
                              .var.name = paste0("msdata$", quantobj_idents_dfname))
      res[[quantobj_idents_dfname]] <- dplyr::semi_join(quantobj_idents_df,
                                                        dplyr::select(mschans_df, !!sym(msrun_idcol)),
                                                        by=msrun_idcol)
    } else {
      warning("msdata$", quantobj_idents_dfname, " not found")
    }
  }
  res[[modelobjs_dfname]] <- modelobjs_df
  res$modelobjects <- dplyr::mutate(modelobjs_df,
                                    object_id = !!sym(paste0(modelobj, "_id")),
                                    object_label = !!sym(paste0(modelobj, "_label")),
                                    chunk = row_number())
  res[[quantobj_intensities_dfname]] <- quantobj_intensities_df
  if (verbose) message("Calculating MS channel (", mschan, ") ", quantobj, " statistics...")
  res[[paste0(mschan, "_", quantobj, "_stats")]] <- mschannel_quantobj_statistics(res)

  checkmate::assert_class(mscalib, "mscalib")
  res[[paste0(quantobj, "_mscalib")]] <- mscalib

  quantobj_labu_shift_name <- paste0(quantobj, "_labu_shift")
  res[[quantobj_labu_shift_name]] <- 0.95*log2(median(quantobj_intensities_df$intensity, na.rm=TRUE))
  if (verbose) message(quantobj_labu_shift_name, "=", res[[quantobj_labu_shift_name]])

  quantobj_labu_min_name <- paste0(quantobj, "_labu_min")
  res[[quantobj_labu_min_name]] <- quantile(quantobj_intensities_df$intensity, min_intensity_quantile, na.rm=TRUE) %>% log2() -
                                            res[[quantobj_labu_shift_name]] + min_intensity_offset
  if (verbose) message(quantobj_labu_min_name, "=", res[[quantobj_labu_min_name]])

  return (res)
}
