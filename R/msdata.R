#' @export
mschannel_statistics <- function(msdata, obj="protgroup") {
    objs_dfname <- str_c(obj, "s")
    obj_idcol = str_c(obj, "_id")
    if (!rlang::has_name(msdata, objs_dfname)) {
        stop("msdata contains no `", objs_dfname, " data")
    }
    objs_df = msdata[[objs_dfname]]
    objs_df$object_id = objs_df[[obj_idcol]]
    obj_tagintensities_dfname <- str_c(obj, "_tagintensities")
    if (!rlang::has_name(msdata, obj_tagintensities_dfname)) {
      stop("msdata contains no `", obj, "_tagintensities` data, msrun_statistics() should be used instead")
    }
    obj_tagintensities_df <- msdata[[obj_tagintensities_dfname]]
    obj_tagintensities_df$object_id <- obj_tagintensities_df[[obj_idcol]]
    res <- dplyr::left_join(tidyr::expand(dplyr::filter(obj_tagintensities_df, !is.na(object_id)),
                                          object_id, msrun, mstag),
                            dplyr::filter(obj_tagintensities_df, !is.na(object_id)) %>%
                            dplyr::select(object_id, msrun, mstag, intensity)) %>%
    dplyr::inner_join(dplyr::select(msdata$mschannels, msrun, mstag, mschannel) %>%
                      dplyr::distinct()) %>%
    #dplyr::group_by(protgroup_id, condition) %>% dplyr::filter(any(!is.na(intensity))) %>%
    dplyr::group_by(mschannel, mstag, msrun) %>%
    summarize(log2_intensity.mean = mean(log2(intensity[!is.na(intensity)])),
              log2_intensity.median = median(log2(intensity[!is.na(intensity)])),
              log2_intensity.sd = sd(log2(intensity[!is.na(intensity)])),
              n = n(),
              n_missing = sum(is.na(intensity))) %>%
    dplyr::ungroup()

    obj_idents_dfname <- str_c(obj, "_idents")
    if (rlang::has_name(msdata, obj_idents_dfname)) {
      obj_idents_df <- msdata[[obj_idents_dfname]]
      obj_idents_df$object_id <- obj_idents_df[[obj_idcol]]
    } else {
      obj_idents_df <- obj_tagintensities_df
    }
    if (!is.null(obj_idents_df)) {
      ident_stats <- dplyr::left_join(dplyr::filter(obj_idents_df, !is.na(object_id)),
                                      dplyr::select(msdata$mschannels, msrun, any_of("msrun_mq")) %>% dplyr::distinct()) %>%
        dplyr::mutate(is_matching = ident_type %in% c("By matching", "MULTI-MATCH", "MULTI-MATCH-MSMS"),
                      is_msms = ident_type %in% c("By MS/MS", "MULTI-MSMS", "ISO-MSMS", "MSMS", "MULTI-SECPEP")) %>%
        dplyr::group_by(msrun) %>%
        summarize(n_matching = n_distinct(object_id[is_matching], na.rm = TRUE),
                  n_msms = n_distinct(object_id[is_msms], na.rm = TRUE)) %>%
        dplyr::ungroup()
      res <- left_join(res, ident_stats)
    } else {
      warning("No protgroup ident_type data found")
    }
    return (res)
}

#' @export
msrun_statistics <- function(msdata, obj="protgroup") {
  objs_dfname <- str_c(obj, "s")
  obj_idcol = str_c(obj, "_id")
  if (!rlang::has_name(msdata, objs_dfname)) {
    stop("msdata contains no `", objs_dfname, " data")
  }
  objs_df = msdata[[objs_dfname]]
  objs_df$object_id = objs_df[[obj_idcol]]
  if (rlang::has_name(msdata, str_c(obj, "_tagintensities"))) {
    warning("msdata contains `", obj, "_tagintensities` data, mschannel_statistics() should be used instead")
  }
  obj_intensities_dfname <- str_c(obj, "_intensities")
  if (!rlang::has_name(msdata, obj_intensities_dfname)) {
    warning("msdata contains no `", obj_intensities_dfname, "` data")
  }
  obj_intensities_df = msdata[[obj_intensities_dfname]]
  obj_intensities_df$object_id = obj_intensities_df[[obj_idcol]]
  res <- dplyr::left_join(tidyr::expand(dplyr::filter(obj_intensities_df, !is.na(object_id)),
                                        object_id, msrun),
                          dplyr::filter(obj_intensities_df, !is.na(object_id)) %>%
                          dplyr::select(object_id, msrun, intensity)) %>%
    dplyr::inner_join(dplyr::select(msdata$msruns, msrun) %>% dplyr::distinct()) %>%
    #dplyr::group_by(protgroup_id, condition) %>% dplyr::filter(any(!is.na(intensity))) %>%
    dplyr::group_by(msrun) %>%
    summarize(log2_intensity.mean = mean(log2(intensity[!is.na(intensity)])),
              log2_intensity.median = median(log2(intensity[!is.na(intensity)])),
              log2_intensity.sd = sd(log2(intensity[!is.na(intensity)])),
              n = n(),
              n_missing = sum(is.na(intensity))) %>%
    dplyr::ungroup()

  obj_idents_dfname <- str_c(obj, "_idents")
  if (rlang::has_name(msdata, obj_idents_dfname)) {
    obj_idents_df <- msdata[[obj_idents_dfname]]
    obj_idents_df$object_id = obj_idents_df[[obj_idcol]]
  } else {
    obj_idents_df <- obj_intensities_df
  }
  if (!is.null(obj_idents_df)) {
    ident_stats <- dplyr::left_join(dplyr::filter(obj_idents_df, !is.na(object_id)),
                                    dplyr::select(msdata$msruns, msrun, any_of("msrun_mq")) %>% dplyr::distinct()) %>%
      dplyr::mutate(is_matching = ident_type %in% c("By matching", "MULTI-MATCH", "MULTI-MATCH-MSMS"),
                    is_msms = ident_type %in% c("By MS/MS", "MULTI-MSMS", "ISO-MSMS", "MSMS", "MULTI-SECPEP")) %>%
      dplyr::group_by(msrun) %>%
      summarize(n_matching = n_distinct(object_id[is_matching], na.rm = TRUE),
                n_msms = n_distinct(object_id[is_msms], na.rm = TRUE)) %>%
      dplyr::ungroup()
    res <- left_join(res, ident_stats)
  } else {
    warning("No protgroup ident_type data found")
  }
  return (res)
}

#' Collect MS data that would be used by MSGLM models.
#'
#' @param msdata list of various MS-related data frames
#' @param model_def *msglm_model* object
#' @param quantobj_mscalib *mscalib* MS intensities noise model
#' @param verbose if TRUE, produces debugging output
#'
#' @return *msglm_data_collection* object with the data for all potential model objects
#' @export
#'
#' @examples
import_msglm_data <- function(msdata, model_def,
                              quantobj_mscalib = get(paste0(model_def$quantobject, "_mscalib")),
                              min_intensity_quantile = 0.001,
                              min_intensity_offset = -5,
                              verbose=model_def$verbose)
{
  checkmate::assert_class(model_def, "msglm_model")
  checkmate::assert_list(msdata)
  res <- list()
  modelobj <- model_def$modelobject
  modelobj_idcol <- paste0(modelobj, "_id")
  quantobj <- model_def$quantobject
  quantobj_idcol <- paste0(quantobj, "_id")

  # FIXME support experiments with mschannels
  if (rlang::has_name(msdata, "mschannels")) {
    if (verbose) message("Using msdata$mschannels for experimental design")
    checkmate::assert_tibble(msdata$mschannels)
    checkmate::assert_names(names(msdata$mschannels),
                            must.include = c("mschannel", "msrun", "mstag", "condition"))
    msexp <- "mschannel"
  } else if (rlang::has_name(msdata, "msruns")) {
    if (verbose) message("Using msdata$msruns for experimental design")
    checkmate::assert_tibble(msdata$msruns)
    checkmate::assert_names(names(msdata$msruns), must.include = c("msrun", "condition"))
    msexp <- "msrun"
  } else {
    stop("Cannot find MS experiments information (msdata$mschannels or msdata$msruns)")
  }
  msexps_dfname <- paste0(msexp, "s")
  msexps_df <- msdata[[msexps_dfname]]
  msexp_idcol <- msexp
  if (rlang::has_name(msexps_df, "is_used")) {
    checkmate::assert_logical(msexps_df$is_used, any.missing=FALSE,
                              .var.name=paste0("msdata$", msexp, "s$is_used"))
    if (verbose) message("Importing ", sum(msexps_df$is_used), " used ", msexp,
                         "(s) of ", nrow(msexps_df))
    msexps_df <- dplyr::filter(msexps_df, is_used)
  } else {
    if (verbose) message("Importing all ", nrow(msexps_df), " ", msexp, "(s)")
  }
  print(setdiff(unique(msexps_df$condition), model_def$conditions$condition))
  checkmate::assert_subset(as.character(msexps_df$condition), as.character(model_def$conditions$condition),
                           .var.name = paste0("msdata$", msexps_dfname, "$condition"))
  res[[msexps_dfname]] <- msexps_df

  quantobj_intensities_dfname <- paste0(quantobj, "_intensities")
  if (rlang::has_name(msdata, quantobj_intensities_dfname)) {
    if (verbose) message("Importing ", quantobj_intensities_dfname, "...")
    quantobj_intensities_df <- msdata[[quantobj_intensities_dfname]]
    checkmate::assert_data_frame(quantobj_intensities_df)
    checkmate::assert_names(colnames(quantobj_intensities_df),
                            must.include = c(msexp_idcol, quantobj_idcol, "intensity"),
                            .var.name = paste0("msdata$", quantobj_intensities_dfname))
    quantobj_intensities_df <- dplyr::semi_join(quantobj_intensities_df,
                                                dplyr::select(msexps_df, !!sym(msexp_idcol)),
                                                by=msexp_idcol)
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
    checkmate::assert_data_frame(msdata[[modelobj_idents_dfname]])
    res[[modelobj_idents_dfname]] <- dplyr::semi_join(msdata[[modelobj_idents_dfname]],
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
  }
  res[[modelobjs_dfname]] <- modelobjs_df
  res$modelobjects <- dplyr::mutate(modelobjs_df,
                                    object_id = !!sym(paste0(modelobj, "_id")),
                                    object_label = !!sym(paste0(modelobj, "_label")),
                                    chunk = row_number())
  res[[quantobj_intensities_dfname]] <- quantobj_intensities_df
  res[[paste0(msexp, "_", quantobj, "_stats")]] <-
    do.call(paste0(msexp, "_statistics"), list(res, obj=quantobj))

  checkmate::assert_class(quantobj_mscalib, "mscalib")
  res[[paste0(quantobj, "_mscalib")]] <- quantobj_mscalib

  quantobj_labu_shift_name <- paste0(quantobj, "_labu_shift")
  res[[quantobj_labu_shift_name]] <- 0.95*log2(median(quantobj_intensities_df$intensity, na.rm=TRUE))
  if (verbose) message(quantobj_labu_shift_name, "=", res[[quantobj_labu_shift_name]])

  quantobj_labu_min_name <- paste0(quantobj, "_labu_min")
  res[[quantobj_labu_min_name]] <- quantile(quantobj_intensities_df$intensity, min_intensity_quantile, na.rm=TRUE) %>% log2() -
                                            res[[quantobj_labu_shift_name]] + min_intensity_offset
  if (verbose) message(quantobj_labu_min_name, "=", res[[quantobj_labu_min_name]])

  return (structure(res, class="msglm_data_collection"))
}
