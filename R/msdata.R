# FIXME
#' @export
mschannel_quantobj_statistics <- function(msdata, verbose=FALSE) {
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
  if (is.na(mschan)) {
    stop("MS channel column not specified")
  } else if (!rlang::has_name(quantobj_intensities_df, mschan)) {
    stop("MS channel column msdata$", quantobj_intensities_dfname, "$", mschan, " not found")
  }
  groupat <- c(mschan)
  if (!is.na(mstag) && rlang::has_name(quantobj_intensities_df, mstag)) {
    groupat <- append(groupat, mstag)
  }
  if (!is.na(msrun) && rlang::has_name(quantobj_intensities_df, msrun)) {
    groupat <- append(groupat, msrun)
  }
  groupat <- unique(groupat)

  quantobj_intensities_df <- dplyr::select_at(quantobj_intensities_df,
      append(groupat, c(quantobj_idcol, "intensity"))) %>%
    dplyr::rename(quantobj_id = !!sym(quantobj_idcol)) %>%
    dplyr::filter(!is.na(quantobj_id))

  expandat <- append(groupat, "quantobj_id")
  intensity_stats_df <- tidyr::expand(quantobj_intensities_df, !!!syms(expandat)) %>%
    dplyr::left_join(quantobj_intensities_df, by=expandat) %>%
    #dplyr::group_by(quantobj_id, condition) %>% dplyr::filter(any(!is.na(intensity))) %>%
    dplyr::group_by_at(groupat) %>%
    summarize(log2_intensity.mean = mean(log2(intensity[!is.na(intensity)])),
              log2_intensity.median = median(log2(intensity[!is.na(intensity)])),
              log2_intensity.sd = sd(log2(intensity[!is.na(intensity)])),
              n = n(),
              n_missing = sum(is.na(intensity)), .groups="drop")

  quantobj_idents_dfname <- str_c(quantobj, "_idents")
  if (!rlang::has_name(msdata, quantobj_idents_dfname)) {
    if (verbose) warning("msdata contains no `", quantobj_idents_dfname,
                         "` data, using intensities instead")
    quantobj_idents_dfname <- quantobj_intensities_dfname
  }
  quantobj_idents_df <- msdata[[quantobj_idents_dfname]]
  if (!rlang::has_name(quantobj_idents_df, "ident_type")) {
    if (verbose) warning("msdata$", quantobj_idents_dfname,
                         " has no ident_type column, skipping identification statistics")
  } else if (!rlang::has_name(quantobj_idents_df, msrun)) {
    if (verbose) warning("msdata$", quantobj_idents_dfname,
                         " has no `", msrun, "` MS run column, skipping identification statistics")
  } else {
    quantobj_idents_df <- dplyr::select_at(quantobj_idents_df,
        c(quantobj_id=quantobj_idcol, msrun, "ident_type"))
    ident_stats_df <- dplyr::mutate(quantobj_idents_df,
        is_matching = ident_type %in% c("By matching", "MULTI-MATCH", "MULTI-MATCH-MSMS"),
        is_msms = ident_type %in% c("By MS/MS", "MULTI-MSMS", "ISO-MSMS", "MSMS", "MULTI-SECPEP")) %>%
      dplyr::group_by(!!sym(msrun)) %>%
      summarize(n_matching = n_distinct(quantobj_id[is_matching], na.rm = TRUE),
                n_msms = n_distinct(quantobj_id[is_msms], na.rm = TRUE),
                .group = "drop")
    intensity_stats_df <- dplyr::left_join(intensity_stats_df, ident_stats_df, by=msrun)
  }
  return (intensity_stats_df)
}

#' Collect MS data that would be used by MSGLM models.
#'
#' @param msdata list of various MS-related data frames
#' @param model_def *msglm_model* object
#' @param object type of the MS object, which abundance would be modeled by MSGLM
#'               (protein group, PTM group etc)
#' @param quantobject type of the object, for which the MS-derived quantifications
#'               are provided. Defaults to the same type as model `object`,
#'               but alternatively could be a MS1 peak intensity (*pepmodstate*)
#' @param mscalib *mscalib* noise model for quantifications (at *quantobject* level)
#' @param verbose if TRUE, produces debugging output
#'
#' @return *msglm_data_collection* object with the data for all potential model objects
#' @export
#'
#' @examples
import_msglm_data <- function(msdata, model_def = NULL,
                              object = c("protgroup", "protregroup", "ptmngroup"),
                              object_cols = NULL,
                              quantobject = object,
                              mscalib = get(paste0(quantobject, "_mscalib")),
                              condition = NA_character_,
                              msexperiment = NA_character_,
                              msprobe = NA_character_,
                              mschannel = NA_character_,
                              msrun = NA_character_,
                              msfraction = NA_character_,
                              mstag = NA_character_,
                              msprotocol = NA_character_,
                              min_intensity_quantile = 0.001,
                              min_intensity_offset = -5,
                              verbose = FALSE)
{
  checkmate::assert_list(msdata)
  res <- structure(list(), class="msglm_data_collection")
  object <- match.arg(object)
  obj_idcol <- paste0(object, "_id")
  checkmate::assert_choice(quantobject, c(object, "pepmodstate"))
  quantobj_idcol <- paste0(quantobject, "_id")
  if (verbose) message("Importing MS data for ", object,
                       "s using ", quantobject, " intensities")

  res$msentities_extra_columns <- list()

  objs_dfname <- paste0(object, "s")
  if (rlang::has_name(msdata, objs_dfname)) {
    objs_df <- msdata[[objs_dfname]]
    if (verbose) message("Found ", nrow(objs_df), " ", object, "(s)")
  } else {
    stop("msdata$", objs_dfname, " not found")
  }
  if (is.null(object_cols)) {
    if (verbose) message('Guessing ', object, ' columns to keep in the report')
    object_cols <- append(object_cols,
        dplyr::select(objs_df, any_of(c("majority_protein_acs", "protein_acs",
                       "gene_names", "protein_names"))) %>% colnames())
    if (object %in% c('ptmngroup')) {
      object_cols <- append(object_cols,
          dplyr::select(objs_df, any_of(c("ptm_type", "nselptms", "ptm_pos", "ptm_AA_seq"))) %>%
          colnames())
    }
    object_cols <- append(object_cols,
        dplyr::select(objs_df, starts_with("is_")) %>% colnames())
    object_cols <- append(object_cols,
        dplyr::select(objs_df, ends_with("_label")) %>% colnames())
  }
  if (verbose) message("Model object (", object, ") columns to use in fit results report: ",
                       paste0(object_cols, collapse=", "))
  res$msentities_extra_columns$object <- object_cols

  if (quantobject != object) {
    quantobjs_dfname <- paste0(quantobject, "s")
    if (rlang::has_name(msdata, quantobjs_dfname)) {
      quantobjs_df <- msdata[[quantobjs_dfname]]
      if (verbose) message("Found ", nrow(quantobjs_df), " ", quantobject, "(s)")
    } else {
      stop("msdata$", quantobjs_dfname, " not found")
    }
  } else {
    quantobjs_dfname <- objs_dfname
    quantobjs_df <- objs_df
  }

  # autodetect/check msglm entities
  if (is.na(condition)) {
    condition <- "condition"
    if (verbose) message("No condition column specified, defaulting to condition=", condition)
  }

  msrun_def <- dplyr::coalesce(msrun, "msrun")
  msruns_def_dfname <- paste0(msrun_def, "s")
  if (!is.na(msexperiment)) {
    msexp <- msexperiment
  } else if (rlang::has_name(msdata, "msexperiments")) {
    if (verbose) message("Detected msdata$msexperiments frame")
    msexp <- "msexperiment"
  } else if (is.na(msfraction) && rlang::has_name(msdata, msruns_def_dfname)
         && !rlang::has_name(msdata[[msruns_def_dfname]], "msfraction")) {
    if (verbose) message("Detected msdata$", msruns_def_dfname,
                         ", no msfractions. Setting msexperiment=", msrun_def)
    msexp <- msrun_def
  } else {
    msexp <- NA_character_
  }

  mschan_def <- dplyr::coalesce(mschannel, "mschannel")
  mschans_def_dfname <- paste0(mschan_def, "s")
  if (!is.na(msprobe)) {
    msprb <- msprobe
  } else if (rlang::has_name(msdata, "msprobes")) {
    if (verbose) message("Detected msdata$msprobes frame")
    msprb <- "msprobe"
  } else if (is.na(mstag) && !is.na(msexp) && rlang::has_name(msdata, paste0(msexp, "s"))) {
    if (verbose) message("No msprobe and mstag columns specified, defaulting to msprobe=", msexp)
    msprb <- msexp
  } else {
    msprb <- NA_character_
  }
  if (!is.na(msprb) && rlang::has_name(msdata, paste0(msprb, "s"))) {
    msprbs_dfname <- paste0(msprb, "s")
    if (verbose) message("Using msdata$", msprbs_dfname, " as MS probes source")
    msprbs_df <- msdata[[msprbs_dfname]]
  } else {
    if (rlang::has_name(msdata, mschans_def_dfname)) {
      msprb_src_dfname <- mschans_def_dfname
      msprb_src <- mschan_def
    } else if (rlang::has_name(msdata, msruns_def_dfname)) {
      msprb_src_dfname <- msruns_def_dfname
      msprb_src <- msrun_def
    } else {
      stop("Cannot autodetect MS probes/experiments/channels/runs frame. ",
           "Make sure msprobe/msexperiment/mschannel/msrun= parameter set correctly and corresponding frames exist.")
    }
    msfrac <- dplyr::coalesce(msfraction, "msfraction")
    if (rlang::has_name(msdata[[msprb_src_dfname]], msfrac)) {
      if (verbose) message("Generating MS probes from msdata$", msprb_src_dfname,
                           " by removing MS fractions")
      msprbs_src_df <- msdata[[msprb_src_dfname]]
      checkmate::assert_data_frame(msprbs_src_df, .var.name = paste0("msdata$", msprb_src_dfname))
      if (is.na(msprb)) {
        if (rlang::has_name(msprbs_src_df, "msprobe")) {
          msprb <- "msprobe"
        } else if (!is.na(msexp) && is.na(mstag) && rlang::has_name(msprbs_src_df, msexp)) {
          msprb <- msexp
        } else if (is.na(msexp) && is.na(mstag) && rlang::has_name(msprbs_src_df, "msexperiment")) {
          msprb <- "msexperiment"
        } else {
          msprb <- if_else(is.na(mstag) && msprb_src_dfname == msruns_def_dfname, "msexperiment", "msprobe")
        }
      }
      if (!rlang::has_name(msprbs_src_df, msprb)) {
        # generate msprobe column in the source data frame
        msprb_key_cols <- intersect(c('condition', 'replicate', mstag), colnames(msprbs_src_df))
        if (verbose) message("Trying to generate ", msprb, " column using: ",
                             paste0(msprb_key_cols, collapse=", "))
        msprbs_df <- dplyr::distinct_at(msprbs_src_df, msprb_key_cols, .keep_all = TRUE) %>%
          dplyr::mutate(!!sym(msprb) := paste(!!!syms(msprb_key_cols), sep='_'))
        # update the source frame with msprb column
        msprb_src_df <- dplyr::inner_join(msdata[[msprb_src_dfname]],
            dplyr::select_at(msprbs_df, c(msprb, msprb_key_cols)), by=msprb_key_cols)
        if (nrow(msprb_src_df) != nrow(msdata[[msprb_src_dfname]])) {
          stop("MSGLM failed to detect columns uniquely identifying ", msprb, ".",
               "Plese add ", msprb, " column to your msdata$", msprb_src_dfname,
               " experimental design data frame")
        }
        msdata[[msprb_src_dfname]] <- msprb_src_df
      } else {
        msprbs_df <- msprbs_src_df
      }
      # try to exclude columns that vary between msfractions
      msprbs_df <- dplyr::select(msprbs_df, -!!sym(msprb_src), -!!sym(msfrac)) %>%
          dplyr::select_if(~n_distinct(paste0(.x, "_", msprbs_df[[msprb]])) == n_distinct(msprbs_df[[msprb]])) %>%
          dplyr::distinct(!!sym(msprb), .keep_all=TRUE)
    } else {
      if (is.na(msprb)) {
        if (verbose) message("Detected msdata$", msprb_src_dfname, ", no msfractions. Setting msprobe=", msprb_src_def)
        msprb <- msprb_src
      }
      msprbs_df <- msdata[[msprb_src_dfname]]
    }
    msprbs_dfname <- paste0(msprb, "s")
  }

  checkmate::assert_data_frame(msprbs_df, .var.name = paste0("msdata$", msprbs_dfname))
  msprb_req_cols <- c(condition, msprb)
  if (!is.na(msexp)) msprb_req_cols <- unique(c(msprb_req_cols, msexp))
  if (!is.na(mstag)) msprb_req_cols <- c(msprb_req_cols, mstag)
  checkmate::assert_names(names(msprbs_df), must.include = msprb_req_cols)
  if (is.na(mstag)) {
    if (rlang::has_name(msprbs_df, "mstag")) {
      if (!is.na(msexp) && (msexp == msprb)) {
        stop("msexperiment=msprobe=", msprb, " so assuming unlabeled MS data, but msdata$",
            msprbs_dfname, "$", mstag, " column detected. Check the correctness of your MS entities")
      } else {
        mstag <- "mstag"
        if (verbose) message("Detected msdata$", msprbs_dfname, "$", mstag, " column, assuming MS tagged data")
      }
    } else if (is.na(msexp)) {
      if (verbose) message("No MS experiment column specified, assuming msexperiment=", msprb)
      msexp <- msprb
    }
  }
  if (is.na(msfraction) && rlang::has_name(msprbs_df, "msfraction")) {
    warning("msdata$", msprbs_dfname, " contains msfraction column. Check the correctness of your MS entities")
  }
  checkmate::assert_character(as.character(msprbs_df[[msprb]]), unique=TRUE, any.missing = FALSE,
                              .var.name = paste0("msdata$", msprbs_dfname, "$", msprb))

  if (msprb == mschan_def) {
    mschan = msprb
    mschans_dfname <- msprbs_dfname
  } else if (rlang::has_name(msdata, mschans_def_dfname)) {
    mschan <- mschan_def
    mschans_dfname <- mschans_def_dfname
  } else if (is.na(mstag) && rlang::has_name(msdata, msruns_def_dfname)) {
    mschan <- msrun_def
    mschans_dfname <- msruns_def_dfname
  } else if (!is.na(msfraction)) {
    stop("MS entities declared to contain fractions (", msfraction,
          " column), but no MS channels/MS runs frame found")
  } else if (is.na(mschannel) && is.na(msfraction) && rlang::has_name(msprbs_df, "raw_file")) {
    mschan <- msprb
    mschans_dfname <- msprbs_dfname
  } else {
    stop("Cannot autodetect MS channels data, specify mschannel= argument")
  }
  checkmate::assert_character(mschans_dfname, any.missing = FALSE)
  if (verbose) message(if_else(mschans_dfname == msprbs_dfname, "Reusing", "Using"),
                       " msdata$", mschans_dfname, " as MS channels source")
  mschans_df <- msdata[[mschans_dfname]]
  checkmate::assert_data_frame(mschans_df, .var.name = paste0("msdata$", mschans_dfname))

  mschan_required_cols <- c(msprb, mschan)
  if (!is.na(msfraction)) {
    if (mschan == msprb) stop("msfraction=", msfraction, " but msprobe=mschannel=", msprb)
    msfrac <- msfaction
  } else if (rlang::has_name(mschans_df, "msfraction")) {
    if (mschan == msprb) stop("MS fraction column (msfraction) detected, but msprobe=mschannel=", msprb)
    msfrac <- "msfraction"
  } else {
    msfrac <- NA_character_
  }
  if (!is.na(msfrac)) {
    mschan_required_cols <- c(mschan_required_cols, msfrac)
  }

  if (!is.na(mstag)) {
    if (!is.na(msrun) && (msrun == mschan)) {
      stop("Specifying mstag column (", mstag,
           ") also requires that MS channel (", mschan, ") and MS run (",
           msrun, ") columns are different")
    }
  } else if (is.na(msrun)) {
    if (verbose) message("Assuming msrun=mschannel=", mschan)
    msrun <- mschan
  }
  if (!is.na(msprotocol)) {
    msproto <- msprotocol
    if (verbose) message("Using MS protocol column ", msproto)
  } else if (rlang::has_name(mschans_df, "msprotocol")) {
    msproto <- "msprotocol"
    if (verbose) message("Detected MS protocol column ", msproto)
  } else {
    msproto <- NA_character_
  }
  if (!is.na(msproto)) {
    mschan_required_cols <- append(mschan_required_cols, msproto)
  }
  checkmate::assert_names(names(mschans_df), must.include = unique(mschan_required_cols),
                          .var.name = paste0("msdata$", mschans_dfname))

  res$msentities <- c(
    object = object,
    quantobject = quantobject,
    condition = condition,
    msexperiment = msexp,
    msprobe = msprb,
    msrun = msrun,
    mschannel = mschan,
    mstag = mstag,
    msfraction = msfrac,
    msprotocol = msproto)

  msprb_idcol <- msprb
  msexp_idcol <- msexp
  mschan_idcol <- mschan
  msrun_idcol <- msrun
  if (rlang::has_name(msprbs_df, "is_used")) {
    checkmate::assert_logical(msprbs_df$is_used, any.missing=FALSE,
                              .var.name=paste0("msdata$", msexp, "s$is_used"))
    if (verbose) message("Importing ", sum(msprbs_df$is_used), " used ", msprb,
                         "(s) of ", nrow(msprbs_df))
    msprbs_df <- dplyr::filter(msprbs_df, is_used)
  } else {
    if (verbose) message("Importing all ", nrow(msprbs_df), " MS probe(s)")
  }
  if (!is.null(model_def)) {
    checkmate::assert_class(model_def, "msglm_model")
    virtual_conditions <- setdiff(unique(msprbs_df[[condition]]), model_def$conditions$condition)
    if (length(virtual_conditions) > 0) {
      warning(length(virtual_conditions), ' condition(s) not covered by MS experiments')
    }
    checkmate::assert_subset(as.character(msprbs_df[[condition]]), as.character(model_def$conditions$condition),
                             .var.name = paste0("msdata$", msprbs_dfname, "$", condition))
  }
  res[[msprbs_dfname]] <- msprbs_df
  if (mschan != msprb) {
    mschans_df <- dplyr::semi_join(msdata[[mschans_dfname]],
                                   dplyr::select(msprbs_df, !!sym(msprb_idcol)),
                                   by=msprb_idcol)
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
    checkmate::assert_set_equal(mschans_df[[msprb_idcol]], msprbs_df[[msprb_idcol]],
                                .var.name = paste0("msdata$", mschans_dfname, "$", msprb_idcol))
  } else {
    mschans_df <- msprbs_df
  }

  quantobj_intensities_dfname <- paste0(quantobject, "_intensities")
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
  obj_idents_dfname <- paste0(object, "_idents")
  if (!rlang::has_name(msprbs_df, msexp_idcol)) {
    if (verbose) message("msdata$", msprbs_dfname, " does not have ", msexp_idcol,
                         " ignoring ", object, "-level identifications")
  } else if (!rlang::has_name(msdata, obj_idents_dfname)) {
    if (verbose) message("msdata$", obj_idents_dfname, " not found")
  } else {
    obj_idents_df <- msdata[[obj_idents_dfname]]
    if (!checkmate::test_data_frame(obj_idents_df)) {
      if (verbose) warning("msdata$", obj_idents_dfname, " is not a data frame, ignoring")
    } else {
      col_check <- checkmate::check_names(colnames(obj_idents_df),
                        must.include = c(msexp_idcol, obj_idcol, "ident_type"),
                        what = "column names")
      if (is.logical(col_check) && col_check) {
        res[[obj_idents_dfname]] <- dplyr::semi_join(obj_idents_df,
                                                     dplyr::select(msprbs_df, !!sym(msexp_idcol)),
                                                     by=msexp_idcol)
      } else {
        if (verbose) warning("Detected msdata$", obj_idents_dfname, " problems, skipping: ", col_check)
      }
    }
  }
  if (object == quantobject) {
    orig_nintensities <- nrow(quantobj_intensities_df)
    quantobj_intensities_df <- dplyr::filter(quantobj_intensities_df, !is.na(intensity))
    if (verbose) message(nrow(quantobj_intensities_df), " of ", orig_nintensities,
                         " ", quantobject, " intensity measurements selected")
    orig_nobjs <- nrow(objs_df)
    objs_df <- dplyr::semi_join(objs_df, quantobj_intensities_df, by=obj_idcol)
    if (verbose) message(nrow(objs_df), " of ", orig_nobjs,
                         " ", object, "s selected")
  } else if (quantobject == "pepmodstate") {
    obj2quantobj_dfname = paste0(object, "2", quantobject)
    if (!rlang::has_name(msdata, obj2quantobj_dfname)) {
      obj2pepmod_dfname = paste0(object, "2pepmod")
      if (verbose) warning("msdata$", obj2quantobj_dfname,
                           " not found, trying msdata$", obj2pepmod_dfname)
      if (!rlang::has_name(msdata, obj2pepmod_dfname)) {
        stop("msdata$", obj2pepmod_dfname, " not found")
      }
      obj2pepmod_df <- msdata[[obj2pepmod_dfname]]
      obj2quantobj_df <- dplyr::inner_join(obj2pepmod_df,
            dplyr::select_at(quantobjs_df, c(quantobj_idcol, "pepmod_id")), by="pepmod_id") %>%
        dplyr::select_at(c(obj_idcol, quantobj_idcol, "pepmod_id", "is_specific")) %>%
        dplyr::distinct()
    } else {
      obj2quantobj_df <- msdata[[obj2quantobj_dfname]]
    }
    checkmate::assert_data_frame(obj2quantobj_df,
                                 .var.name=paste0("msdata$", obj2quantobj_dfname))
    checkmate::assert_names(colnames(obj2quantobj_df),
                            must.include=c(obj_idcol, quantobj_idcol, "is_specific"),
                            .var.name = paste0("msdata$", obj2quantobj_dfname))
    if (verbose) message("Filtering for specific and quantified ", quantobject,
                         "s and corresponding ", object, "s...")
    quantobjs_orig_df <- quantobjs_df
    quantobj_intensities_orig_df <- quantobj_intensities_df
    objs_orig_df <- objs_df
    obj2quantobj_df <- dplyr::filter(obj2quantobj_df, is_specific)
    # TODO require identification in a selected set of conditions
    # (e.g. for AP-MS data should be identified in AP-MS conditions,
    #  whereas identification in the background is not required)
    quantobj_intensities_df <- dplyr::filter(quantobj_intensities_df, !is.na(intensity)) %>%
        dplyr::semi_join(obj2quantobj_df, by=quantobj_idcol)
    if (verbose) message("  * ", nrow(quantobj_intensities_df), " of ", nrow(quantobj_intensities_orig_df),
                         " specific ", quantobject, " intensity measurements selected")
    quantobjs_df <- dplyr::semi_join(quantobjs_df, quantobj_intensities_df, by=quantobj_idcol)
    if (verbose) message("  * ", nrow(quantobjs_df), " of ", nrow(quantobjs_orig_df),
                         " specific ", quantobject, "s selected")
    obj2quantobj_df <- dplyr::semi_join(obj2quantobj_df, quantobjs_df, by=quantobj_idcol)
    objs_df <- dplyr::semi_join(objs_df, obj2quantobj_df, by=obj_idcol)
    if (verbose) message("  * ", nrow(objs_df), " of ", nrow(objs_orig_df),
                         " ", object, "s have specific ", quantobject, " quantitations")
    res[[obj2quantobj_dfname]] <- obj2quantobj_df
    res[[quantobjs_dfname]] <- quantobjs_df

    quantobj_idents_dfname <- paste0(quantobject, "_idents")
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
      if (verbose) warning("msdata$", quantobj_idents_dfname, " not found")
    }
  }
  res[[objs_dfname]] <- objs_df
  res$objects <- dplyr::mutate(objs_df,
                               object_id = !!sym(obj_idcol),
                               object_label = !!sym(paste0(object, "_label")),
                               chunk = row_number())
  res[[quantobj_intensities_dfname]] <- quantobj_intensities_df
  if (verbose) message("Calculating MS channel (", mschan, ") ", quantobject, " statistics...")
  res[[paste0(mschan, "_", quantobject, "_stats")]] <- mschannel_quantobj_statistics(res, verbose=verbose)

  checkmate::assert_class(mscalib, "mscalib")
  res[[paste0(quantobject, "_mscalib")]] <- mscalib

  quantobj_labu_min_name <- paste0(quantobject, "_labu_min")
  res[[quantobj_labu_min_name]] <- quantile(quantobj_intensities_df$intensity, min_intensity_quantile, na.rm=TRUE) %>% log2() -
                                            mscalib$zShift + min_intensity_offset
  if (verbose) message(quantobj_labu_min_name, "=", res[[quantobj_labu_min_name]])

  return (res)
}
