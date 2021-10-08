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
