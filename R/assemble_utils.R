#' @export
load_fit_chunk <- function(chunk, chunks = fit_files.df, chunks_path = fit_path,
                           load_vars = c("msglm_results", "results_info"),
                           process.f = NULL) {
  chunk_file <- chunks[chunk, 'filename']
  message('Loading ', chunk_file, '...')
  tmp.env <- new.env(parent = baseenv())
  load(file.path(chunks_path, chunk_file), envir = tmp.env)
  if (!is.null(process.f)) {
    process.f(envir = tmp.env, chunk_file)
  }
  gc()
  res <- mget(load_vars, envir=tmp.env)
  return ( res )
}

#' @export
combine_fit_chunk_section <- function(section, chunks, type, objs_df, report_var="msglm_results") {
  message('Combining ', type, ' reports for ', section, '...')
  report_df <- purrr::map2_df(seq_along(chunks), chunks, function(chunk_ix, chunk) {
    chunk <- chunks[[chunk_ix]]
    frame <- chunk[[report_var]][[section]][[type]]
    if (!is.null(frame)) {
      frame$index_chunk <- chunk_ix
    }
    return (frame)
  })
  # attach extra object information columns to the report
  if (length(report_df) > 0) {
    if (!is.null(objs_df)) {
      extra_info_cols <- setdiff(colnames(objs_df), c(colnames(report_df), "index_chunk"))
      if (length(extra_info_cols) > 0) {
        message('  attaching object info: ', paste0(extra_info_cols, collapse=", "))
        report_df <- dplyr::left_join(report_df,
                                      dplyr::select_at(objs_df, c(extra_info_cols, "index_chunk")),
                                      by="index_chunk")
      }
    }
    report_df$index_chunk <- NULL
  } else {
    report_df <- NULL
  }
  return(report_df)
}

#' @export
combine_fit_chunks <- function(chunks, type, report_var="msglm_results") {
  hasobj_mask <- purrr::map_lgl(chunks, ~is.data.frame(.[[report_var]]$objects$stats))
  if (any(hasobj_mask)) {
    report1 <- chunks[[which(hasobj_mask)[[1]]]][[report_var]]
    # detect columns with object metadata
    obj_info_cols <- str_subset(setdiff(colnames(report1$objects$stats), colnames(report1$global$stats)),
                                "^is_|_(id|label)$")
    # detect object ID columns
    obj_id_cols <- str_subset(obj_info_cols, "_id$")
    #message('  object ID columns: ', paste0(obj_id_cols, collapse=", "))

    objs_df <- purrr::map2_df(chunks[hasobj_mask], which(hasobj_mask),
                              ~ dplyr::select_at(.x[[report_var]]$objects$stats, obj_info_cols) %>%
                                dplyr::distinct() %>% dplyr::mutate(index_chunk = .y))
    if (nrow(objs_df) != sum(hasobj_mask)) {
      warning("Problems extracting objects information")
      objs_df <- NULL
    }
  } else {
    report1 <- chunks[[1]][[report_var]]
    objs_df <- NULL
  }

  sections <- names(report1)
  res <- rlang::set_names(lapply(sections, combine_fit_chunk_section,
                                 chunks, type, objs_df, report_var=report_var),
                          sections)
  res <- res[!sapply(res, is.null)]
  return(res)
}
