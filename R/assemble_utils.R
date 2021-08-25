#' @export
load_fit_chunk <- function(chunk, chunks = fit_files.df, chunks_path = fit_path,
                           load_vars = c("msglm_results", "results_info", "dims_info"),
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
combine_fit_chunk_section <- function(section, chunks, type, dims_info, report_var="msglm_results") {
  message('Combining ', type, ' reports for ', section, '...')
  obj_id_cols <- str_subset(colnames(dims_info$object), "_id$")
  #message('  object ID columns: ', paste0(obj_id_cols, collapse=", "))
  report_df <- dplyr::bind_rows(lapply(chunks, function(chunk) {
    frame <- chunk[[report_var]][[section]][[type]]
    if (!is.null(frame)) {
      obj_info = chunk$dims_info$object
      for (col in setdiff(obj_id_cols, colnames(frame))) {
        frame[[col]] <- rep_len(obj_info[[col]][[1]], nrow(frame))
      }
    }
    return (frame)
  }))
  # attach extra object information columns to the report
  if (length(report_df) > 0) {
    if (length(obj_id_cols) > 0) {
      objs_df <- bind_rows(lapply(chunks[sapply(chunks, function(chunk) is.data.frame(chunk$dims_info$object))],
                                  function(chunk) chunk$dims_info$object))
      obj_info_cols <- setdiff(colnames(objs_df), colnames(report_df))
      if (length(obj_info_cols) > 0) {
        message('  attaching object info: ', paste0(obj_info_cols, collapse=", "))
        report_df <- dplyr::left_join(report_df, dplyr::select_at(objs_df, c(obj_id_cols, obj_info_cols)), by=obj_id_cols)
      }
    } else {
      warning("No object id columns for the ", type, " reports for ", section)
    }
  } else {
    report_df <- NULL
  }
  return(report_df)
}

#' @export
combine_fit_chunks <- function(reports, type, report_var="msglm_results") {
  dims_info <- reports[[1]]$dims_info
  sections <- names(reports[[1]][[report_var]])
  res <- rlang::set_names(lapply(sections, combine_fit_chunk_section,
                                 reports, type, dims_info, report_var=report_var),
                          sections)
  res <- res[!sapply(res, is.null)]
  return(res)
}
