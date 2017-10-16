process_msglm_chunk <- function(file_ix,
                                strip_samples=FALSE, strip_stats=FALSE,
                                postprocess.f = NULL) {
  fit_file <- fit_files.df[file_ix, 'filename']
  message( 'Loading ', fit_file, '...' )
  tmp.env <- new.env(parent = baseenv())
  load(file.path(fit_path, fit_file), envir = tmp.env)
  if (!is.null(postprocess.f)) {
    postprocess.f(envir = tmp.env, fit_file)
  }
  if (strip_samples || strip_stats) {
    tmp.env$msglm_results <- lapply( tmp.env$msglm_results, function( var_results ) {
      if (strip_samples) { var_results$samples <- NULL }
      if (strip_stats) { var_results$stats <- NULL }
      return ( var_results )
    } )
  }
  gc()
  res <- mget(ls(envir=tmp.env), envir=tmp.env)
  return ( res )
}

join_report_frames <- function(reports, frame_extractor=function(report) stop('no frame extractor'),
                               global_vars = c('model_dataset', 'version', 'chunk'),
                               prot_info=NULL)
{
  if (is.null(names(reports))) {
    report_names <- seq_along(reports)
  } else {
    report_names <- names(reports)
  }
  joined_res.df <- bind_rows(lapply(report_names, function(report_name) {
    report <- reports[[report_name]]
    frame <- frame_extractor(report)
    if (is.null(frame)) return ( NULL )
    if (nrow(frame) > 0L) {
      frame$report_name <- report_name
    } else {
      frame$report_name <- character()
    }
    for (var in global_vars) {
      frame[[var]] <- report[[var]]
    }
    return ( frame )
  } ) )
  # FIXME other object types?
  if ( !is.null(prot_info) && 'protein_ac_noiso' %in% colnames( joined_res.df ) ) {
    joined_res.df <- dplyr::inner_join( joined_res.df, prot_info[,c('protein_ac_noiso','protein_label','description')] )
  }
  return ( joined_res.df )
}

join_msglm_reports <- function(section, reports, type, results_tag="msglm_results") {
  message( 'Assembling joint ', type, ' report for ', section, '...' )
  res <- join_report_frames(reports, frame_extractor = function(report) {
    res <- report[[results_tag]][[section]][[type]]
    # add object id columns
    # FIXME sites support?
    # FIXME model_data$objects support?
    if (!is.null(res) && !('majority_protein_acs' %in% colnames(res))
        && "majority_protein_acs" %in% colnames(report$model_data$objects)) {
      res$majority_protein_acs <- rep_len(report$model_data$objects$majority_protein_acs[1], nrow(res))
      res$protgroup_id <- rep_len(report$model_data$objects$protgroup_id[1], nrow(res))
    }
    if (!is.null(res) && !('pepmod_id' %in% colnames(res))
        && "pepmod_id" %in% colnames(report$model_data$objects)) {
      res$pepmod_id <- rep_len(report$model_data$objects$pepmod_id[1], nrow(res))
    }
    res
  },
  NULL, global_vars = c())
  return(res)
}

join_msglm_reports_allsections <- function(reports, type, results_tag="msglm_results") {
  sections <- names(reports[[1]][[results_tag]])
  res <- lapply(sections, join_msglm_reports, reports, type, results_tag=results_tag)
  names(res) <- sections
  return(res)
}
