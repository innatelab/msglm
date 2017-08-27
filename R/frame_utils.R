# Utilities for assembling Stan report frames
#
# Author: Alexey Stukalov
###############################################################################

calc_contrasts_subset <- function(vars_results, vars_info, dims_info,
                                  contrasts = NULL, condition_shifts = NULL,
                                  val_trans = NULL,
                                  condition.reported = "lhs",
                                  condition.quantiles_lhs = c(0, 1), condition.quantiles_rhs = c(0, 1),
                                  condition_agg_col = "condition")
{
    # FIXME do once per assembly
    sel_contrastXmetacondition <- contrastXmetacondition[contrasts, , drop=FALSE]
    sel_contrastXmetacondition <- sel_contrastXmetacondition[, colSums(abs(sel_contrastXmetacondition)) != 0, drop=FALSE]
    sel_conditionXmetacondition.df <- dplyr::filter(conditionXmetacondition.df,
                                                    metacondition %in% colnames(sel_contrastXmetacondition))
    if (!is.null(condition_shifts)) {
      sel_conditionXmetacondition.df <- dplyr::inner_join(sel_conditionXmetacondition.df, condition_shifts)
    }
    sel_contrastXcondition.df <- dplyr::semi_join(contrastXcondition.df, sel_conditionXmetacondition.df) %>%
        dplyr::filter(contrast %in% contrasts)
    calc_contrasts(vars_results, vars_info, dims_info,
                   sel_contrastXmetacondition, sel_conditionXmetacondition.df,
                   sel_contrastXcondition.df, val_trans = val_trans,
                   condition_agg_col = condition_agg_col,
                   condition.reported = condition.reported,
                   condition.quantiles_lhs = condition.quantiles_lhs,
                   condition.quantiles_rhs = condition.quantiles_rhs)
}

process_msglm_chunk <- function(file_ix,
                                strip_samples=FALSE, strip_stats=FALSE,
                                postprocess.f = NULL) {
  fit_model.file <- fit_model.files.df[ file_ix, 'filename' ]
  message( 'Loading ', fit_model.file, '...' )
  tmp.env <- new.env(parent = baseenv())
  load(file.path(fit_model_path, fit_model.file ), envir = tmp.env)
  if (!is.null(postprocess.f)) {
    postprocess.f(envir = tmp.env, fit_model.file)
  }
  if (strip_samples || strip_stats) {
    tmp.env$vars_results <- lapply( tmp.env$vars_results, function( var_results ) {
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
  if ( !is.null(prot_info) && 'protein_ac_noiso' %in% colnames( joined_res.df ) ) {
    joined_res.df <- dplyr::inner_join( joined_res.df, prot_info[,c('protein_ac_noiso','protein_label','description')] )
  }
  return ( joined_res.df )
}

join_msglm_reports <- function(section, reports, type) {
  message( 'Assembling joint ', type, ' report for ', section, '...' )
  res <- join_report_frames(reports, frame_extractor = function(report) {
    res <- report$vars_results[[section]][[type]]
    if (!is.null(res) && !('majority_protein_acs' %in% colnames(res))
        && "majority_protein_acs" %in% colnames(report$model_data$protgroups)) {
      res$majority_protein_acs <- rep_len(report$model_data$protgroups$majority_protein_acs[1], nrow(res))
      res$protgroup_id <- rep_len(report$model_data$protgroups$protgroup_id[1], nrow(res))
    }
    if (!is.null(res) && !('pepmod_id' %in% colnames(res))
        && "pepmod_id" %in% colnames(report$model_data$pepmods)) {
      res$pepmod_id <- rep_len(report$model_data$pepmods$pepmod_id[1], nrow(res))
    }
    res
  },
  NULL, global_vars = c())
  return(res)
}

rbind_all_frames <- function(frames_coll_list, frame_names = NULL, link_col = NULL)
{
  if ( is.null(frame_names) ) {
    # select names that appear in all frames
    all_frame_names <- unlist( lapply( frames_coll_list, function( frames_coll ) {
      res <- names(frames_coll)
      res[ sapply( res, function(fname) is.data.frame(frames_coll[[fname]]) ) ]
    } ) )
    frame_names_freq <- table(all_frame_names)
    frame_names <- names(frame_names_freq)[ frame_names_freq == max(frame_names_freq) ]
  }
  res <- lapply( frame_names, function( fname ) {
    bind_rows( lapply( seq_along(frames_coll_list), function( coll_ix ) {
      res <- frames_coll_list[[coll_ix]][[fname]]
      if ( !is.null(link_col) && nrow(res)>0 ) res[[link_col]] <- if (!is.null(names(frames_coll_list))) names(frames_coll_list)[[coll_ix]] else coll_ix
      return ( res )
    } ) )
  } )
  names(res) <- frame_names
  return ( res )
}

expand_collapsed <- function(df, collapsed_col, separated_col, extra_cols=NULL, sep=";") {
  exp_list <- strsplit(df[[collapsed_col]], sep, fixed=TRUE)
  exp_lengths <- sapply(exp_list, length)
  res <- df[rep.int(1:nrow(df), exp_lengths), c(collapsed_col, extra_cols)]
  res[[separated_col]] <- unlist(exp_list)
  res[,c(separated_col, collapsed_col, extra_cols)]
}
