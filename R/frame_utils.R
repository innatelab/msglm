# Utilities for pcp fitting
#
# Author: Alexey Stukalov
###############################################################################

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
    frame$report_name <- report_name
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
