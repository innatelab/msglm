# Utilities for pcp fitting
#
# Author: Alexey Stukalov
###############################################################################

join_report_frames <- function(reports, frame_extractor=function(report) stop('no frame extractor'),
                               global_vars = c('model_dataset', 'version', 'chunk'),
                               prot_info=NULL)
{
  joined_res.df <- bind_rows(lapply( reports, function(report) {
    frame <- frame_extractor(report)
    if (is.null(frame)) return ( NULL )
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
