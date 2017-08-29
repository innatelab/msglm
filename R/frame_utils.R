# Utilities for assembling Stan report frames
#
# Author: Alexey Stukalov
###############################################################################

rbind_all_frames <- function(frames_coll_list, frame_names = NULL, link_col = NULL)
{
  if (is.null(frame_names)) {
    # select names that appear in all frames
    all_frame_names <- unlist(lapply(frames_coll_list, function(frames_coll) {
      res <- names(frames_coll)
      res[sapply(res, function(fname) is.data.frame(frames_coll[[fname]]))]
    }))
    frame_names_freq <- table(all_frame_names)
    all_frames_mask <- frame_names_freq == max(frame_names_freq)
    if (!all(all_frames_mask)) {
      warning("Frames ", paste0(names(frame_names_freq)[!all_frames_mask], collapse=", "),
              " are missing in some of the lists, ignored")
    }
    frame_names <- names(frame_names_freq)[all_frame_names]
  }
  res <- lapply(frame_names, function(fname) {
    bind_rows(lapply(seq_along(frames_coll_list), function(coll_ix) {
      res <- frames_coll_list[[coll_ix]][[fname]]
      if (!is.null(link_col) && nrow(res)>0) res[[link_col]] <- if (!is.null(names(frames_coll_list))) names(frames_coll_list)[[coll_ix]] else coll_ix
      return (res)
    }))
  })
  names(res) <- frame_names
  return (res)
}

expand_collapsed <- function(df, collapsed_col, separated_col, extra_cols=NULL, sep=";") {
  exp_list <- strsplit(df[[collapsed_col]], sep, fixed=TRUE)
  exp_lengths <- sapply(exp_list, length)
  res <- df[rep.int(1:nrow(df), exp_lengths), c(collapsed_col, extra_cols)]
  res[[separated_col]] <- unlist(exp_list)
  res[,c(separated_col, collapsed_col, extra_cols)]
}
