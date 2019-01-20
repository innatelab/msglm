# Utilities for assembling Stan report frames
#
# Author: Alexey Stukalov
###############################################################################

rbind_all_frames <- function(frames_coll_list, frame_names = NULL, link_col = NULL, verbose=FALSE)
{
  if (is.null(frame_names)) {
    if (verbose) message("Determining which frames to bind...")
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
    frame_names <- names(frame_names_freq)[all_frames_mask]
  }
  if (verbose) message("Binding rows of the following frame(s): ", paste0(frame_names, collapse=", "))
  res <- lapply(frame_names, function(fname) {
    if (verbose) message("Binding rows of ", fname)
    bind_rows(lapply(seq_along(frames_coll_list), function(coll_ix) {
      res <- frames_coll_list[[coll_ix]][[fname]]
      if (!is.null(link_col) && nrow(res)>0) res[[link_col]] <- if (!is.null(names(frames_coll_list))) names(frames_coll_list)[[coll_ix]] else coll_ix
      return (res)
    }))
  })
  names(res) <- frame_names
  return (res)
}

# Expands a frame by the column that contains sep-separated list of entities.
# In the new frame there is a row for each entity.
expand_collapsed <- function(df, collapsed_col, separated_col,
                             extra_cols=NULL, sep=fixed(";")) {
  exp_list <- str_split(df[[collapsed_col]], sep)
  exp_lengths <- sapply(exp_list, length)
  res <- df[rep.int(1:nrow(df), exp_lengths), c(collapsed_col, extra_cols)]
  res[[separated_col]] <- unlist(exp_list)
  res[,c(separated_col, collapsed_col, extra_cols)]
}

# converts data.frame df (long format) into a matrix
# using row_col and col_col as its rows and columns and val_col as its values
frame2matrix <- function(df, row_col, col_col, val_col="w", cols=NULL, rows=NULL) {
  mtx_dims <- list(if (!is.null(rows) && length(rows) > 0) rows else if (is.factor(df[[row_col]])) levels(df[[row_col]]) else as.character(unique(df[[row_col]])),
                   if (!is.null(cols) && length(cols) > 0) cols else if (is.factor(df[[col_col]])) levels(df[[col_col]]) else as.character(unique(df[[col_col]])))
  names(mtx_dims) <- c(row_col, col_col)
  mtx <- do.call(zero_matrix, mtx_dims)
  row_vals <- df[[row_col]]
  col_vals <- df[[col_col]]
  w_vals <- df[[val_col]]
  if (nrow(df) > 0L) {
    for (i in 1:nrow(df)) {
      mtx[row_vals[[i]], col_vals[[i]]] <- w_vals[[i]]
    }
  }
  return(mtx)
}
