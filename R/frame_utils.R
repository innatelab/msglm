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

# Expands a frame by the column that contains sep-separated list of entities.
# In the new frame there is a row for each entity.
expand_collapsed <- function(df, collapsed_col, separated_col, extra_cols=NULL, sep=";") {
  exp_list <- strsplit(df[[collapsed_col]], sep, fixed=TRUE)
  exp_lengths <- sapply(exp_list, length)
  res <- df[rep.int(1:nrow(df), exp_lengths), c(collapsed_col, extra_cols)]
  res[[separated_col]] <- unlist(exp_list)
  res[,c(separated_col, collapsed_col, extra_cols)]
}

# converts data.frame df (long format) into a matrix
# using row_col and col_col as its rows and columns and val_col as its values
frame2matrix <- function(df, row_col, col_col, val_col, cols=NULL, rows=NULL) {
  df_wide <- reshape(dplyr::select_(df, .dots=c(row_col, col_col, val_col)),
        direction = "wide", timevar=col_col, idvar=row_col,
        v.names = val_col) %>%
  dplyr::mutate_at(vars(starts_with(val_col)), funs(if_else(is.na(.), 0.0, .)))
  nz_mtx <- as.matrix(dplyr::select(df_wide, -one_of(row_col)))
  rownames(nz_mtx) <- df_wide[[row_col]]
  colnames(nz_mtx) <- str_replace(colnames(nz_mtx), fixed(paste0(val_col, ".")), "")
  mtx_dims <- list(rows %||% rownames(nz_mtx),
                   cols %||% colnames(nz_mtx))
  names(mtx_dims) <- c(row_col, col_col)
  mtx <- do.call(zero_matrix, mtx_dims)
  mtx[rownames(nz_mtx), colnames(nz_mtx)] <- nz_mtx
  return(mtx)
}
