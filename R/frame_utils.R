# Utilities for assembling Stan report frames
#
# Author: Alexey Stukalov
###############################################################################

#' Binds the rows of the data frames from the list
#' @param frames_coll_list list of lists of data frames
#'
#' @param frame_names character vector of data frame names (subelements of `frames_coll_list` elements) to bind
#' @param collection_idcol name of the additional column that will be set to the name of the `frames_coll_list`
#'        element where the rows are coming from
#' @param verbose produce verbose output
#'
#' @export
rbind_all_frames <- function(frames_coll_list, frame_names = NULL, collection_idcol = NULL, verbose=FALSE)
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
    dplyr::bind_rows(lapply(seq_along(frames_coll_list), function(coll_ix) {
      res <- frames_coll_list[[coll_ix]][[fname]]
      if (!is.null(collection_idcol) && nrow(res)>0) {
        res[[collection_idcol]] <- if (!is.null(names(frames_coll_list))) {
          names(frames_coll_list)[[coll_ix]]
        } else coll_ix
      }
      return (res)
    }))
  })
  names(res) <- frame_names
  return (res)
}

maybe_rename <- function(df, cols, verbose=FALSE) {
  for (i in seq_along(cols)) {
    old_col <- cols[[i]]
    if (rlang::has_name(df, old_col)) {
      new_col <- names(cols)[[i]]
      if (is.na(new_col) || (new_col == old_col)) next
      if (rlang::has_name(df, new_col)) {
        message(new_col, " already exists, ", old_col, " not renamed")
      } else {
        message("renaming ", old_col, " to ", new_col)
        df <- dplyr::rename_at(df, old_col, ~new_col)
      }
    }
  }
  return(df)
}

#' Creates a matrix filled with a given value.
#'
#' @param val value for the matrix elements
#' @param ... matrix dimnames specification
#'
#' @export
constant_matrix <- function(val, dimnames, .var.name = varname(val))
{
  checkmate::assert_scalar(val, .var.name = .var.name)
  checkmate::assert_list(dimnames, len=2, names="unique")
  matrix(val, ncol = length(dimnames[[2]]), nrow = length(dimnames[[1]]),
         dimnames = dimnames)
}

#' Converts data.frame into a matrix
#' using row_col and col_col as its rows and columns and val_col as its values
#'
#' @param df data.frame (in a long format) to convert
#' @param row_col the Id of the row in the resulting matrix
#' @param col_col the Id of the column in the resulting matrix
#' @param val_col the values of the matrix elements
#' @param val_default the default value if no corresponding row in `df` exist
#' @param cols (optional) column ids of the resulting matrix (in a specified order)
#' @param rows (optional) row ids of the resulting matrix (in a specified order)
#'
#' @export
frame2matrix <- function(df, row_col, col_col, val_col="w", val_default=0, cols=NULL, rows=NULL) {
  # fix the rows/cols order
  if (is.null(rows)) {
    df_rows <- df[[row_col]]
    rows <- if (is.factor(df_rows)) {
      levels(df_rows)
    } else {
      sort(unique(df_rows))
    }
  }
  if (is.null(cols)) {
    df_cols <- df[[col_col]]
    cols <- if (is.factor(df_cols)) {
      levels(df_cols)
    } else {
      sort(unique(df_cols))
    }
  }
  df_expanded <- tidyr::expand_grid(!!sym(col_col) := cols, !!sym(row_col) := rows) %>%
    dplyr::left_join(df, by=c(col_col, row_col)) %>%
    dplyr::mutate(!!sym(val_col) := coalesce(!!sym(val_col), val_default),
                  !!sym(row_col) := factor(!!sym(row_col), levels=rows),
                  !!sym(col_col) := factor(!!sym(col_col), levels=cols)) %>%
    dplyr::arrange_at(c(col_col, row_col))
  mtx <- stats::xtabs(as.formula(paste0(val_col, " ~ ", row_col, " + ", col_col)), data=df_expanded)
  checkmate::assert_set_equal(names(dimnames(mtx)), c(row_col, col_col), ordered=TRUE)
  if (nrow(mtx) == 0 || is.integer(rows) && vctrs::vec_equal(rows, seq_len(nrow(mtx)))) {
    rownames(mtx) <- NULL
  } else {
    checkmate::assert_set_equal(rownames(mtx), rows, ordered=TRUE)
  }
  if (ncol(mtx) == 0 || is.integer(cols) && vctrs::vec_equal(cols, seq_len(ncol(mtx)))) {
    colnames(mtx) <- NULL
  } else {
    checkmate::assert_set_equal(colnames(mtx), cols, ordered=TRUE)
  }
  return(mtx)
}

#' Converts a matrix into a 3-column frame (matrix row, matrix column, cell value).
#'
#' @param mtx matrix to convert
#' @param row_col the column of the result containing matrix row names
#' @param col_col the column of the result containing matrix column names
#' @param val_col the column of the result containing matrix element values
#' @param skip_val if not null, the element of the matrix that have this value would not be exported into data frame
#'
#' @export
matrix2frame <- function(mtx, row_col = NULL, col_col = NULL, val_col = "w", skip_val = 0) {
  dnn <- names(dimnames(mtx))
  if (is.null(row_col)) {
    row_col <- dnn[[1]]
    checkmate::assert_string(row_col, na.ok=FALSE, null.ok=FALSE)
  } else {
    dnn[[1]] <- row_col # overriding the existing name
  }
  if (is.null(col_col)) {
    col_col <- dnn[[2]]
    checkmate::assert_string(col_col, na.ok=FALSE, null.ok=FALSE)
  } else {
    dnn[[2]] <- col_col # overriding the existing name
  }
  df <- as.data.frame.table(mtx, dnn=dnn, stringsAsFactors=TRUE, responseName=val_col)
  if (!is.null(skip_val)) {
    df <- dplyr::filter(df, !!sym(val_col) != skip_val)
  }

  # fix that for empty matrices, rows/cols columns are not added
  if (nrow(mtx) == 0) {
    # add missing row column
    df[[row_col]] <- if (is.na(names(dimnames(mtx))[[1]])) integer(0) else character(0)
  } else if (any(rownames(mtx) == as.character(seq_len(nrow(mtx))))) {
    # don't keep indices as factors
    df[[row_col]] <- as.integer(df[[row_col]])
  }
  if (ncol(mtx) == 0) {
    # add missing col column
    df[[col_col]] <- if (is.na(names(dimnames(mtx))[[2]])) integer(0) else character(0)
  } else if (any(colnames(mtx) == as.character(seq_len(ncol(mtx))))) {
    # don't keep indices as factors
    df[[col_col]] <- as.integer(df[[col_col]])
  }
  return(df)
}
