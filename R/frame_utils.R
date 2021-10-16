# Utilities for assembling Stan report frames
#
# Author: Alexey Stukalov
###############################################################################

#' Binds the rows of the data frames from the list.
#'
#' @param frames_coll_list list of lists of data frames
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

# maybe rename the columns of the data.frame is they exist
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

#' Create a matrix filled with a given value.
#'
#' Create a constant matrix with all elements equal to *val*.
#'
#' @param val value for the matrix elements
#' @param ... matrix dimnames specification (as in [base::matrix()])
#' @param .var.name name of the matrix to use in diagnostic messages
#'
#' @export
constant_matrix <- function(val, dimnames, .var.name = checkmate::varname(val))
{
  checkmate::assert_scalar(val, .var.name = .var.name)
  checkmate::assert_list(dimnames, len=2, names="unique")
  matrix(val, ncol = length(dimnames[[2]]), nrow = length(dimnames[[1]]),
         dimnames = dimnames)
}

#' Convert "long" format data frame into a matrix.
#'
#' Convert *data.frame* into a matrix. The data frame is assumed to be
#' in a "long" format, i.e. each row of *df* should correspond to a single
#' element of the resulting matrix. The row and column of the element are
#' specified by the values in *row_col* and *col_col*, correspondingly,
#' and the value is taken from the *val_col* column.
#'
#' @param df *data.frame* (in a long format) to convert
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
    checkmate::assert_set_equal(rownames(mtx), as.character(rows), ordered=TRUE)
  }
  if (ncol(mtx) == 0 || is.integer(cols) && vctrs::vec_equal(cols, seq_len(ncol(mtx)))) {
    colnames(mtx) <- NULL
  } else {
    checkmate::assert_set_equal(colnames(mtx), as.character(cols), ordered=TRUE)
  }
  return(mtx)
}

#' Convert a matrix into a data frame.
#'
#' Convert a matrix into a "long" format data frame.
#' The resulting *data.frame* object will have three columns
#' (*row_col*, *col_col*, *val_col*) and for each matrix element there would
#' be exactly one row with values id of matrix row, id of matrix column
#' and the cell value, correspondingly.
#'
#' @param mtx matrix to convert
#' @param row_col the column of the result containing matrix row names
#' @param col_col the column of the result containing matrix column names
#' @param val_col the column of the result containing matrix element values
#' @param skip_val if not null, the element of the matrix that have this value would not be exported into data frame
#'
#' @export
matrix2frame <- function(mtx, row_col = NULL, col_col = NULL, val_col = "w", skip_val = 0) {
  df <- as.data.frame.table(mtx, stringsAsFactors=TRUE, responseName=val_col)
  dnn <- colnames(df)
  if (!is.null(row_col)) {
    # overriding the existing name
    df <- dplyr::rename(df, !!sym(row_col) := !!sym(dnn[[1]]))
  } else {
    row_col <- dnn[[1]]
  }
  # workaround df column being null (bug in as.data.frame.table()?)
  if (is.null(df[[1]])) {
    df[[1]] <- character(0)
  }
  if (is.null(df[[2]])) {
    df[[2]] <- character(0)
  }
  if (!is.null(col_col)) {
    # overriding the existing name
    df <- dplyr::rename(df, !!sym(col_col) := !!sym(dnn[[2]]))
  } else {
    col_col <- dnn[[2]]
  }
  if (!rlang::has_name(df, val_col)) {
    df[[val_col]] <- rlang::exec(typeof(skip_val), nrow(df))
  }
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

# check if the index column exists and equals to 1:nrow
# or create it otherwise
ensure_primary_index_column <- function(df, index_col, id_col=NA_character_, ids_ordered=NULL,
                                        create=FALSE, .var.name=checkmate::vname(df)) {
  checkmate::assert_data_frame(df, .var.name=.var.name)
  if (!is.na(id_col)) {
    colname <- paste0(.var.name, "$", id_col)
    if (rlang::has_name(df, id_col)) {
      if (is.factor(df[[id_col]])) df[[id_col]] <- as.character(df[[id_col]])
      if (!is.null(ids_ordered)) {
        checkmate::assert_character(ids_ordered, any.missing=FALSE, names="unnamed", unique=TRUE)
        checkmate::assert_set_equal(df[[id_col]], ids_ordered, ordered=FALSE, .var.name = colname)
      } else {
        checkmate::assert_character(as.character(df[[id_col]]), any.missing=FALSE, names="unnamed", unique=TRUE)
      }
    } else {
      stop("ID column '", id_col, "' not found in ", .var.name)
    }
  }
  if (rlang::has_name(df, index_col)) {
    colname <- paste0(.var.name, "$", index_col)
    checkmate::assert_integer(df[[index_col]], .var.name=colname)
    checkmate::assert_set_equal(df[[index_col]], seq_len(nrow(df)),
                                ordered=FALSE, .var.name=colname)
    if (!is.na(id_col)) {
      # check the order of the ids
      df_ids_ordered <- df[[id_col]]
      df_ids_ordered[df[[index_col]]] <- df_ids_ordered
      checkmate::assert_set_equal(df_ids_ordered, ids_ordered)
      df <- dplyr::arrange_at(df, index_col)
    }
  } else {
    if (create) {
      if (!is.na(id_col)) {
        if (!is.null(ids_ordered)) {
          df[[index_col]] <- match(df[[id_col]], ids_ordered)
          df <- dplyr::arrange_at(df, index_col)
        } else {
          df <- dplyr::arrange_at(df, id_col)
          df[[index_col]] <- seq_len(nrow(df))
        }
      } else {
        df[[index_col]] <- seq_len(nrow(df))
      }
    } else {
      stop('data.frame ', .var.name, ' doesn\'t have ', index_col, ' primary index')
    }
  }
  return(df)
}
