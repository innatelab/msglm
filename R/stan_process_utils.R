# Misc utilities for STAN output processing
#
# Author: astukalov
###############################################################################

# converts matrix to Stan-compatible compressed sparse row storage (CSR) and
# return a list of <mtx_name>_Nw, <mtx_name>_u, <mtx_name>_v, <mtx_name>_w
matrix2stancsr <- function(mtx, mtx_name) to_standata(as(Matrix::Matrix(mtx), "RsparseMatrix"), mtx_name)

#' @export
to_standata.RsparseMatrix <- function(spmtx, mtx_name) {
    mtx_csr <- list(
      Nw = length(spmtx@x),
      u = as.array(spmtx@p + 1L),
      v = as.array(spmtx@j + 1L),
      w = as.array(spmtx@x))
    names(mtx_csr) <- paste0(mtx_name, "_", names(mtx_csr))
    return(mtx_csr)
}

extract_var <- function(varspecs) {
  stringr::str_remove(varspecs, '\\[(?:\\d+\\,)*\\d+\\]$')
}

process_varspecs <- function(varspecs, vars_info, dims_info) {
  all_vars.df <- bind_rows(lapply(names(vars_info), function(cat) tibble(category = cat,
                                                                         var = vars_info[[cat]]$names)))
  all_varspecs.df <- tibble(varspec = varspecs,
                            index_varspec = seq_along(varspecs)) %>%
    tidyr::extract(varspec, c("var", "var_index"),
                   '^(\\w[^[]+)(?:\\[((?:\\d+\\,)*\\d+)\\])?$', remove=FALSE) %>%
    dplyr::mutate(var_index = ifelse(var_index == '', NA_character_, var_index)) %>%
    dplyr::inner_join(all_vars.df, by = 'var')
  cat_specs.df <- dplyr::distinct(all_varspecs.df, category, var_index) %>%
    dplyr::group_by(category)
  cats_info <- dplyr::group_map(cat_specs.df, function(cat_df, cat_row) {
      cat <- cat_row$category
      cat_info <- vars_info[[cat]]

      dims.df <- tibble(ix = seq_along(cat_info$dims),
                        name = if (is.null(cat_info$dims)) {character(0)} else {cat_info$dims}) %>%
        dplyr::mutate(col = paste0('index_', name)) %>%
        dplyr::group_by(name) %>%
        dplyr::mutate(local_ix = row_number(),
                      col = if (n() > 1) { paste0(col, '.', local_ix) } else { col }) %>%
        dplyr::ungroup()
      if (any(dims.df$local_ix > 1)) {
        warning('Category ', cat, ' dimensions ',
                dims.df$name[dims.df$local_ix > 1], ' not unique, using suffixes')
      }
      cat_df <- tidyr::separate(cat_df, var_index, c(dims.df$col), remove=FALSE, convert=TRUE, sep=",")
      # append info of each dimension
      for (dim_ix in seq_len(nrow(dims.df))) {
        dim_name = dims.df$name
        if (!(dim_name %in% names(dims_info))) {
          warning('No information for dimension #', dim_ix, " (", dim_name, ")")
        } else {
          dim_info = dplyr::select(dims_info[[dim_name]], -all_of(paste0('index_', dim_name)))
          ixs = cat_df[[dims.df$col[[dim_ix]]]]
          avail_ixs <- unique(ixs)
          def_ixs = 1:nrow(dim_info)
          missing_ixs = setdiff(def_ixs, avail_ixs)
          if (length(missing_ixs) > 0) {
            stop('Dimension #', dim_ix, '(', dim_name,'): indices ',
                 paste0(missing_ixs, collapse=' '), ' not found in the results')
          }
          extra_ixs = setdiff(avail_ixs, def_ixs)
          if (length(extra_ixs) > 0) {
            stop('Dimension #', dim_ix, '(', dim_name,'): indices ',
                paste0(extra_ixs, collapse=' '), ' have no information')
          }
          cat_df <- bind_cols(cat_df, dim_info[ixs, , drop=FALSE])
        }
      }
      return(cat_df)
    })
  names(cats_info) <- group_keys(cat_specs.df)$category
  return (list(spec_info = all_varspecs.df, cats_info = cats_info))
}

# compresses -log10(p_value) so that very significant p-values (mlog10(pvalue) >= threshold) are constrained
mlog10pvalue_compress_scalar <- function(x, threshold = 10) {
  if (x < threshold) {
    return (x)
  } else if (is.finite(x)) {
    t <- x - threshold
    return (threshold + sqrt(t))
  } else {
    return (3*threshold)
  }
}

mlog10pvalue_compress <- function(x, threshold = 10) {
  if (rlang::is_scalar_vector(x)) {
    mlog10pvalue_compress_scalar(x, threshold=threshold)
  } else {
    vapply(x, mlog10pvalue_compress_scalar, 0.0, threshold=threshold)
  }
}

pvalue_compare <- function(xsamples, y=0, tail = c("both", "negative", "positive"),
                           nsteps = 100, bandwidth = NA,
                           mlog10_threshold = 10,
                           mlog10_hard_threshold_factor = 3)
{
  tail = match.arg(tail)
  if (length(xsamples) == 0L) {
    warning("No samples provided, returning P-value=NA")
    return(NA_real_)
  } else if (tail == "negative") {
    res <- ProbabilityLessSmoothed(xsamples, y, nsteps = nsteps, bandwidth = bandwidth)
  } else if (tail == "positive") {
    res <- ProbabilityLessSmoothed(-xsamples, -y, nsteps = nsteps, bandwidth = bandwidth)
  } else if (tail == "both") {
    # 2x correction as both tails are tested
    res <- 2 * min(c(0.5,
                     ProbabilityLessSmoothed(xsamples, y, nsteps = nsteps, bandwidth = bandwidth,
                                             mlog10_threshold = mlog10_threshold,
                                             mlog10_hard_threshold_factor = mlog10_hard_threshold_factor),
                     ProbabilityLessSmoothed(-xsamples, -y, nsteps = nsteps, bandwidth = bandwidth,
                                             mlog10_threshold = mlog10_threshold,
                                             mlog10_hard_threshold_factor = mlog10_hard_threshold_factor)))
  }
  return(res)
}

# quantiles for symmetric 50% and 95% credible intervals
#' @export
quantiles_ci <- function(x) { posterior::quantile2(x, probs=c(0.025, 0.25, 0.75, 0.975)) }

# the default metrics the posterior::summary() function should calculate
posterior_summary_metrics <- c("mean", "median", "sd", "mad", "quantiles_ci",
                               "rhat", "ess_bulk", "ess_tail")
