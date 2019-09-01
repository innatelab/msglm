# Misc utilities for STAN output processing
#
# Author: astukalov
###############################################################################

# converts matrix to compressed row storage (CSR)
# wrapper for extract_sparse_parts that uses the specified matrix name
matrix2csr <- function(mtx_name, mtx) {
    mtx_sparse <- rstan::extract_sparse_parts(mtx)
    mtx_sparse$Nw <- length(mtx_sparse$w)
    names(mtx_sparse) <- paste0(mtx_name, "_", names(mtx_sparse))
    return(mtx_sparse)
}

stan.iterations_frame <- function(stan_result)
{
  n_iterations <- dim(stan_result)[[1]]
  n_chains <- dim(stan_result)[[2]]
  n_thin <- stan_result@sim$thin

  tidyr::crossing(chain = seq_len(n_chains),
                  sample_ix = seq_len(n_iterations)) %>%
    mutate(iteration = sample_ix * n_thin,
           unpermuted_ix = row_number())
}

stan.extract_samples <- function(stan_result, pars, chains = NA, min.iteration = NA, permuted = FALSE)
{
  res <- rstan::extract(stan_result, pars = pars, permuted = permuted, inc_warmup = !is.na(min.iteration))
  if (permuted) {
    iter_info.df <- stan.iterations_frame(stan_result)
  } else {
    iter_info.df <- NULL
  }
  if (!is.na(min.iteration) || !is.na(chains)) {
    n_thin <- stan_result@sim$thin
    if (permuted) {
      if (!is.na(min.iteration)) {
        iter_info.df <- dplyr::filter(iter_info.df, iteration >= min.iteration - stan_result@sim$warmup)
      }
      if (!is.na(chains)) {
        iter_info.df <- dplyr::filter(iter_info.df, chain %in% chains)
      }
      sample_ixs <- sort(iter_info.df$unpermuted_ix)
      n_samples_per_chain <- n_distinct(sample_ixs)
      res <- lapply(res, function(var_samples) {
        new_dims <- dim(var_samples)
        new_dims[[1]] <- n_samples_per_chain
        array(var_samples[slice.index(var_samples, 1) %in% sample_ixs], dim = new_dims) })
    } else {
      if (!is.na(min.iteration)) {
        res <- res[seq(from=as.integer(ceiling(min.iteration/n_thin)), dim(res)[[1]]), , ]
      }
      if (!is.na(chains)) {
        res <- res[, chains, , drop=FALSE]
      }
    }
  }
  if (!is.null(iter_info.df)) {
    attr(res, 'iter_info') <- iter_info.df
  }
  return(res)
}

extract_index <- function(var_names) {
  var_ndims <- str_count(var_names, stringr::fixed(",")) + 1L
  res <- str_remove(var_names, '\\]$') %>% str_remove('^[^\\[]+\\[') %>%
      str_split_fixed(stringr::fixed(','), var_ndims) %>%
      as.integer() %>% matrix(ncol=var_ndims)
  return(res)
}
extract_var <- function(var_names) {
  stringr::str_remove(var_names, '\\[(?:\\d+\\,)*\\d+\\]$')
}

pvalue_not_zero <- function(samples, tail = c("both", "negative", "positive"))
{
  tail = match.arg(tail)
  if (tail == "negative") {
    return(ProbabilityLessZeroSmoothed(samples, nsteps = 100, bandwidth = NA))
  } else if (tail == "positive") {
    return(ProbabilityLessZeroSmoothed(-samples, nsteps = 100, bandwidth = NA))
  } else if (tail == "both") {
    # 2x correction as both tails are tested
    2 * min(c(0.5,ProbabilityLessZeroSmoothed(samples, nsteps = 100, bandwidth = NA),
                  ProbabilityLessZeroSmoothed(-samples, nsteps = 100, bandwidth = NA)))
  }
}
