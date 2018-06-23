# Misc utilities for STAN output processing
#
# Author: astukalov
###############################################################################

append_sparse <- function(stan_data, mtx_name, mtx) {
    mtx_sparse <- extract_sparse_parts(mtx)
    stan_data[[paste0(mtx_name, "_Nw")]] <- length(mtx_sparse$w)
    stan_data[[paste0(mtx_name, "_w")]] <- mtx_sparse$w
    stan_data[[paste0(mtx_name, "_u")]] <- mtx_sparse$u
    stan_data[[paste0(mtx_name, "_v")]] <- mtx_sparse$v
    return(stan_data)
}

stan.iterations_frame <- function(stan_result)
{
  n_iterations <- dim( stan_result )[[1]]
  n_chains <- dim( stan_result )[[2]]
  n_thin <- stan_result@sim$thin

  expand.grid( sample_ix = seq_len(n_iterations),
               chain = seq_len(n_chains) ) %>%
    mutate( iteration = sample_ix * n_thin,
            unpermuted_ix = row_number() )
}

stan.extract_samples <- function(stan_result, pars, min.iteration = NA, permuted = FALSE)
{
  res <- rstan::extract(stan_result, pars = pars, permuted = permuted, inc_warmup = !is.na(min.iteration))
  if ( permuted ) {
    attr(res, 'iter_info') <- stan.iterations_frame(stan_result)
  }
  if ( !is.na(min.iteration) ) {
    n_thin <- stan_result@sim$thin
    if ( permuted ) {
      iter_info <- attr(res, 'iter_info') %>% dplyr::filter(iteration >= min.iteration - stan_result@sim$warmup)
      sample_ixs <- sort(iter_info$unpermuted_ix)
      res <- lapply( res, function(var_samples) {
        new_dims <- dim(var_samples)
        new_dims[[1]] <- n_distinct(sample_ixs)
        array( var_samples[slice.index(var_samples,1) %in% sample_ixs], dim = new_dims ) } )
      attr(res, 'iter_info') <- iter_info
      res
    } else {
      res[seq(from=as.integer(ceiling(min.iteration/n_thin)), dim(res)[[1]] ), , ]
    }
  } else {
    res
  }
}

extract_index <- function(var_names) {
  res <- apply(do.call(cbind, strsplit(gsub(']$', '', gsub('^[^[]+\\[', '', var_names)), ',')),
         1, as.integer )
  if (is.vector(res)) res <- matrix(res, nrow = length(var_names))
  return (res)
}
extract_var <- function(var_names) {
  gsub('\\[(\\d+\\,?)+]$', '', var_names)
}

pvalue_not_zero <- function(samples, tail = c("both", "negative", "positive"))
{
  tail = match.arg(tail)
  if (tail == "negative") {
    return(insilicoMop:::ProbabilityLessZeroSmoothed(samples, nsteps = 100, bandwidth = NA))
  } else if (tail == "positive") {
    return(insilicoMop:::ProbabilityLessZeroSmoothed(-samples, nsteps = 100, bandwidth = NA))
  } else if (tail == "both") {
    # 2x correction as both tails are tested
    2 * min(c(0.5,insilicoMop:::ProbabilityLessZeroSmoothed(samples, nsteps = 100, bandwidth = NA),
                  insilicoMop:::ProbabilityLessZeroSmoothed(-samples, nsteps = 100, bandwidth = NA)))
  }
}
