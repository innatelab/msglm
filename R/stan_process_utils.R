# Misc utilities for STAN output processing
# 
# Author: astukalov
###############################################################################

stan.iterations_frame <- function( stan_result )
{
  n_iterations <- dim( stan_result )[[1]]
  n_chains <- dim( stan_result )[[2]]
  n_thin <- stan_result@sim$thin

  expand.grid( sample_ix = seq_len(n_iterations),
               chain = seq_len(n_chains) ) %>%
    mutate( iteration = sample_ix * n_thin,
            unpermuted_ix = row_number() )
}

stan.extract_samples <- function( stan_result, pars, min.iteration = NA, permuted = FALSE )
{
  res <- rstan::extract( stan_result, pars = pars, permuted = permuted, inc_warmup = !is.na(min.iteration) )
  if ( permuted ) {
    attr(res, 'iter_info') <- stan.iterations_frame( stan_result )
  }
  if ( !is.na(min.iteration) ) {
    n_thin <- stan_result@sim$thin
    if ( permuted ) {
      iter_info <- attr(res, 'iter_info') %>% dplyr::filter( iteration >= min.iteration - stan_result@sim$warmup )
      sample_ixs <- sort(iter_info$unpermuted_ix)
      res <- lapply( res, function(var_samples) {
        new_dims <- dim(var_samples)
        new_dims[[1]] <- n_distinct(sample_ixs)
        array( var_samples[slice.index(var_samples,1) %in% sample_ixs], dim = new_dims ) } )
      attr(res, 'iter_info') <- iter_info
      res
    } else {
      res[ seq(from=as.integer(ceiling(min.iteration/n_thin)), dim(res)[[1]] ), , ]    
    }
  } else {
    res
  }
}

extract_index <- function( var_names ) {
  res <- apply( do.call( cbind, strsplit( gsub( ']$', '', gsub( '^[^[]+\\[', '', var_names ) ), ',' ) ),
         1, as.integer )
  if ( is.vector( res ) ) res <- matrix( res, nrow = length(var_names) )
  return ( res )
}
extract_var <- function( var_names ) {
  gsub( '\\[(\\d+\\,?)+]$', '', var_names )
}

pvalue_not_zero <- function( samples )
{
  2 * min( c( 0.5, insilicoMop:::ProbabilityLessZeroSmoothed( 
         samples, nsteps = 100, bandwidth = NA ),
       insilicoMop:::ProbabilityLessZeroSmoothed( 
         -samples, nsteps = 100, bandwidth = NA ) ) )
}

vars_pvalue_not_zero <- function( samples.df, vars_cat_info, dim_info )
{
   group_cols <- paste0( 'index_', vars_cat_info$dims )
   #print( str(samples.df) )
   #print( str(group_cols ) )
   p_value_all_samples <- function( samples ) {
       do.call( rbind, lapply( vars_cat_info$names, function(col) {
           data.frame( var = col,
                       p_value = pvalue_not_zero( samples[[col]] ),
                       stringsAsFactors = FALSE )
           } ) )
   }
   p_value.df <- samples.df %>%
        group_by_( .dots = group_cols ) %>%
        do( p_value_all_samples(.) )
   if ( 'fraction' %in% group_cols ) {
     p_value.df$fraction <- as.integer( p_value.df$fraction )
   }
   return ( p_value.df )
}

vars_contrast_stats <- function(samples.df, var_names, group_cols,
                                condition_col, contrastXcondition,
                                experiment_col = condition_col, condition2experiments.df = NULL,
                                nsteps = 100, maxBandwidth = NA,
                                quant.probs = c(0.025, 0.25, 0.50, 0.75, 0.975),
                                val_trans = NULL)
{
  n_samples <- n_distinct(samples.df$unpermuted_ix)
  if (n_samples == 0L) stop("No samples for ", paste0(var_names, collapse=" "), " ", condition_col)
  if (experiment_col != condition_col) {
    if (condition_col %in% colnames(samples.df)) {
      if (!is.null(condition2experiments.df)) {
        stop("Ambiguous experiment design specification: condition2experiment.df provided, but ", condition_col, " is also present in the samples.df")
      }
      condition2experiments.df <- dplyr::distinct(dplyr::select_(samples.df, .dots=c("condition" = condition_col, "experiment" = experiment_col)))
    } else {
      if (is.null(condition2experiments.df)) {
        stop("No experiment design: provide condition2experiment.df or ", condition_col, " in the samples.df")
      }
      condition2experiments.df <- dplyr::select_(condition2experiments.df, .dots=c("condition" = condition_col, "experiment" = experiment_col))
    }
  } else {
    condition2experiments.df <- dplyr::distinct(dplyr::select_(samples.df, .dots=c("condition" = condition_col))) %>%
      dplyr::mutate(experiment = condition)
  }

  contrast_stats_all_samples <- function( samples ) {
    #print(str(samples))
    # recode experiment indices to match what is in the samples
    cur_experiments = unique(samples[[experiment_col]])
    cur_cond2expr.df <- dplyr::filter(condition2experiments.df, experiment %in% cur_experiments) %>%
      dplyr::mutate(experiment = factor(as.character(experiment), levels=cur_experiments))
    # recode condition indices
    cur_contrastXcondition <- contrastXcondition[,colnames(contrastXcondition) %in% unique(cur_cond2expr.df$condition),drop=FALSE]
    cur_cond2expr.df <- dplyr::mutate(cur_cond2expr.df,
      condition = factor(as.character(condition), levels=colnames(cur_contrastXcondition)))

    bind_rows(lapply(var_names, function(var_col) {
      sample_vals <- samples[[var_col]]
      if (!is.null(val_trans)) {
        if (val_trans == "exp") {
          # FIXME that doesn't work for assembling per-protein reports
          sample_vals <- exp(sample_vals-median(sample_vals))
        }
        else {
          stop("Unknown val_trans=\"", val_trans, "\". Supported transforms: exp")
        }
      }
      samples.arr <- matrix(sample_vals, nrow=n_samples)

      #print(str(samples.arr))
      res <- insilicoMop:::ContrastStatistics( 
        samples.arr, as.integer(cur_cond2expr.df$experiment), as.integer(cur_cond2expr.df$condition),
        cur_contrastXcondition,
        nsteps = nsteps, maxBandwidth = maxBandwidth,
        quant_probs = quant.probs ) %>%
        as.data.frame(stringsAsFactors = FALSE, optional = TRUE) %>%
        mutate(var = var_col)
      # TODO se_mean, n_eff, Rhat
      return(res)
    } ) )
  }
  contrast_stats.df <- dplyr::arrange_(samples.df,
                                       .dots = c(group_cols, experiment_col, 'iteration', 'chain') ) %>%
    group_by_( .dots = group_cols ) %>%
    do( contrast_stats_all_samples(.) )
  if ( 'fraction' %in% group_cols ) {
    contrast_stats.df$fraction <- as.integer( contrast_stats.df$fraction )
  }
  if ( !is.null(dimnames(contrastXcondition)[[1]]) ) {
    contrast_stats.df$contrast <- dimnames(contrastXcondition)[[1]][contrast_stats.df$index_contrast]
  }
  return ( contrast_stats.df )
}

.attach_dim_info.by_array_index <- function( samples.df, stan_samples, dim_info )
{
    indexes.df <- do.call( cbind, lapply( seq_along(dim_info), function(dim_ix) {
      slice.index( stan_samples[[ colnames(samples.df)[[ncol(samples.df)]] ]], dim_ix+1 ) %>% as.vector()
    } ) ) %>% as.data.frame()
    colnames(indexes.df) <- paste0( 'index_', names(dim_info) )
    samples.df <- cbind( samples.df, indexes.df )
    for ( dim_ix in seq_along(dim_info) ) {
      dim_name <- names(dim_info)[[dim_ix]]
      if ( !is.null(dim_info[[dim_ix]]) ) {
        message( "Adding ", dim_name, " information" )
        samples.df <- cbind( samples.df, dim_info[[dim_ix]][ indexes.df[[dim_ix]], , drop=FALSE ] )
      } else {
        warning( "No ", dim_name, " information provided" )
      }
    }
    return ( samples.df )
}

.attach_dim_info.by_var_name <- function( data.df, dim_info )
{
    res.df <- data.df
    res_index <- extract_index( res.df$var_name )
    dims.df <- data.frame( ix = seq_along(dim_info),
                           name = names(dim_info),
                           stringsAsFactors = FALSE ) %>%
      group_by( name ) %>%
      mutate( local_ix = row_number(), suffix = if ( n() > 1 ) paste0( '.', local_ix ) else '' ) %>%
      ungroup()
    indexes.df <- do.call( cbind, lapply( seq_along(dim_info), function(dim_ix) {
      res_index[,dim_ix]
    } ) )
    colnames(indexes.df) <- paste0( 'index_', names(dim_info) )

    for ( dim_ix in seq_len(nrow(dims.df)) ) {
        if ( dim_ix > ncol(res_index) ) {
            stop( 'Data contain less dimensions than the dim_info' )
        }
        dim_name <- dims.df$name[[dim_ix]]
        if ( dims.df$local_ix[[dim_ix]] > 1 ) {
          warning( 'Dimension ', dim_name, ' is not unique, using suffixes' )
        }
        if ( !is.null(dim_info[[dim_ix]]) ) {
            res.df[[dim_name]] <- NULL
            if ( nrow(dim_info[[dim_ix]]) != max(res_index[,dim_ix]) ) {
                stop( 'Dimension ', dim_ix, ' info contains ', nrow(dim_info[[dim_ix]]),
                      ' elements, the data contain ', max(res_index[,dim_ix]) )
            }
            res_dim_info.df <- dim_info[[dim_ix]][ res_index[,dim_ix], , drop=FALSE ]
            if ( dims.df$suffix[[dim_ix]] != '' ) {
              colnames(res_dim_info.df) <- paste0( colnames(res_dim_info.df), dims.df$suffix[dim_ix] )
            }
            res.df <- cbind( res.df, res_dim_info.df )
        }
    }
    return ( cbind( indexes.df, res.df ) )
}

stan_samples_frame <- function(stan_samples, var_names, var_dims) {
    # extract only samples of selected variables
    # convert arrays of samples of selected variables into data frames,
    # add back single copy of dimensions and real iteration information
    samples.df <- do.call(data.frame, lapply(var_names, function(var_name) as.vector(stan_samples[[var_name]])))
    if (ncol(samples.df) > 0) {
        colnames(samples.df) <- var_names
        samples.df <- cbind(attr(stan_samples, 'iter_info'), samples.df)

        # add additional dimension information
        if (length(var_dims) > 0) {
            samples.df <- .attach_dim_info.by_array_index(samples.df, stan_samples, var_dims)
        }
    }
}

vars_statistics <- function( vars_category, stan_stats, stan_samples, vars_info, dim_info ) {
    message( 'Extracting ', vars_category, ' variables...' )
    vars_cat_info <- vars_info[[vars_category]]
    samples.df <- stan_samples_frame(stan_samples, vars_cat_info$names, dim_info[vars_cat_info$dims])
    # process convergence information
    stats.df <- subset( stan_stats, grepl( paste0( '^(', paste0( vars_cat_info$names, collapse='|' ), ')(\\[|$)' ), var_name ) )
    stats.df$var <- extract_var( stats.df$var_name )
    if ( nrow(stats.df) > 0 && length(vars_cat_info$dims) > 0 ) {
        stats.df <- .attach_dim_info.by_var_name( stats.df, dim_info[vars_cat_info$dims] )
    }
    return ( list( samples=samples.df, stats=stats.df ) )
}

vars_opt_convert <- function( vars_category, opt_results, vars_info, dim_info ) {
    message( 'Extracting ', vars_category, ' variables...' )
    vars_cat_info <- vars_info[[vars_category]]

    # extract only samples of selected variables
    # convert arrays of samples of selected variables into data frames,
    # add back single copy of dimensions and real iteration information
    res_mask <- grepl( paste0( '^(', paste0( vars_cat_info$names, collapse='|' ), ')(\\[|$)' ), names(opt_results$par) )
    if ( !any(res_mask) ) {
        # checkf if the variable is degenerated
        if ( prod( sapply( dim_info[vars_cat_info$dims], nrow ) ) > 0 ) {
          stop( 'Variables ', paste0( "'", vars_cat_info$names, "'", collapse = ", " ), ' not found in the opt. results' ) 
        } else {
          # skip the degenerated variable
          return ( data.frame() )
        }
    }
    res.df <- data.frame(
            var_name = names(opt_results$par)[res_mask],
            mean = opt_results$par[ res_mask ],
            stringsAsFactors = FALSE ) %>%
      mutate( var = extract_var( var_name ) )

    # add additional dimension information
    if ( length(vars_cat_info$dims) > 0 ) {
        miss_mask <- !(vars_cat_info$dims %in% names(dim_info))
        if ( any(miss_mask) ) {
            stop( sum(miss_mask), ' dimension info are missing: ',
                  paste0( vars_cat_info$dims[miss_mask], collapse = ' ' ) )
        }
        res.df <- .attach_dim_info.by_var_name( res.df, dim_info[vars_cat_info$dims] )
    }
    return ( res.df )
}
