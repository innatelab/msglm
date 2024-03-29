# Normalization of MS experiments
#
# Author: astukalov
###############################################################################

#' Calculate normalization shifts within one MS condition group.
#' @export
norm_shifts.condgroup <- function(stan_norm_model, quantobj_mscalib,
        msdata_df, mschan_preshifts, condgroup_id, cond_col,
        max_objs=1000L, quant_ratio.max=NA,
        stan_method = c("mcmc", "variational", "optimizing"),
        mcmc.iter=2000L, mcmc.chains=getOption("mc.cores", 4), mcmc.thin=1L, mcmc.adapt_delta=0.9,
        vb.iter=100000L,
        Rhat.max = 1.1, neff_ratio.min = 0.25,
        shifts_constraint=c("none", "mean=0", "median=0"),
        verbose=FALSE)
{
    objs <- msdata_df %>%
      dplyr::group_by(obj, n_max_mschannels, n_max_conditions) %>%
      dplyr::summarise(quant_med = median(norm_quant, na.rm=TRUE),
                       quant_min = min(norm_quant, na.rm=TRUE),
                       quant_max = max(norm_quant, na.rm=TRUE),
                       sd = sd(norm_quant, na.rm=TRUE)) %>% dplyr::ungroup() %>%
      dplyr::mutate(sd_rel = sd/quant_med,
                    sd_quantile = percent_rank(replace_na(sd_rel, Inf)),
                    is_valid = is.na(quant_ratio.max) |
                               (quant_min*quant_ratio.max > quant_med &
                                quant_max < quant_ratio.max*quant_med))
    # max sure at least 25% of group objects or max_objs (if specified) are valid for normalization
    if (sum(objs$is_valid) < pmin(0.25*nrow(objs),
                                  ifelse(max_objs > 0L, max_objs, Inf))) {
        objs$is_valid <- objs$sd_quantile <= 0.25
    }
    valid_objs <- dplyr::filter(objs, is_valid)
    sel_objs <- if (max_objs > 0L && nrow(valid_objs) > max_objs) {
        dplyr::sample_n(valid_objs, max_objs)
    } else {
        valid_objs
    }
    msdata_df_orig <- msdata_df
    msdata_df <- dplyr::semi_join(msdata_df, dplyr::select(sel_objs, obj), by="obj")
    stan_input <- list(
      Nobjects = nrow(sel_objs)
    )

    stan_input <- modifyList(stan_input, to_standata(quantobj_mscalib, convert_base=TRUE, silent=!verbose))

    if (nrow(msdata_df) == 0L) {
        # degenerated case, no data
        warning("No valid observations in group ", condgroup_id)
        res <- tibble(condition = unique(msdata_df$condition),
                      shift = 0.0) %>% rlang::set_names(nm=c(cond_col, "shift"))
        return (res)
    } else if (n_distinct(msdata_df$condition) == 1L) {
        # another degenerated case, single condition
        res <- tibble(condition = msdata_df$condition[1],
                      shift = 0.0) %>% rlang::set_names(nm=c(cond_col, "shift"))
        return (res)
    }
    #if (any(is.na(msdata_df$quant))) {
    #  stop("Non-quanted observations detected")
    #}
    mschan_df <- dplyr::select(msdata_df, mschannel, condition, sumgroup) %>%
      dplyr::distinct() %>% dplyr::arrange(sumgroup, condition, mschannel) %>%
      dplyr::mutate(mschannel = as.character(mschannel),
                    condition = as.character(condition),
                    sumgroup = as.character(sumgroup)) %>%
      dplyr::inner_join(dplyr::select(mschan_preshifts, mschannel, preshift), by="mschannel") %>%
      dplyr::mutate(mschannel = factor(mschannel, levels=mschannel),
                    condition = factor(condition, levels=unique(condition)),
                    sumgroup = factor(sumgroup, levels=unique(sumgroup))) %>%
      dplyr::arrange(as.integer(mschannel))
    if (any(is.na(mschan_df$preshift))) {
        warning("No pre-shifts for mschannels: ",
                paste0(mschan_df$mschannel[is.na(mschan_df$preshift)], collapse=" "))
    }
    stan_input$Nmschannels <- nrow(mschan_df)
    stan_input$Nshifts <- n_distinct(mschan_df$condition)
    stan_input$Nsumgroups <- n_distinct(mschan_df$sumgroup)
    message("Normalizing group '", condgroup_id, "', ",
            stan_input$Nobjects, " object(s), ",
            stan_input$Nshifts, " shift(s) to identify, ",
            stan_input$Nmschannels, " mschannel(s), ",
            stan_input$Nsumgroups, " averaging mschannel group(s)")
    stan_input$mschannel2shift <- as.array(as.integer(mschan_df$condition))
    stan_input$mschannel2sumgroup <- as.array(as.integer(mschan_df$sumgroup))
    stan_input$mschannel_preshift <- as.array(mschan_df$preshift)
    msdata_df <- dplyr::mutate(msdata_df,
                               obj = factor(as.character(obj)),
                               mschannel = factor(as.character(mschannel),
                                                  levels=levels(mschan_df$mschannel)))
    msdata_df <- dplyr::left_join(tidyr::expand(msdata_df, obj, mschannel), msdata_df,
                                  by = c("obj", "mschannel")) %>%
      dplyr::arrange(as.integer(obj), as.integer(mschannel)) %>%
      dplyr::mutate(safe_quant = replace_na(quant, 1.0))

    stan_input$qData <- matrix(msdata_df$safe_quant,
                               ncol = nrow(mschan_df),
                               byrow=TRUE) # nrow=nrow(valid_objs)
    message("Running Stan optimization...")
    out_params <- c("data_sigma", "shift_sigma", "shift")
    if (stan_method == 'optimizing') {
        norm_fit <- stan_norm_model$optimize(stan_input, algorithm="LBFGS",
                               init=list(sigma=1.0, shift0=as.array(rep.int(0.0, stan_input$Nshifts-1L))),
                               history_size=10L)
        shift_pars <- norm_fit$par[str_detect(names(norm_fit$par), "^shift\\[\\d+\\]$")]
        shift_ixs <- as.integer(str_match(names(shift_pars), "\\[(\\d+)\\]$")[,2])
        res <- tibble(condition = levels(mschan_df$condition)[shift_ixs],
                      shift = as.numeric(shift_pars))
    } else if (stan_method == 'mcmc') {
        norm_fit <- stan_norm_model$sample(stan_input, chains=mcmc.chains,
                                           iter_warmup=0.5*mcmc.iter, iter_sampling=0.5*mcmc.iter,
                                           thin=mcmc.thin,
                                    #pars=out_params, include=TRUE,
                                    #init=function() list(shift_sigma=1.0, shift_shift0=as.array(rep.int(0.0, stan_input$Nshifts-1L))),
                                    adapt_delta=mcmc.adapt_delta, show_messages=verbose)
        neff_min <- neff_ratio.min * mcmc.iter
    } else if (stan_method == 'variational') {
      for (i in 1:10) {
        norm_fit <- stan_norm_model$variational(stan_input, iter=vb.iter, tol_rel_obj=0.01,
                                                algorithm="meanfield", eta=0.1,
                                                init=function() list(data_sigma_a=1.0, data_sigma_t=1.0,
                                                                     shift_sigma_a=1.0, shift_sigma_t=1.0,
                                                                     shift0_unscaled=as.array(rnorm(stan_input$Nshifts-1L, sd=0.1))))
        if (!rlang::inherits_any(try(norm_fit$metadata()), "try-error")) break
      }
      neff_min <- neff_ratio.min * vb.iter
    } else {
      stop('Unknown method ', stan_method)
    }
    res <- norm_fit$summary(variables = out_params, posterior_summary_metrics) %>%
      tidyr::extract("variable", "shift_ix", "^shift\\[(\\d+)\\]$", remove=FALSE, convert=TRUE) %>%
      dplyr::filter(!is.na(shift_ix)) %>%
      dplyr::mutate(condition = levels(mschan_df$condition)[shift_ix],
                    converged = (rhat <= Rhat.max) & (ess_bulk >= neff_min)) %>%
      dplyr::select(condition, shift=mean, shift_sd=sd, shift_median=median, shift_mad=mad,
                    shift_q25=q25, shift_q75=q75, shift_q2.5=q2.5, shift_q97.5=q97.5,
                    rhat, ess_bulk, ess_tail, converged)
    if (!all(res$converged)) {
        warning("Convergence problems for ", sum(!res$converged), " shift(s)")
    }
    if (shifts_constraint == "median=0") {
        res <- dplyr::mutate(res, shift = shift - median(shift))
    } else if (shifts_constraint == "mean=0") {
        res <- dplyr::mutate(res, shift = shift - mean(shift))
    }
    # fix back the condition col into logical
    if (is.logical(msdata_df_orig$condition)) {
      res <- dplyr::mutate(res, condition = as.logical(condition))
    }

    res <- dplyr::rename(res, !!sym(cond_col) := condition) %>%
      dplyr::mutate(
        stan_method = stan_method,
        n_objects = stan_input$Nobjects,
        n_mschannels = stan_input$Nmschannels)
    return (res)
}

#' @export
normalize_experiments <- function(quantobj_mscalib, msdata_df,
                                  quant_col = "intensity", obj_col = "protgroup_id",
                                  mschan_col = "mschannel", cond_col="condition", condgroup_col = NULL, sumgroup_col = NULL,
                                  mschan_preshifts = NULL, preshift_col="shift",
                                  stan_method = c("mcmc", "variational", "optimizing"),
                                  max_objs=1000L,
                                  quant_ratio.max=NA, mschan_shift.min = -0.5,
                                  nmschan_ratio.min=0.9, ncond_ratio.min=if (mschan_col==cond_col) nmschan_ratio.min else 1.0,
                                  mcmc.iter=2000L, mcmc.chains=getOption("mc.cores", 4L), mcmc.thin=1L, mcmc.adapt_delta=0.9,
                                  vb.iter=100000L,
                                  Rhat.max = 1.1, neff_ratio.min = 0.25,
                                  shifts_constraint=c("none", "mean=0", "median=0"),
                                  verbose=FALSE)
{
    stan_method <- match.arg(stan_method)
    shifts_constraint <- match.arg(shifts_constraint)
    stan_norm_model <- msglm_stan_model("msglm_normalize")
    if (is.null(condgroup_col)) {
      if (verbose) message("No condgroup_col specified, putting all conditions in a single group")
      msdata_df$`__fake_condgroup__` <- TRUE
      condgroup_col <- "__fake_condgroup__"
    }
    if (is.null(sumgroup_col)) {
      if (verbose) message("No sumgroup_col specified, putting all mschannels in a single averaging group")
      msdata_df$`__fake_sumgroup__` <- TRUE
      sumgroup_col <- "__fake_sumgroup__"
    }
    # prepare "standartized" data frame (note that some original columns may be duplicated)
    msdata_df_std <- dplyr::select_at(msdata_df, c(obj=obj_col, mschannel=mschan_col,
                                                   condition=cond_col, condgroup=condgroup_col,
                                                   sumgroup=sumgroup_col, quant=quant_col))
    if (is.null(mschan_preshifts)) {
      # no shifts by default
      mschan_preshifts <- tibble(mschannel = sort(unique(msdata_df_std$mschannel)),
                                 preshift = rep_len(0.0, n_distinct(msdata_df_std$mschannel)))
    } else {
      mschan_preshifts <- dplyr::select_at(mschan_preshifts, c(mschannel = mschan_col, preshift = preshift_col))
    }
    if (!is.na(mschan_shift.min) && is.finite(mschan_shift.min)) {
      mschan_preshifts <- dplyr::inner_join(mschan_preshifts,
                                            dplyr::distinct(dplyr::select(msdata_df_std, mschannel, condgroup, sumgroup)),
                                            by = "mschannel") %>%
          dplyr::group_by(condgroup) %>%
          dplyr::mutate(is_used = preshift >= mean(preshift) + mschan_shift.min) %>%
          dplyr::ungroup() %>%
          dplyr::select(-condgroup, -sumgroup)
      if (!all(mschan_preshifts$is_used)) {
        warning(sum(!mschan_preshifts$is_used), " mschannel(s) with low pre-shifts excluded from ", cond_col, " normalization: ",
                paste0(mschan_preshifts$mschannel[!mschan_preshifts$is_used], collapse=" "))
        mschan_preshifts <- dplyr::filter(mschan_preshifts, is_used)
      }
      msdata_df_std <- dplyr::inner_join(msdata_df_std,
                                         dplyr::select(mschan_preshifts, mschannel, preshift),
                                         by = "mschannel") %>%
        dplyr::mutate(norm_quant = quant*exp(-preshift)) %>% dplyr::select(-preshift)
    } else {
      mschan_preshifts <- dplyr::mutate(mschan_preshifts, is_used = TRUE)
      msdata_df_std$norm_quant <- msdata_df_std$quant
    }
    mschan_preshifts$mschannel <- as.character(mschan_preshifts$mschannel)
    obj_condgroup_stats <- dplyr::filter(msdata_df_std, !is.na(quant)) %>%
      dplyr::group_by(obj, condgroup) %>%
      dplyr::summarise(n_mschannels = n_distinct(mschannel),
                       n_conditions = n_distinct(condition)) %>%
      dplyr::ungroup()
    condgroup_stats <- dplyr::group_by(obj_condgroup_stats, condgroup) %>%
      dplyr::summarise(n_max_mschannels = max(n_mschannels),
                       n_max_conditions = max(n_conditions))
    message(n_distinct(condgroup_stats$condgroup), " group(s)")
    valid_objs <- dplyr::inner_join(obj_condgroup_stats, condgroup_stats, by = "condgroup") %>%
      dplyr::filter(n_mschannels > 1L & n_mschannels >= nmschan_ratio.min*n_max_mschannels &
                    n_conditions >= ncond_ratio.min*n_max_conditions)
    res <- dplyr::group_by(valid_objs, condgroup) %>%
           dplyr::group_modify(.keep=TRUE, ~ norm_shifts.condgroup(stan_norm_model, quantobj_mscalib,
                              dplyr::inner_join(msdata_df_std, .x, by = c('obj', 'condgroup')), .y$condgroup,
                              cond_col=cond_col,
                              mschan_preshifts=mschan_preshifts,
                              stan_method=stan_method,
                              max_objs=max_objs,
                              quant_ratio.max=quant_ratio.max,
                              mcmc.iter=mcmc.iter, mcmc.chains=mcmc.chains,
                              mcmc.thin=mcmc.thin, mcmc.adapt_delta=mcmc.adapt_delta,
                              vb.iter=vb.iter,
                              Rhat.max=Rhat.max, neff_ratio.min=neff_ratio.min,
                              shifts_constraint=shifts_constraint,
                              verbose=verbose)
    ) %>% dplyr::ungroup()
    # back to user-specified name of the group column
    if (condgroup_col == "__fake_condgroup__") {
      res$condgroup <- NULL
    } else {
      res <- dplyr::rename(res, !!sym(condgroup_col) := condgroup)
    }
    # remove sumgroup column
    res$sumgroup <- NULL
    return(res)
}

#' @export
multilevel_normalize_experiments <- function(quantobj_mscalib,
                                             mschannels_df, msdata_df,
                                  quant_col = "intensity", obj_col = "protgroup_id",
                                  mschan_col = "mschannel",
                                  mschan_preshifts_df = NULL, mschan_preshift_col = "shift",
                                  norm_levels,
                                  max_objs=1000L,
                                  nmschan_ratio.min=0.9, ncond_ratio.min=1.0,
                                  quant_ratio.max=NA, mschan_shift.min = -0.5, Rhat.max = 1.1, neff_ratio.min = 0.25,
                                  stan_method = c("mcmc", "optimizing", "vb"),
                                  mcmc.iter=2000L, mcmc.chains=getOption("mc.cores", 4), mcmc.thin=4, mcmc.adapt_delta=0.9,
                                  vb.iter=100000L,
                                  shifts_constraint=c("none", "median=0", "mean=0"),
                                  verbose=FALSE)
{
  stan_method <- match.arg(stan_method)
  shifts_constraint <- match.arg(shifts_constraint)

  # FIXME compose msrun info
  lev_cols <- unique(c(unlist(lapply(norm_levels, function(lev) c(lev$cond_col, lev$condgroup_col, lev$sumgroup_col))), mschan_col))
  lev_cols <- lev_cols[!is.na(lev_cols)]
  mschan_df <- dplyr::select_at(mschannels_df, lev_cols) %>% dplyr::distinct()
  msdata_df <- dplyr::select_at(msdata_df, c(obj_col, mschan_col, quant_col)) %>%
    dplyr::left_join(mschan_df, by=mschan_col)
  total_shift_col <- paste0("total_",mschan_col,"_shift")
  # initialize mschannel shifts
  if (is.null(mschan_preshifts_df)) {
    mschan_shifts_df <- mschan_df
    mschan_shifts_df[[total_shift_col]] <- 0.0
  } else {
    # use existing shifts
    mschan_shifts_df <- dplyr::select_at(mschan_preshifts_df, c(mschan_col, mschan_preshift_col))
    mschan_shifts_df <- dplyr::left_join(mschan_df, mschan_shifts_df, by = mschan_col)
    mschan_shifts_df[[total_shift_col]] <- replace_na(mschan_shifts_df[[mschan_preshift_col]])
  }
  lev_norm_res <- list()
  for (i in seq_along(norm_levels)) {
    lev_info <- norm_levels[[i]]
    next_lev_info <- if (i < length(norm_levels)) norm_levels[[i+1]] else list()
    lev_name <- names(norm_levels)[[i]]
    if (verbose) message("Normalizing ", lev_name, " (level #", i, ")...")
    lev_shifts_df <- normalize_experiments(quantobj_mscalib, msdata_df,
                                           obj_col = obj_col, quant_col = quant_col,
                                           mschan_col = mschan_col,
                                           cond_col = lev_info$cond_col,
                                           condgroup_col = if (rlang::has_name(lev_info, "condgroup_col")) {lev_info$condgroup_col} else {next_lev_info$cond_col %||% NULL}, # NULL in lev_info is possible
                                           sumgroup_col = lev_info$sumgroup_col %||% NULL,
                                           mschan_preshifts = mschan_shifts_df,
                                           preshift_col = if (is.null(mschan_shifts_df)) NA_character_ else total_shift_col,
                                           shifts_constraint = lev_info$shifts_contraint %||% shifts_constraint,
                                           stan_method = lev_info$stan_method %||% stan_method,
                                           max_objs = lev_info$max_objs %||% max_objs,
                                           Rhat.max = Rhat.max, neff_ratio.min = neff_ratio.min,
                                           nmschan_ratio.min = lev_info$nmschan_ratio.min %||% nmschan_ratio.min,
                                           ncond_ratio.min = lev_info$ncond_ratio.min %||% ncond_ratio.min,
                                           quant_ratio.max = lev_info$quant_ratio.max %||% quant_ratio.max,
                                           mschan_shift.min = lev_info$mschan_shift.min %||% mschan_shift.min,
                                           mcmc.iter=mcmc.iter, mcmc.chains=mcmc.chains, mcmc.thin=mcmc.thin, mcmc.adapt_delta=mcmc.adapt_delta,
                                           vb.iter=vb.iter, verbose=verbose)
    lev_shift_col <- paste0(lev_name, "_shift")
    mschan_shifts_df <- dplyr::left_join(mschan_shifts_df,
                                         dplyr::select_at(lev_shifts_df, c(lev_info$cond_col, "shift")),
                                         by = lev_info$cond_col)
    colnames(mschan_shifts_df)[colnames(mschan_shifts_df)=="shift"] <- lev_shift_col
    mschan_shifts_df[[total_shift_col]] <- mschan_shifts_df[[total_shift_col]] +
          tidyr::replace_na(mschan_shifts_df[[lev_shift_col]], 0.0)
    lev_norm_res[[lev_name]] <- list(level_shifts = lev_shifts_df,
                                     mschannel_shifts = mschan_shifts_df)
  }
  return(list(levels=lev_norm_res,
              mschannel_shifts=mschan_shifts_df))
}
