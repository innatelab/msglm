# Normalization of MS experiments
#
# Author: astukalov
###############################################################################

normalize_experiments <- function(stan_norm_model, stan_input_base, msdata_df,
                                  quant_col = "intensity", obj_col = "protgroup_id",
                                  mschan_col = "mschannel", cond_col = "condition",
                                  cond_group_col = NULL, mschan_shifts = NULL, shift_col="shift",
                                  stan_method = c("optimizing", "mcmc", "vb"),
                                  max_objs=1000L,
                                  quant_ratio.max=NA, mschan_shift.min = -0.5,
                                  nmschan_ratio.min=0.9, ncond_ratio.min=if (mschan_col==cond_col) nmschan_ratio.min else 1.0,
                                  mcmc.iter=2000L, mcmc.chains=4, mcmc.thin=4, mcmc.adapt_delta=0.9,
                                  vb.iter=100000L,
                                  Rhat.max = 1.1, neff_ratio.min = 0.25,
                                  shifts_constraint=c("none", "mean=0", "median=0"),
                                  verbose=FALSE)
{
    stan_method <- match.arg(stan_method)
    shifts_constraint <- match.arg(shifts_constraint)
    if (is.null(cond_group_col)) {
      msdata_df$`__all__` <- TRUE
      cond_group_col <- "__all__"
    }
    # prepare "standartized" data frame
    msdata_df_std <- msdata_df[,c(obj_col, mschan_col, cond_col, cond_group_col, quant_col)]
    names(msdata_df_std) <-c("obj", "mschannel", "condition", "cond_group", "quant")
    if (is.null(mschan_shifts)) {
      # no shifts by default
      mschan_shifts <- data.frame(mschannel = sort(unique(msdata_df_std$mschannel)),
                                  shift = rep.int(0.0, n_distinct(msdata_df_std$mschannel)),
                                  stringsAsFactors = TRUE)
    } else {
      mschan_shifts <- mschan_shifts[,c(mschan_col, shift_col)]
      colnames(mschan_shifts) <- c("mschannel", "shift")
    }
    if (!is.na(mschan_shift.min) && is.finite(mschan_shift.min)) {
      mschan_shifts <- dplyr::mutate(mschan_shifts, is_used = shift >= mschan_shift.min)
      if (!all(mschan_shifts$is_used)) {
        warning(sum(!mschan_shifts$is_used), " mschannel(s) with low shifts excluded from ", cond_col, " normalization: ",
                paste0(mschan_shifts$mschannel[!mschan_shifts$is_used], collapse=" "))
        mschan_shifts <- dplyr::filter(mschan_shifts, is_used)
      }
      msdata_df_std <- dplyr::inner_join(msdata_df_std, dplyr::select(mschan_shifts, mschannel, shift)) %>%
        dplyr::mutate(norm_quant = quant*exp(-shift)) %>% dplyr::select(-shift)
    } else {
      mschan_shifts <- dplyr::mutate(mschan_shifts, is_used = TRUE)
      msdata_df_std$norm_quant <- msdata_df_std$quant
    }
    mschan_shifts$mschannel <- as.character(mschan_shifts$mschannel)
    obj_condgrp_stats <- msdata_df_std %>% dplyr::filter(!is.na(quant)) %>%
      dplyr::group_by(obj, cond_group) %>%
      dplyr::summarise(n_mschannels = n_distinct(mschannel),
                       n_conditions = n_distinct(condition)) %>%
      dplyr::ungroup()
    condgrp_stats <- dplyr::group_by(obj_condgrp_stats, cond_group) %>%
      dplyr::summarise(n_max_mschannels = max(n_mschannels),
                       n_max_conditions = max(n_conditions))
    message(n_distinct(condgrp_stats$cond_group), " condition group(s)")
    valid_objs <- dplyr::inner_join(obj_condgrp_stats, condgrp_stats) %>%
      dplyr::filter(n_mschannels > 1L & n_mschannels >= nmschan_ratio.min*n_max_mschannels &
                    n_conditions >= ncond_ratio.min*n_max_conditions)
    res <- dplyr::group_by(valid_objs, cond_group) %>% do({
        condgrp_msdata_df <- dplyr::inner_join(msdata_df_std, .)
        condgrp_objs <- condgrp_msdata_df %>%
          dplyr::group_by(obj, n_max_mschannels, n_max_conditions, cond_group) %>%
          dplyr::summarise(quant_med = median(norm_quant, na.rm=TRUE),
                           quant_min = min(norm_quant, na.rm=TRUE),
                           quant_max = max(norm_quant, na.rm=TRUE),
                           sd = sd(norm_quant, na.rm=TRUE)) %>% dplyr::ungroup() %>%
          dplyr::mutate(sd_rel = sd/quant_med,
                        sd_quantile = percent_rank(if_else(!is.na(sd_rel), sd_rel, Inf)),
                        is_valid = is.na(quant_ratio.max) |
                                   (quant_min*quant_ratio.max > quant_med &
                                    quant_max < quant_ratio.max*quant_med))
        # max sure at least 25% of group objects or max_objs (if specified) are valid for normalization
        if (sum(condgrp_objs$is_valid) < pmin(0.25*nrow(condgrp_objs),
                                            ifelse(max_objs > 0L, max_objs, Inf))) {
            condgrp_objs$is_valid <- condgrp_objs$sd_quantile <= 0.25
        }
        valid_condgrp_objs <- dplyr::filter(condgrp_objs, is_valid)
        sel_condgrp_objs <- if (max_objs > 0L && nrow(valid_condgrp_objs) > max_objs) {
          sample_n(valid_condgrp_objs, max_objs)
        } else {
          valid_condgrp_objs
        }
        condgrp_msdata_df <- dplyr::semi_join(condgrp_msdata_df, sel_condgrp_objs)
        stan_input <- stan_input_base
        stan_input$Nobjects <- nrow(sel_condgrp_objs)
        message("Normalizing group '", sel_condgrp_objs$cond_group[1], "', ",
                stan_input$Nobjects, " object(s), ",
                stan_input$Nconditions, " condition(s), ",
                stan_input$Nmschannels, " mschannel(s)")
        if (nrow(condgrp_msdata_df) == 0L) {
            # degenerated case
            warning("No valid observations in group ", .$cond_group[1])
            res <- data.frame(condition = as.character(condgrp_msdata_df$condition[1]),
                              shift = 0.0,
                              stringsAsFactors = FALSE)
            colnames(res) <- c(cond_col, "shift")
        } else if (n_distinct(condgrp_msdata_df$condition) == 1L) {
            # another degenerated case
            res <- data.frame(condition = as.character(condgrp_msdata_df$condition[1]),
                              shift = 0.0,
                              stringsAsFactors = FALSE)
            colnames(res) <- c(cond_col, "shift")
        } else {
        #if (any(is.na(condgrp_msdata_df$quant))) {
        #  stop("Non-quanted observations detected")
        #}
        mschan_df <- dplyr::select(condgrp_msdata_df, mschannel, condition) %>%
          dplyr::distinct() %>% dplyr::arrange(condition, mschannel) %>%
          dplyr::mutate(mschannel = as.character(mschannel),
                        condition= as.character(condition)) %>%
          dplyr::inner_join(mschan_shifts) %>%
          dplyr::mutate(mschannel = factor(mschannel, levels=mschannel),
                        condition = factor(condition, levels=unique(condition))) %>%
          dplyr::arrange(as.integer(mschannel))
        if (any(is.na(mschan_df$shift))) {
          warning("No shifts for mschannels: ",
                  paste0(mschan_df$mschannel[is.na(mschan_df$shift)], collapse=" "))
        }
        stan_input$Nmschannels <- nrow(mschan_df)
        stan_input$Nconditions <- n_distinct(mschan_df$condition)
        stan_input$mschannel2condition <- as.array(as.integer(mschan_df$condition))
        stan_input$mschannel_shift <- as.array(mschan_df$shift)
        condgrp_msdata_df <- dplyr::mutate(condgrp_msdata_df,
                                  mschannel = factor(as.character(mschannel),
                                                     levels=levels(mschan_df$mschannel)))
        condgrp_msdata_df <- dplyr::left_join(tidyr::expand(condgrp_msdata_df, obj, mschannel), condgrp_msdata_df) %>%
            dplyr::arrange(as.integer(obj), as.integer(mschannel)) %>%
            dplyr::mutate(safe_quant = if_else(is.na(quant), 0.0, quant))

        stan_input$qData <- matrix(condgrp_msdata_df$safe_quant, ncol = nrow(mschan_df), byrow=TRUE) #nrow=nrow(valid_condgrp_objs)
        message("Running Stan optimization...")
        out_params <- c("data_sigma", "condition_sigma", "condition_shift")
        if (stan_method == 'optimizing') {
          norm_fit <- optimizing(stan_norm_model, stan_input, algorithm="LBFGS",
                                 init=list(condition_sigma=1.0, condition_shift0=as.array(rep.int(0.0, stan_input$Nconditions-1L))),
                                 history_size=10L)
          cond_shift_pars <- norm_fit$par[str_detect(names(norm_fit$par), "^condition_shift\\[\\d+\\]$")]
          cond_shift_ixs <- as.integer(str_match(names(cond_shift_pars), "\\[(\\d+)\\]$")[,2])
          res <- data.frame(condition = levels(mschan_df$condition)[cond_shift_ixs],
                            shift = as.numeric(cond_shift_pars),
                            stringsAsFactors=FALSE)
        } else if (stan_method == 'mcmc') {
          norm_fit <- sampling(stan_norm_model, stan_input, chains=mcmc.chains, iter=mcmc.iter, thin=mcmc.thin,
                               pars=out_params, include=TRUE,
                               init=function() list(condition_sigma=1.0, condition_shift0=as.array(rep.int(0.0, stan_input$Nconditions-1L))),
                               control = list(adapt_delta=mcmc.adapt_delta))
          norm_fit_stat <- monitor(norm_fit)
          cond_shift_mask <- str_detect(rownames(norm_fit_stat), "^condition_shift\\[\\d+\\]$")
          nonconv_mask <- norm_fit_stat[cond_shift_mask, 'Rhat'] > Rhat.max
          if (any(nonconv_mask)) {
            warning("Rhat>", Rhat.max, " for ", sum(nonconv_mask), " shift(s)")
          }
          neff_min <- neff_ratio.min*mcmc.iter
          nonconv_mask <- norm_fit_stat[cond_shift_mask, 'n_eff'] < neff_min
          if (any(nonconv_mask)) {
            warning("n_eff<", neff_min, " for ", sum(nonconv_mask), " shift(s)")
          }
          cond_shift_pars <- norm_fit_stat[cond_shift_mask, 'mean']
          cond_shift_ixs <- as.integer(str_match(rownames(norm_fit_stat)[cond_shift_mask], "\\[(\\d+)\\]$")[,2])
          res <- data.frame(condition = levels(mschan_df$condition)[cond_shift_ixs],
                            shift = as.numeric(cond_shift_pars),
                            Rhat = norm_fit_stat[cond_shift_mask, 'Rhat'],
                            n_eff = norm_fit_stat[cond_shift_mask, 'n_eff'],
                            stringsAsFactors=FALSE)
        } else if (stan_method == 'vb') {
          norm_fit <- vb(stan_norm_model, stan_input, iter=vb.iter,
                         pars=out_params, include=TRUE,
                         init=function() list(condition_sigma=1.0, condition_shift0=as.array(rep.int(0.0, stan_input$Nconditions-1L))) )
          norm_fit_stat <- monitor(norm_fit)
          cond_shift_mask <- str_detect(rownames(norm_fit_stat), "^condition_shift\\[\\d+\\]$")
          cond_shift_pars <- norm_fit_stat[cond_shift_mask, 'mean']
          cond_shift_ixs <- as.integer(str_match(rownames(norm_fit_stat)[cond_shift_mask], "\\[(\\d+)\\]$")[,2])
          res <- data.frame(condition = levels(mschan_df$condition)[cond_shift_ixs],
                            shift = as.numeric(cond_shift_pars),
                            stringsAsFactors=FALSE)
        } else {
          stop('Unknown method ', stan_method)
        }
        if (shifts_constraint == "median=0") {
          res$shift <- res$shift - median(res$shift)
        } else if (shifts_constraint == "mean=0") {
          res$shift <- res$shift - mean(res$shift)
        }
        col_renames <- "condition"
        names(col_renames) <- cond_col
        res <- dplyr::rename_(res, .dots=col_renames) %>%
          dplyr::mutate(
            stan_method = stan_method,
            n_objects = stan_input$Nobjects,
            n_mschannels = stan_input$Nmschannels)
        }
        res
    }) %>% dplyr::ungroup()
    # back to user-specified name of the condition group column
    if (cond_group_col == "__all__") {
      res$cond_group <- NULL
    } else {
      new_colnames <- colnames(res)
      new_colnames[new_colnames == "cond_group"] <- cond_group_col
      colnames(res) <- new_colnames
    }
    res
}

multilevel_normalize_experiments <- function(stan_norm_model, instr_calib,
                                             mschannels_df, msdata_df,
                                  quant_col = "intensity", obj_col = "protgroup_id",
                                  mschan_col = "mschannel",
                                  norm_levels,
                                  max_objs=1000L,
                                  nmschan_ratio.min=0.9, ncond_ratio.min=1.0,
                                  quant_ratio.max=NA, mschan_shift.min = -0.5, Rhat.max = 1.1, neff_ratio.min = 0.25,
                                  stan_method = c("mcmc", "optimizing", "vb"),
                                  mcmc.iter=2000L, mcmc.chains=4, mcmc.thin=4, mcmc.adapt_delta=0.9,
                                  vb.iter=100000L,
                                  shifts_constraint=c("none", "median=0", "mean=0"),
                                  verbose=FALSE)
{
  stan_method <- match.arg(stan_method)
  shifts_constraint <- match.arg(shifts_constraint)
  stan_input_base <- instr_calib[c('zDetectionFactor', 'zDetectionIntercept',
                               'detectionMax', 'sigmaScaleHi', 'sigmaScaleLo',
                               'sigmaOffset', 'sigmaBend', 'sigmaSmooth',
                               'zShift', 'zScale')]
  # FIXME compose msrun info
  lev_cols <- unique(c(unlist(lapply(norm_levels, function(lev) c(lev$cond_col, lev$cond_group_col))), mschan_col))
  lev_cols <- lev_cols[!is.na(lev_cols)]
  mschan_df <- dplyr::select_(mschannels_df, .dots=lev_cols) %>% dplyr::distinct()
  msdata_df <- dplyr::select_(msdata_df, .dots=c(obj_col, mschan_col, quant_col)) %>%
    dplyr::left_join(mschan_df)

  mschan_shifts_df <- NULL
  lev_norm_res <- list()
  for (i in seq_along(norm_levels)) {
    lev_info <- norm_levels[[i]]
    next_lev_info <- if (i < length(norm_levels)) norm_levels[[i+1]] else list()
    lev_name <- names(norm_levels)[[i]]
    if (verbose) {
      message("Normalizing ", lev_name, " (level #", i, ")...")
    }
    lev_shifts_df <- normalize_experiments(stan_norm_model, stan_input_base, msdata_df,
                                           obj_col = obj_col, quant_col = quant_col,
                                           mschan_col = mschan_col, cond_col = lev_info$cond_col,
                                           cond_group_col = lev_info$cond_group_col %||% (next_lev_info$cond_col %||% NULL),
                                           mschan_shifts = mschan_shifts_df,
                                           shift_col = if (is.null(mschan_shifts_df)) NA_character_ else "total_mschannel_shift",
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
    if (is.null(mschan_shifts_df)) {
      # initialize mschannel shifts
      mschan_shifts_df <- dplyr::mutate(dplyr::inner_join(lev_shifts_df, mschan_df), total_mschannel_shift = 0.0)
    } else {
      mschan_shifts_df <- dplyr::inner_join(mschan_shifts_df, dplyr::select_(lev_shifts_df, .dots=c(lev_info$cond_col, "shift")))
    }
    lev_shift_col <- paste0(lev_name, "_shift")
    colnames(mschan_shifts_df)[colnames(mschan_shifts_df)=="shift"] <- lev_shift_col
    mschan_shifts_df$total_mschannel_shift <- mschan_shifts_df$total_mschannel_shift + mschan_shifts_df[[lev_shift_col]]
    lev_norm_res[[lev_name]] <- list(level_shifts = lev_shifts_df,
                                     mschannel_shifts = mschan_shifts_df)
  }
  return(list(levels=lev_norm_res,
              mschannel_shifts=mschan_shifts_df))
}
