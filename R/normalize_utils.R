# Normalization of MS experiments
#
# Author: astukalov
###############################################################################

normalize_experiments <- function(norm_model, def_norm_data, df,
                                  quant_col = "intensity", obj_col = "protgroup_id",
                                  exp_col = "mschannel_ix", cond_col = "bait_orf",
                                  cond_group_col = NULL, exp_shifts = NULL, shift_col="shift",
                                  missing_exp.ratio=0.1, missing_cond.ratio=if (exp_col==cond_col) missing_exp.ratio else 0.0,
                                  method = c("optimizing", "mcmc", "vb"), center_shifts=TRUE,
                                  max_quant.ratio=NA, exp_shift.min = -0.5, Rhat_max = 1.1,
                                  max_objs=1000L,
                                  mcmc.iter=2000L, mcmc.chains=4, mcmc.thin=4, mcmc.adapt_delta=0.9,
                                  vb.iter=100000L,
                                  verbose=FALSE)
{
    method <- match.arg(method)
    if (is.null(cond_group_col)) {
      df$`__all__` <- TRUE
      cond_group_col <- "__all__"
    }
    # prepare "standartized" data frame
    df_std <- df[,c(obj_col, exp_col, cond_col, cond_group_col, quant_col)]
    names(df_std) <-c("obj", "experiment", "condition", "cond_group", "quant")
    if (is.null(exp_shifts)) {
      # no shifts by default
      exp_shifts <- data.frame(experiment = sort(unique(df_std$experiment)),
                               shift = rep.int(0.0, n_distinct(df_std$experiment)),
                               stringsAsFactors = TRUE)
    } else {
      exp_shifts <- exp_shifts[,c(exp_col, shift_col)]
      colnames(exp_shifts) <- c("experiment", "shift")
    }
    if (!is.na(exp_shift.min) && is.finite(exp_shift.min)) {
      exp_shifts <- dplyr::mutate(exp_shifts, is_used = shift >= exp_shift.min)
      if (!all(exp_shifts$is_used)) {
        warning(sum(!exp_shifts$is_used), " experiment(s) with low shifts excluded from ", cond_col, " normalization: ",
                paste0(exp_shifts$experiment[!exp_shifts$is_used], collapse=" "))
        exp_shifts <- dplyr::filter(exp_shifts, is_used)
      }
      df_std <- dplyr::inner_join(df_std, exp_shifts %>% dplyr::select(experiment, shift)) %>%
        dplyr::mutate(norm_quant = quant*exp(-shift)) %>% dplyr::select(-shift)
    } else {
      exp_shifts <- dplyr::mutate(exp_shifts, is_used = TRUE)
      df_std$norm_quant <- df_std$quant
    }
    exp_shifts$experiment <- as.character(exp_shifts$experiment)
    obj_group_stats <- df_std %>% dplyr::filter(!is.na(quant)) %>%
      dplyr::group_by(obj, cond_group) %>%
      dplyr::summarise(n_experiments = n_distinct(experiment),
                       n_conditions = n_distinct(condition)) %>%
      dplyr::ungroup()
    group_stats <- dplyr::group_by(obj_group_stats, cond_group) %>%
      dplyr::summarise(n_max_experiments = max(n_experiments),
                       n_max_conditions = max(n_conditions))
    valid_objs <- dplyr::inner_join(obj_group_stats, group_stats) %>%
      dplyr::filter(n_experiments > 1L & n_experiments >= (1-missing_exp.ratio)*n_max_experiments &
                    n_conditions >= (1-missing_cond.ratio)*n_max_conditions)
    message(n_distinct(group_stats$cond_group), " condition group(s)")
    res <- dplyr::group_by(valid_objs, cond_group) %>% do({
        df_group <- dplyr::inner_join(df_std, .)
        group_objs <- df_group %>%
          dplyr::group_by(obj, n_max_experiments, n_max_conditions, cond_group) %>%
          dplyr::summarise(quant_med = median(norm_quant, na.rm=TRUE),
                           quant_min = min(norm_quant, na.rm=TRUE),
                           quant_max = max(norm_quant, na.rm=TRUE),
                           sd = sd(norm_quant, na.rm=TRUE)) %>% dplyr::ungroup() %>%
          dplyr::mutate(sd_rel = sd/quant_med,
                        sd_quantile = percent_rank(if_else(!is.na(sd_rel), sd_rel, Inf)),
                        is_valid = is.na(max_quant.ratio) |
                                   (quant_min*max_quant.ratio > quant_med &
                                    quant_max < max_quant.ratio*quant_med))
        # max sure at least 25% of group objects or max_objs (if specified) are valid for normalization
        if (sum(group_objs$is_valid) < pmin(0.25*nrow(group_objs),
                                            ifelse(max_objs > 0L, max_objs, Inf))) {
            group_objs$is_valid <- group_objs$sd_quantile <= 0.25
        }
        valid_group_objs <- dplyr::filter(group_objs, is_valid)
        sel_group_objs <- if (max_objs > 0L && nrow(valid_group_objs) > max_objs) {
          sample_n(valid_group_objs, max_objs)
        } else {
          valid_group_objs
        }
        df_group <- dplyr::semi_join(df_group, sel_group_objs)
        norm_data <- def_norm_data
        norm_data$Nexperiments <- sel_group_objs$n_max_experiments[1]
        norm_data$Nconditions <- sel_group_objs$n_max_conditions[1]
        norm_data$Nobjects <- nrow(sel_group_objs)
        message("Normalizing group '", sel_group_objs$cond_group[1], "', ",
                norm_data$Nobjects, " object(s), ",
                norm_data$Nconditions, " condition(s), ",
                norm_data$Nexperiments, " experiment(s)")
        if (nrow(df_group) == 0L) {
            # degenerated case
            warning("No valid observations in group ", .$cond_group[1])
            res <- data.frame(condition = as.character(df_group$condition[1]),
                              shift = 0.0,
                              stringsAsFactors = FALSE)
            colnames(res) <- c(cond_col, "shift")
            return (res)
        } else if (norm_data$Nconditions == 1L) {
            # another degenerated case
            res <- data.frame(condition = as.character(df_group$condition[1]),
                              shift = 0.0,
                              stringsAsFactors = FALSE)
            colnames(res) <- c(cond_col, "shift")
            return (res)
        }
        #if (any(is.na(df_group$quant))) {
        #  stop("Non-quanted observations detected")
        #}
        exp_df <- dplyr::select(df_group, experiment, condition) %>%
          dplyr::distinct() %>% dplyr::arrange(condition, experiment) %>%
          dplyr::mutate(experiment = as.character(experiment),
                        condition= as.character(condition)) %>%
          dplyr::left_join(exp_shifts) %>%
          dplyr::mutate(experiment = factor(experiment, levels=experiment),
                        condition = factor(condition, levels=unique(condition))) %>%
          dplyr::arrange(as.integer(experiment))
        if (any(is.na(exp_df$shift))) {
          warning("No shifts for experiments: ",
                  paste0(exp_df$experiment[is.na(exp_df$shift)], collapse=" "))
        }
        norm_data$experiment2condition <- as.array(as.integer(exp_df$condition))
        norm_data$experiment_shift <- as.array(exp_df$shift)
        df_group <- dplyr::mutate(df_group,
                                  experiment = factor(as.character(experiment),
                                                      levels=levels(exp_df$experiment)))
        df_group <- dplyr::left_join(tidyr::expand(df_group, obj, experiment), df_group) %>%
            dplyr::arrange(as.integer(obj), as.integer(experiment)) %>%
            dplyr::mutate(safe_quant = if_else(is.na(quant), 0.0, quant))

        norm_data$qData <- matrix(df_group$safe_quant, ncol = nrow(exp_df), byrow=TRUE) #nrow=nrow(valid_group_objs)
        message("Running Stan optimization...")
        out_params <- c("data_sigma", "condition_sigma", "condition_shift")
        if (method == 'optimizing') {
          norm_fit <- optimizing(norm_model, norm_data, algorithm="LBFGS",
                                 init=list(condition_sigma=1.0, condition_shift0=as.array(rep.int(0.0, norm_data$Nconditions-1L))),
                                 history_size=10L)
          cond_shift_pars <- norm_fit$par[str_detect(names(norm_fit$par), "^condition_shift\\[\\d+\\]$")]
          cond_shift_ixs <- as.integer(str_match(names(cond_shift_pars), "\\[(\\d+)\\]$")[,2])
          res <- data.frame(condition = levels(exp_df$condition)[cond_shift_ixs],
                            shift = as.numeric(cond_shift_pars),
                            stringsAsFactors=FALSE)
        } else if (method == 'mcmc') {
          norm_fit <- sampling(norm_model, norm_data, chains=mcmc.chains, iter=mcmc.iter, thin=mcmc.thin,
                               pars=out_params, include=TRUE,
                               init=function() list(condition_sigma=1.0, condition_shift0=as.array(rep.int(0.0, norm_data$Nconditions-1L))),
                               control = list(adapt_delta=mcmc.adapt_delta))
          norm_fit_stat <- monitor(norm_fit)
          cond_shift_mask <- str_detect(rownames(norm_fit_stat), "^condition_shift\\[\\d+\\]$")
          nonconv_mask <- norm_fit_stat[cond_shift_mask, 'Rhat'] > Rhat_max
          if (any(nonconv_mask)) {
            warning("Rhat>", Rhat_max, " for ", sum(nonconv_mask), " shift(s)")
          }
          n_eff_min <- 0.1*mcmc.iter
          nonconv_mask <- norm_fit_stat[cond_shift_mask, 'n_eff'] < n_eff_min
          if (any(nonconv_mask)) {
            warning("n_eff<", n_eff_min, " for ", sum(nonconv_mask), " shift(s)")
          }
          cond_shift_pars <- norm_fit_stat[cond_shift_mask, 'mean']
          cond_shift_ixs <- as.integer(str_match(rownames(norm_fit_stat)[cond_shift_mask], "\\[(\\d+)\\]$")[,2])
          res <- data.frame(condition = levels(exp_df$condition)[cond_shift_ixs],
                            shift = as.numeric(cond_shift_pars),
                            Rhat = norm_fit_stat[cond_shift_mask, 'Rhat'],
                            n_eff = norm_fit_stat[cond_shift_mask, 'n_eff'],
                            stringsAsFactors=FALSE)
        } else if (method == 'vb') {
          norm_fit <- vb(norm_model, norm_data, iter=vb.iter,
                         pars=out_params, include=TRUE,
                         init=function() list(condition_sigma=1.0, condition_shift0=as.array(rep.int(0.0, norm_data$Nconditions-1L))) )
          norm_fit_stat <- monitor(norm_fit)
          cond_shift_mask <- str_detect(rownames(norm_fit_stat), "^condition_shift\\[\\d+\\]$")
          cond_shift_pars <- norm_fit_stat[cond_shift_mask, 'mean']
          cond_shift_ixs <- as.integer(str_match(rownames(norm_fit_stat)[cond_shift_mask], "\\[(\\d+)\\]$")[,2])
          res <- data.frame(condition = levels(exp_df$condition)[cond_shift_ixs],
                            shift = as.numeric(cond_shift_pars),
                            stringsAsFactors=FALSE)
        } else {
          stop('Unknown method ', method)
        }
        if (center_shifts) {
          res$shift <- res$shift - median(res$shift)
        }
        col_renames <- "condition"
        names(col_renames) <- cond_col
        res <- dplyr::rename_(res, .dots=col_renames) %>%
          dplyr::mutate(
            method = method,
            n_objects = norm_data$Nobjects,
            n_experiments = norm_data$Nexperiments)
    }) %>% dplyr::ungroup()
    if (cond_group_col == "__all__") {
      res$cond_group <- NULL
    } else {
      new_colnames <- colnames(res)
      new_colnames[new_colnames == "cond_group"] <- cond_group_col
      colnames(res) <- new_colnames
    }
    res
}
