require(tibble)
require(dplyr)
require(stringr)

gen_objects <- function(modelobject, n, quantobject = NA_character_, avg_nquantobj = 5,
                        labu_mean = 10, labu_sigma = 5, quantobj_labu_shift_shape = 3) {
    modelobj_id <- dplyr::sym(paste0(modelobject, "_id"))
    modelobj_label <- dplyr::sym(paste0(modelobject, "_label"))
    modelobjs_df = tibble(modelobj_id = seq_len(n),
                          !!modelobj_id := modelobj_id,
                          !!modelobj_label := paste0(modelobject, "_", !!modelobj_id),
                          labu_base = rnorm(n, mean=labu_mean, sd=labu_sigma))
    if (!is.na(quantobject)) {
        modelobjs_df <- dplyr::mutate(modelobjs_df, nquantobjs = 1 + rpois(length(labu_base),  avg_nquantobj - 1))
        quantobj_id <- dplyr::sym(paste0(quantobject, "_id"))
        quantobj_label <- dplyr::sym(paste0(quantobject, "_label"))
        quantobjs_df <- dplyr::group_by(modelobjs_df, !!modelobj_id, modelobj_id) %>%
            dplyr::group_modify(~{
                tibble(quantobj_lbl = paste0(.x[[modelobj_label]], "_", seq_len(.x$nquantobjs)),
                       labu_shift = c(0, -rgamma(.x$nquantobjs - 1L, quantobj_labu_shift_shape)))
            }) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(quantobj_id = row_number(),
                          !!quantobj_id := quantobj_id)
    } else {
        quantobjs_df <- NULL
    }
    return(list(modelobjects = modelobjs_df, quantobjects = quantobjs_df))
}

gen_effects <- function(effects_df, effect_range = c(-10, 10)) {
    dplyr::group_by(effects_df, effect) %>%
    dplyr::summarise(value = prior_mean + rnorm(1) * prior_tau * abs(rnorm(1) * rnorm(1))/
                             (rgamma(1, prior_df1/2, prior_df1/2)*rgamma(1, prior_df2/2, prior_df2/2)),
                     .groups = "drop") %>%
    dplyr::mutate(value = pmax(effect_range[[1]], pmin(value, effect_range[[2]])))
}

gen_modelobj_effects <- function(modelobjs_df, effects_df, effect_range = c(-10, 10)) {
    dplyr::group_by(modelobjs_df, modelobj_id) %>%
    dplyr::group_modify(~{
        gen_effects(effects_df = effects_df, effect_range = effect_range) %>%
        dplyr::mutate(labu_base= .x$labu_base)
    }) %>% dplyr::ungroup()
}

gen_interactions <- function(model_def, modelobjXeffects_df) {
    dplyr::group_by(modelobjXeffects_df, modelobj_id) %>%
    dplyr::group_modify(~{
        effect_values <- rlang::set_names(.x$value, .x$effect)
        iaction_values <- model_def$conditionXeffect %*% effect_values[colnames(model_def$conditionXeffect)]
        tibble(condition = rownames(model_def$conditionXeffect),
               value = as.vector(iaction_values) + .x$labu_base[[1]])
    }) %>% dplyr::ungroup()
}

gen_modelobj_intensities <- function(iactions_df, msprobes_df, labu_base_shift = 0.0, replicate_sigma = 0.1) {
    dplyr::inner_join(iactions_df, dplyr::select(msprobes_df, msprobe, condition), by="condition") %>%
    dplyr::mutate(replicate_shift = rnorm(n(), 0, replicate_sigma),
                  intensity = 2^(value + replicate_shift + labu_base_shift))
}

gen_quantobj_intensities <- function(iactions_df, quantobjs_df, mschannels_df,
                                     labu_base_shift = 0.0, replicate_sigma = 0.1) {
    dplyr::inner_join(iactions_df, dplyr::select(mschannels_df, msprobe, any_of("mschannel"), condition) %>% dplyr::distinct(),
                      by="condition") %>%
    dplyr::mutate(replicate_shift = rnorm(n(), 0, replicate_sigma)) %>%
    dplyr::inner_join(quantobjs_df) %>%
    dplyr::mutate(intensity = 2^(value + replicate_shift + labu_base_shift + labu_shift))
}

gen_msdata <- function(model_def, mschannels_df,
                       msprobe="msprobe", mschannel=msprobe,
                       msfraction = "msfraction", mstag = "mstag",
                       modelobject = "protgroup", quantobject = NA_character_,
                       nmodelobjects = 10) {
    res <- list()
    objects <- gen_objects(modelobject, n = nmodelobjects, quantobject = quantobject)
    res[[paste0(modelobject, 's')]] <- objects$modelobjects
    res$object_effects <- gen_modelobj_effects(objects$modelobjects, model_def$effects)
    res$interactions <- gen_interactions(model_def, res$object_effects)
    intensities_cols <- c("intensity")
    #if (rlang::has_name(mschannels_df, "msfraction")) {
    #    intensities_cols[[msfraction]] <- "msfraction"
    #}
    #if (rlang::has_name(mschannels_df, "mstag")) {
    #    intensities_cols[[mstag]] <- "mstag"
    #s}
    if (!is.na(quantobject)) {
        if (rlang::has_name(mschannels_df, "msfraction")) {
            objects$quantobjects <- dplyr::mutate(objects$quantobjects,
                                                  msfraction := sample(unique(mschannels_df$msfraction), n(), replace=TRUE))
        }
        res[[paste0(quantobject, 's')]] <- objects$quantobjects
        if (rlang::has_name(mschannels_df, "msfraction") && msfraction != "msfraction") {
            res[[paste0(quantobject, 's')]] <- dplyr::rename(res[[paste0(quantobject, 's')]],
                                                             !!sym(msfraction) := msfraction)
        }
        intensities_cols[[paste0(quantobject, "_id")]] <- "quantobj_id"
        intensities_cols[[mschannel]] <- "mschannel"
        res[[paste0(quantobject, "_intensities")]] <- gen_quantobj_intensities(res$interactions, objects$quantobjects, mschannels_df) %>%
            dplyr::select_at(intensities_cols)
        res[[paste0(modelobject, "2", quantobject)]] <- dplyr::select_at(objects$quantobjects, paste0(c(modelobject, quantobject), "_id")) %>%
            dplyr::mutate(is_specific = TRUE)
    } else {
        intensities_cols[[paste0(modelobject, "_id")]] <- "modelobj_id"
        intensities_cols[[msprobe]] <- "msprobe"
        res[[paste0(modelobject, "_intensities")]] <- gen_modelobj_intensities(res$interactions, mschannels_df) %>%
            dplyr::select_at(intensities_cols)
    }
    mschan_cols <- c("condition")
    mschan_cols[[msprobe]] <- "msprobe"
    if ((msprobe != mschannel) && rlang::has_name(mschannels_df, "mschannel")) {
        mschan_cols[[mschannel]] <- "mschannel"
    }
    if (rlang::has_name(mschannels_df, "mstag")) { mschan_cols[[mstag]] <- "mstag" }
    if (rlang::has_name(mschannels_df, "msfraction")) { mschan_cols[[msfraction]] <- "msfraction" }
    if (rlang::has_name(mschannels_df, "replicate")) { mschan_cols <- c(mschan_cols, "replicate") }
    res[[paste0(mschannel, 's')]] <- dplyr::select_at(mschannels_df, mschan_cols)
    return(res)
}
