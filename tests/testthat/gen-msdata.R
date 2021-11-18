require(tibble)
require(dplyr)
require(stringr)

gen_objects <- function(object, n, quantobject = NA_character_,
                        avg_npepmods = 5, avg_ncharges = 1.5,
                        labu_mean = 10, labu_sigma = 5, quantobj_labu_shift_shape = 3) {
    obj_id <- dplyr::sym(paste0(object, "_id"))
    obj_label <- dplyr::sym(paste0(object, "_label"))
    objs_df = tibble(obj_id = seq_len(n),
                     !!obj_id := obj_id,
                     !!obj_label := paste0(object, "_", !!obj_id),
                     labu_base = rnorm(n, mean=labu_mean, sd=labu_sigma))
    if (!is.na(quantobject)) {
        objs_df <- dplyr::mutate(objs_df, npepmods = 1 + rpois(length(labu_base), avg_npepmods - 1))
        quantobj_id <- dplyr::sym(paste0(quantobject, "_id"))
        quantobj_label <- dplyr::sym(paste0(quantobject, "_label"))
        pepmods_df <- dplyr::group_by(objs_df, !!obj_id, obj_id) %>%
            dplyr::group_modify(~{
                tibble(pepmod_label = paste0(.x[[obj_label]], "_", seq_len(.x$npepmods)),
                       labu_shift = c(0, -rgamma(.x$npepmods - 1L, quantobj_labu_shift_shape)),
                       ncharges = 1L + rpois(.x$npepmods, avg_ncharges - 1))
            }) %>% dplyr::ungroup() %>%
            dplyr::mutate(pepmod_id = row_number())
        quantobjs_df <- dplyr::group_by(pepmods_df, !!obj_id, obj_id,
                                        pepmod_id, pepmod_label) %>%
            dplyr::group_modify(~{
                tibble(charge = 2L + seq_len(.x$ncharges),
                       labu_shift = .x$labu_shift + c(0, -rgamma(.x$ncharges - 1L, quantobj_labu_shift_shape)))
            }) %>% dplyr::ungroup() %>%
            dplyr::mutate(quantobj_id = row_number(),
                          !!quantobj_id := quantobj_id,
                          !!quantobj_label := paste0(pepmod_label, ".", charge),
                          pepmod_label = NULL)
    } else {
        quantobjs_df <- NULL
        pepmods_df <- NULL
    }
    return(list(objects = objs_df, pepmods = pepmods_df, quantobjects = quantobjs_df))
}

gen_effects <- function(effects_df, effect_range = c(-10, 10)) {
    dplyr::group_by(effects_df, effect) %>%
    dplyr::summarise(value = prior_mean + rnorm(1) * prior_tau * abs(rnorm(1) * rnorm(1))/
                             (rgamma(1, prior_df1/2, prior_df1/2)*rgamma(1, prior_df2/2, prior_df2/2)),
                     .groups = "drop") %>%
    dplyr::mutate(value = pmax(effect_range[[1]], pmin(value, effect_range[[2]])))
}

gen_obj_effects <- function(objs_df, effects_df, effect_range = c(-10, 10)) {
    dplyr::group_by(objs_df, obj_id) %>%
    dplyr::group_modify(~{
        gen_effects(effects_df = effects_df, effect_range = effect_range) %>%
        dplyr::mutate(labu_base= .x$labu_base)
    }) %>% dplyr::ungroup()
}

gen_obj_conditions <- function(model_def, objXeffects_df) {
    dplyr::group_by(objXeffects_df, obj_id) %>%
    dplyr::group_modify(~{
        effect_values <- rlang::set_names(.x$value, .x$effect)
        objcond_values <- model_def$conditionXeffect %*% effect_values[colnames(model_def$conditionXeffect)]
        tibble(condition = rownames(model_def$conditionXeffect),
               value = as.vector(objcond_values) + .x$labu_base[[1]])
    }) %>% dplyr::ungroup()
}

gen_obj_intensities <- function(objconds_df, mschannels_df, labu_base_shift = 0.0,
                                replicate_sigma = 0.1, tech_replicate_sigma = 0.1) {
    dplyr::inner_join(objconds_df, dplyr::select(mschannels_df, msprobe, condition) %>% dplyr::distinct(), by="condition") %>%
    dplyr::mutate(replicate_shift = rnorm(n(), 0, replicate_sigma),
                  intensity = 2^(value + replicate_shift + labu_base_shift)) %>%
    dplyr::inner_join(dplyr::select(mschannels_df, mschannel, msprobe), by="msprobe") %>%
    dplyr::group_by(!!!syms(names(objconds_df)), msprobe) %>%
    dplyr::mutate(tech_replicate_shift = if(n() > 1) rnorm(n(), 0, tech_replicate_sigma) else 0) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(intensity = intensity * 2^tech_replicate_shift)
}

gen_quantobj_intensities <- function(objconds_df, quantobjs_df, mschannels_df,
                                     labu_base_shift = 0.0, replicate_sigma = 0.1) {
    dplyr::inner_join(objconds_df, dplyr::select(mschannels_df, msprobe, any_of("mschannel"), condition) %>% dplyr::distinct(),
                      by="condition") %>%
    dplyr::mutate(replicate_shift = rnorm(n(), 0, replicate_sigma)) %>%
    dplyr::inner_join(quantobjs_df) %>%
    dplyr::mutate(intensity = 2^(value + replicate_shift + labu_base_shift + labu_shift))
}

gen_msdata <- function(model_def, mschannels_df,
                       msprobe="msprobe", mschannel=msprobe,
                       msfraction = "msfraction", mstag = "mstag",
                       object = "protgroup", quantobject = NA_character_,
                       nobjects = 10) {
    res <- list()
    objects <- gen_objects(object, n = nobjects, quantobject = quantobject)
    res[[paste0(object, 's')]] <- objects$objects
    res$object_effects <- gen_obj_effects(objects$objects, model_def$effects)
    res$object_conditions <- gen_obj_conditions(model_def, res$object_effects)
    intensities_cols <- c("intensity")
    intensities_cols[[mschannel]] <- "mschannel"
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
        res$pepmods <- objects$pepmods
        if (rlang::has_name(mschannels_df, "msfraction") && msfraction != "msfraction") {
            res[[paste0(quantobject, 's')]] <- dplyr::rename(res[[paste0(quantobject, 's')]],
                                                             !!sym(msfraction) := msfraction)
        }
        intensities_cols[[paste0(quantobject, "_id")]] <- "quantobj_id"
        res[[paste0(quantobject, "_intensities")]] <- gen_quantobj_intensities(res$object_conditions, objects$quantobjects, mschannels_df) %>%
            dplyr::select_at(intensities_cols)
        res[[paste0(object, "2pepmod")]] <- dplyr::select_at(objects$pepmods, paste0(c(object, "pepmod"), "_id")) %>%
            dplyr::mutate(is_specific = TRUE)
        res[[paste0(object, "2", quantobject)]] <- dplyr::select_at(objects$quantobjects, paste0(c(object, quantobject), "_id")) %>%
            dplyr::mutate(is_specific = TRUE)
    } else {
        intensities_cols[[paste0(object, "_id")]] <- "obj_id"
        res[[paste0(object, "_intensities")]] <- gen_obj_intensities(res$object_conditions, mschannels_df) %>%
            dplyr::select_at(intensities_cols)
    }
    mschan_cols <- character()
    for (col in c('condition', 'replicate', 'tech_replicate', 'raw_file', 'rawfile')) {
        if (rlang::has_name(mschannels_df, col)) {
            mschan_cols <- append(mschan_cols, col)
        }
    }
    renamed_cols <- c("mstag", "msfraction", "msprobe")
    if (msprobe != mschannel) renamed_cols <- append(renamed_cols, "mschannel")
    for (col in renamed_cols) {
        colval <- get(col, inherits = FALSE)
        if (rlang::has_name(mschannels_df, col)) {
            mschan_cols[[colval]] <- col
        }
    }
    res[[paste0(mschannel, 's')]] <- dplyr::select_at(mschannels_df, mschan_cols)
    return(res)
}
