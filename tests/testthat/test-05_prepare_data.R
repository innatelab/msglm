context("data preration pipleine")

library(checkmate)

source(test_path("gen-msdata.R"))

# generate simplest model
conditions_df = tibble(infection = factor(c("control", "virus", "control", "virus"), levels=c("control", "virus")),
                       treatment = factor(c("mock", "mock", "drug", "drug"), levels=c("mock", "drug")),
                       condition = paste0(infection, "_", treatment))

conditionXeffect <- model.matrix(~ 1 + infection * treatment, data = conditions_df)
conditionXeffect <- conditionXeffect[, -1, drop=FALSE]
rownames(conditionXeffect) <- conditions_df$condition

effects_df = tibble(effect = colnames(conditionXeffect)) %>%
    dplyr::filter(effect != "(Intercept)") %>%
    dplyr::mutate(is_interaction = str_detect(effect, ':'),
                  prior_mean = 0,
                  prior_tau = if_else(is_interaction, 0.5, 1.0),
                  prior_df1 = if_else(is_interaction, 4, 2),
                  prior_df2 = if_else(is_interaction, 4, 2))

model_def <- msglm_model(conditionXeffect, conditions_df, effects_df)

mscalib <- read_mscalib_json(test_path("test_mscalib.json")) %>% msglm::convert_logintensityBase(2)

for (obj in c("protgroup", "ptmngroup")) {
for (mstag in c(NA_character_, "mstag")) {

msprobes_df <- tidyr::expand_grid(condition = conditions_df$condition,
                                  replicate = 1:3) %>%
               dplyr::mutate(msprobe = paste0(condition, "_", replicate))
if (!is.na(mstag)) {
    msprobes_df <- dplyr::mutate(msprobes_df,
                                 !!sym(mstag) := if_else(str_detect(condition, "drug"), "H", "L"),
                                 msexperiment = msprobe,
                                 msprobe = paste0(msexperiment, "_", mstag))
}

for (msprobe in if (is.na(mstag)) c("msprobe", "msexperiment", "msrun", "mschannel") else c("msprobe", "mschannel")) {

msprobe_shifts_df <- tibble(!!sym(msprobe) := msprobes_df$msprobe,
                            !!sym(paste0("total_", msprobe, "_shift")) := rnorm(nrow(msprobes_df), 0.0, 0.1))

test_that(paste0(obj, "/", obj, " model, no msfractions, ",
          if_else(!is.na(mstag), "", "no "), "mstags, specifying ", msprobe, "s"), {
    orig_msdata <- gen_msdata(model_def, dplyr::mutate(msprobes_df, mschannel=msprobe),
                              msprobe = msprobe,
                              object = obj, nobjects = 3)
    msprobes_dfname <- paste0(msprobe, "s")
    if (!(msprobe %in% c("msrun", "mschannel"))) {
        # add rawfile column to facilitate autodetection
        orig_msdata[[msprobes_dfname]] <- dplyr::mutate(orig_msdata[[msprobes_dfname]],
                                                        rawfile = paste0(!!sym(msprobe), ".raw"))
    }

    test_that(paste0("import fails if no ", msprobe, " data provided"), {
        bad_msdata <- orig_msdata
        bad_msdata[[msprobes_dfname]] <- NULL
        expect_error(import_msglm_data(bad_msdata, model_def, mscalib, object=obj),
                     "Cannot autodetect MS probes")
        bad_msdata <- orig_msdata
        if (!(msprobe %in% c("msrun", "mschannel"))) {
            bad_msdata[[msprobes_dfname]] <- dplyr::mutate(orig_msdata[[msprobes_dfname]],
                                                           rawfile = NULL)
            expect_error(import_msglm_data(bad_msdata, model_def, mscalib, object=obj),
                         "Cannot autodetect MS channels")
        }
    })

    if (obj != "protgroup") {
        test_that("expecting protgroups by default", {
            expect_error(import_msglm_data(orig_msdata, model_def, mscalib=mscalib),
                        "msdata\\$protgroups not found")
        })
    }

    test_that(paste0("bad ", obj, "_idents frame (", obj, " X ", msprobe, ") is skipped without error"), {
        bad_idents_msdata <- orig_msdata
        bad_idents_msdata[[paste0(obj, "_idents")]] <- tibble()
        msdata <- import_msglm_data(bad_idents_msdata, model_def, object=obj, mscalib=mscalib)
        expect_s3_class(msdata, "msglm_data_collection")
        expect_names(names(msdata), disjunct.from = paste0(obj, "_idents"))
    })

    msdata <- import_msglm_data(orig_msdata, model_def, object=obj, mscalib=mscalib)
    expect_s3_class(msdata, "msglm_data_collection")
    expect_names(names(msdata), must.include = c("msentities", msprobes_dfname,
                 paste0(obj, "_intensities"), paste0(obj, "s")))
    expect_tibble(msdata[[msprobes_dfname]])
    expect_names(names(msdata[[msprobes_dfname]]), must.include = c("condition", msprobe))
    if (msprobe %in% c("msprobe", "msexperiment")) {
        # don't include mschannel/msrun frames
        expect_names(names(msdata), disjunct.from = c("msruns", "mschannels"))
    }
    obj_intens_df <- msdata[[paste0(obj, "_intensities")]]
    expect_tibble(obj_intens_df)
    expect_names(colnames(obj_intens_df),
                 must.include = c(paste0(obj, "_id"), msprobe, "intensity"))

    expect_equal(msdata$msentities, c(object=obj, quantobject=obj,
        condition="condition",
        msexperiment=if_else(is.na(mstag), msprobe, NA_character_), msprobe=msprobe,
        msrun=if_else(is.na(mstag), msprobe, NA_character_), mschannel=msprobe,
        mstag=mstag, msfraction=NA_character_, msprotocol=NA_character_))

    msdata[[paste0(msprobe, "_shifts")]] <- msprobe_shifts_df

    model_data1 <- msglm_data(model_def, msdata, 1L)
    expect_s3_class(model_data1, "msglm_model_data")
    expect_names(names(model_data1), must.include = c("model_def", "object_id", "msentities",
                 "objects", "msprobes", "mschannels"))
    expect_tibble(model_data1$objects)
    expect_equal(model_data1$object_id, 1L)
    expect_equal(unique(model_data1$objects$object_id), 1L)
    stan_data1 <- to_standata(model_data1)
    expect_list(stan_data1)
    expect_names(names(stan_data1), must.include=c("Nconditions", "Nobjects",
        "Nprobes", "NobjProbes",
        "Nmschannels", "mschannel_shift",
        "Neffects", "NobjEffects", "NbatchEffects", "NobjBatchEffects",
        "NobjConditions", "Nquanted", "Nmissed",
        "quant2mschannel", "miss2mschannel",
        "quant2obj_probe", "miss2obj_probe", "qData",
        "obj_condXeff_Nw", "obj_probeXbatcheff_Nw"))

    model_data2 <- msglm_data(model_def, msdata, 2L)
    expect_s3_class(model_data2, "msglm_model_data")
    expect_equal(model_data2$object_id, 2L)
    expect_equal(unique(model_data2$objects$object_id), 2L)

    expect_error(msglm_data(model_def, msdata, 5L), "Objects not found")
})

test_that(paste0(obj, "/pepmodstate model, no msfractions, ",
          if_else(!is.na(mstag), "", "no "), "mstags, specifying ", msprobe, "s"), {
    mschannels_df <- dplyr::mutate(msprobes_df, mschannel = msprobe)
    orig_msdata <- gen_msdata(model_def, mschannels_df, msprobe = msprobe, mschannel = msprobe,
                              object = obj, quantobject = "pepmodstate", nobjects=3)
    msprobes_dfname <- paste0(msprobe, "s")
    orig_msdata[[msprobes_dfname]] <- dplyr::rename(msprobes_df, !!sym(msprobe) := msprobe) %>%
        dplyr::mutate(rawfile = paste0(!!sym(msprobe), ".raw"))
    if (msprobe != "msprobe") {
        orig_msdata[[msprobes_dfname]]$msprobe <- NULL
    }
    msdata <- import_msglm_data(orig_msdata, model_def,
                                object = obj, quantobject = "pepmodstate",
                                mscalib = mscalib)
    expect_s3_class(msdata, "msglm_data_collection")
    expect_equal(msdata$msentities, c(object=obj, quantobject="pepmodstate",
        condition="condition",
        msexperiment=if_else(is.na(mstag), msprobe, NA_character_), msprobe=msprobe,
        msrun=if_else(is.na(mstag), msprobe, NA_character_), mschannel=msprobe,
        mstag=mstag, msfraction=NA_character_, msprotocol=NA_character_))

    msdata[[paste0(msprobe, "_shifts")]] <- msprobe_shifts_df
    model_data1 <- msglm_data(model_def, msdata, 1L)
    expect_s3_class(model_data1, "msglm_model_data")
    expect_names(names(model_data1), must.include = c("model_def", "object_id", "msentities",
                 "objects", "quantobjects", "msprobes", "mschannels"))
    expect_tibble(model_data1$objects)
    expect_equal(model_data1$object_id, 1L)
    expect_equal(unique(model_data1$objects$object_id), 1L)
    expect_equal(unique(model_data1$quantobjects$index_object), 1L)
})

}

for (msprobe in if(is.na(mstag)) c("msprobe", "msexperiment") else "msprobe") {

test_that(paste0(obj, "/pepmodstate model, msfractions, ",
          if_else(!is.na(mstag), "", "no "), "mstags, specifying ", msprobe, "s"), {

    mschannel <- c(msprobe = "mschannel", msexperiment = "msrun")[[msprobe]]
    mschannels_df <- dplyr::full_join(msprobes_df, tibble(msfraction = c(1L, 2L)), by=character()) %>%
            dplyr::mutate(mschannel = paste0(condition, "_", replicate, "_F", msfraction))
    if (!is.na(mstag)) {
        mschannels_df <- dplyr::mutate(mschannels_df,
                                       msrun = mschannel,
                                       mschannel = paste0(mschannel, "_", mstag))
    }
    mschannels_df <- dplyr::mutate(mschannels_df,
                                   rawfile = paste0(condition, '_', replicate,
                                                    '_F', msfraction, ".raw"))
    orig_msdata <- gen_msdata(model_def, mschannels_df, msprobe = msprobe, mschannel = mschannel,
                              object = obj, quantobject = "pepmodstate", nobjects=3)
    msprobes_dfname <- paste0(msprobe, "s")
    mschannels_dfname <- paste0(mschannel, "s")
    # without specifying mschannels, by with hint that msfractions exist
    bad_msdata <- orig_msdata
    bad_msdata[[paste0(mschannel, "s")]] <- NULL
    expect_error(msglm::import_msglm_data(bad_msdata, model_def, mscalib,
                                          object=obj, msfraction="msfraction"),
                "Cannot autodetect MS probes")

    # without specifying msprobes
    orig_msdata1 <- orig_msdata
    msdata <- import_msglm_data(orig_msdata1, model_def,
                                object=obj, quantobject="pepmodstate",
                                mscalib=mscalib)
    expect_s3_class(msdata, "msglm_data_collection")
    expect_names(names(msdata), must.include = c("msentities", msprobes_dfname, mschannels_dfname,
                 "pepmodstate_intensities", paste0(obj, "s"), "pepmodstates",
                 paste0(obj, "2pepmodstate")))
    expect_equal(msdata$msentities, c(object=obj, quantobject="pepmodstate",
        condition="condition",
        msexperiment=if_else(is.na(mstag), msprobe, NA_character_), msprobe=msprobe,
        msrun=if_else(is.na(mstag), mschannel, NA_character_), mschannel=mschannel,
        mstag=mstag, msfraction="msfraction", msprotocol = NA_character_))
    expect_tibble(msdata[[msprobes_dfname]], nrows=nrow(msprobes_df))
    expect_names(names(msdata[[msprobes_dfname]]), must.include = c("condition", msprobe))
    expect_tibble(msdata[[mschannels_dfname]], nrows=nrow(mschannels_df))
    expect_names(names(msdata[[mschannels_dfname]]), must.include = c(msprobe, mschannel))
    expect_tibble(msdata$pepmodstate_intensities)
    expect_names(names(msdata$pepmodstate_intensities),
                 must.include = c("pepmodstate_id", mschannel, "intensity"))

    # without specifying msprobes and without msprobe column in mschannels
    orig_msdata1b <- orig_msdata
    orig_msdata1b[[mschannels_dfname]] <- dplyr::mutate(orig_msdata1b[[mschannels_dfname]], !!sym(msprobe) := NULL)
    msdata <- msglm::import_msglm_data(orig_msdata1b, model_def,
                                       object = obj, quantobject = "pepmodstate",
                                       mscalib = mscalib)
    expect_s3_class(msdata, "msglm_data_collection")
    expect_names(names(msdata), must.include = c("msentities", msprobes_dfname, mschannels_dfname,
                 "pepmodstate_intensities",
                 paste0(obj, "s"), "pepmodstates", paste0(obj, "2pepmodstate")))
    expect_equal(msdata$msentities, c(object=obj, quantobject="pepmodstate",
        condition="condition",
        msexperiment=if_else(is.na(mstag), msprobe, NA_character_), msprobe=msprobe,
        msrun=if_else(is.na(mstag), mschannel, NA_character_), mschannel=mschannel,
        mstag=mstag, msfraction="msfraction", msprotocol=NA_character_))
    expect_tibble(msdata[[msprobes_dfname]], nrows=nrow(msprobes_df))
    expect_tibble(msdata[[mschannels_dfname]], nrows=nrow(mschannels_df))

    # with obj2pepmod but without obj2pepmodstate
    orig_msdata1c <- orig_msdata
    obj2pms_dfname <- paste0(obj, "2pepmodstate")
    orig_msdata1c[[obj2pms_dfname]] <- NULL
    msdata <- msglm::import_msglm_data(orig_msdata1c, model_def,
                                       object = obj, quantobject = "pepmodstate",
                                       mscalib = mscalib)
    expect_s3_class(msdata, "msglm_data_collection")
    # check that obj2pepmodstate is automatically generated
    expect_names(names(msdata), must.include = obj2pms_dfname)
    expect_tibble(msdata[[obj2pms_dfname]], nrows=nrow(orig_msdata[[obj2pms_dfname]]))
    expect_names(names(msdata[[obj2pms_dfname]]),
                 must.include = c(paste0(c(obj, "pepmodstate"), "_id"), "is_specific"))

    # without specifying mschannels
    orig_msdata2 <- orig_msdata
    orig_msdata2[[msprobes_dfname]] <- dplyr::rename(msprobes_df, !!sym(msprobe) := msprobe)
    orig_msdata2[[paste0(mschannel, "s")]] <- NULL
    expect_error(msglm::import_msglm_data(orig_msdata2, model_def,
                                          object = obj, quantobject = "pepmodstate",
                                          mscalib = mscalib),
                 "Cannot autodetect MS channels data")

    # with both
    orig_msdata3 <- orig_msdata
    orig_msdata3[[msprobes_dfname]] <- dplyr::rename(msprobes_df, !!sym(msprobe) := msprobe)
    msdata <- import_msglm_data(orig_msdata3, model_def,
                                object = obj, quantobject = "pepmodstate",
                                mscalib = mscalib)
    expect_s3_class(msdata, "msglm_data_collection")
    expect_equal(msdata$msentities, c(object=obj, quantobject="pepmodstate",
        condition="condition",
        msexperiment=if_else(is.na(mstag), msprobe, NA_character_), msprobe=msprobe,
        msrun=if_else(is.na(mstag), mschannel, NA_character_), mschannel=mschannel,
        mstag=mstag, msfraction="msfraction", msprotocol=NA_character_))

    if (mschannel != msprobe) {
        msdata[[paste0(mschannel, "_shifts")]] <- tibble(
            !!sym(mschannel) := mschannels_df$mschannel,
            !!sym(paste0("total_", mschannel, "_shift")) := rnorm(nrow(mschannels_df), 0.0, 0.1))
    } else {
        msdata[[paste0(msprobe, "_shifts")]] <- msprobe_shifts_df
    }
    model_data1 <- msglm_data(model_def, msdata, 1L)
    expect_s3_class(model_data1, "msglm_model_data")
    expect_names(names(model_data1), must.include = c("model_def", "object_id", "msentities",
                 "objects", "quantobjects", "msprobes", "mschannels"))
    expect_tibble(model_data1$objects)
    expect_equal(model_data1$object_id, 1L)
    expect_equal(unique(model_data1$objects$object_id), 1L)
    expect_equal(unique(model_data1$quantobjects$object_id), 1L)

    stan_data1 <- to_standata(model_data1)
    expect_list(stan_data1)
    expect_names(names(stan_data1), must.include=c(
        "Nconditions", "Nobjects", "Nquantobjects",
        "Nprobes", "NobjProbes", "NqobjChannels",
        "Nmschannels", "mschannel_shift",
        "Neffects", "NobjEffects", "NbatchEffects", "NobjBatchEffects",
        "NquantBatchEffects", "NqobjBatchEffects",
        "NobjConditions", "Nquanted", "Nmissed",
        "quant2qobj_channel", "miss2qobj_channel", "qData",
        "obj_condXeff_Nw", "obj_probeXbatcheff_Nw", "qobj_channelXqbatcheff_Nw"))
})

}}}

test_that("Technical MS replicates are supported", {
    msprobes_df <- tidyr::expand_grid(condition = conditions_df$condition,
                                      replicate = 1:4) %>%
               dplyr::mutate(msprobe = paste0(condition, "_", replicate))
    mschannels_df <- tidyr::expand_grid(msprobe = msprobes_df$msprobe,
                                        tech_replicate = 1:2) %>%
               dplyr::left_join(dplyr::select(msprobes_df, msprobe, condition), by="msprobe") %>%
               dplyr::mutate(mschannel = paste0(msprobe, "_", tech_replicate))
    mschannel_shifts_df <- tibble(mschannel = mschannels_df$mschannel,
                                  total_mschannel_shift = rnorm(nrow(mschannels_df), 0.0, 0.1))

    test_that("Technical MS replicates are supported on pepmodstate level", {
        orig_msdata <- gen_msdata(model_def, mschannels_df, msprobe = "msprobe", mschannel = "mschannel",
                                  object = "protgroup", quantobject = "pepmodstate", nobjects=3)
        msdata <- import_msglm_data(orig_msdata, model_def,
                                    object="protgroup", quantobject="pepmodstate",
                                    mscalib=mscalib)
        msdata$mschannel_shifts <- mschannel_shifts_df
        model_data <- msglm_data(model_def, msdata, 1L)
        expect_s3_class(model_data, "msglm_model_data")
        expect_equal(nrow(model_data$msprobes), nrow(msprobes_df))
        expect_equal(nrow(model_data$mschannels), nrow(mschannels_df))
        stan_data <- to_standata(model_data)
        expect_list(stan_data)

        test_that("Technical MS replicates + quant batch effects work", {
            mschannels_be_df <- dplyr::mutate(dplyr::left_join(mschannels_df, dplyr::select(msprobes_df, msprobe, replicate),
                                                               by="msprobe"),
                                              batch = LETTERS[pmin(3L, 1L + (replicate + tech_replicate) %% 4)])
            batch_effect_mtx <- model.matrix(~ 1 + batch, data = mschannels_be_df)
            batch_effect_mtx <- batch_effect_mtx[, -1, drop=FALSE]
            dimnames(batch_effect_mtx) <- list(mschannel = mschannels_be_df$mschannel,
                                               quant_batch_effect = colnames(batch_effect_mtx))

            model_def_be <- set_batch_effects(model_def, batch_effect_mtx, applies_to = "quantobject")
            expect_matrix(model_def_be$mschannelXquantBatchEffect, any.missing = FALSE,
                          nrows=nrow(mschannels_be_df), ncols=n_distinct(mschannels_be_df$batch) - 1L)
            expect_data_frame(model_def_be$quant_batch_effects, nrows=n_distinct(mschannels_be_df$batch) - 1L)
            model_data_be <- msglm_data(model_def_be, msdata, 1L)
            expect_data_frame(model_data_be$quantobject_batch_effects)
            expect_names(names(model_data_be$quantobject_batch_effects),
                         must.include = c("quantobject_batch_effect", "quant_batch_effect", "index_quantobject"))
            expect_matrix(model_data_be$quantobject_mschannelXquant_batch_effect, any.missing = FALSE,
                          nrows = nrow(model_data_be$msdata),
                          ncols = nrow(model_data_be$quantobject_batch_effects))
            stan_data_be <- to_standata(model_data_be)
            expect_equal(stan_data_be$NquantBatchEffects, n_distinct(mschannels_be_df$batch) - 1L)
            expect_equal(stan_data_be$NqobjBatchEffects, nrow(model_data_be$quantobject_batch_effects))
        })
    })

    test_that("Technical MS replicates are supported on protgroup level", {
        orig_msdata <- gen_msdata(model_def, mschannels_df, msprobe = "msprobe", mschannel = "mschannel",
                                  object = "protgroup", quantobject = NA_character_, nobjects=3)
        msdata <- import_msglm_data(orig_msdata, model_def,
                                    object="protgroup", quantobject="protgroup",
                                    mscalib=mscalib, verbose=TRUE)
        msdata$mschannel_shifts <- mschannel_shifts_df
        model_data <- msglm_data(model_def, msdata, 1L)
        expect_s3_class(model_data, "msglm_model_data")
        expect_data_frame(model_data$msprobes, nrows=nrow(msprobes_df))
        expect_data_frame(model_data$mschannels, nrows=nrow(mschannels_df))
        stan_data <- to_standata(model_data)
        expect_list(stan_data)
        expect_equal(stan_data$Nprobes, nrow(msprobes_df))
        expect_equal(stan_data$Nmschannels, nrow(mschannels_df))
        expect_integer(stan_data$quant2mschannel,
                       lower=1, upper=stan_data$Nmschannels,
                       len=sum(model_data$msdata$is_observed))
        expect_integer(stan_data$miss2mschannel,
                       lower=1, upper=stan_data$Nmschannels,
                       len=sum(!model_data$msdata$is_observed))
    })
})
