context("data preration pipleine")

library(checkmate)

source("gen-msdata.R")

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

for (modelobj in c("protgroup", "ptmngroup")) {
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

test_that(paste0(modelobj, "/", modelobj, " model, no msfractions, ",
          if_else(!is.na(mstag), "", "no "), "mstags, specifying ", msprobe, "s"), {
    orig_msdata <- gen_msdata(model_def, msprobes_df, msprobe = msprobe,
                              modelobject = modelobj, nmodelobjects = 3)
    msprobes_dfname <- paste0(msprobe, "s")
    bad_msdata <- orig_msdata
    bad_msdata[[msprobes_dfname]] <- NULL
    expect_error(import_msglm_data(bad_msdata, model_def, mscalib, modelobject=modelobj),
                 "Cannot autodetect MS probes")
    if (!(msprobe %in% c("msrun", "mschannel"))) {
        expect_error(import_msglm_data(orig_msdata, model_def, mscalib, modelobject=modelobj),
                    "Cannot autodetect MS channels")
        # add raw_file column to facilitate autodetection
        orig_msdata[[msprobes_dfname]] <- dplyr::mutate(orig_msdata[[msprobes_dfname]],
                                                        raw_file = paste0(!!sym(msprobe), ".raw"))
    }

    if (modelobj != "protgroup") {
        # expecting protgroups by default
        expect_error(import_msglm_data(orig_msdata, model_def, mscalib),
                     "msdata\\$protgroups not found")
    }

    msdata <- import_msglm_data(orig_msdata, model_def, mscalib, modelobject=modelobj)
    expect_s3_class(msdata, "msglm_data_collection")
    expect_names(names(msdata), must.include = c("msentities", msprobes_dfname,
                 paste0(modelobj, "_intensities"), paste0(modelobj, "s")))
    expect_tibble(msdata[[msprobes_dfname]])
    expect_names(names(msdata[[msprobes_dfname]]), must.include = c("condition", msprobe))
    if (msprobe %in% c("msprobe", "msexperiment")) {
        # don't include mschannel/msrun frames
        expect_names(names(msdata), disjunct.from = c("msruns", "mschannels"))
    }
    modelobj_intens_df <- msdata[[paste0(modelobj, "_intensities")]]
    expect_tibble(modelobj_intens_df)
    expect_names(colnames(modelobj_intens_df),
                 must.include = c(paste0(modelobj, "_id"), msprobe, "intensity"))

    expect_equal(msdata$msentities, c(quantobject=modelobj, modelobject=modelobj,
        condition="condition",
        msexperiment=if_else(is.na(mstag), msprobe, NA_character_), msprobe=msprobe,
        msrun=if_else(is.na(mstag), msprobe, NA_character_), mschannel=msprobe,
        mstag=mstag, msfraction=NA_character_))

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
        "Nprobes", "Nobservations",
        "Nmschannels", "mschannel_shift",
        "Neffects", "NobjEffects", "NbatchEffects", "NobjBatchEffects",
        "Niactions", "Nquanted", "Nmissed",
        "iactXobjeff_Nw", "obsXobjbatcheff_Nw"))

    model_data2 <- msglm_data(model_def, msdata, 2L)
    expect_s3_class(model_data2, "msglm_model_data")
    expect_equal(model_data2$object_id, 2L)
    expect_equal(unique(model_data2$objects$object_id), 2L)

    expect_error(msglm_data(model_def, msdata, 5L), "Objects not found")
})

test_that(paste0(modelobj, "/pepmodstate model, no msfractions, ",
          if_else(!is.na(mstag), "", "no "), "mstags, specifying ", msprobe, "s"), {
    mschannels_df <- dplyr::mutate(msprobes_df, mschannel = msprobe)
    orig_msdata <- gen_msdata(model_def, mschannels_df, msprobe = msprobe, mschannel = msprobe,
                              modelobject = modelobj, quantobject = "pepmodstate", nmodelobjects=3)
    msprobes_dfname <- paste0(msprobe, "s")
    orig_msdata[[msprobes_dfname]] <- dplyr::rename(msprobes_df, !!sym(msprobe) := msprobe) %>%
        dplyr::mutate(raw_file = paste0(!!sym(msprobe), ".raw"))
    msdata <- import_msglm_data(orig_msdata, model_def, mscalib,
                                modelobject=modelobj, quantobject = "pepmodstate")
    expect_s3_class(msdata, "msglm_data_collection")
    expect_equal(msdata$msentities, c(quantobject="pepmodstate", modelobject=modelobj,
        condition="condition",
        msexperiment=if_else(is.na(mstag), msprobe, NA_character_), msprobe=msprobe,
        msrun=if_else(is.na(mstag), msprobe, NA_character_), mschannel=msprobe,
        mstag=mstag, msfraction=NA_character_))

    msdata[[paste0(msprobe, "_shifts")]] <- msprobe_shifts_df
    model_data1 <- msglm_data(model_def, msdata, 1L)
    expect_s3_class(model_data1, "msglm_model_data")
    expect_names(names(model_data1), must.include = c("model_def", "object_id", "msentities",
                 "objects", "subobjects", "msprobes", "mschannels"))
    expect_tibble(model_data1$objects)
    expect_equal(model_data1$object_id, 1L)
    expect_equal(unique(model_data1$objects$object_id), 1L)
    expect_equal(unique(model_data1$subobjects$index_object), 1L)
})

}

for (msprobe in if(is.na(mstag)) c("msprobe", "msexperiment") else "msprobe") {

test_that(paste0(modelobj, "/pepmodstate model, msfractions, ",
          if_else(!is.na(mstag), "", "no "), "mstags, specifying ", msprobe, "s"), {

    mschannel <- c(msprobe = "mschannel", msexperiment = "msrun")[[msprobe]]
    mschannels_df <- dplyr::full_join(msprobes_df, tibble(msfraction = c(1L, 2L)), by=character()) %>%
                     dplyr::mutate(mschannel = paste0(msprobe, "_F", msfraction))
    if (!is.na(mstag)) {
        mschannels_df <- dplyr::mutate(mschannels_df,
                                       msrun = mschannel,
                                       mschannel = paste0(mschannel, "_", mstag))
    }
    orig_msdata <- gen_msdata(model_def, mschannels_df, msprobe = msprobe, mschannel = mschannel,
                              modelobject = modelobj, quantobject = "pepmodstate", nmodelobjects=3)
    msprobes_dfname <- paste0(msprobe, "s")
    mschannels_dfname <- paste0(mschannel, "s")
    # without specifying mschannels, by with hint that msfractions exist
    bad_msdata <- orig_msdata
    bad_msdata[[paste0(mschannel, "s")]] <- NULL
    expect_error(msglm::import_msglm_data(bad_msdata, model_def, mscalib,
                                          modelobject=modelobj, msfraction="msfraction"),
                "Cannot autodetect MS probes")

    # without specifying msprobes
    orig_msdata1 <- orig_msdata
    msdata <- import_msglm_data(orig_msdata1, model_def, mscalib,
                                modelobject=modelobj, quantobject="pepmodstate")
    expect_s3_class(msdata, "msglm_data_collection")
    expect_names(names(msdata), must.include = c("msentities", msprobes_dfname, mschannels_dfname,
                 "pepmodstate_intensities", paste0(modelobj, "s"), "pepmodstates",
                 paste0(modelobj, "2pepmodstate")))
    expect_equal(msdata$msentities, c(quantobject="pepmodstate", modelobject=modelobj,
        condition="condition",
        msexperiment=if_else(is.na(mstag), msprobe, NA_character_), msprobe=msprobe,
        msrun=if_else(is.na(mstag), mschannel, NA_character_), mschannel=mschannel,
        mstag=mstag, msfraction="msfraction"))
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
    msdata <- msglm::import_msglm_data(orig_msdata1b, model_def, mscalib,
                                       modelobject=modelobj, quantobject="pepmodstate")
    expect_s3_class(msdata, "msglm_data_collection")
    expect_names(names(msdata), must.include = c("msentities", msprobes_dfname, mschannels_dfname,
                 "pepmodstate_intensities",
                 paste0(modelobj, "s"), "pepmodstates", paste0(modelobj, "2pepmodstate")))
    expect_equal(msdata$msentities, c(quantobject="pepmodstate", modelobject=modelobj,
        condition="condition",
        msexperiment=if_else(is.na(mstag), msprobe, NA_character_), msprobe=msprobe,
        msrun=if_else(is.na(mstag), mschannel, NA_character_), mschannel=mschannel,
        mstag=mstag, msfraction="msfraction"))
    expect_tibble(msdata[[msprobes_dfname]], nrows=nrow(msprobes_df))
    expect_tibble(msdata[[mschannels_dfname]], nrows=nrow(mschannels_df))

    # without specifying mschannels
    orig_msdata2 <- orig_msdata
    orig_msdata2[[msprobes_dfname]] <- dplyr::rename(msprobes_df, !!sym(msprobe) := msprobe)
    orig_msdata2[[paste0(mschannel, "s")]] <- NULL
    expect_error(msglm::import_msglm_data(orig_msdata2, model_def, mscalib,
                                          modelobject = modelobj, quantobject = "pepmodstate"),
                 "Cannot autodetect MS channels data")

    # with both
    orig_msdata3 <- orig_msdata
    orig_msdata3[[msprobes_dfname]] <- dplyr::rename(msprobes_df, !!sym(msprobe) := msprobe)
    msdata <- import_msglm_data(orig_msdata3, model_def, mscalib,
                                modelobject=modelobj, quantobject="pepmodstate")
    expect_s3_class(msdata, "msglm_data_collection")
    expect_equal(msdata$msentities, c(quantobject="pepmodstate", modelobject= modelobj,
        condition="condition",
        msexperiment=if_else(is.na(mstag), msprobe, NA_character_), msprobe=msprobe,
        msrun=if_else(is.na(mstag), mschannel, NA_character_), mschannel=mschannel,
        mstag=mstag, msfraction="msfraction"))

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
                 "objects", "subobjects", "msprobes", "mschannels"))
    expect_tibble(model_data1$objects)
    expect_equal(model_data1$object_id, 1L)
    expect_equal(unique(model_data1$objects$object_id), 1L)
    expect_equal(unique(model_data1$subobjects$object_id), 1L)

    stan_data1 <- to_standata(model_data1)
    expect_list(stan_data1)
    expect_names(names(stan_data1), must.include=c(
        "Nconditions", "Nobjects", "Nsubobjects",
        "Nprobes", "Nobservations", "Nsubobservations",
        "Nmschannels", "mschannel_shift",
        "Neffects", "NobjEffects", "NbatchEffects", "NobjBatchEffects",
        "NquantBatchEffects", "NsubobjBatchEffects",
        "Niactions", "Nquanted", "Nmissed",
        "iactXobjeff_Nw", "obsXobjbatcheff_Nw", "subobsXsubobjbatcheff_Nw"))
})

}}}
