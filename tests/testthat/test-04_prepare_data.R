context("msdata import")

source("gen-msdata.R")

# generate simplest model
conditions_df = tibble(infection = factor(c("control", "virus", "control", "virus"), levels=c("control", "virus")),
                       treatment = factor(c("mock", "mock", "drug", "drug"), levels=c("mock", "drug")),
                       condition = paste0(infection, "_", treatment))

msprobes_df <- tidyr::expand_grid(condition = conditions_df$condition,
                                  replicate = 1:3) %>%
                dplyr::mutate(msprobe = paste0(condition, "_", replicate))

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

model_def <- msglm::msglm_model(conditionXeffect, conditions_df, effects_df)

mscalib <- msglm::read_mscalib_json(test_path("test_mscalib.json"))

for (msprobe in c("msprobe", "msexperiment", "msrun", "mschannel")) {
test_that(paste0("protgroup/protgroup model, no msfractions, no mstags, specifying ", msprobe, "s"), {
    orig_msdata <- gen_msdata(model_def, msprobes_df, msprobe = msprobe,
                              modelobject = "protgroup", nmodelobjects = 3)
    msprobes_dfname <- paste0(msprobe, "s")
    bad_msdata <- orig_msdata
    bad_msdata[[msprobes_dfname]] <- NULL
    expect_error(msglm::import_msglm_data(bad_msdata, model_def),
                 "Cannot autodetect MS probes")
    if (!(msprobe %in% c("msrun", "mschannel"))) {
        expect_error(msglm::import_msglm_data(orig_msdata, model_def),
                    "Cannot autodetect MS channels")
        # add raw_file column to facilitate autodetection
        orig_msdata[[msprobes_dfname]] <- dplyr::mutate(orig_msdata[[msprobes_dfname]],
                                                        raw_file = paste0(!!sym(msprobe), ".raw"))
    }
    msdata <- msglm::import_msglm_data(orig_msdata, model_def, mscalib)
    testthat::expect_s3_class(msdata, "msglm_data_collection")
    testthat::expect_equal(msdata$msentities, c(quantobject="protgroup", modelobject="protgroup",
        condition="condition", msexperiment=msprobe, msprobe=msprobe,
        msrun=msprobe, mschannel=msprobe, mstag=NA_character_, msfraction=NA_character_))
})

test_that(paste0("protgroup/pepmodstate model, no msfractions, no mstags, specifying ", msprobe, "s"), {
    model_def <- msglm::msglm_model(conditionXeffect, conditions_df, effects_df)
    mschannels_df <- dplyr::mutate(msprobes_df, mschannel = msprobe)
    orig_msdata <- gen_msdata(model_def, mschannels_df, msprobe = msprobe, mschannel = msprobe,
                              modelobject = "protgroup", quantobject = "pepmodstate", nmodelobjects=3)
    msprobes_dfname <- paste0(msprobe, "s")
    orig_msdata[[msprobes_dfname]] <- dplyr::rename(msprobes_df, !!sym(msprobe) := msprobe) %>%
        dplyr::mutate(raw_file = paste0(!!sym(msprobe), ".raw"))
    msdata <- msglm::import_msglm_data(orig_msdata, model_def, quantobject = "pepmodstate", mscalib)
    testthat::expect_s3_class(msdata, "msglm_data_collection")
    testthat::expect_equal(msdata$msentities, c(quantobject="pepmodstate", modelobject="protgroup",
        condition="condition", msexperiment=msprobe, msprobe=msprobe,
        msrun=msprobe, mschannel=msprobe, mstag=NA_character_, msfraction=NA_character_))
})

}

for (msprobe in c("msprobe", "msexperiment")) {
test_that(paste0("protgroup/pepmodstate model, msfractions, no mstags, specifying ", msprobe, "s"), {
    mschannel <- c(msprobe = "mschannel", msexperiment = "msrun")[[msprobe]]
    model_def <- msglm::msglm_model(conditionXeffect, conditions_df, effects_df)
    mschannels_df <- dplyr::full_join(msprobes_df, tibble(msfraction = c(1L, 2L)), by=character()) %>%
                     dplyr::mutate(mschannel = paste0(msprobe, "_F", msfraction))
    orig_msdata <- gen_msdata(model_def, mschannels_df, msprobe = msprobe, mschannel = mschannel,
                              modelobject = "protgroup", quantobject = "pepmodstate", nmodelobjects=3)
    msprobes_dfname <- paste0(msprobe, "s")
    mschannels_dfname <- paste0(mschannel, "s")
    # without specifying mschannels, by with hint that msfractions exist
    bad_msdata <- orig_msdata
    bad_msdata[[paste0(mschannel, "s")]] <- NULL
    expect_error(msglm::import_msglm_data(bad_msdata, model_def, msfraction="msfraction"),
                "Cannot autodetect MS probes")

    # without specifying msprobes
    orig_msdata1 <- orig_msdata
    msdata <- msglm::import_msglm_data(orig_msdata1, model_def, quantobject = "pepmodstate", mscalib)
    testthat::expect_s3_class(msdata, "msglm_data_collection")
    checkmate::test_names(names(msdata), must.include = c("msentities", msprobes_dfname, mschannels_dfname,
                          "pepmodstate_intensities", "protgroups", "pepmodstates", "protgroup2pepmodstate"))
    testthat::expect_equal(msdata$msentities, c(quantobject="pepmodstate", modelobject="protgroup",
        condition="condition", msexperiment=msprobe, msprobe=msprobe,
        msrun=mschannel, mschannel=mschannel, mstag=NA_character_, msfraction="msfraction"))
    checkmate::test_tibble(msdata[[msprobes_dfname]], nrows=nrow(msprobes_df))
    checkmate::test_tibble(msdata[[mschannels_dfname]], nrows=nrow(mschannels_df))

    # without specifying msprobes and without msprobe column in mschannels
    orig_msdata1b <- orig_msdata
    orig_msdata1b[[mschannels_dfname]] <- dplyr::mutate(orig_msdata1b[[mschannels_dfname]], !!sym(msprobe) := NULL)
    msdata <- msglm::import_msglm_data(orig_msdata1b, model_def, quantobject = "pepmodstate", mscalib)
    testthat::expect_s3_class(msdata, "msglm_data_collection")
    checkmate::test_names(names(msdata), must.include = c("msentities", msprobes_dfname, mschannels_dfname,
                          "pepmodstate_intensities", "protgroups", "pepmodstates", "protgroup2pepmodstate"))
    testthat::expect_equal(msdata$msentities, c(quantobject="pepmodstate", modelobject="protgroup",
        condition="condition", msexperiment=msprobe, msprobe=msprobe,
        msrun=mschannel, mschannel=mschannel, mstag=NA_character_, msfraction="msfraction"))
    checkmate::test_tibble(msdata[[msprobes_dfname]], nrows=nrow(msprobes_df))
    checkmate::test_tibble(msdata[[mschannels_dfname]], nrows=nrow(mschannels_df))

    # without specifying mschannels
    orig_msdata2 <- orig_msdata
    orig_msdata2[[msprobes_dfname]] <- dplyr::rename(msprobes_df, !!sym(msprobe) := msprobe)
    orig_msdata2[[paste0(mschannel, "s")]] <- NULL
    testthat::expect_error(msglm::import_msglm_data(orig_msdata2, model_def, quantobject = "pepmodstate", mscalib),
                           "Cannot autodetect MS channels data")

    # with both
    orig_msdata3 <- orig_msdata
    orig_msdata3[[msprobes_dfname]] <- dplyr::rename(msprobes_df, !!sym(msprobe) := msprobe)
    msdata <- msglm::import_msglm_data(orig_msdata3, model_def, quantobject = "pepmodstate", mscalib)
    testthat::expect_s3_class(msdata, "msglm_data_collection")
    testthat::expect_equal(msdata$msentities, c(quantobject="pepmodstate", modelobject="protgroup",
        condition="condition", msexperiment=msprobe, msprobe=msprobe,
        msrun=mschannel, mschannel=mschannel, mstag=NA_character_, msfraction="msfraction"))
})

}
