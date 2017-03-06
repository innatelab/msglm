job.args <- c("dhaas_vacv", "vacv_apms_glm", "20160229", "0" )
if ( !exists( 'job.args' ) ) {
  job.args <- commandArgs( trailingOnly = TRUE )
}

project_id <- job.args[[1]]
message( 'Project ID=', project_id )

job_name <- as.character( job.args[[2]] )
message( 'Job name is ', job_name )

job_version <- job.args[[3]]
message( 'Job dataset version is ', job_version )

job_id <- as.integer( job.args[[4]] )
message( 'Job ID is ', job_id )

source("~/R/config.R")
source( file.path( base_scripts_path, 'pcp/setup_base_paths.R' ) )
source( file.path( base_scripts_path, 'pcp/reshape.R' ) )
source( file.path( base_scripts_path, 'pcp/stan_process_utils.R') )
source( file.path( base_scripts_path, 'misc/plot_themes.R' ) )
source( file.path( base_scripts_path, 'pcp/setup_project_paths.R' ) )
source( file.path( project_scripts_path, 'apms_glm_fit_utils.R' ) )

mop.max_nprocesses <- 8
mop.nprocesses <- 8
source( file.path( pipeline_scripts_path, 'init_cluster.R' ) )

require(dplyr)
require(ggplot2)
require(Biostrings)
require(RGraphML)
require(rstan)

read_innate_fasta <- function(fasta_filename) {
  ens.fasta <- readAAStringSet(fasta_filename)
  ens.matches <- regexec( "sp\\|([A-Za-z0-9_]+)\\s(.+)\\sOS=([^=]+)\\sGN=([^=]+)(?:\\sPE=([^=]+)\\sSV=([^=]+))?",
                          names(ens.fasta) )
  pre_df <- regmatches(names(ens.fasta), ens.matches)
  max_cols <- max(sapply(pre_df, length))
  pre_df <- lapply(pre_df, function(matches) if (length(matches)<max_cols) c(matches, rep.int(NA, max_cols - length(matches))) else matches )
  res.df <- do.call(rbind, pre_df) %>% as.data.frame(stringsAsFactors=FALSE)
  colnames(res.df) <- c("src", "protein_id", "description", "organism", "gene_name", "PE", "SV")[1:ncol(res.df)]
  res.df$seqlen <- nchar(ens.fasta)
  res.df$seq <- as.character(ens.fasta)
  res.df
}

require(rjson)
require(rstan)
apms_glm.stan_model <- stan_model(file.path(project_scripts_path, "apms_glm.stan"), "apms_glm")

data_path <- file.path(base_analysis_path, project_id)

instr_calib_json_filename <- "instr_orbi2_calib_laudenbach_orbi2calib_20161005_borg.json"
message( 'Loading MS instrument calibration data from ', instr_calib_json_filename, '...' )
instr_calib <- fromJSON(file = file.path(data_path, instr_calib_json_filename))$instr_calib

vacv_baits.df <- read_innate_fasta(file.path(data_path, "VACWR_selected.fasta"))
vacv_orf.df <- read.table(file.path(data_path, "vacv_orf.txt"), stringsAsFactors = FALSE, header = TRUE, sep = '\t')

ms_data.wide <- read.table(file.path(data_path, "proteinGroups.txt"), sep='\t',
                           header=TRUE, check.names = FALSE, stringsAsFactors = FALSE) %>%
  dplyr::select(starts_with("Intensity "), contains("protein"), contains("Gene")) %>%
  mutate( `Protein names` = ifelse(`Protein names` == '' & grepl("_VACCW", `Protein IDs`),
                                   gsub(" OS=.+$", "", gsub("^.+_VACCW ", "", `Protein IDs`)),
                                   `Protein names`),
          `Protein IDs` = gsub(">sp\\|", "", gsub("(?:_VACCW).+", "", `Protein IDs` ) ),
          `Majority protein IDs` = gsub(">sp\\|", "", gsub("(?:_VACCW).+", "", `Majority protein IDs` ) ),
          `Gene names` = ifelse(`Gene names` != "", `Gene names`, `Protein IDs` )) %>%
  mutate_each(funs(ifelse(.==0.0, NA, .))) %>%
  dplyr::rename(majority_protein_ids = `Majority protein IDs`)

ms_data.mtx <- as.matrix(ms_data.wide %>% dplyr::select(starts_with("Intensity ")))
msruns.df <- data.frame(msrun_id = 1:ncol(ms_data.mtx),
                        msrun = gsub("Intensity ", "", colnames(ms_data.mtx)[grepl("Intensity", colnames(ms_data.mtx))]),
                        stringsAsFactors=FALSE) %>%
  dplyr::mutate(bait_orf = gsub( "-4", "", gsub("_.+$", "", msrun) ),
                replicate = as.integer(gsub("^.+_", "", msrun)))
baits.df <- msruns.df %>% dplyr::select(bait_orf) %>% dplyr::distinct() %>%
  left_join(vacv_orf.df) %>%
  dplyr::arrange(desc(bait_orf == 'ctrl'), bait_orf) %>%
  dplyr::mutate(bait_orf = factor(bait_orf, levels=bait_orf),
                bait_ix = as.integer(bait_orf))
msruns.df <- dplyr::rename(msruns.df, bait_orf_str = bait_orf) %>%
  dplyr::inner_join(dplyr::mutate(baits.df, bait_orf_str = as.character(bait_orf)) %>%
                    dplyr::select(bait_orf_str, bait_orf, bait_ix)) %>%
  dplyr::select(-bait_orf_str) %>%
  dplyr::arrange(bait_ix, msrun_id) %>%
  dplyr::mutate(msrun_ix = as.integer(factor(msrun_id)))

# setup experimental design

all_conditions <- as.character(baits.df$bait_orf)
control_conditions <- c('ctrl')
noncontrol_conditions <- setdiff(all_conditions, control_conditions)
all_effects <- setdiff(all_conditions, control_conditions)# dplyr::filter(baits.df, bait_orf != 'ctrl') %>% .$bait_orf

effectXcondition <- zero_matrix(effect = all_effects, condition = all_conditions)
for (i in 1:nrow(effectXcondition)) {
  effectXcondition[all_effects[[i]], all_effects[[i]]] <- 1.0
  #effectXcondition['background', all_conditions[[i]]] <- 1.0
}
effectXcondition.df <- as.data.frame(as.table(effectXcondition)) %>%
  dplyr::filter(Freq != 0) %>% dplyr::select(-Freq)

compound_metaconditions <- paste0("all_minus_", noncontrol_conditions)
all_metaconditions <- c(all_conditions, compound_metaconditions)
conditionXmetacondition <- false_matrix(condition = all_conditions,
                                        metacondition = all_metaconditions)
for (cname in all_conditions) {
  conditionXmetacondition[cname, cname] <- TRUE
}
for (cname in noncontrol_conditions) {
  conditionXmetacondition[, paste0("all_minus_", cname)] <- rep.int(TRUE, length(all_conditions))
  conditionXmetacondition[cname, paste0("all_minus_", cname)] <- FALSE
}
vs_ctrl_contrasts <- paste0(noncontrol_conditions, "_vs_ctrl")
vs_others_contrasts <- paste0(noncontrol_conditions, "_vs_others")
all_contrasts <- c(vs_ctrl_contrasts, vs_others_contrasts)
contrastXmetacondition <- zero_matrix(contrast = all_contrasts, metacondition = all_metaconditions)
for (cname in noncontrol_conditions) {
  contrastXmetacondition[paste0(cname, "_vs_ctrl"), cname] <- 1.0
  contrastXmetacondition[paste0(cname, "_vs_ctrl"), control_conditions] <- -1.0
  contrastXmetacondition[paste0(cname, "_vs_others"), cname] <- 1.0
  contrastXmetacondition[paste0(cname, "_vs_others"), paste0("all_minus_", cname)] <- -1.0
}
contrastXmetacondition.df <- as.data.frame(as.table(contrastXmetacondition)) %>% dplyr::filter(Freq != 0) %>%
  dplyr::rename(weight = Freq) %>%
  dplyr::mutate(contrast_type = 'filtering',
                condition_role = if_else(contrast_type == 'filtering',
                                         if_else(weight > 0, 'signal', 'background' ),
                                         'signal') )
conditionXmetacondition.df <- as.data.frame(as.table(conditionXmetacondition)) %>%
  dplyr::filter(Freq) %>% dplyr::select(- Freq)

contrastXcondition.df <- as.data.frame(as.table(conditionXmetacondition)) %>% dplyr::filter(Freq) %>%
  dplyr::select(-Freq) %>%
  dplyr::inner_join(contrastXmetacondition.df) %>%
  dplyr::mutate(bait_orf = condition) %>%
  dplyr::inner_join(baits.df %>% dplyr::select(bait_orf, bait_ix, gene_name)) %>%
  dplyr::arrange(contrast, metacondition, condition)

batch_effects.df <- dplyr::select(msruns.df, bait_ix, bait_orf) %>%
  dplyr::group_by(bait_orf, bait_ix) %>%
  dplyr::mutate(n_replicates = n()) %>%
  dplyr::filter(row_number() > 1L) %>%
  dplyr::mutate(in_condition_ix = row_number()) %>% dplyr::ungroup() %>%
  dplyr::arrange(bait_ix) %>%
  dplyr::mutate(batch_effect = paste0(bait_orf, "_", in_condition_ix)) %>%
  dplyr::mutate(batch_effect = factor(batch_effect, levels=batch_effect)) %>% dplyr::ungroup()

batch_effects_coeffs.df <- inner_join(msruns.df, batch_effects.df) %>%
  mutate(coef = (ifelse(in_condition_ix == replicate, n_replicates, 0.0) - 1.0)/n_replicates)

all_batch_effects <- batch_effects.df$batch_effect
batchEffectXexperiment <- zero_matrix(batch_effect = all_batch_effects, experiment = msruns.df$msrun)
for (i in 1L:nrow(batch_effects_coeffs.df)) {
  batchEffectXexperiment[batch_effects_coeffs.df$batch_effect[i], batch_effects_coeffs.df$msrun_id[i]] <-
    batch_effects_coeffs.df$coef[i]
}

sel_proteins <- ms_data.wide %>% dplyr::filter(`Gene names` == 'SPI2') %>% .$majority_protein_ids %>% unique()
sel_proteins <- c(sample(ms_data.wide$majority_protein_ids, 100),
                  ms_data.wide %>% dplyr::filter(`Gene names` %in% c('F17', 'A57', 'A47', 'B2', 'KRB2', 'SPI2', 'ORF50', "VRK1", "ASNS", "MRPS18C")) %>%
                    .$majority_protein_ids) %>% unique()
sel_proteins <- ms_data.wide$majority_protein_ids
message(length(sel_proteins), " protein group(s) selected")

# all MS data including non-observed interactions
ms_data_all.df <- reshape(ms_data.wide, direction="long", idvar=c("majority_protein_ids"),
                      drop = c("Protein IDs", "Protein names", "Number of proteins", "Gene names"),
                      timevar = "msrun", sep=" ", v.names = "Intensity") %>%
  dplyr::rename(intensity = Intensity) %>%
  dplyr::inner_join(msruns.df) %>%
  dplyr::mutate(iaction_id = paste(bait_orf, majority_protein_ids))

def_norm_data <- instr_calib[c('zDetectionFactor', 'zDetectionIntercept',
                               'detectionMax', 'sigmaScaleHi', 'sigmaScaleLo',
                               'sigmaOffset', 'sigmaBend', 'sigmaSmooth',
                               'zShift', 'zScale')]

exp_shifts.df <- normalize_experiments(apms_normalize.stan_model, def_norm_data,
                                       ms_data_all.df, # %>% dplyr::semi_join(model_data$proteins %>% dplyr::filter(`Gene names` %in% c('RPS13', 'RPS21'))),
                                       exp_col = "msrun_ix", cond_col = "msrun_ix", cond_group = "bait_ix") %>%
  dplyr::mutate(msrun_ix = as.integer(msrun_ix))
bait_shifts.df <- normalize_experiments(apms_normalize.stan_model, def_norm_data,
                                        ms_data_all.df, # %>% dplyr::semi_join(model_data$proteins %>% dplyr::filter(`Gene names` == 'RPS13')),
                                        exp_col = "msrun_ix", cond_col = "bait_ix", cond_group = NULL,
                                    exp_shifts = exp_shifts.df) %>%
  dplyr::mutate(bait_ix = as.integer(bait_ix)) %>% dplyr::inner_join(baits.df)

total_exp_shifts.df <- dplyr::select(exp_shifts.df, bait_ix, msrun_ix, exp_shift = shift) %>%
  dplyr::inner_join(dplyr::select(bait_shifts.df, bait_ix, bait_shift = shift)) %>%
  dplyr::mutate(total_shift = exp_shift + bait_shift,
                msrun_ix = as.integer(msrun_ix)) %>% dplyr::arrange(msrun_ix)

ms_data.df <- dplyr::filter(ms_data_all.df, is.finite(intensity) & intensity > 0 & majority_protein_ids %in% sel_proteins)

model_data <- list()
model_data$conditions <- contrastXcondition.df %>% dplyr::select(-contrast, -weight) %>%
  dplyr::distinct() %>% dplyr::arrange(bait_ix) %>%
  dplyr::mutate(condition_ix = as.integer(factor(condition)))
model_data$msruns <- dplyr::select(msruns.df, bait_ix, bait_orf, msrun_ix, msrun) %>%
  dplyr::inner_join(model_data$conditions) %>%
  dplyr::select(bait_ix, bait_orf, condition, msrun, msrun_ix) %>% dplyr::distinct() %>%
  dplyr::inner_join(total_exp_shifts.df) %>%
  dplyr::arrange(msrun_ix)

# for interactions use all non-control proteins and all their observations in the contrasts as well
model_data$interactions <- ms_data.df %>% dplyr::filter(majority_protein_ids %in% sel_proteins) %>%
  dplyr::select(majority_protein_ids, bait_ix) %>% dplyr::distinct() %>%
  dplyr::inner_join(model_data$conditions %>% dplyr::select(bait_ix, condition)) %>%
  dplyr::inner_join(dplyr::filter(contrastXcondition.df, condition_role != 'background') %>% dplyr::select(contrast, condition)) %>%
  dplyr::select(-bait_ix, -condition) %>% dplyr::distinct() %>%
  dplyr::inner_join(contrastXcondition.df) %>% dplyr::inner_join(ms_data_all.df) %>%
  dplyr::inner_join(model_data$conditions) %>%
  dplyr::group_by(iaction_id, majority_protein_ids, condition, condition_role, condition_ix, bait_ix, bait_orf) %>%
  dplyr::summarize(is_virtual = all(is.na(intensity))) %>% dplyr::ungroup() %>%
  dplyr::group_by(majority_protein_ids, condition_role) %>%
  dplyr::filter(!is_virtual | (condition_role != 'background') | (condition_role == 'background' & condition %in% control_conditions)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-condition_role) %>% dplyr::distinct() %>%
  dplyr::arrange(bait_ix, majority_protein_ids) %>%
  dplyr::mutate(iaction_ix = as.integer(factor(iaction_id, levels=iaction_id)),
                protgroup_ix = as.integer(factor(majority_protein_ids)))

model_data$proteins <- dplyr::select(model_data$interactions, protgroup_ix, majority_protein_ids) %>%
  dplyr::inner_join( dplyr::select(ms_data.wide, majority_protein_ids, `Protein names`, `Protein IDs`, `Gene names`) %>% dplyr::distinct() ) %>%
  dplyr::distinct() %>%
  dplyr::arrange(protgroup_ix)

# entries for an interaction in all replicate experiments
model_data$ms_data <- dplyr::inner_join(model_data$interactions, msruns.df) %>% dplyr::left_join(ms_data.df) %>%
  dplyr::arrange(msrun_ix, protgroup_ix) %>%
  dplyr::mutate(observation_ix = seq_len(n()),
                qdata_ix = if_else(!is.na(intensity), cumsum(!is.na(intensity)), NA_integer_),
                mdata_ix = if_else(is.na(intensity), cumsum(is.na(intensity)), NA_integer_))

model_data$protein_effects <- obj2effect(effectXcondition, "protgroup_ix", model_data$proteins$protgroup_ix,
                                         model_data$interactions$protgroup_ix,
                                         model_data$interactions$condition_ix) %>%
  dplyr::mutate(protgroup_ix = as.integer(protgroup_ix))

model_data$protein_batch_effects <- obj2effect(batchEffectXexperiment, "protgroup_ix",
                                               model_data$proteins$protgroup_ix,
                                               model_data$ms_data$protgroup_ix,
                                               model_data$ms_data$msrun_ix)

# detect proteins that have no quantifications for estimating protein_shift
pe2iact.df <- model_data$protein_effects %>%
  dplyr::inner_join(effectXcondition.df) %>%
  dplyr::inner_join(model_data$interactions)
effectless_iactions.df <- model_data$interactions %>% dplyr::anti_join(pe2iact.df)
underdefined_proteins.df <- dplyr::anti_join(model_data$proteins, effectless_iactions.df %>% dplyr::filter(!is_virtual)) %>% dplyr::select(protgroup_ix)

model_data$proteins <- dplyr::mutate(model_data$proteins,
                                     is_underdefined = protgroup_ix %in% underdefined_proteins.df$protgroup_ix)

apms_glm.stan_data <- list(
  Niactions = nrow(model_data$interactions),
  Nobservations = nrow(model_data$ms_data),
  Nexperiments = n_distinct(model_data$msruns$msrun_ix),
  Nconditions = length(all_conditions),
  Nproteins = n_distinct(model_data$ms_data$protgroup_ix),
  experiment2condition = model_data$msruns$bait_ix,
  effectXcondition = effectXcondition,
  experiment_shift = model_data$msruns$total_shift,
  batchEffectXexperiment = batchEffectXexperiment,
  observation2experiment = model_data$ms_data$msrun_ix,
  observation2iaction = model_data$ms_data$iaction_ix,
  iaction2protein = model_data$interactions$protgroup_ix,
  iaction2condition = model_data$interactions$bait_ix,
  Neffects = length(all_effects),
  NproteinEffects = nrow(model_data$protein_effects),
  NeffectsPerProteinCumsum = effects_cumsum(model_data$protein_effects, 'protgroup_ix'),
  protein_effect2protein = model_data$protein_effects$protgroup_ix,
  protein_effect2effect = as.integer(model_data$protein_effects$effect),
  NbatchEffects = length(all_batch_effects),
  NproteinBatchEffects = nrow(model_data$protein_batch_effects),
  NbatchEffectsPerProteinCumsum = effects_cumsum(model_data$protein_batch_effects, 'protgroup_ix'),
  protein_batch_effect2batch_effect = as.integer(model_data$protein_batch_effects$batch_effect),
  NunderdefProteins = sum(model_data$proteins$is_underdefined),
  underdef_proteins = dplyr::filter(model_data$proteins, is_underdefined) %>% .$protgroup_ix,
  Nquanted = sum(!is.na(model_data$ms_data$intensity)),
  Nmissed = sum(is.na(model_data$ms_data$intensity)),
  quant2observation = model_data$ms_data$observation_ix[!is.na(model_data$ms_data$qdata_ix)],
  miss2observation = model_data$ms_data$observation_ix[!is.na(model_data$ms_data$mdata_ix)],
  qData = model_data$ms_data$intensity[!is.na(model_data$ms_data$intensity)]
) %>%
  c(instr_calib[c('zDetectionFactor', 'zDetectionIntercept',
                  'detectionMax', 'sigmaScaleHi', 'sigmaScaleLo',
                  'sigmaOffset', 'sigmaBend', 'sigmaSmooth',
                  'zShift', 'zScale')])

apms_glm.stan_data$protein_base <- 0.95*median(log(apms_glm.stan_data$qData))

apms_glm.stan_model <- stan_model(file.path(project_scripts_path, "apms_glm.stan"), "apms_glm", save_dso = TRUE)
message( 'Running STAN in HMC mode...' )
if ( !is.null( mop.cluster ) ) {
clusterEvalQ( mop.cluster, library(rstan) )
clusterEvalQ( mop.cluster, library(dplyr) )
stan_sampling_seed = sample.int(.Machine$integer.max, 1)
clusterExport(mop.cluster, varlist=c('stan_sampling_seed', 'apms_glm.stan_data', 'apms_glm.stan_model'))
apms_glm.stan_mcmc <- clusterApplyLB(mop.cluster, seq_len(mop.nprocesses), function( chain_ix ) {
         sampling(apms_glm.stan_model,
                  data = apms_glm.stan_data,
                  pars = c('condition_shift0', 'experiment_shift0',
                           #'protein_effect_lambda',
                           'protein_effect_unscaled'), include=FALSE,
                  #init = function() { pcp_peaks_glm.generate_init_params(pcp_peaks_glm.model_data) },
                  iter = 5000L, chains = 1L, thin = 4L, chain_id=chain_ix, seed=stan_sampling_seed)
})
apms_glm.stan_mcmc <- sflist2stanfit( apms_glm.stan_mcmc )
apms_glm.stan_fit <- apms_glm.stan_mcmc
} else {
apms_glm.stan_mcmc <- sampling(apms_glm.stan_model,
                               pars = c('condition_shift0', 'experiment_shift0',
                                        #'protein_effect_lambda',
                                        'protein_effect_unscaled'), include=FALSE,
                               data = apms_glm.stan_data,
                               #init = function() { pcp_peaks_glm.generate_init_params(pcp_peaks_glm.model_data) },
                               iter = 2000, chains = 4, thin = 4 )
apms_glm.stan_fit <- apms_glm.stan_mcmc
}

apms_glm.stan_vb <- vb(apms_glm.stan_model,
                       pars = c('condition_shift0', 'experiment_shift0', 'protein_shift0',
                                #'protein_effect_lambda', 'protein_effect_tau',
                                'protein_effect_unscaled'), include=FALSE,
                       data = apms_glm.stan_data,
                       #init = function() { pcp_peaks_glm.generate_init_params(pcp_peaks_glm.model_data) },
                       adapt_iter = 2000, iter = 50000, algorithm = "meanfield", grad_samples = 3)
apms_glm.stan_fit <- apms_glm.stan_vb
options(shiny.launch.browser = FALSE)
launch_shinystan(as.shinystan(apms_glm.stan_fit), port=9080, host="0.0.0.0")

vars_info <- list(
  global = list( names = c('protein_shift_sigma', 'protein_effect_tau', 'underdef_protein_shift'), dims = c() ),
  #effects = list( names = c('protein_effect_tau'), dims = c('effect') ),
  batch_effects = list(names = c('protein_batch_effect_sigma'), dims = c('batch_effect')),
  iactions = list( names = c('iaction_shift'), dims = c('iaction')),
  observations = list( names = c('iaction_repl_shift'), dims = c('observation')),
  proteins = list( names = c('protein_shift'), dims = c('protein') ),
  protein_effects = list( names = c('protein_effect_lambda',
                                    'protein_effect'), dims = c('protein_effect') ),
  protein_batch_effects = list( names = c('protein_batch_effect'), dims = c('protein_batch_effect') )
)

min.iteration <- as.integer(1.25 * apms_glm.stan_fit@sim$warmup)

message( 'Computing parameters statistics...' )
apms_glm.stan_stats <- apms_glm.stan_fit %>%
  stan.extract_samples( pars = unlist( sapply( vars_info, function( vi ) vi$names ) ),
                        min.iteration = min.iteration ) %>%
  monitor( print = FALSE ) %>% as.data.frame
apms_glm.stan_stats$var_name <- rownames( apms_glm.stan_stats )

require( insilicoMop )

message( 'Extracting MCMC samples...' )
apms_glm.stan_samples <- stan.extract_samples( apms_glm.stan_fit,
                                               pars = unlist( sapply( vars_info, function( vi ) vi$names ) ),
                                               min.iteration = min.iteration,
                                               permuted = TRUE )

dim_info <- list( iteration = NULL,
                  effect = data.frame( param = all_effects,
                                       param_type = 'effect', stringsAsFactors = FALSE ),
                  batch_effect = data.frame( param = all_batch_effects,
                                             param_type = 'effect', stringsAsFactors = FALSE ),
                  condition = data.frame( param = all_conditions, index_condition = seq_along(all_conditions),
                                          param_type = 'condition', stringsAsFactors = FALSE ),
                  experiment = model_data$msruns,
                  iaction = dplyr::select(model_data$interactions, iaction_ix, protgroup_ix, condition_ix, bait_ix, condition, bait_orf) %>%
                            dplyr::inner_join(model_data$proteins),
                  observation = dplyr::select(model_data$ms_data, iaction_ix, protgroup_ix, condition_ix, bait_ix, condition, bait_orf, msrun, msrun_ix),
                  protein = model_data$proteins,
                  protein_effect = model_data$protein_effects %>%
                                    dplyr::mutate(protgroup_ix = as.integer(protgroup_ix)) %>%
                                    dplyr::inner_join(model_data$proteins),
                  protein_batch_effect = model_data$protein_batch_effects %>%
                                    dplyr::mutate(protgroup_ix = as.integer(protgroup_ix))
                )

vars_results <- lapply( names(vars_info), vars_statistics, apms_glm.stan_stats, apms_glm.stan_samples, vars_info, dim_info )
names(vars_results) <- names(vars_info)

message( 'Calculating P-values' )
local({
  for ( vars_category in names(vars_results) ) {
    vars_cat_subset_info <- vars_info[[vars_category]]
    vars_cat_subset_info$names <- intersect(vars_cat_subset_info$names, c('protein_effect') )

    if ( length(vars_cat_subset_info$names) > 0 ) {
      message( 'Calculating P-values for ', vars_category, ' variables...' )
      p_value.df <- vars_pvalue_not_zero( vars_results[[vars_category]]$samples, vars_cat_subset_info, dim_info )
      #print(str(p_value.df))
      #print(str(stats.df))
      vars_results[[vars_category]]$stats <<- left_join( vars_results[[vars_category]]$stats %>% mutate( p_value = NULL ), p_value.df )
    }
}})

message( 'Calculating contrasts' )
local({
  for ( vars_category in names(vars_results) ) {
    vars_cat_subset_info <- vars_info[[vars_category]]
    vars_cat_subset_info$names <- intersect(vars_cat_subset_info$names, c('iaction_shift', 'iaction_repl_shift'))

    if ( length(vars_cat_subset_info$names) > 0 ) {
      message( 'Calculating contrasts for ', vars_category, ' variables...' )
      samples.df <- dplyr::inner_join(vars_results[[vars_category]]$samples,
                                      conditionXmetacondition.df)
      #print(vars_cat_subset_info)
      #print(str(samples.df))
      obj_dims <- 'protgroup_ix'
      contrast_stats.df <- vars_contrast_stats( samples.df, vars_cat_subset_info$names,
                                                'protgroup_ix', 'metacondition',
                                                if ('msrun' %in% colnames(samples.df)) 'msrun'
                                                else if ('condition' %in% colnames(samples.df)) 'condition'
                                                else 'condition',
                                                contrastXmetacondition, transform = 'exp' )
      var_info_cols <- 'protgroup_ix'
      var_info.df <- vars_results[[vars_category]]$stats %>%
        dplyr::select_(.dots = var_info_cols ) %>% dplyr::distinct()
      # inject contrast statistics into vars_results
      vars_results[[vars_category]]$contrast_stats <<- left_join( var_info.df,
                                                               contrast_stats.df ) %>%
        dplyr::mutate(prob_nonpos_fake = pmax(if_else(mean > 0, 1E-300, 1E-100), rgamma(n(), shape=1E-2, scale=1E-150), prob_nonpos),
                      prob_nonneg_fake = pmax(if_else(mean > 0, 1E-300, 1E-100), rgamma(n(), shape=1E-2, scale=1E-150), prob_nonneg)) %>%
      left_join(contrastXcondition.df %>% dplyr::filter(weight > 0)) %>%
        dplyr::inner_join(dim_info$iaction)
    }
}})

ggplot(vars_results$iactions$contrast_stats %>%
       dplyr::filter(bw > 0 & !is.na(iaction_ix)),
       aes(x=mean, y=2*pmin(prob_nonpos_fake, prob_nonneg_fake))) +
  geom_point(alpha=0.1) +
  geom_text(data=dplyr::filter(vars_results$iactions$contrast_stats,
                               mean != 0 & bw > 0 & !is.na(iaction_ix) & pmin(prob_nonpos, prob_nonneg)<1E-0) %>%
              dplyr::left_join(model_data$proteins),
            aes(label=if_else(is.na(`Gene names`), `Protein IDs`, `Gene names`)),
            cex=3, vjust=-1.1) +
  scale_y_log10() + facet_wrap( ~ contrast, scales = "free", ncol = 4)

ggplot(vars_results$protein_batch_effects$stats) +
  geom_density(aes(x=mean)) +
  geom_rug(data=dplyr::filter(model_data$proteins, `Gene names` == 'KRB2') %>%
             dplyr::inner_join(vars_results$protein_batch_effects$stats),
           aes(x=mean, color=gsub("_\\d+$", "", batch_effect))) +
  facet_wrap(~ var)# + scale_x_continuous(limits=c(-10,10))

ggplot(vars_results$protein_effects$stats) +
  geom_density(aes(x=mean)) +
  geom_rug(data=dplyr::filter(model_data$proteins, `Gene names` == 'B2') %>%
             dplyr::inner_join(vars_results$protein_effects$stats),
           aes(x=mean, color=effect)) +
  facet_wrap(~ var)# + scale_x_continuous(limits=c(-10,10))

ggplot(vars_results$proteins$stats) +
  geom_density(aes(x=mean)) +
  geom_rug(data=dplyr::filter(model_data$proteins, `Gene names` == 'B2') %>%
             dplyr::inner_join(vars_results$proteins$stats)) +
  facet_wrap(~ var, scales = "free")# + scale_x_continuous(limits=c(-10,10))

ggplot(vars_results$effects$stats) +
  geom_density(aes(x=mean)) +
  facet_wrap(~ var)# + scale_x_continuous(limits=c(-10,10))

ggplot(vars_results$experiments$samples %>% dplyr::inner_join(dim_info$condition %>% dplyr::rename(condition=param))) +
  geom_density(aes(x=experiment_shift, fill=factor(index_experiment)), alpha=0.3) +
  facet_wrap(~ condition, scales="free")# + scale_x_continuous(limits=c(-10,10))

ggplot(vars_results$observations$contrast_stats,
       aes(x=mean, y=pmax(1E-300,2*pmin(prob_nonpos, prob_nonneg)))) +
  geom_point() +
  geom_text(data=dplyr::filter(vars_results$observations$contrast_stats,
                               prob_nonpos<1E-10) %>%
              dplyr::left_join(model_data$proteins),
            aes(label=if_else(is.na(`Gene names`), `Protein IDs`, `Gene names`)),
            cex=3, vjust=-1.1) +
  scale_y_log10() + facet_wrap( ~ contrast, ncol = 4)

if (!dir.exists(file.path(scratch_path, "apms_glm"))) {
  dir.create(file.path(scratch_path, "apms_glm"))
}
rdata_filepath <- file.path(scratch_path, "apms_glm", paste0( 'apms_glm_', job_version, '.RData' ) )
message( 'Saving STAN results to ', rdata_filepath, '...' )
save( job_version, # apms_glm.stan_opt,
      apms_glm.stan_data, apms_glm.stan_fit,
      model_data, vars_results, file = rdata_filepath )
message( 'Done.' )
