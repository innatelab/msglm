functions {
    real intensity_log_std(real z, real scaleHi, real scaleLo, real offset, real bend, real smooth) {
        return 0.5*(scaleHi+scaleLo)*(z-bend) + 0.5*(scaleHi-scaleLo)*sqrt((z-bend)*(z-bend)+smooth) + offset;
    }
}

data {
  int<lower=1> Nexperiments;    # number of experiments
  int<lower=1> Nconditions;     # number of experimental conditions
  int<lower=0> Nproteins;       # number of proteins
  int<lower=0> Niactions;       # number of interactions
  int<lower=0> Nobservations;   # number of potential observations of interactions
  int<lower=0> Neffects;        # number of effects (that define conditions)
  int<lower=0> NbatchEffects;   # number of batch effects (that define experimental variation)
  int<lower=0> NunderdefProteins; # number of virtual interactions (the ones not detected but required for comparison)
  int<lower=1,upper=Nproteins> iaction2protein[Niactions];
  int<lower=1,upper=Nconditions> iaction2condition[Niactions];
  int<lower=1,upper=Nproteins> underdef_proteins[NunderdefProteins];

  matrix[Neffects, Nconditions] effectXcondition;
  matrix[NbatchEffects, Nexperiments] batchEffectXexperiment;
  int<lower=1,upper=Nconditions> experiment2condition[Nexperiments];
  vector[Nexperiments] experiment_shift;

  int<lower=1,upper=Nexperiments> observation2experiment[Nobservations];
  int<lower=1,upper=Niactions> observation2iaction[Nobservations];

  # map from labelXreplicateXprotein to observed/missed data
  int<lower=0> Nquanted;       # total number of quantified proteinsXexperiments
  int<lower=1,upper=Nobservations>  quant2observation[Nquanted];
  int<lower=0> Nmissed;         # total number of missed proteinsXexperiments
  int<lower=1,upper=Nobservations> miss2observation[Nmissed];

  int<lower=0> NproteinEffects;
  int<lower=1,upper=Neffects> protein_effect2effect[NproteinEffects];
  int<lower=1,upper=Nproteins> protein_effect2protein[NproteinEffects];
  int<lower=0> NeffectsPerProteinCumsum[Nproteins+1];

  int<lower=0> NproteinBatchEffects;
  int<lower=1,upper=NbatchEffects> protein_batch_effect2batch_effect[NproteinBatchEffects];
  int<lower=0> NbatchEffectsPerProteinCumsum[Nproteins+1];

  vector<lower=0>[Nquanted] qData; # quanted data

  real protein_base;

  # instrument calibrated parameters 
  real<lower=0> zDetectionFactor;
  real zDetectionIntercept;
  real<lower=0, upper=1> detectionMax;

  real<lower=0> sigmaScaleHi;
  real<lower=0> sigmaScaleLo;
  real sigmaOffset;
  real sigmaBend;
  real sigmaSmooth;

  real zShift;
  real zScale;
}

transformed data {
  int<lower=0> NreplEffects; # number of replicate effects
  real mzShift; # zShift for the missing observation intensity (zShift shifted by protein_base)
  vector[Nquanted] zScore;
  vector[Nquanted] qLogStd;
  vector<lower=0>[Nquanted] qDataNorm;
  real<lower=0> min_intensity;
  real log_min_intensity;
  matrix[Nexperiments, Nconditions] experimentXcondition;

  int<lower=1> NconditionExperiments[Nconditions];
  matrix[Nexperiments, Nexperiments-Nconditions] experimentXrepl_effect;

  int<lower=1,upper=Niactions> quant2iaction[Nquanted];
  int<lower=1,upper=Nexperiments> quant2experiment[Nquanted];
  int<lower=1,upper=Niactions> miss2iaction[Nmissed];
  int<lower=1,upper=Nexperiments> miss2experiment[Nmissed];

  NreplEffects = Nexperiments-Nconditions;
  zScore = (log(qData) - zShift) * zScale;
  mzShift = zShift - protein_base;

  # process the intensity data to optimize likelihood calculation
  for (i in 1:Nquanted) {
    qLogStd[i] = intensity_log_std(zScore[i], sigmaScaleHi, sigmaScaleLo, sigmaOffset, sigmaBend, sigmaSmooth);
    qDataNorm[i] = qData[i] * exp(-qLogStd[i]);
    qLogStd[i] = qLogStd[i] - protein_base;
  }
  quant2experiment = observation2experiment[quant2observation];
  quant2iaction = observation2iaction[quant2observation];
  miss2experiment = observation2experiment[miss2observation];
  miss2iaction = observation2iaction[miss2observation];

  experimentXcondition = rep_matrix(0.0, Nexperiments, Nconditions);
  for (i in 1:Nexperiments) {
    experimentXcondition[i, experiment2condition[i]] = 1.0;
  }

  {
    int experiment_index[Nexperiments];
    int condition_shift0_offset[Nconditions+1];

    NconditionExperiments = rep_array(0, Nconditions);
    for (i in 1:Nexperiments) {
      int cond_ix;
      cond_ix = experiment2condition[i];
      NconditionExperiments[cond_ix] = NconditionExperiments[cond_ix] + 1;
      experiment_index[i] = NconditionExperiments[cond_ix];
    }

    condition_shift0_offset[1] = 0;
    for (i in 1:Nconditions) {
      condition_shift0_offset[i+1] = condition_shift0_offset[i] + NconditionExperiments[i]-1;
    }

    experimentXrepl_effect = rep_matrix(0.0, Nexperiments, Nexperiments-Nconditions);
    for (i in 1:Nexperiments-Nconditions) {
      for (j in 1:Nexperiments) {
        int cond_ix;
        cond_ix = experiment2condition[j];
        if (condition_shift0_offset[cond_ix]<i && i<=condition_shift0_offset[cond_ix+1]) {
          int param_ix;
          int n_cond_exps;
          param_ix = i - condition_shift0_offset[cond_ix];
          n_cond_exps = NconditionExperiments[cond_ix];
          if (experiment_index[j] == param_ix) {
            experimentXrepl_effect[j, i] = (n_cond_exps-1.0)/n_cond_exps;
          } else {
            experimentXrepl_effect[j, i] = -1.0/n_cond_exps;
          }
        } else {
          experimentXrepl_effect[j, i] = 0.0;
        }
      }
    }
  }

  log_min_intensity = -100;
  min_intensity = exp(log_min_intensity);

  print("experimentXrepl_effect=", experimentXrepl_effect);
  print("min(qDataNorm)=", min(qDataNorm), " max(qDataNorm)=", max(qDataNorm) );
  print("min(qLogStd)=", min(qLogStd), " max(qLogStd)=", max(qLogStd) );

  {
    matrix[Nproteins, Neffects] proteinXeffect;
    matrix[Nproteins, Nconditions] proteinXcondition;
    matrix[Nproteins, NbatchEffects] proteinXbatchEffect;
    matrix[Nproteins, Nexperiments] proteinXexperiment;

    proteinXeffect = csr_to_dense_matrix(Nproteins, Neffects, cumulative_sum(rep_vector(1.0, NproteinEffects)),
                                         protein_effect2effect, NeffectsPerProteinCumsum);
    proteinXcondition = proteinXeffect * effectXcondition;
    proteinXbatchEffect = csr_to_dense_matrix(Nproteins, NbatchEffects, cumulative_sum(rep_vector(1.0, NproteinBatchEffects)),
                                              protein_batch_effect2batch_effect, NbatchEffectsPerProteinCumsum);
    proteinXexperiment = proteinXbatchEffect * batchEffectXexperiment;
    print("sample proteinXeffect=", proteinXeffect);
    print("sample proteinXcondition=", proteinXcondition);
    print("sample proteinXbatchEffect=", proteinXbatchEffect);
    print("sample proteinXexperiment=", proteinXexperiment);
  }
}

parameters {
  #real protein_base;
  real<lower=-20, upper=-2> underdef_protein_shift;
  real<lower=0, upper=5.0> protein_shift_sigma;
  #vector<lower=0>[Nconditions] condition_repl_effect_sigma;

  row_vector<lower=-10,upper=10>[Nproteins] protein_shift0;

  real<lower=0.01,upper=5.0> protein_effect_tau;
  vector<lower=0.01,upper=5.0>[NproteinEffects] protein_effect_lambda;
  vector<lower=-10.0,upper=10.0>[NproteinEffects] protein_effect_unscaled;

  #vector<lower=0, upper=3.0>[NbatchEffects] protein_batch_effect_sigma;
  #vector<lower=-5.0, upper=5.0>[NproteinBatchEffects] protein_batch_effect;
}

transformed parameters {
  row_vector[Nproteins] protein_shift;
  vector[NproteinEffects] protein_effect;
  vector[Niactions] iaction_shift;
  #vector[Nobservations] iaction_repl_shift;

  protein_shift = protein_shift0;
  for (i in 1:NunderdefProteins) {
    protein_shift[underdef_proteins[i]] = protein_shift[underdef_proteins[i]] + underdef_protein_shift;
  }
  protein_effect = protein_effect_unscaled .* protein_effect_lambda * protein_effect_tau;
  # calculate iaction_shift
  {
    matrix[Nproteins, Nconditions] proteinXcondition;

    proteinXcondition = csr_to_dense_matrix(Nproteins, Neffects, protein_effect,
                                            protein_effect2effect, NeffectsPerProteinCumsum) * effectXcondition;
    for (i in 1:Niactions) {
      iaction_shift[i] = protein_shift[iaction2protein[i]] + proteinXcondition[iaction2protein[i], iaction2condition[i]];
    }
  }
  # calculate log_observation
  // {
  //   matrix[Nproteins, Nexperiments] proteinXexperiment;
  // 
  //   proteinXexperiment = csr_to_dense_matrix(Nproteins, NbatchEffects, protein_batch_effect,
  //                                            protein_batch_effect2batch_effect, NbatchEffectsPerProteinCumsum) * batchEffectXexperiment;
  //   for (i in 1:Nobservations) {
  //     iaction_repl_shift[i] = iaction_shift[observation2iaction[i]]
  //                        + proteinXexperiment[iaction2protein[observation2iaction[i]], observation2experiment[i]];
  //   }
  // }

  #print("protein_base=", protein_base, " base_sigma=", base_sigma,
  #      " condition_sigma=", condition_sigma, " experiment_sigma=", experiment_sigma);
  #print("min(effect_tau)=", min(effect_tau), " max(effect_tau)=", max(effect_tau));
  #print("min(repl_effect_sigma)=", min(repl_effect_sigma), " max(repl_effect_sigma)=", max(repl_effect_sigma));
  #print("min(protein_effect_lambda)=", min(protein_effect_lambda), " max(protein_effect_lambda)=", max(protein_effect_lambda));
  #print("min(condition_shift_unscaled)=", min(condition_shift_unscaled), " max(condition_shift_unscaled)=", max(condition_shift_unscaled));
  #print("min(experiment_shift_unscaled)=", min(experiment_shift_unscaled), " max(experiment_shift_unscaled)=", max(experiment_shift_unscaled));
  #print("min(log_protein_base_unscaled)=", min(log_protein_base_unscaled), " max(log_protein_base_unscaled)=", max(log_protein_base_unscaled));
  #print("min(log_protein_effect)=", min(log_protein_effect), " max(log_protein_effect)=", max(log_protein_effect));
  #print("min(log_protein_condition)=", min(log_protein_condition), " max(log_protein_condition)=", max(log_protein_condition));
  #print("min(log_protein_repl_effect_unscaled)=", min(log_protein_repl_effect_unscaled),
  #      " max(log_protein_repl_effect_unscaled)=", max(log_protein_repl_effect_unscaled));
}

model {
    # abundance distribution
    #print("peak_param=", peak_param)
    #print("peak_param_sigma=", param_sigma)
    #protein_base ~ normal(zShift, 1.0);
    protein_shift_sigma ~ inv_gamma(2.0, 0.33/zScale); # mode is 1/zScale
    protein_shift0 ~ normal(0, protein_shift_sigma);
    # treatment effect parameters, horseshoe prior
    protein_effect_tau ~ student_t(3, 0.0, 1.0);
    protein_effect_lambda ~ student_t(3, 0.0, 1.0);
    protein_effect_unscaled ~ normal(0.0, 1.0);
    # batch effect parameters, cauchy prior on sigma
    #condition_repl_effect_sigma ~ inv_gamma(1.5, 1.0);
    #print("after params: ", get_lp())

    underdef_protein_shift ~ normal(0.0, 10.0);

    #protein_batch_effect_sigma ~ cauchy(0.0, 0.1);
    #protein_batch_effect ~ normal(0.0, protein_batch_effect_sigma[protein_batch_effect2batch_effect]);

    # calculate the likelihood
    {
        vector[Nquanted] qLogAbu;
        vector[Nmissed] mLogAbu;

        #qLogAbu = iaction_repl_shift[quant2observation] + experiment_shift[observation2experiment[quant2observation]];
        #mLogAbu = iaction_repl_shift[miss2observation] + experiment_shift[observation2experiment[miss2observation]];
        qLogAbu = iaction_shift[quant2iaction] + experiment_shift[observation2experiment[quant2observation]];
        mLogAbu = iaction_shift[miss2iaction] + experiment_shift[observation2experiment[miss2observation]];
        #print( dot_self(to_vector(log_oData) - oAbu - oAdjust) );
        #print("oAbu=", oAbu);
        #oAbuNorm = oAbu .* oInvDataSigma - oDataNorm;
        #print("oAbuNorm=", oAbuNorm );
        #print("min(oAbuNorm)=", min(oAbuNorm), " max(oAbuNorm)=", max(oAbuNorm) );
        exp(qLogAbu - qLogStd) - qDataNorm ~ double_exponential(0, 1);
        #print("after oAbu: ", get_lp())
        # model missing data
        0 ~ bernoulli_logit(mLogAbu * (zScale * zDetectionFactor) - mzShift * zScale * zDetectionFactor + zDetectionIntercept);
        //print("after mZscore: ", get_lp())
    }
}
