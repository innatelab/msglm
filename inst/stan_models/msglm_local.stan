functions {
    real intensity_log_std(real z, real scaleHi, real scaleLo, real offset, real bend, real smooth) {
        return 0.5*(scaleHi+scaleLo)*(z-bend) + 0.5*(scaleHi-scaleLo)*sqrt((z-bend)*(z-bend)+smooth) + offset;
    }
}

data {
  int<lower=1> Nexperiments;    # number of experiments
  int<lower=1> Nconditions;     # number of experimental conditions
  int<lower=0> Nobjects;        # number of objects (proteins/peptides/sites etc)
  int<lower=0> Niactions;       # number of interactions (observed objectXcondition pairs)
  int<lower=0> Nobservations;   # number of observations of interactions (objectXexperiment pairs for all iactions and experiments of its condition)
  int<lower=0> Neffects;        # number of effects (that define conditions)
  int<lower=0> NreplEffects;    # number of repl effects (that define biological experimental variation, iaction_repl_shift)
  int<lower=0> NbatchEffects;   # number of batch effects (that define assay experimental variation, but not biology)
  int<lower=0> NunderdefObjs;   # number of virtual interactions (the ones not detected but required for comparison)
  int<lower=1,upper=Nobjects> iaction2obj[Niactions];
  int<lower=1,upper=Nconditions> iaction2condition[Niactions];
  int<lower=1,upper=Nobjects> underdef_objs[NunderdefObjs];

  matrix[Neffects, Nconditions] effectXcondition;
  matrix[NreplEffects, Nexperiments] replEffectXexperiment;
  matrix[NbatchEffects, Nexperiments] batchEffectXexperiment;
  int<lower=1,upper=Nconditions> experiment2condition[Nexperiments];
  vector[Nexperiments] experiment_shift;

  int<lower=1,upper=Nexperiments> observation2experiment[Nobservations];
  int<lower=1,upper=Niactions> observation2iaction[Nobservations];

  # map from labelXreplicateXobject to observed/missed data
  int<lower=0> Nquanted;        # total number of quantified objectsXexperiments
  int<lower=1,upper=Nobservations>  quant2observation[Nquanted];
  int<lower=0> Nmissed;         # total number of missed objectsXexperiments
  int<lower=1,upper=Nobservations> miss2observation[Nmissed];

  int<lower=0> NobjEffects;
  int<lower=1,upper=Neffects> obj_effect2effect[NobjEffects];
  int<lower=1,upper=Nobjects> obj_effect2obj[NobjEffects];
  int<lower=0> NeffectsPerObjCumsum[Nobjects+1];

  int<lower=0> NobjReplEffects;
  int<lower=1,upper=NreplEffects> obj_repl_effect2repl_effect[NobjReplEffects];
  int<lower=0> NreplEffectsPerObjCumsum[Nobjects+1];

  int<lower=0> NobjBatchEffects;
  int<lower=1,upper=NbatchEffects> obj_batch_effect2batch_effect[NobjBatchEffects];
  int<lower=0> NbatchEffectsPerObjCumsum[NobjBatchEffects > 0 ? Nobjects+1 : 1];

  vector<lower=0>[Nquanted] qData; # quanted data

  # global model constants
  real global_labu_shift;   # shift to be applied to all XXX_labu variables to get the real log intensity
  real<lower=0> obj_effect_tau;
  real<lower=0> obj_repl_shift_tau;
  real<lower=0> obj_batch_effect_tau;
  real<lower=0> obj_base_labu_sigma;
  real<upper=0> underdef_obj_shift;

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
  real mzShift; # zShift for the missing observation intensity (zShift shifted by obj_base)
  vector[Nquanted] zScore; # log(qData) transformed in zScore
  vector[Nquanted] qLogStd; # log(sd(qData))-obj_base
  vector<lower=0>[Nquanted] qDataNorm; # qData/sd(qData)
  matrix[Nexperiments, Nconditions] experimentXcondition;

  int<lower=1> NconditionExperiments[Nconditions];
  matrix[Nexperiments, Nexperiments-Nconditions] experimentXrepl_effect;

  int<lower=1,upper=Niactions> quant2iaction[Nquanted];
  int<lower=1,upper=Nexperiments> quant2experiment[Nquanted];
  int<lower=1,upper=Niactions> miss2iaction[Nmissed];
  int<lower=1,upper=Nexperiments> miss2experiment[Nmissed];

  {
    vector[Nquanted] qLogData;
    qLogData = log(qData);
    zScore = (qLogData - zShift) * zScale;
    mzShift = zShift - global_labu_shift;
  
    # process the intensity data to optimize likelihood calculation
    for (i in 1:Nquanted) {
      qLogStd[i] = intensity_log_std(zScore[i], sigmaScaleHi, sigmaScaleLo, sigmaOffset, sigmaBend, sigmaSmooth);
      qDataNorm[i] = exp(qLogData[i] - qLogStd[i]);
      qLogStd[i] = qLogStd[i] - global_labu_shift; # obs_labu is modeled without obj_base
    }
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
}

parameters {
  #real obj_base;
  #real<lower=-20, upper=-2> underdef_obj_shift;
  #real<lower=0, upper=5.0> obj_shift_sigma;
  #vector<lower=0>[Nconditions] condition_repl_effect_sigma;

  row_vector[Nobjects] obj_base_labu0; # baseline object abundance without underdefinedness adjustment

  #real<lower=0.0> obj_effect_tau;
  vector<lower=0.0>[NobjEffects] obj_effect_lambda_t;
  vector<lower=0.0>[NobjEffects] obj_effect_lambda_a;
  vector[NobjEffects] obj_effect_unscaled;

  #real<lower=0> obj_repl_effect_sigma;
  #vector<lower=0>[Nobjects*Nexperiments] repl_shift_lambda;
  vector<lower=0>[Nconditions] obj_repl_shift_lambda_t;
  vector<lower=0>[Nconditions] obj_repl_shift_lambda_a;
  vector[NobjReplEffects] obj_repl_effect_unscaled;

  #real<lower=0> obj_batch_effect_sigma;
  vector<lower=0>[NobjBatchEffects] obj_batch_effect_lambda_t;
  vector<lower=0>[NobjBatchEffects] obj_batch_effect_lambda_a;
  vector[NobjBatchEffects] obj_batch_effect_unscaled;
}

transformed parameters {
  row_vector[Nobjects] obj_base_labu;
  vector[NobjEffects] obj_effect;
  vector<lower=0>[NobjEffects] obj_effect_lambda;

  vector[Niactions] iaction_labu;

  vector[Nconditions] obj_repl_shift_sigma;
  matrix[Nobjects, Nexperiments] objXexp_repl_shift_unscaled; # unscaled replicate shifts for all potential observations (including unobserved)
  vector[Nobservations] obs_labu; # iaction_labu + objXexp_repl_shift * obj_repl_shift_sigma

  matrix[Nobjects, Nexperiments] objXexp_batch_shift; # batch shifts for all potential observations (including unobserved)

  # correct baseline abundances of underdefined objects
  obj_base_labu = obj_base_labu0;
  for (i in 1:NunderdefObjs) {
    obj_base_labu[underdef_objs[i]] = obj_base_labu[underdef_objs[i]] + underdef_obj_shift;
  }

  # calculate effects lambdas and scale effects
  obj_effect_lambda = obj_effect_lambda_a ./ sqrt(obj_effect_lambda_t);
  obj_effect = obj_effect_unscaled .* obj_effect_lambda * obj_effect_tau;
  obj_repl_shift_sigma = obj_repl_shift_lambda_a ./ sqrt(obj_repl_shift_lambda_t) * obj_repl_shift_tau;

  # calculate iaction_labu
  {
    matrix[Nobjects, Nconditions] objXcondition;

    objXcondition = csr_to_dense_matrix(Nobjects, Neffects, obj_effect,
                                        obj_effect2effect, NeffectsPerObjCumsum) * effectXcondition;
    for (i in 1:Niactions) {
      iaction_labu[i] = obj_base_labu[iaction2obj[i]] + objXcondition[iaction2obj[i], iaction2condition[i]];
    }
  }
  # calculate objXexp_repl_shift and obs_labu
  {
    vector[NobjReplEffects] obj_repl_effect;
    vector[NobjReplEffects] obj_repl_effect_lambda;

    # note: linear transform of obj_repl_effect_unscaled, Jacobian is zero
    objXexp_repl_shift_unscaled = csr_to_dense_matrix(Nobjects, NreplEffects, obj_repl_effect_unscaled,
                                        obj_repl_effect2repl_effect, NreplEffectsPerObjCumsum) * replEffectXexperiment;
    for (i in 1:Nobservations) {
      int iaction_ix;
      iaction_ix = observation2iaction[i];
      obs_labu[i] = iaction_labu[iaction_ix]
                  + objXexp_repl_shift_unscaled[iaction2obj[iaction_ix], observation2experiment[i]] *
                    obj_repl_shift_sigma[iaction2condition[iaction_ix]];
    }
  }
  # calculate objXexp_batch_shift (doesn't make sense to add to obs_labu)
  if (NbatchEffects > 0) {
    vector[NobjBatchEffects] obj_batch_effect;
    vector[NobjBatchEffects] obj_batch_effect_lambda;

    obj_batch_effect_lambda = obj_batch_effect_lambda_a ./ sqrt(obj_batch_effect_lambda_t);
    obj_batch_effect = obj_batch_effect_unscaled .* obj_batch_effect_lambda * obj_batch_effect_tau;

    objXexp_batch_shift = csr_to_dense_matrix(Nobjects, NbatchEffects, obj_batch_effect,
                                              obj_batch_effect2batch_effect, NbatchEffectsPerObjCumsum) * batchEffectXexperiment;
  } else {
    objXexp_batch_shift = rep_matrix(0.0, Nobjects, Nexperiments);
  }
}

model {
    # abundance distribution
    #obj_base ~ normal(zShift, 1.0);
    #obj_shift_sigma ~ inv_gamma(2.0, 0.33/zScale); # mode is 1/zScale
    obj_base_labu0 ~ normal(0, obj_base_labu_sigma);
    # treatment effect parameters, horseshoe prior
    #obj_effect_tau ~ student_t(2, 0.0, 1.0);
    #obj_effect_lambda ~ student_t(2, 0.0, obj_effect_tau);
    obj_effect_lambda_t ~ gamma(1.0, 1.0); // 1.0 = 2/2
    obj_effect_lambda_a ~ normal(0.0, 1.0); // 1.0 = 2/2
    obj_effect_unscaled ~ normal(0.0, 1.0);
    # batch effect parameters, cauchy prior on sigma
    #condition_repl_effect_sigma ~ inv_gamma(1.5, 1.0);

    #underdef_obj_shift ~ normal(0.0, 10.0);

    #repl_shift_lambda ~ student_t(2, 0.0, repl_shift_tau);
    #obj_repl_effect_lambda ~ student_t(2, 0.0, obj_repl_effect_tau);
    obj_repl_shift_lambda_t ~ gamma(1.0, 1.0);
    obj_repl_shift_lambda_a ~ normal(0.0, 1.0);
    #obj_repl_effect ~ normal(0.0, obj_repl_effect_lambda);
    to_vector(objXexp_repl_shift_unscaled) ~ normal(0.0, 1.0);
    #to_vector(repl_shift) ~ normal(0.0, repl_shift_lambda);

    #obj_batch_effect_lambda ~ student_t(2, 0.0, obj_batch_effect_tau);
    obj_batch_effect_lambda_t ~ gamma(1.0, 1.0);
    obj_batch_effect_lambda_a ~ normal(0.0, 1.0);
    #obj_batch_effect ~ normal(0.0, obj_batch_effect_lambda);
    obj_batch_effect_unscaled ~ normal(0.0, 1.0);

    # calculate the likelihood
    {
        vector[Nquanted] q_labu;
        vector[Nmissed] m_labu;

        q_labu = obs_labu[quant2observation] + experiment_shift[observation2experiment[quant2observation]];
        m_labu = obs_labu[miss2observation] + experiment_shift[observation2experiment[miss2observation]];
        #qLogAbu = iaction_shift[quant2iaction] + experiment_shift[observation2experiment[quant2observation]];
        #mLogAbu = iaction_shift[miss2iaction] + experiment_shift[observation2experiment[miss2observation]];
        if (NbatchEffects > 0) {
          # adjust by objXexp_batch_shift
          vector[Nobservations] obs_batch_shift;
          for (i in 1:Nobservations) {
              obs_batch_shift[i] = objXexp_batch_shift[iaction2obj[observation2iaction[i]], observation2experiment[i]];
          }
          q_labu = q_labu + obs_batch_shift[quant2observation];
          m_labu = m_labu + obs_batch_shift[miss2observation];
        }

        qDataNorm ~ double_exponential(exp(q_labu - qLogStd), 1);
        # model quantitations and missing data
        1 ~ bernoulli_logit(q_labu * (zScale * zDetectionFactor) + (-mzShift * zScale * zDetectionFactor + zDetectionIntercept));
        0 ~ bernoulli_logit(m_labu * (zScale * zDetectionFactor) + (-mzShift * zScale * zDetectionFactor + zDetectionIntercept));
    }
}
