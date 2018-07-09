functions {
    real intensity_log_std(real z, real scaleHi, real scaleLo, real offset, real bend, real smooth) {
        return 0.5*(scaleHi+scaleLo)*(z-bend) + 0.5*(scaleHi-scaleLo)*sqrt((z-bend)*(z-bend)+smooth) + offset;
    }
}

data {
  int<lower=1> Nexperiments;    // number of experiments
  int<lower=1> Nconditions;     // number of experimental conditions
  int<lower=0> Nobjects;        // number of objects (proteins/peptides/sites etc)
  int<lower=0> Nsubobjects;     // number of objects subcomponents (peptides of proteins etc), 0 if not supported
  int<lower=0> Nmsprotocols;    // number of MS protocols used
  int<lower=0> Niactions;       // number of interactions (observed objectXcondition pairs)
  int<lower=0> Nobservations;   // number of observations of interactions (objectXexperiment pairs for all iactions and experiments of its condition)
  int<lower=0> Neffects;        // number of effects (that define conditions)
  int<lower=0> NreplEffects;    // number of repl effects (that define biological experimental variation, iaction_repl_shift)
  int<lower=0> NbatchEffects;   // number of batch effects (that define assay experimental variation, but not biology)
  int<lower=0> NunderdefObjs;   // number of virtual interactions (the ones not detected but required for comparison)
  int<lower=1,upper=Nobjects> suo2obj[Nsubobjects];
  int<lower=1,upper=Nobjects> iaction2obj[Niactions];
  int<lower=1,upper=Nconditions> iaction2condition[Niactions];
  int<lower=1,upper=Nobjects> underdef_objs[NunderdefObjs];

  matrix[Neffects, Nconditions] effectXcondition;
  matrix<lower=0>[Neffects, Nconditions] inv_effectXcondition;
  matrix[NreplEffects, Nexperiments] replEffectXexperiment;
  matrix[NbatchEffects, Nexperiments] batchEffectXexperiment;
  int<lower=1,upper=Nconditions> experiment2condition[Nexperiments];
  vector[Nexperiments] experiment_shift;

  int<lower=1,upper=Nexperiments> observation2experiment[Nobservations];
  int<lower=1,upper=Niactions> observation2iaction[Nobservations];
  int<lower=1,upper=Nmsprotocols> experiment2msproto[Nmsprotocols > 0 ? Nexperiments : 0];

  // map from labelXreplicateXobject to observed/missed data
  int<lower=0> Nquanted;        // total number of quantified subobjectsXexperiments
  int<lower=1,upper=Nobservations>  quant2observation[Nquanted];
  int<lower=1,upper=Nsubobjects> quant2suo[Nsubobjects > 0 ? Nquanted : 0];
  int<lower=0> Nmissed;         // total number of missed subobjectsXexperiments
  int<lower=1,upper=Nobservations> miss2observation[Nmissed];
  int<lower=1,upper=Nsubobjects> miss2suo[Nsubobjects > 0 ? Nmissed : 0];

  int<lower=0> NobjEffects;
  int<lower=1,upper=Neffects> obj_effect2effect[NobjEffects];
  int<lower=1,upper=Nobjects> obj_effect2obj[NobjEffects];
  int<lower=0,upper=1> effect_is_positive[Neffects];
  int<lower=0> NeffectsPerObjCumsum[Nobjects+1];

  int<lower=0> NobjReplEffects;
  int<lower=1,upper=NreplEffects> obj_repl_effect2repl_effect[NobjReplEffects];
  int<lower=0> NreplEffectsPerObjCumsum[NobjReplEffects > 0 ? Nobjects+1 : 1];

  int<lower=0> NobjBatchEffects;
  int<lower=1,upper=NbatchEffects> obj_batch_effect2batch_effect[NobjBatchEffects];
  int<lower=0,upper=1> batch_effect_is_positive[NbatchEffects];
  int<lower=0> NbatchEffectsPerObjCumsum[NobjBatchEffects > 0 ? Nobjects+1 : 1];

  vector<lower=0>[Nquanted] qData; // quanted data

  // global model constants
  real global_labu_shift;   // shift to be applied to all XXX_labu variables to get the real log intensity
  vector<lower=0>[Neffects] effect_tau;
  real<lower=0> obj_repl_shift_tau;
  real<lower=0> obj_batch_effect_tau;
  real<lower=0> obj_base_labu_sigma;
  real<upper=0> underdef_obj_shift;

  // instrument calibrated parameters
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
  real mzShift; // zShift for the missing observation intensity (zShift shifted by obj_base)
  vector[Nquanted] zScore; // log(qData) transformed in zScore
  vector[Nquanted] qLogStd; // log(sd(qData))-obj_base
  vector<lower=0>[Nquanted] qDataNorm; // qData/sd(qData)
  matrix[Nexperiments, Nconditions] experimentXcondition;

  int<lower=0> NconditionExperiments[Nconditions];
  matrix[Nexperiments, Nexperiments-Nconditions] experimentXrepl_effect;

  int<lower=0> NobjectSubs[Nobjects];
  matrix[Nsubobjects, Nsubobjects > 0 ? Nsubobjects - Nobjects : 0] suoXsuo_effect;
  vector[Nsubobjects] zero_suo_shifts;

  int<lower=1,upper=Niactions> quant2iaction[Nquanted];
  int<lower=1,upper=Nexperiments> quant2experiment[Nquanted];
  int<lower=1,upper=Nmsprotocols*Nsubobjects> quant2msprotoXsuo[((Nmsprotocols > 1) && (Nsubobjects > 0)) ? Nquanted : 0];
  int<lower=1,upper=Niactions> miss2iaction[Nmissed];
  int<lower=1,upper=Nexperiments> miss2experiment[Nmissed];
  int<lower=1,upper=Nmsprotocols*Nsubobjects> miss2msprotoXsuo[((Nmsprotocols > 1) && (Nsubobjects > 0)) ? Nmissed : 0];
  int<lower=0,upper=NobjEffects> NobjEffectsPos;
  int<lower=0,upper=NobjEffects> NobjEffectsOther;
  int<lower=1,upper=NobjEffects> obj_effect_reshuffle[NobjEffects];
  vector<lower=0>[NobjEffects] obj_effect_tau;

  int<lower=0,upper=NobjBatchEffects> NobjBatchEffectsPos;
  int<lower=0,upper=NobjBatchEffects> NobjBatchEffectsOther;
  int<lower=1,upper=NobjBatchEffects> obj_batch_effect_reshuffle[NobjBatchEffects];

  NobjEffectsPos = sum(effect_is_positive[obj_effect2effect]);
  NobjEffectsOther = NobjEffects - NobjEffectsPos;
  {
    int cur_pos_eff;
    int cur_other_eff;
    cur_pos_eff = 0;
    cur_other_eff = NobjEffectsPos;
    for (i in 1:NobjEffects) {
      if (effect_is_positive[obj_effect2effect[i]]) {
        cur_pos_eff = cur_pos_eff + 1;
        obj_effect_reshuffle[i] = cur_pos_eff;
      } else {
        cur_other_eff = cur_other_eff + 1;
        obj_effect_reshuffle[i] = cur_other_eff;
      }
    }
  }
  obj_effect_tau = effect_tau[obj_effect2effect];

  NobjBatchEffectsPos = sum(batch_effect_is_positive[obj_batch_effect2batch_effect]);
  NobjBatchEffectsOther = NobjBatchEffects - NobjBatchEffectsPos;
  {
    int cur_pos_eff;
    int cur_other_eff;
    cur_pos_eff = 0;
    cur_other_eff = NobjBatchEffectsPos;
    for (i in 1:NobjBatchEffects) {
      if (batch_effect_is_positive[obj_batch_effect2batch_effect[i]]) {
        cur_pos_eff = cur_pos_eff + 1;
        obj_batch_effect_reshuffle[i] = cur_pos_eff;
      } else {
        cur_other_eff = cur_other_eff + 1;
        obj_batch_effect_reshuffle[i] = cur_other_eff;
      }
    }
  }

  {
    vector[Nquanted] qLogData;
    qLogData = log(qData);
    zScore = (qLogData - zShift) * zScale;
    mzShift = zShift - global_labu_shift;

    // process the intensity data to optimize likelihood calculation
    for (i in 1:Nquanted) {
      qLogStd[i] = intensity_log_std(zScore[i], sigmaScaleHi, sigmaScaleLo, sigmaOffset, sigmaBend, sigmaSmooth);
      qDataNorm[i] = exp(qLogData[i] - qLogStd[i]);
      qLogStd[i] = qLogStd[i] - global_labu_shift; // obs_labu is modeled without obj_base
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

  if (Nsubobjects > 0) {
    int suo_index[Nsubobjects];
    int object_shift0_offset[Nobjects+1];

    NobjectSubs = rep_array(0, Nobjects);
    for (i in 1:Nsubobjects) {
      int obj_ix;
      obj_ix = sub2obj[i];
      NobjectSubs[obj_ix] = NobjectSubs[obj_ix] + 1;
      suo_index[i] = NobjectSubs[obj_ix];
    }

    object_shift0_offset[1] = 0;
    for (i in 1:Nobjects) {
      object_shift0_offset[i+1] = object_shift0_offset[i] + NobjectSubs[i]-1;
    }

    suoXsuo_effect = rep_matrix(0.0, Nsubobjects, Nsubobjects-Nobjects);
    for (i in 1:Nsubobjects-Nobjects) {
      for (j in 1:Nsubobjects) {
        int obj_ix;
        obj_ix = sub2obj[j];
        if (object_shift0_offset[obj_ix]<i && i<=object_shift0_offset[obj_ix+1]) {
          int param_ix;
          int n_obj_subs;
          param_ix = i - object_shift0_offset[obj_ix];
          n_obj_subs = NobjectSubs[obj_ix];
          if (suo_index[j] == param_ix) {
            suoXsuo_effect[j, i] = (n_obj_subs-1.0)/n_obj_subs;
          } else {
            suoXsuo_effect[j, i] = -1.0/n_obj_subs;
          }
        } else {
          suoXsuo_effect[j, i] = 0.0;
        }
      }
    }
    if (Nmsprotocols > 1) {
      for (i in 1:Nquanted) {
        quant2msprotoXsub[i] = (experiment2msproto[quant2experiment[i]]-1)*Nsubobjects + quant2sub[i];
      }
      for (i in 1:Nmissed) {
        miss2msprotoXsub[i] = (experiment2msproto[miss2experiment[i]]-1)*Nsubobjects + miss2sub[i];
      }
      zero_suo_shifts = rep_vector(0.0, Nsubobjects);
    }
  }
}

parameters {
  //real obj_base;
  //real<lower=-20, upper=-2> underdef_obj_shift;
  //real<lower=0, upper=5.0> obj_shift_sigma;
  //vector<lower=0>[Nconditions] condition_repl_effect_sigma;

  row_vector[Nobjects] obj_base_labu0; // baseline object abundance without underdefinedness adjustment
  real<lower=0.0> suo_shift_sigma;
  vector[Nsubobjects > 0 ? Nsubobjects-Nobjects : 0] suo_labu0_shift_unscaled; // subcomponent shift within object
  real<lower=0.0> suo_msproto_shift_sigma;
  matrix[Nsubobjects > 0 ? Nsubobjects : 0, Nmsprotocols > 1 ? Nmsprotocols-1 : 0] suo_msproto_shift_unscaled;

  //real<lower=0.0> obj_effect_tau;
  vector<lower=0.0>[NobjEffects] obj_effect_lambda_t;
  vector<lower=0.0>[NobjEffects] obj_effect_lambda_a;
  vector<lower=0.0>[NobjEffectsPos] obj_effect_unscaled_pos;
  vector[NobjEffectsOther] obj_effect_unscaled_other;

  //real<lower=0> obj_repl_effect_sigma;
  //vector<lower=0>[Nobjects*Nexperiments] repl_shift_lambda;
  vector<lower=0>[NreplEffects > 0 ? Nconditions : 0] obj_repl_shift_lambda_t;
  vector<lower=0>[NreplEffects > 0 ? Nconditions : 0] obj_repl_shift_lambda_a;
  vector[NobjReplEffects] obj_repl_effect_unscaled;

  //real<lower=0> obj_batch_effect_sigma;
  vector<lower=0>[NobjBatchEffects] obj_batch_effect_lambda_t;
  vector<lower=0>[NobjBatchEffects] obj_batch_effect_lambda_a;
  vector<lower=0.0>[NobjBatchEffectsPos] obj_batch_effect_unscaled_pos;
  vector[NobjBatchEffectsOther] obj_batch_effect_unscaled_other;
}

transformed parameters {
  row_vector[Nobjects] obj_base_labu;
  vector[Nsubobjects] suo_shift_unscaled; // subcomponent shift within object
  vector[NobjEffects] obj_effect;
  vector<lower=0>[NobjEffects] obj_effect_lambda;
  vector[NobjBatchEffects] obj_batch_effect;

  vector[Niactions] iaction_labu;

  vector[NreplEffects > 0 ? Nconditions : 0] obj_repl_shift_sigma;
  matrix[Nobjects, Nexperiments] objXexp_repl_shift_unscaled; // unscaled replicate shifts for all potential observations (including unobserved)
  vector[Nobservations] obs_labu; // iaction_labu + objXexp_repl_shift * obj_repl_shift_sigma

  matrix[Nobjects, Nexperiments] objXexp_batch_shift; // batch shifts for all potential observations (including unobserved)

  // correct baseline abundances of underdefined objects
  obj_base_labu = obj_base_labu0;
  for (i in 1:NunderdefObjs) {
    obj_base_labu[underdef_objs[i]] = obj_base_labu[underdef_objs[i]] + underdef_obj_shift;
  }

  // calculate effects lambdas and scale effects
  obj_effect_lambda = obj_effect_lambda_a ./ sqrt(obj_effect_lambda_t);
  obj_effect = append_row(obj_effect_unscaled_pos, obj_effect_unscaled_other)[obj_effect_reshuffle] .* obj_effect_lambda .* obj_effect_tau;

  // calculate iaction_labu
  {
    matrix[Nobjects, Nconditions] objXcondition;

    objXcondition = csr_to_dense_matrix(Nobjects, Neffects, obj_effect,
                                        obj_effect2effect, NeffectsPerObjCumsum) * effectXcondition;
    for (i in 1:Niactions) {
      iaction_labu[i] = obj_base_labu[iaction2obj[i]] + objXcondition[iaction2obj[i], iaction2condition[i]];
    }
  }
  // calculate objXexp_repl_shift and obs_labu
  if (NreplEffects > 0) {
    vector[NobjReplEffects] obj_repl_effect;
    vector[NobjReplEffects] obj_repl_effect_lambda;

    // note: linear transform of obj_repl_effect_unscaled, Jacobian is zero
    obj_repl_shift_sigma = obj_repl_shift_lambda_a ./ sqrt(obj_repl_shift_lambda_t) * obj_repl_shift_tau;
    objXexp_repl_shift_unscaled = csr_to_dense_matrix(Nobjects, NreplEffects, obj_repl_effect_unscaled,
                                        obj_repl_effect2repl_effect, NreplEffectsPerObjCumsum) * replEffectXexperiment;
    for (i in 1:Nobservations) {
      int iaction_ix;
      iaction_ix = observation2iaction[i];
      obs_labu[i] = iaction_labu[iaction_ix]
                  + objXexp_repl_shift_unscaled[iaction2obj[iaction_ix], observation2experiment[i]] *
                    obj_repl_shift_sigma[iaction2condition[iaction_ix]];
    }
  } else {
    obs_labu = iaction_labu[observation2iaction];
    objXexp_repl_shift_unscaled = rep_matrix(0.0, Nobjects, Nexperiments);
  }
  // calculate objXexp_batch_shift (doesn't make sense to add to obs_labu)
  if (NbatchEffects > 0) {
    vector[NobjBatchEffects] obj_batch_effect_lambda;

    obj_batch_effect_lambda = obj_batch_effect_lambda_a ./ sqrt(obj_batch_effect_lambda_t);
    obj_batch_effect = append_row(obj_batch_effect_unscaled_pos, obj_batch_effect_unscaled_other)[obj_batch_effect_reshuffle] .* obj_batch_effect_lambda * obj_batch_effect_tau;

    objXexp_batch_shift = csr_to_dense_matrix(Nobjects, NbatchEffects, obj_batch_effect,
                                              obj_batch_effect2batch_effect, NbatchEffectsPerObjCumsum) * batchEffectXexperiment;
  } else {
    objXexp_batch_shift = rep_matrix(0.0, Nobjects, Nexperiments);
  }
  // calculate suo_shift
  if (Nsubobjects > 1) {
    suo_shift_unscaled = suoXsuo_effect * suo_labu0_shift_unscaled;
  } else if (Nsubobjects == 1) {
    suo_shift_unscaled = rep_vector(0.0, Nsubobjects);
  }
}

model {
    // abundance distribution
    //obj_base ~ normal(zShift, 1.0);
    //obj_shift_sigma ~ inv_gamma(2.0, 0.33/zScale); // mode is 1/zScale
    obj_base_labu0 ~ normal(0, obj_base_labu_sigma);
    // treatment effect parameters, horseshoe prior
    //obj_effect_tau ~ student_t(2, 0.0, 1.0);
    //obj_effect_lambda ~ student_t(2, 0.0, obj_effect_tau);
    obj_effect_lambda_t ~ chi_square(2.0);
    obj_effect_lambda_a ~ normal(0.0, 1.0); // 1.0 = 2/2
    obj_effect_unscaled_pos ~ normal(0.0, 1.0);
    obj_effect_unscaled_other ~ normal(0.0, 1.0);
    // batch effect parameters, cauchy prior on sigma
    //condition_repl_effect_sigma ~ inv_gamma(1.5, 1.0);

    //underdef_obj_shift ~ normal(0.0, 10.0);

    //repl_shift_lambda ~ student_t(2, 0.0, repl_shift_tau);
    //obj_repl_effect ~ normal(0.0, obj_repl_effect_lambda);
    if (NreplEffects > 0) {
      obj_repl_shift_lambda_t ~ chi_square(4.0);
      obj_repl_shift_lambda_a ~ normal(0.0, 1.0);
      to_vector(objXexp_repl_shift_unscaled) ~ normal(0.0, 1.0);
    }
    //to_vector(repl_shift) ~ normal(0.0, repl_shift_lambda);

    //obj_batch_effect_lambda ~ student_t(2, 0.0, obj_batch_effect_tau);
    if (NbatchEffects > 0) {
      obj_batch_effect_lambda_t ~ chi_square(3.0);
      obj_batch_effect_lambda_a ~ normal(0.0, 1.0);
      //obj_batch_effect ~ normal(0.0, obj_batch_effect_lambda);
      obj_batch_effect_unscaled_pos ~ normal(0.0, 1.0);
      obj_batch_effect_unscaled_other ~ normal(0.0, 1.0);
    }
    if (Nsubobjects > 0) {
      suo_shift_sigma ~ cauchy(0.0, 1.0);
      suo_shift_unscaled ~ normal(0.0, 1.0);
      if (Nmsprotocols > 1) {
        suo_msproto_shift_sigma ~ cauchy(0.0, 1.0);
        to_vector(suo_msproto_shift_unscaled) ~ normal(0.0, 1.0);
      }
    }

    // calculate the likelihood
    {
        vector[Nquanted] q_labu;
        vector[Nmissed] m_labu;

        q_labu = obs_labu[quant2observation] + experiment_shift[quant2experiment];
        m_labu = obs_labu[miss2observation] + experiment_shift[miss2experiment];
        //qLogAbu = iaction_shift[quant2iaction] + experiment_shift[observation2experiment[quant2observation]];
        //mLogAbu = iaction_shift[miss2iaction] + experiment_shift[observation2experiment[miss2observation]];
        if (Nsubobjects > 0) {
          // adjust by subcomponent shift
          if (Nmsprotocols > 1) {
            vector[Nmsprotocols*Nsubobjects] suo_shift;
            suo_shift = to_vector(rep_matrix(suo_shift_unscaled * suo_shift_sigma, Nmsprotocols) +
              append_col(zero_suo_shifts, suo_msproto_shift_unscaled * suo_msproto_shift_sigma));
            q_labu = q_labu + suo_shift[quant2msprotoXsub];
            m_labu = m_labu + suo_shift[miss2msprotoXsub];
          } else {
            vector[Nsubobjects] suo_shift;
            suo_shift = suo_shift_unscaled * suo_shift_sigma;
            q_labu = q_labu + suo_shift[quant2sub];
            m_labu = m_labu + suo_shift[miss2sub];
          }
        }
        if (NbatchEffects > 0) {
          // adjust by objXexp_batch_shift
          vector[Nobservations] obs_batch_shift;
          for (i in 1:Nobservations) {
              obs_batch_shift[i] = objXexp_batch_shift[iaction2obj[observation2iaction[i]], observation2experiment[i]];
          }
          q_labu = q_labu + obs_batch_shift[quant2observation];
          m_labu = m_labu + obs_batch_shift[miss2observation];
        }

        qDataNorm ~ double_exponential(exp(q_labu - qLogStd), 1);
        // model quantitations and missing data
        1 ~ bernoulli_logit(q_labu * (zScale * zDetectionFactor) + (-mzShift * zScale * zDetectionFactor + zDetectionIntercept));
        0 ~ bernoulli_logit(m_labu * (zScale * zDetectionFactor) + (-mzShift * zScale * zDetectionFactor + zDetectionIntercept));
    }
}

generated quantities {
  vector[Neffects] effect_shift_replCI_sigma;
  vector[Niactions] iaction_labu_replCI;
  vector[NobjEffects] obj_effect_replCI;

  for (i in 1:Niactions) {
    iaction_labu_replCI[i] = normal_rng(iaction_labu[i], obj_repl_shift_sigma[iaction2condition[i]]);
  }
  effect_shift_replCI_sigma = inv_effectXcondition * obj_repl_shift_sigma;
  for (i in 1:NobjEffects) {
    obj_effect_replCI[i] = normal_rng(obj_effect[i], effect_shift_replCI_sigma[obj_effect2effect[i]]);
  }
}
