functions {
    real intensity_log_std(real z, real scaleHi, real scaleLo, real offset, real bend, real smooth) {
        return 0.5*(scaleHi+scaleLo)*(z-bend) + 0.5*(scaleHi-scaleLo)*sqrt((z-bend)*(z-bend)+smooth) + offset;
    }
}

data {
  int<lower=1> Nexperiments;    // number of experiments
  int<lower=1> Nconditions;     // number of experimental conditions
  int<lower=0> Nobjects;        // number of objects (proteins/peptides/sites etc)
  int<lower=0> Niactions;       // number of interactions (observed objectXcondition pairs)
  int<lower=0> Nobservations;   // number of observations of interactions (objectXexperiment pairs for all iactions and experiments of its condition)
  int<lower=0> Neffects;        // number of effects (that define conditions)
  int<lower=0> NbatchEffects;   // number of batch effects (that define assay experimental variation, but not biology)
  int<lower=0> NunderdefObjs;   // number of virtual interactions (the ones not detected but required for comparison)
  int<lower=1,upper=Nobjects> iaction2obj[Niactions];
  int<lower=1,upper=Nconditions> iaction2condition[Niactions];
  int<lower=1,upper=Nobjects> underdef_objs[NunderdefObjs];

  vector[Nexperiments] experiment_shift;

  int<lower=1,upper=Nexperiments> observation2experiment[Nobservations];
  int<lower=1,upper=Niactions> observation2iaction[Nobservations];

  // map from labelXreplicateXobject to observed/missed data
  int<lower=0> Nquanted;        // total number of quantified objectsXexperiments
  int<lower=1,upper=Nobservations>  quant2observation[Nquanted];
  int<lower=0> Nmissed;         // total number of missed objectsXexperiments
  int<lower=1,upper=Nobservations> miss2observation[Nmissed];

  int<lower=0> NobjEffects;
  int<lower=1,upper=Neffects> obj_effect2effect[NobjEffects];
  int<lower=0,upper=1> effect_is_positive[Neffects];

  int<lower=0> NobjBatchEffects;
  int<lower=1,upper=NbatchEffects> obj_batch_effect2batch_effect[NobjBatchEffects];
  int<lower=0,upper=1> batch_effect_is_positive[NbatchEffects];

  // iactXobjeff (interaction X object_effect) sparse matrix
  int<lower=0> iactXobjeff_Nw;
  vector[iactXobjeff_Nw] iactXobjeff_w;
  int<lower=0, upper=iactXobjeff_Nw+1> iactXobjeff_u[Niactions+1];
  int<lower=0, upper=NobjEffects> iactXobjeff_v[iactXobjeff_Nw];

  // obsXobj_batcheff (observation X batch_effect) sparse matrix
  int<lower=0> obsXobjbatcheff_Nw;
  vector[obsXobjbatcheff_Nw] obsXobjbatcheff_w;
  int<lower=0, upper=obsXobjbatcheff_Nw+1> obsXobjbatcheff_u[Nobservations + 1];
  int<lower=0, upper=NobjBatchEffects> obsXobjbatcheff_v[obsXobjbatcheff_Nw];

  vector<lower=0>[Nquanted] qData; // quanted data

  // global model constants
  real global_labu_shift;   // shift to be applied to all XXX_labu variables to get the real log intensity
  vector<lower=0>[Neffects] effect_tau;
  real<lower=0> obj_base_repl_shift_tau;
  real<lower=0> obj_effect_repl_shift_tau;
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

  int<lower=1,upper=Niactions> quant2iaction[Nquanted];
  int<lower=1,upper=Nexperiments> quant2experiment[Nquanted];
  int<lower=1,upper=Niactions> miss2iaction[Nmissed];
  int<lower=1,upper=Nexperiments> miss2experiment[Nmissed];
  int<lower=0,upper=NobjEffects> NobjEffectsPos;
  int<lower=0,upper=NobjEffects> NobjEffectsOther;
  int<lower=1,upper=NobjEffects> obj_effect_reshuffle[NobjEffects];
  vector<lower=0>[NobjEffects] obj_effect_tau;

  int<lower=0,upper=NobjBatchEffects> NobjBatchEffectsPos;
  int<lower=0,upper=NobjBatchEffects> NobjBatchEffectsOther;
  int<lower=1,upper=NobjBatchEffects> obj_batch_effect_reshuffle[NobjBatchEffects];

  vector[iactXobjeff_Nw] iactXobjeff4sigma_w;

  vector[Niactions] iactXobjbase_w;
  int<lower=0> iactXobjbase_u[Niactions + 1];

  vector[Nobservations] obsXiact_w;
  int<lower=0> obsXiact_u[Nobservations + 1];

  int<lower=0> Nobservations0;  // number of observations degrees of freedom ()

  vector[2*(Nobservations - Niactions)] obsXobs_shift0_w;
  int<lower=1, upper=2*(Nobservations - Niactions) + 1> obsXobs_shift0_u[Nobservations + 1];
  int<lower=1, upper=Nobservations - Niactions> obsXobs_shift0_v[2*(Nobservations - Niactions)];

  // prepare reshuffling of positive/other obj effects
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

  // prepare reshuffling of positive/other batch effects
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

  // preprocess signals (MS noise)
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

  // prepare obsXobs_shift0
  {
    int iaction2nobs[Niactions];
    int iaction2maxobsix[Niactions];

    int iaction2nobs_2ndpass[Niactions];
    int iaction2obseff_offset[Niactions];
    int iaction2obsXobseff_offset[Niactions];

    int obsXobs_shift0_offset;
    int obs_shift0_offset;

    iaction2nobs = rep_array(0, Niactions);
    iaction2maxobsix = rep_array(0, Niactions);
    for (i in 1:Nobservations) {
      int iact_ix;
      iact_ix = observation2iaction[i];
      iaction2nobs[iact_ix] += 1;
      if (iaction2maxobsix[iact_ix] < i) {
        iaction2maxobsix[iact_ix] = i;
      }
    }
    //print("iaction2nobs=", iaction2nobs);
    //print("iaction2maxobsix=", iaction2maxobsix);
    obs_shift0_offset = 0;
    obsXobs_shift0_offset = 0;
    obsXobs_shift0_u[1] = 1;
    iaction2nobs_2ndpass = rep_array(0, Niactions);
    iaction2obseff_offset = rep_array(0, Niactions);
    iaction2obsXobseff_offset = rep_array(0, Niactions);
    for (i in 1:Nobservations) {
      int iact_ix;
      int iact_nobs;
      iact_ix = observation2iaction[i];
      //print("iact_ix[", i, "]=", iact_ix);
      iact_nobs = iaction2nobs[iact_ix];
      if (iact_nobs > 1) {
        if (iaction2maxobsix[iact_ix] == i) {
          //print("last(", iact_nobs, ") of iact_ix=", iact_ix);
          // last observation
          obsXobs_shift0_u[i+1] = obsXobs_shift0_u[i] + (iact_nobs-1);
          //print("obsXobs_shift0_u[", i+1, "]=", obsXobs_shift0_u[i+1]);
          for (j in 1:(iaction2nobs[iact_ix]-1)) {
            obsXobs_shift0_v[iaction2obsXobseff_offset[iact_ix] + (iact_nobs-1) + j] = iaction2obseff_offset[iact_ix] + j;
            obsXobs_shift0_w[iaction2obsXobseff_offset[iact_ix] + (iact_nobs-1) + j] = -1.0;
          }
        } else {
          int obs_ix;
          if (iaction2nobs_2ndpass[iact_ix] == 0) {
            // reserve nobs-1 obs_shift0 variables
            iaction2obseff_offset[iact_ix] = obs_shift0_offset;
            obs_shift0_offset += iact_nobs-1;
            //print("obs_shift0_offset=", obs_shift0_offset);

            // reserve 2(nobs-1) obsXobs_shift0 variables
            iaction2obsXobseff_offset[iact_ix] = obsXobs_shift0_offset;
            obsXobs_shift0_offset += 2*(iact_nobs-1);
            //print("obsXobs_shift0_offset=", obsXobs_shift0_offset);
          }
          iaction2nobs_2ndpass[iact_ix] += 1;
          obs_ix = iaction2nobs_2ndpass[iact_ix];
          //print(obs_ix, "-th of iact_ix=", iact_ix);
          obsXobs_shift0_u[i+1] = obsXobs_shift0_u[i] + 1;
          obsXobs_shift0_w[iaction2obsXobseff_offset[iact_ix] + obs_ix] = 1.0;
          obsXobs_shift0_v[iaction2obsXobseff_offset[iact_ix] + obs_ix] = iaction2obseff_offset[iact_ix] + obs_ix;
        }
      }
    }
  }
  Nobservations0 = Nobservations - Niactions;
  //print("obsXobs_shift0=", csr_to_dense_matrix(Nobservations, Nobservations0, obsXobs_shift0_w, obsXobs_shift0_v, obsXobs_shift0_u));

  iactXobjeff4sigma_w = square(iactXobjeff_w);
  iactXobjbase_w = rep_vector(1.0, Niactions);
  for (i in 1:(Niactions+1)) iactXobjbase_u[i] = i;
  obsXiact_w = rep_vector(1.0, Nobservations);
  for (i in 1:(Nobservations+1)) obsXiact_u[i] = i;
}

parameters {
  //real obj_base;
  //real<lower=-20, upper=-2> underdef_obj_shift;
  //real<lower=0, upper=5.0> obj_shift_sigma;
  //vector<lower=0>[Nconditions] condition_repl_effect_sigma;

  vector[Nobjects] obj_base_labu0; // baseline object abundance without underdefinedness adjustment
  vector<lower=0.01>[Nobservations0 > 0 ? Nobjects : 0] obj_base_repl_shift_sigma;

  //real<lower=0.0> obj_effect_tau;
  vector<lower=0.0>[NobjEffects] obj_effect_lambda_t;
  vector<lower=0.0>[NobjEffects] obj_effect_lambda_a;
  vector<lower=0.0>[NobjEffectsPos] obj_effect_unscaled_pos;
  vector[NobjEffectsOther] obj_effect_unscaled_other;

  //real<lower=0> obj_repl_effect_sigma;
  //vector<lower=0>[Nobjects*Nexperiments] repl_shift_lambda;
  vector<lower=0.01>[Nobservations0 > 0 ? NobjEffects : 0] obj_effect_repl_shift_sigma;
  vector[Nobservations0] obs_shift0;

  //real<lower=0> obj_batch_effect_sigma;
  vector<lower=0>[NobjBatchEffects] obj_batch_effect_lambda_t;
  vector<lower=0>[NobjBatchEffects] obj_batch_effect_lambda_a;
  vector<lower=0.0>[NobjBatchEffectsPos] obj_batch_effect_unscaled_pos;
  vector[NobjBatchEffectsOther] obj_batch_effect_unscaled_other;
}

transformed parameters {
  vector[Nobjects] obj_base_labu;
  vector[NobjEffects] obj_effect;
  vector<lower=0>[NobjEffects] obj_effect_sigma;
  vector[NobjBatchEffects] obj_batch_effect;

  vector[Niactions] iaction_labu;

  vector[Nobservations0 > 0 ? Niactions : 0] iact_repl_shift_sigma;
  vector[Nobservations] obs_labu; // iaction_labu + objXexp_repl_shift * obj_repl_shift_sigma
  vector[Nobservations0 > 0 ? Nobservations : 0] obs_repl_shift; // replicate shifts for all potential observations (including missing)
  vector[NobjBatchEffects > 0 ? Nobservations : 0] obs_batch_shift;

  // correct baseline abundances of underdefined objects
  obj_base_labu = obj_base_labu0;
  for (i in 1:NunderdefObjs) {
    obj_base_labu[underdef_objs[i]] = obj_base_labu[underdef_objs[i]] + underdef_obj_shift;
  }

  // calculate effects lambdas and scale effects
  obj_effect_sigma = obj_effect_lambda_a ./ sqrt(obj_effect_lambda_t) .* obj_effect_tau;
  obj_effect = append_row(obj_effect_unscaled_pos, obj_effect_unscaled_other)[obj_effect_reshuffle] .* obj_effect_sigma;

  // calculate iaction_labu
  iaction_labu = csr_matrix_times_vector(Niactions, Nobjects, iactXobjbase_w, iaction2obj, iactXobjbase_u, obj_base_labu) +
                 csr_matrix_times_vector(Niactions, NobjEffects, iactXobjeff_w, iactXobjeff_v, iactXobjeff_u, obj_effect);

  obs_labu = csr_matrix_times_vector(Nobservations, Niactions, obsXiact_w, observation2iaction, obsXiact_u, iaction_labu);
  // calculate obs_shift and obs_labu
  if (Nobservations0 > 0) {
    // FIXME: non-linear transform of obj_effect_repl_shift_sigma, Jacobian is not zero
    iact_repl_shift_sigma = sqrt(csr_matrix_times_vector(Niactions, NobjEffects, iactXobjeff4sigma_w, iactXobjeff_v, iactXobjeff_u,
                                                         square(obj_effect_repl_shift_sigma)) +
                                 csr_matrix_times_vector(Niactions, Nobjects, iactXobjbase_w, iaction2obj, iactXobjbase_u,
                                                         square(obj_base_repl_shift_sigma)));
    obs_repl_shift = csr_matrix_times_vector(Nobservations, Nobservations0, obsXobs_shift0_w, obsXobs_shift0_v, obsXobs_shift0_u, obs_shift0);
    obs_labu += obs_repl_shift;
  }
  // calculate objXexp_batch_shift (doesn't make sense to add to obs_labu)
  if (NbatchEffects > 0) {
    vector[NobjBatchEffects] obj_batch_effect_sigma;

    obj_batch_effect_sigma = obj_batch_effect_lambda_a ./ sqrt(obj_batch_effect_lambda_t) * obj_batch_effect_tau;
    obj_batch_effect = append_row(obj_batch_effect_unscaled_pos, obj_batch_effect_unscaled_other)[obj_batch_effect_reshuffle] .* obj_batch_effect_sigma;
    obs_batch_shift = csr_matrix_times_vector(Nobservations, NobjBatchEffects, obsXobjbatcheff_w, obsXobjbatcheff_v, obsXobjbatcheff_u, obj_batch_effect);
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
    if (Nobservations0 > 0) {
      vector[Nobservations] obs_repl_shift_sigma;
      obj_base_repl_shift_sigma ~ student_t(4, 0.0, obj_base_repl_shift_tau);
      obj_effect_repl_shift_sigma ~ student_t(4, 0.0, obj_effect_repl_shift_tau);
      obs_repl_shift_sigma = csr_matrix_times_vector(Nobservations, Niactions,
              obsXiact_w, observation2iaction, obsXiact_u, iact_repl_shift_sigma);
      //print("iact_repl_shift_sigma=", iact_repl_shift_sigma);
      //print("obsXiact=", csr_to_dense_matrix(Nobservations, Niactions,
      //          obsXiact_w, observation2iaction, obsXiact_u));
      obs_repl_shift ~ normal(0.0, obs_repl_shift_sigma);
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

    // calculate the likelihood
    {
        vector[Nquanted] q_labu;
        vector[Nmissed] m_labu;

        q_labu = obs_labu[quant2observation] + experiment_shift[observation2experiment[quant2observation]];
        m_labu = obs_labu[miss2observation] + experiment_shift[observation2experiment[miss2observation]];
        //qLogAbu = iaction_shift[quant2iaction] + experiment_shift[observation2experiment[quant2observation]];
        //mLogAbu = iaction_shift[miss2iaction] + experiment_shift[observation2experiment[miss2observation]];
        if (NbatchEffects > 0) {
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
    vector[Nobjects] obj_base_labu_replCI;
    vector[NobjEffects] obj_effect_replCI;
    vector[Niactions] iaction_labu_replCI;

    for (i in 1:Nobjects) {
        obj_base_labu_replCI[i] = normal_rng(obj_base_labu[i], obj_base_repl_shift_sigma[i]);
    }
    for (i in 1:NobjEffects) {
        obj_effect_replCI[i] = normal_rng(obj_effect[i], obj_effect_repl_shift_sigma[i]);
    }
    iaction_labu_replCI = csr_matrix_times_vector(Niactions, Nobjects, iactXobjbase_w, iaction2obj, iactXobjbase_u, obj_base_labu_replCI) +
                          csr_matrix_times_vector(Niactions, NobjEffects, iactXobjeff_w, iactXobjeff_v, iactXobjeff_u, obj_effect_replCI);
}
