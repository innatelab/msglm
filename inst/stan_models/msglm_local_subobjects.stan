functions {
    real intensity_log_std(real z, real scaleHi, real scaleLo, real offset, real bend, real smooth) {
        return 0.5*(scaleHi+scaleLo)*(z-bend) + 0.5*(scaleHi-scaleLo)*sqrt((z-bend)*(z-bend)+smooth) + offset;
    }

    // reimplementation of R contr.poly::make.poly()
    matrix contr_poly(int n) {
        vector[n] scores;
        matrix[n, n] x;
        matrix[n, n] r;

        for (i in 1:n) scores[i] = i;
        scores -= mean(scores);

        for (i in 1:n) {
            for (j in 1:n) {
                x[i, j] = scores[i]^(j-1);
            }
        }
        r = qr_Q(x) * diag_matrix(diagonal(qr_R(x)));
        r ./= rep_matrix(sqrt(columns_dot_self(r)), n);
        return block(r, 1, 2, n, n-1);
    }

    // count nonzero elements in effXeff0 matrix
    // (block-diagonal matrix with block generated by contr_poly())
    int effXeff0_Nw(int ngroups, int[] effect2group) {
        int neffs[ngroups];
        int nw;
        for (i in 1:ngroups) neffs[i] = 0;
        for (i in 1:num_elements(effect2group)) neffs[effect2group[i]] += 1;
        nw = 0;
        for (i in 1:ngroups) nw += (neffs[i]-1)*neffs[i];
        return nw;
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
  int<lower=0> NbatchEffects;   // number of batch effects (that define assay experimental variation, but not biology)
  int<lower=0> NunderdefObjs;   // number of virtual interactions (the ones not detected but required for comparison)
  int<lower=1,upper=Nobjects> suo2obj[Nsubobjects];
  int<lower=1,upper=Nobjects> iaction2obj[Niactions];
  int<lower=1,upper=Nobjects> underdef_objs[NunderdefObjs];

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
  int<lower=0,upper=1> effect_is_positive[Neffects];
  vector[Neffects] effect_mean;

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
  int<lower=1,upper=Nmsprotocols*Nsubobjects> quant2msprotoXsuo[((Nmsprotocols > 1) && (Nsubobjects > 0)) ? Nquanted : 0];
  int<lower=1,upper=Niactions> miss2iaction[Nmissed];
  int<lower=1,upper=Nexperiments> miss2experiment[Nmissed];
  int<lower=1,upper=Nmsprotocols*Nsubobjects> miss2msprotoXsuo[((Nmsprotocols > 1) && (Nsubobjects > 0)) ? Nmissed : 0];
  int<lower=0,upper=NobjEffects> NobjEffectsPos;
  int<lower=0,upper=NobjEffects> NobjEffectsOther;
  int<lower=1,upper=NobjEffects> obj_effect_reshuffle[NobjEffects];
  vector<lower=0>[NobjEffects] obj_effect_tau;
  vector[NobjEffects] obj_effect_mean;

  int<lower=0,upper=NobjBatchEffects> NobjBatchEffectsPos;
  int<lower=0,upper=NobjBatchEffects> NobjBatchEffectsOther;
  int<lower=1,upper=NobjBatchEffects> obj_batch_effect_reshuffle[NobjBatchEffects];

  vector[iactXobjeff_Nw] iactXobjeff4sigma_w;

  vector[Niactions] iactXobjbase_w;
  int<lower=0> iactXobjbase_u[Niactions + 1];

  vector[Nobservations] obsXiact_w;
  int<lower=0> obsXiact_u[Nobservations + 1];

  int<lower=0> Nobservations0;  // number of observations degrees of freedom ()

  int<lower=0> obsXobs0_Nw;
  vector[effXeff0_Nw(Niactions, observation2iaction)] obsXobs_shift0_w;
  int<lower=1, upper=effXeff0_Nw(Niactions, observation2iaction) + 1> obsXobs_shift0_u[Nobservations + 1];
  int<lower=1, upper=Nobservations - Niactions> obsXobs_shift0_v[effXeff0_Nw(Niactions, observation2iaction)];

  int<lower=0> suoXsuo0_Nw;
  vector[effXeff0_Nw(Nobjects, suo2obj)] suoXsuo_shift0_w;
  int<lower=0, upper=effXeff0_Nw(Nobjects, suo2obj) + 1> suoXsuo_shift0_u[Nsubobjects + 1];
  int<lower=0, upper=Nsubobjects - Nobjects> suoXsuo_shift0_v[effXeff0_Nw(Nobjects, suo2obj)];

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
        cur_pos_eff += 1;
        obj_effect_reshuffle[i] = cur_pos_eff;
      } else {
        cur_other_eff += 1;
        obj_effect_reshuffle[i] = cur_other_eff;
      }
    }
  }
  obj_effect_tau = effect_tau[obj_effect2effect];
  obj_effect_mean = effect_mean[obj_effect2effect];

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
        cur_pos_eff += 1;
        obj_batch_effect_reshuffle[i] = cur_pos_eff;
      } else {
        cur_other_eff += 1;
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
      qLogStd[i] -= global_labu_shift; // obs_labu is modeled without obj_base
    }
  }
  quant2experiment = observation2experiment[quant2observation];
  quant2iaction = observation2iaction[quant2observation];
  miss2experiment = observation2experiment[miss2observation];
  miss2iaction = observation2iaction[miss2observation];

  // prepare obsXobs_shift0
  obsXobs0_Nw = effXeff0_Nw(Niactions, observation2iaction);
  {
    int iaction2nobs[Niactions];

    int iaction2nobs_2ndpass[Niactions];
    int iaction2obs_shift0_offset[Niactions];

    int obsXobs_shift0_offset;
    int obs_shift0_offset;

    iaction2nobs = rep_array(0, Niactions);
    for (i in 1:Nobservations) {
      int iact_ix;
      iact_ix = observation2iaction[i];
      iaction2nobs[iact_ix] += 1;
    }
    //print("iaction2nobs=", iaction2nobs);
    obs_shift0_offset = 0;
    obsXobs_shift0_offset = 0;
    obsXobs_shift0_u[1] = 1;
    iaction2nobs_2ndpass = rep_array(0, Niactions);
    iaction2obs_shift0_offset = rep_array(0, Niactions);
    for (i in 1:Nobservations) {
        int iact_ix;
        int iact_nobs;
        iact_ix = observation2iaction[i];
        //print("iact_ix[", i, "]=", iact_ix);
        iact_nobs = iaction2nobs[iact_ix];
        if (iact_nobs > 1) {
            // (re)generate contr_poly for interaction FIXME pre-build contr_poly for 2..max_nobs
            matrix[iact_nobs, iact_nobs-1] iact_obsXobs0 = contr_poly(iact_nobs);

            if (iaction2nobs_2ndpass[iact_ix] == 0) {
                // reserve (nobs-1) obs_shift0 variables
                iaction2obs_shift0_offset[iact_ix] = obs_shift0_offset;
                obs_shift0_offset += iact_nobs-1;
                //print("obs_shift0_offset=", obs_shift0_offset);
            }
            iaction2nobs_2ndpass[iact_ix] += 1;
            // add 2npass-th row of iact_obsXobs0 to the obsXobs0
            for (j in 1:cols(iact_obsXobs0)) {
                obsXobs_shift0_v[obsXobs_shift0_offset + j] = iaction2obs_shift0_offset[iact_ix] + j;
                obsXobs_shift0_w[obsXobs_shift0_offset + j] = iact_obsXobs0[iaction2nobs_2ndpass[iact_ix], j];
            }
            obsXobs_shift0_u[i+1] = obsXobs_shift0_u[i] + cols(iact_obsXobs0);
            obsXobs_shift0_offset += cols(iact_obsXobs0);
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

  suoXsuo0_Nw = effXeff0_Nw(Nobjects, suo2obj);
  suoXsuo_shift0_u = rep_array(0, Nsubobjects+1);
  if (Nsubobjects > 0) {
    int obj2nsuo[Nobjects];

    int obj2nsuo_2ndpass[Nobjects];
    int obj2suo_shift0_offset[Nobjects];

    int suoXsuo_shift0_offset;
    int suo_shift0_offset;

    obj2nsuo = rep_array(0, Nobjects);
    for (i in 1:Nsubobjects) {
      int obj_ix;
      obj_ix = suo2obj[i];
      obj2nsuo[obj_ix] += 1;
    }
    //print("object2nsuo=", object2nsuo);
    suo_shift0_offset = 0;
    suoXsuo_shift0_offset = 0;
    suoXsuo_shift0_u[1] = 1;
    obj2nsuo_2ndpass = rep_array(0, Nobjects);
    obj2suo_shift0_offset = rep_array(0, Nobjects);
    for (i in 1:Nsubobjects) {
        int obj_ix;
        int obj_nsuo;
        obj_ix = suo2obj[i];
        //print("obj_ix[", i, "]=", obj_ix);
        obj_nsuo = obj2nsuo[obj_ix];
        if (obj_nsuo > 1) {
            // (re)generate contr_poly for interaction FIXME pre-build contr_poly for 2..max_nsuo
            matrix[obj_nsuo, obj_nsuo-1] obj_suoXsuo0 = contr_poly(obj_nsuo);

            if (obj2nsuo_2ndpass[obj_ix] == 0) {
                // reserve (nsuo-1) suo_shift0 variables
                obj2suo_shift0_offset[obj_ix] = suo_shift0_offset;
                suo_shift0_offset += obj_nsuo-1;
                //print("suo_shift0_offset=", suo_shift0_offset);
            }
            obj2nsuo_2ndpass[obj_ix] += 1;
            // add 2npass-th row of obj_suoXsuo0 to the suoXsuo0
            for (j in 1:cols(obj_suoXsuo0)) {
                suoXsuo_shift0_v[suoXsuo_shift0_offset + j] = obj2suo_shift0_offset[obj_ix] + j;
                suoXsuo_shift0_w[suoXsuo_shift0_offset + j] = obj_suoXsuo0[obj2nsuo_2ndpass[obj_ix], j];
            }
            suoXsuo_shift0_u[i+1] = suoXsuo_shift0_u[i] + cols(obj_suoXsuo0);
            suoXsuo_shift0_offset += cols(obj_suoXsuo0);
        }
    }

    if (Nmsprotocols > 1) {
      // all references to the 1st protocol are redirected to index 1 (this shift is always 0)
      for (i in 1:Nquanted) {
        int msproto;
        msproto = experiment2msproto[quant2experiment[i]];
        quant2msprotoXsuo[i] = msproto > 1 ? (msproto-1)*Nsubobjects + quant2suo[i] + 1 : 1;
      }
      for (i in 1:Nmissed) {
        int msproto;
        msproto = experiment2msproto[quant2experiment[i]];
        miss2msprotoXsuo[i] = msproto > 1 ? (msproto-1)*Nsubobjects + miss2suo[i] + 1 : 1;
      }
    }
  }
}

parameters {
  //real obj_base;
  //real<lower=-20, upper=-2> underdef_obj_shift;
  //real<lower=0, upper=5.0> obj_shift_sigma;
  //vector<lower=0>[Nconditions] condition_repl_effect_sigma;

  vector[Nobjects] obj_base_labu0; // baseline object abundance without underdefinedness adjustment
  vector<lower=0.01>[Nobservations0 > 0 ? Nobjects : 0] obj_base_repl_shift_sigma;

  real<lower=0.0> suo_shift_sigma;
  vector[Nsubobjects > 0 ? Nsubobjects-Nobjects : 0] suo_shift0_unscaled; // subobject shift within object
  real<lower=0.0> suo_msproto_shift_sigma;
  vector[Nsubobjects > 0 && Nmsprotocols > 1 ? Nsubobjects * (Nmsprotocols-1) : 0] suo_msproto_shift_unscaled;

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

  vector[Nsubobjects] suo_shift_unscaled; // subcomponent shift within object

  // correct baseline abundances of underdefined objects
  obj_base_labu = obj_base_labu0;
  for (i in 1:NunderdefObjs) {
    obj_base_labu[underdef_objs[i]] = obj_base_labu[underdef_objs[i]] + underdef_obj_shift;
  }

  // calculate effects lambdas and scale effects
  obj_effect_sigma = obj_effect_lambda_a ./ sqrt(obj_effect_lambda_t) .* obj_effect_tau;
  obj_effect = obj_effect_mean + append_row(obj_effect_unscaled_pos, obj_effect_unscaled_other)[obj_effect_reshuffle] .* obj_effect_sigma;

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
  // calculate suo_labu_shift
  if (Nsubobjects > 1) {
    suo_shift_unscaled = csr_matrix_times_vector(Nsubobjects, Nsubobjects - Nobjects, suoXsuo_shift0_w, suoXsuo_shift0_v, suoXsuo_shift0_u, suo_shift0_unscaled);
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
    if (Nsubobjects > 0) {
      suo_shift_sigma ~ inv_gamma(1.0, 1.0);
      suo_shift_unscaled ~ normal(0.0, 1.0);
      suo_msproto_shift_sigma ~ inv_gamma(2.0, 1.0); // if Nmsprotocols==1, fake something far from 0 and +Inf
      if (Nmsprotocols > 1) {
        to_vector(suo_msproto_shift_unscaled) ~ normal(0.0, 1.0);
      }
    }

    // calculate the likelihood
    {
        vector[Nquanted] q_labu;
        vector[Nmissed] m_labu;

        q_labu = obs_labu[quant2observation] + experiment_shift[quant2experiment];
        m_labu = obs_labu[miss2observation] + experiment_shift[miss2experiment];
        //qLogAbu = iaction_shift[quant2iaction] + experiment_shift[quant2experiment];
        //mLogAbu = iaction_shift[miss2iaction] + experiment_shift[miss2experiment];
        if (Nsubobjects > 0) {
            // adjust by subcomponent shift
            vector[Nsubobjects] suo_shift;
            suo_shift = suo_shift_unscaled * suo_shift_sigma;
            q_labu += suo_shift[quant2suo];
            m_labu += suo_shift[miss2suo];

            if (Nmsprotocols > 1) {
                vector[(Nmsprotocols-1)*Nsubobjects+1] suo_msproto_shift;
                suo_msproto_shift[1] = 0.0; // the 1st protocol has no shift (taken care by suo_shift)
                suo_msproto_shift[2:((Nmsprotocols-1)*Nsubobjects+1)] = suo_msproto_shift_unscaled * suo_msproto_shift_sigma;

                q_labu += suo_msproto_shift[quant2msprotoXsuo];
                m_labu += suo_msproto_shift[miss2msprotoXsuo];
            }
        }
        if (NbatchEffects > 0) {
          q_labu = q_labu + obs_batch_shift[quant2observation];
          m_labu = m_labu + obs_batch_shift[miss2observation];
        }

        // model quantitations and missing data
        qDataNorm ~ double_exponential(exp(q_labu - qLogStd), 1);
        1 ~ bernoulli_logit(q_labu * (zScale * zDetectionFactor) + (-mzShift * zScale * zDetectionFactor + zDetectionIntercept));
        0 ~ bernoulli_logit(m_labu * (zScale * zDetectionFactor) + (-mzShift * zScale * zDetectionFactor + zDetectionIntercept));
    }
}

generated quantities {
    vector[Nobjects] obj_base_labu_replCI;
    vector[NobjEffects] obj_effect_replCI;
    vector[Niactions] iaction_labu_replCI;
    vector[Nsubobjects] suo_llh;

    for (i in 1:Nobjects) {
        obj_base_labu_replCI[i] = normal_rng(obj_base_labu[i], obj_base_repl_shift_sigma[i]);
    }
    for (i in 1:NobjEffects) {
        obj_effect_replCI[i] = normal_rng(obj_effect[i], obj_effect_repl_shift_sigma[i]);
    }
    iaction_labu_replCI = csr_matrix_times_vector(Niactions, Nobjects, iactXobjbase_w, iaction2obj, iactXobjbase_u, obj_base_labu_replCI) +
                          csr_matrix_times_vector(Niactions, NobjEffects, iactXobjeff_w, iactXobjeff_v, iactXobjeff_u, obj_effect_replCI);
    // per-subobject loglikelihood (the code copied from "model" section)
    if (Nsubobjects > 0) {
        vector[Nquanted] q_labu;
        vector[Nmissed] m_labu;
        vector[Nsubobjects] suo_shift;

        suo_shift = suo_shift_unscaled * suo_shift_sigma;
        // prepare predicted abundances
        q_labu = obs_labu[quant2observation] + experiment_shift[quant2experiment] + suo_shift[quant2suo];
        m_labu = obs_labu[miss2observation] + experiment_shift[miss2experiment] + suo_shift[miss2suo];

        if (Nmsprotocols > 1) {
            vector[(Nmsprotocols-1)*Nsubobjects+1] suo_msproto_shift;
            suo_msproto_shift[1] = 0.0; // the 1st protocol has no shift (taken care by suo_shift)
            suo_msproto_shift[2:((Nmsprotocols-1)*Nsubobjects+1)] = suo_msproto_shift_unscaled * suo_msproto_shift_sigma;

            q_labu += suo_msproto_shift[quant2msprotoXsuo];
            m_labu += suo_msproto_shift[miss2msprotoXsuo];
        }
        if (NbatchEffects > 0) {
          q_labu = q_labu + obs_batch_shift[quant2observation];
          m_labu = m_labu + obs_batch_shift[miss2observation];
        }

        // calculate log-likelihood per subobject
        suo_llh = rep_vector(0.0, Nsubobjects);
        for (i in 1:Nquanted) {
          suo_llh[quant2suo[i]] += double_exponential_lpdf(qDataNorm[i] | exp(q_labu[i] - qLogStd[i]), 1) +
              bernoulli_logit_lpmf(1 | q_labu[i] * (zScale * zDetectionFactor) + (-mzShift * zScale * zDetectionFactor + zDetectionIntercept));
        }
        for (i in 1:Nmissed) {
          suo_llh[miss2suo[i]] += bernoulli_logit_lpmf(0 | m_labu[i] * (zScale * zDetectionFactor) + (-mzShift * zScale * zDetectionFactor + zDetectionIntercept));
        }
    }
}
