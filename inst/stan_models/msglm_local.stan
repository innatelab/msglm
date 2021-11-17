functions {
#include include/array_utils.stan
#include include/contrast_matrices.stan
#include include/msnoise.stan
}

data {
  int<lower=1> Nconditions;     // number of experimental conditions
  int<lower=0> Nobjects;        // number of objects (proteins/peptides/sites etc)

  int<lower=0> NobjConditions;  // number of object-in-condition pairs
  int<lower=1,upper=Nobjects> obj_cond2obj[NobjConditions];

  int<lower=1> Nprobes;         // number of MS probes (= MS channels)
  vector[Nprobes] mschannel_shift;

  int<lower=0> NobjProbes;      // number of object-in-msprobe pairs with MS data
  int<lower=1,upper=Nprobes> obj_probe2probe[NobjProbes];
  int<lower=1,upper=NobjConditions> obj_probe2obj_cond[NobjProbes];

  // map from quantitations/missings to object-in-msprobe
  int<lower=0> Nquanted;        // total number of quantified object-in-msprobes
  int<lower=1,upper=NobjProbes>  quant2obj_probe[Nquanted];
  int<lower=0,upper=1> quant_isreliable[Nquanted];
  int<lower=0> Nmissed;         // total number of missed object-in-msprobes
  int<lower=1,upper=NobjProbes> miss2obj_probe[Nmissed];
  vector<lower=0>[Nquanted] qData; // quanted data
  vector<lower=0, upper=1>[Nmissed] missing_sigmoid_scale; // sigmoid scales for indiv. object-in-msprobes (<1 for higher uncertainty)

  // linear model specification
  int<lower=0> Neffects;        // number of effects (that define conditions)
  int<lower=0> NobjEffects;
  int<lower=1,upper=Neffects> obj_effect2effect[NobjEffects];
  int<lower=0,upper=1> effect_is_positive[Neffects];
  vector[Neffects] effect_mean;

  int<lower=0> NbatchEffects;   // number of batch effects (that define assay experimental variation, but not biology)
  int<lower=0> NobjBatchEffects;
  int<lower=1,upper=NbatchEffects> obj_batch_effect2batch_effect[NobjBatchEffects];
  int<lower=0,upper=1> batch_effect_is_positive[NbatchEffects];

  // obj_condXeff (object-in-condition X object-effect) sparse matrix
  int<lower=0> obj_condXeff_Nw;
  vector[obj_condXeff_Nw] obj_condXeff_w;
  int<lower=1, upper=obj_condXeff_Nw+1> obj_condXeff_u[NobjConditions+1];
  int<lower=1, upper=NobjEffects> obj_condXeff_v[obj_condXeff_Nw];

  // obj_probeXeff (object-in-msprobe X object-effect) sparse matrix
  int<lower=0> obj_probeXeff_Nw;
  vector[obj_probeXeff_Nw] obj_probeXeff_w;
  int<lower=1, upper=obj_probeXeff_Nw+1> obj_probeXeff_u[NobjProbes+1];
  int<lower=1, upper=NobjEffects> obj_probeXeff_v[obj_probeXeff_Nw];

  // obj_probeXbatcheff (object-in-msprobe X object-batch_effect) sparse matrix
  int<lower=0> obj_probeXbatcheff_Nw;
  vector[obj_probeXbatcheff_Nw] obj_probeXbatcheff_w;
  int<lower=1, upper=obj_probeXbatcheff_Nw+1> obj_probeXbatcheff_u[NbatchEffects > 0 ? NobjProbes + 1 : 1];
  int<lower=1, upper=NobjBatchEffects> obj_probeXbatcheff_v[obj_probeXbatcheff_Nw];

  // global model constants
  real obj_labu_min; // minimal average abundance of an object
  real<lower=0> obj_labu_min_scale; // scale that defines the softness of lower abundance limit

  real<lower=0.0> hsprior_lambda_a_offset; // prevent lambda_a/t being too close to zero, because that negatively affects MCMC convergence
  real<lower=0.0> hsprior_lambda_t_offset;
  vector<lower=0>[Neffects] effect_tau;
  vector<lower=0>[Neffects] effect_df;
  vector<lower=0>[Neffects] effect_df2;
  real<lower=0> effect_slab_df;
  real<lower=0> effect_slab_scale;

  real<lower=0> obj_base_labu_sigma; // sigma of average abundance distribution
  real<lower=0> obj_probe_shift_tau; // tau for the object-in-msprobe shift (relative to object-in-condition) prior
  real<lower=0> obj_probe_shift_df;  // degrees of freedom for object-in-msprobe shift
  real<lower=0> batch_effect_sigma;

  // instrument calibrated parameters (FIXME: msproto-dependent)
  real<lower=0> zDetectionFactor;
  real zDetectionIntercept;

  real<lower=0> sigmaScaleHi;
  real<lower=0> sigmaScaleLo;
  real sigmaOffset;
  real sigmaBend;
  real sigmaSmooth;
  real outlierProb;

  real zShift;
  real zScale;
}

transformed data {
  real compress_a = cauchy_compress_a(outlierProb);
  vector[Nquanted] zScore = (log2(qData) - zShift) * zScale;
  vector[Nquanted] qLog2Std; // log2(sd(qData))-zShift
  vector[Nquanted] qLogShift;
  vector<lower=0>[Nquanted] qDataNorm; // qData/sd(qData)
  int<lower=0,upper=Nquanted> NreliableQuants = sum(quant_isreliable);
  int<lower=1,upper=Nquanted> reliable_quants[NreliableQuants];

  int<lower=1> Nmschannels = Nprobes;
  int<lower=1,upper=NobjConditions> quant2obj_cond[Nquanted] = obj_probe2obj_cond[quant2obj_probe];
  int<lower=1,upper=Nmschannels> quant2mschannel[Nquanted] = obj_probe2probe[quant2obj_probe];
  int<lower=1,upper=NobjConditions> miss2obj_cond[Nmissed] = obj_probe2obj_cond[miss2obj_probe];
  int<lower=1,upper=Nmschannels> miss2mschannel[Nmissed] = obj_probe2probe[miss2obj_probe];

  int<lower=0,upper=NobjEffects> NobjEffectsPos = sum(effect_is_positive[obj_effect2effect]);
  int<lower=0,upper=NobjEffects> NobjEffectsOther = NobjEffects - NobjEffectsPos;
  int<lower=1,upper=NobjEffects> obj_effect_reshuffle[NobjEffects] =
      nonzeros_first(effect_is_positive[obj_effect2effect]);
  vector<lower=0>[NobjEffects] obj_effect_tau = effect_tau[obj_effect2effect];
  vector[NobjEffects] obj_effect_mean = effect_mean[obj_effect2effect];
  vector<lower=0>[NobjEffects] obj_effect_df = effect_df[obj_effect2effect];
  vector<lower=0>[NobjEffects] obj_effect_df2 = effect_df2[obj_effect2effect];

  int<lower=0,upper=NobjBatchEffects> NobjBatchEffectsPos = sum(batch_effect_is_positive[obj_batch_effect2batch_effect]);
  int<lower=0,upper=NobjBatchEffects> NobjBatchEffectsOther = NobjBatchEffects - NobjBatchEffectsPos;
  int<lower=1,upper=NobjBatchEffects> obj_batch_effect_reshuffle[NobjBatchEffects] =
      nonzeros_first(batch_effect_is_positive[obj_batch_effect2batch_effect]);

  vector[NobjConditions] obj_condXbase_w = rep_vector(1.0, NobjConditions);
  int<lower=0> obj_condXbase_u[NobjConditions + 1] = one_to(NobjConditions + 1);

  vector[NobjProbes] obj_probeXcond_w = rep_vector(1.0, NobjProbes);
  int<lower=0> obj_probeXcond_u[NobjProbes + 1] = one_to(NobjProbes + 1);

  vector[NobjProbes] obj_probeXbase_w = rep_vector(1.0, NobjProbes);
  int<lower=1> obj_probeXbase_u[NobjProbes + 1] = one_to(NobjProbes + 1);
  int<lower=1, upper=Nobjects> obj_probe2obj[NobjProbes] = obj_cond2obj[obj_probe2obj_cond];

  int<lower=0> NobjConditionsObserved = ndistinct(obj_probe2obj_cond, NobjConditions);
  int<lower=0> NobjProbes0 = NobjProbes - NobjConditionsObserved; // number of object-in-msprobe degrees of freedom ()
  int<lower=0> obj_probe_shiftXshift0_Nw = NobjProbes0 > 0 ? contr_poly_Nw(NobjConditions, obj_probe2obj_cond) : 0;
  vector[obj_probe_shiftXshift0_Nw] obj_probe_shiftXshift0_w;
  int<lower=1, upper=obj_probe_shiftXshift0_Nw + 1> obj_probe_shiftXshift0_u[NobjProbes0 > 0 ? NobjProbes + 1 : 0];
  int<lower=1, upper=NobjProbes0> obj_probe_shiftXshift0_v[NobjProbes0 > 0 ? obj_probe_shiftXshift0_Nw : 0];

  matrix[Nobjects + NobjEffects, NobjConditions] obj_baseffXcond; // OLS operator converting obj_cond to obj_effect and obj_base
  //matrix[NobjConditions, Nobjects] obj_condXbase = csr_to_dense_matrix(NobjConditions, Nobjects, obj_condXbase_w, obj_cond2obj, obj_condXbase_u);
  //matrix[NobjConditions, NobjEffects] obj_condXeff = csr_to_dense_matrix(NobjConditions, NobjEffects, obj_condXeff_w, obj_condXeff_v, obj_condXeff_u);

  // process the intensity data to optimize likelihood calculation
  {
    for (i in 1:Nquanted) {
      qLog2Std[i] = intensity_log2_std(zScore[i], sigmaScaleHi, sigmaScaleLo, sigmaOffset, sigmaBend, sigmaSmooth);
      qDataNorm[i] = exp2(log2(qData[i]) - qLog2Std[i]);
      qLog2Std[i] -= zShift; // obj_probe_labu is normalized to zShift
      qLogShift[i] = -qLog2Std[i] * log(2);
    }
  }

  // collect indices of reliable quantifications
  {
    int j = 1;
    for (i in 1:Nquanted) {
      if (quant_isreliable[i]) {
        reliable_quants[j] = i;
        j += 1;
      }
    }
  }

  // prepare obj_probe_shiftXshift0
  if (NobjProbes0 > 0) {
    matrix[NobjProbes, NobjProbes0] obj_probe_shiftXshift0 = block_contr_poly(NobjConditions, obj_probe2obj_cond, positive_infinity());
    obj_probe_shiftXshift0_w = csr_extract_w_0(obj_probe_shiftXshift0, positive_infinity());
    obj_probe_shiftXshift0_u = csr_extract_u(obj_probe_shiftXshift0);
    obj_probe_shiftXshift0_v = csr_extract_v(obj_probe_shiftXshift0);
    //print("obj_probe_shiftXshift0=", obs_shiftXobs_shift0);
  }

  {
    matrix[NobjConditions, Nobjects+NobjEffects] obj_condXbaseff;
    // = append_col(obj_condXbase, obj_condXeff);
    // wierd way of convertion
    for (i in 1:Nobjects) {
      vector[Nobjects] obj = rep_vector(0.0, Nobjects);
      vector[NobjConditions] objcond;
      obj[i] = 1.0;
      objcond = csr_matrix_times_vector(NobjConditions, Nobjects, obj_condXbase_w, obj_cond2obj, obj_condXbase_u, obj);
      for (j in 1:NobjConditions) {
        obj_condXbaseff[j, i] = objcond[j];
      }
    }
    for (i in 1:NobjEffects) {
      vector[NobjEffects] objeff = rep_vector(0.0, NobjEffects);
      vector[NobjConditions] objcond;
      objeff[i] = 1.0;
      objcond = csr_matrix_times_vector(NobjConditions, NobjEffects, obj_condXeff_w, obj_condXeff_v, obj_condXeff_u, objeff);
      for (j in 1:NobjConditions) {
        obj_condXbaseff[j, i + Nobjects] = objcond[j];
      }
    }
    // OLS operator
    obj_baseffXcond = crossprod(obj_condXbaseff)\(obj_condXbaseff');
    //print("obj_baseffXcond=", obj_baseffXcond);
  }
}

parameters {
  //real obj_base;
  //real<lower=0, upper=5.0> obj_shift_sigma;
  //vector<lower=0>[Nconditions] condition_repl_effect_sigma;

  vector[Nobjects] obj_base_labu0; // baseline object abundance

  //real<lower=0.0> obj_effect_tau;
  real<lower=0.0> effect_slab_c_t;
  vector<lower=0.0>[NobjEffects] obj_effect_lambda_t;
  vector<lower=0.0>[NobjEffects] obj_effect_lambda_a;
  vector<lower=0.0>[NobjEffects] obj_effect_eta_t;
  vector<lower=0.0>[NobjEffects] obj_effect_eta_a;
  vector<lower=0.0>[NobjEffectsPos] obj_effect_unscaled_pos;
  vector[NobjEffectsOther] obj_effect_unscaled_other;

  //real<lower=0> obj_repl_effect_sigma;
  //vector<lower=0>[Nobjects*Nmschannels] repl_shift_lambda;
  vector<lower=0>[NobjProbes0 > 0 ? NobjConditions : 0] obj_probe_shift_lambda_t;
  vector<lower=0>[NobjProbes0 > 0 ? NobjConditions : 0] obj_probe_shift_lambda_a;
  vector[NobjProbes0] obj_probe_shift0;

  //real<lower=0> obj_batch_effect_sigma;
  vector<lower=0.0>[NobjBatchEffectsPos] obj_batch_effect_unscaled_pos;
  vector[NobjBatchEffectsOther] obj_batch_effect_unscaled_other;
}

transformed parameters {
  vector[Nobjects] obj_base_labu = obj_base_labu0 * obj_base_labu_sigma;
  vector[NobjEffects] obj_effect;
  real<lower=0> effect_slab_c;
  vector<lower=0>[NobjEffects] obj_effect_sigma; // AKA lambda_tilde*tau in rstanarm
  vector[NobjBatchEffects] obj_batch_effect;
  //vector<lower=0>[NobjBatchEffects] obj_batch_effect_sigma;

  vector<lower=0>[NobjProbes0 > 0 ? NobjConditions : 0] obj_probe_shift_sigma;

  vector[NobjProbes] obj_probe_labu; // obj_cond_labu + obj_probe_shift * obj_probe_shift_sigma
  vector[NobjProbes0 > 0 ? NobjProbes : 0] obj_probe_shift; // replicate shifts for all potential object-in-msprobes (including missing)
  vector[NobjBatchEffects > 0 ? NobjProbes : 0] obj_probe_batch_shift;

  // calculate effects lambdas and scale effects
  {
      vector[NobjEffects] obj_effect_sigma_pre; // AKA lambda_eta2 in rstanarm
      real effect_slab_c2 = square(effect_slab_scale) * effect_slab_c_t; // AKA c2 in rstanarm
      effect_slab_c = effect_slab_scale * sqrt(effect_slab_c_t);
      obj_effect_sigma_pre = square(obj_effect_lambda_a .* obj_effect_eta_a) .* obj_effect_lambda_t .* obj_effect_eta_t;
      obj_effect_sigma = sqrt(effect_slab_c2 * obj_effect_sigma_pre ./ (effect_slab_c2 + square(obj_effect_tau) .* obj_effect_sigma_pre)) .* obj_effect_tau;
  }
  obj_effect = obj_effect_mean + append_row(obj_effect_unscaled_pos, obj_effect_unscaled_other)[obj_effect_reshuffle] .* obj_effect_sigma;

  // calculate object-in-msprobe log abundance
  obj_probe_labu = csr_matrix_times_vector(NobjProbes, Nobjects, obj_probeXbase_w, obj_probe2obj, obj_probeXbase_u, obj_base_labu) +
             csr_matrix_times_vector(NobjProbes, NobjEffects, obj_probeXeff_w, obj_probeXeff_v, obj_probeXeff_u, obj_effect);

  // calculate obj_probe_shift and obj_probe_labu
  if (NobjProbes0 > 0) {
    obj_probe_shift_sigma = obj_probe_shift_lambda_a .* sqrt(obj_probe_shift_lambda_t) * obj_probe_shift_tau;
    obj_probe_shift = csr_matrix_times_vector(NobjProbes, NobjProbes0, obj_probe_shiftXshift0_w, obj_probe_shiftXshift0_v, obj_probe_shiftXshift0_u, obj_probe_shift0) .*
        csr_matrix_times_vector(NobjProbes, NobjConditions, obj_probeXcond_w, obj_probe2obj_cond, obj_probeXcond_u, obj_probe_shift_sigma);
    obj_probe_labu += obj_probe_shift;
  }
  // calculate obj_probe_batch_shift (doesn't make sense to add to obj_probe_labu)
  if (NbatchEffects > 0) {
    obj_batch_effect = append_row(obj_batch_effect_unscaled_pos, obj_batch_effect_unscaled_other)[obj_batch_effect_reshuffle] * batch_effect_sigma;
    obj_probe_batch_shift = csr_matrix_times_vector(NobjProbes, NobjBatchEffects, obj_probeXbatcheff_w, obj_probeXbatcheff_v, obj_probeXbatcheff_u, obj_batch_effect);
  }
}

model {
    // abundance distribution
    //obj_base ~ normal(zShift, 1.0);
    //obj_shift_sigma ~ inv_gamma(2.0, 0.33/zScale); // mode is 1/zScale
    obj_base_labu0 ~ std_normal();
    // treatment effect parameters, horseshoe prior
    //obj_effect_tau ~ student_t(2, 0.0, 1.0);
    //obj_effect_lambda ~ student_t(2, 0.0, obj_effect_tau);
    obj_effect_lambda_t - hsprior_lambda_t_offset ~ inv_gamma(0.5 * obj_effect_df, 0.5 * obj_effect_df);
    obj_effect_lambda_a - hsprior_lambda_a_offset ~ std_normal(); // 1.0 = 2/2
    obj_effect_eta_t - hsprior_lambda_t_offset ~ inv_gamma(0.5 * obj_effect_df2, 0.5 * obj_effect_df2);
    obj_effect_eta_a - hsprior_lambda_a_offset ~ std_normal(); // 1.0 = 2/2
    obj_effect_unscaled_pos ~ std_normal();
    obj_effect_unscaled_other ~ std_normal();
    effect_slab_c_t - hsprior_lambda_t_offset ~ inv_gamma(0.5 * effect_slab_df, 0.5 * effect_slab_df);
    // batch effect parameters, cauchy prior on sigma
    //condition_repl_effect_sigma ~ inv_gamma(1.5, 1.0);

    //repl_shift_lambda ~ student_t(2, 0.0, repl_shift_tau);
    //obj_repl_effect ~ normal(0.0, obj_repl_effect_lambda);
    if (NobjProbes0 > 0) {
      vector[NobjProbes] obj_probe_shift_unscaled; // unscaled object-in-msprobe shifts

      obj_probe_shift_lambda_t - hsprior_lambda_t_offset ~ inv_gamma(0.5 * obj_probe_shift_df, 0.5 * obj_probe_shift_df);
      obj_probe_shift_lambda_a - hsprior_lambda_a_offset ~ std_normal();

      //print("obj_probe_shift_sigma=", obj_probe_shift_sigma);
      //print("obj_probeXcond=", csr_to_dense_matrix(NobjProbes, NobjConditions,
      //      obj_probeXcond_w, obj_probe2obj_cond, obj_probeXcond_u));
      obj_probe_shift_unscaled = csr_matrix_times_vector(NobjProbes, NobjProbes0, obj_probe_shiftXshift0_w, obj_probe_shiftXshift0_v, obj_probe_shiftXshift0_u, obj_probe_shift0);
      obj_probe_shift_unscaled ~ std_normal();
    }
    //to_vector(repl_shift) ~ normal(0.0, repl_shift_lambda);

    //obj_batch_effect_lambda ~ student_t(2, 0.0, obj_batch_effect_tau);
    if (NbatchEffects > 0) {
      obj_batch_effect_unscaled_pos ~ std_normal();
      obj_batch_effect_unscaled_other ~ std_normal();
    }

    // soft lower limit of protein abundance for each object-in-msprobe
    1 ~ bernoulli_logit((obj_probe_labu - obj_labu_min) * obj_labu_min_scale);

    // calculate the likelihood
    {
        vector[Nquanted] q_labu;
        vector[Nmissed] m_labu;

        q_labu = obj_probe_labu[quant2obj_probe] + mschannel_shift[quant2mschannel];
        m_labu = obj_probe_labu[miss2obj_probe] + mschannel_shift[miss2mschannel];
        if (NbatchEffects > 0) {
          q_labu += obj_probe_batch_shift[quant2obj_probe];
          m_labu += obj_probe_batch_shift[miss2obj_probe];
        }

        // model quantitations and missing data
        //{ // ~10% slower version with explicit normal/cauchy mixture
        //  vector[Nquanted] delta = exp2(q_labu - qLog2Std) - qDataNorm;
        //  for (i in 1:Nquanted) {
        //    target += log_mix(outlierProb, cauchy_lpdf(delta[i] | 0, 1), std_normal_lpdf(delta[i]));
        //  }
        //}
        // 10% faster version with cauchy_compressv() transform
        cauchy_compressv(exp2(q_labu - qLog2Std) - qDataNorm, compress_a, 4.0) ~ std_normal();
        // soft-lower-limit for object intensities of reliable quantifications
        1 ~ bernoulli_logit(q_labu[reliable_quants] * (zScale * zDetectionFactor) + zDetectionIntercept);
        0 ~ bernoulli_logit(missing_sigmoid_scale .* (m_labu * (zScale * zDetectionFactor) + zDetectionIntercept));
    }
}

generated quantities {
    // calculate object-in-condition log abundance
    vector[NobjConditions] obj_cond_labu =
        csr_matrix_times_vector(NobjConditions, Nobjects, obj_condXbase_w, obj_cond2obj, obj_condXbase_u, obj_base_labu) +
        csr_matrix_times_vector(NobjConditions, NobjEffects, obj_condXeff_w, obj_condXeff_v, obj_condXeff_u, obj_effect);
        //obj_base_labu[obj_cond2obj] + obj_condXeff * obj_effect;

    // object-in-condition log-abundance with object-in-msprobe variation (replCI)
    vector[NobjConditions] obj_cond_labu_replCI = NobjProbes0 > 0 ? to_vector(normal_rng(obj_cond_labu, obj_probe_shift_sigma)) : obj_cond_labu;
    vector[Nobjects] obj_base_labu_replCI;
    vector[NobjEffects] obj_effect_replCI;
    {
      vector[Nobjects + NobjEffects] obj_baseff = obj_baseffXcond * obj_cond_labu_replCI;
      for (i in 1:Nobjects) {
        obj_base_labu_replCI[i] = obj_baseff[i];
      }
      for (i in 1:NobjEffects) {
        obj_effect_replCI[i] = obj_baseff[i+Nobjects];
      }
    }
}
