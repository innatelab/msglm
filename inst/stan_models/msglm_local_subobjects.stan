functions {
    int ndistinct(int[] ixs, int maxix) {
      int ix2used[maxix] = rep_array(0, maxix);
      for (ix in ixs) {
        ix2used[ix] = 1;
      }
      return sum(ix2used);
    }

    int[] one_to(int n) {
      int res[n];
      for (i in 1:n) res[i] = i;
      return res;
    }

    real intensity_log_std(real z, real scaleHi, real scaleLo, real offs, real bend, real smooth) {
        return 0.5*(scaleHi+scaleLo)*(z-bend) + 0.5*(scaleHi-scaleLo)*sqrt((z-bend)*(z-bend)+smooth) + offs;
    }

    // compresses x: x~0 -> logcompress(x)~x, abs(x)>>0 -> logcompress(x)~sign(x)*log(abs(x))
    real logcompress(real x, data real s) {
      return x * (1 + log1p(fabs(s*x))) / (1 + fabs(s*x));
    }

    vector logcompressv(vector x, data real s) {
      return x .* (1 + log1p(fabs(s*x))) ./ (1 + fabs(s*x));
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

    matrix contr_treatment(int n) {
        matrix[n, n-1] r;

        r = rep_matrix(0.0, n, n-1);
        for (i in 1:(n-1)) r[i+1, i] = 1;
        return r;
    }

    // count nonzero elements in effXeff0 matrix
    // (block-diagonal matrix with block generated by contr_poly())
    int contr_poly_Nw(int ngroups, int[] obj2group) {
        int nobjs[ngroups];
        int nw;
        for (i in 1:ngroups) nobjs[i] = 0;
        for (i in 1:num_elements(obj2group)) nobjs[obj2group[i]] += 1;
        nw = 0;
        for (i in 1:ngroups) nw += (nobjs[i]-1)*nobjs[i];
        return nw;
    }

    // (block-diagonal matrix with block generated by contr_treatment())
    int contr_treatment_Nw(int ngroups, int[] obj2group) {
        int nobjs[ngroups];
        int nw;
        for (i in 1:ngroups) nobjs[i] = 0;
        for (i in 1:num_elements(obj2group)) nobjs[obj2group[i]] += 1;
        nw = 0;
        for (i in 1:ngroups) nw += (nobjs[i]-1);
        return nw;
    }

    int[] objeffects_reshuffle(int[] objeff2eff, int[] iseffpositive) {
        int cur_pos_eff = 0;
        int cur_other_eff = sum(iseffpositive[objeff2eff]);
        int res[size(objeff2eff)];
        for (i in 1:size(objeff2eff)) {
            if (iseffpositive[objeff2eff[i]]) {
                cur_pos_eff += 1;
                res[i] = cur_pos_eff;
            } else {
                cur_other_eff += 1;
                res[i] = cur_other_eff;
            }
        }
        return res;
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
  int<lower=0> NunderdefObjs;   // number of virtual interactions (the ones not detected but required for comparison)
  int<lower=1,upper=Nobjects> suo2obj[Nsubobjects];
  int<lower=1,upper=Nobjects> iaction2obj[Niactions];
  int<lower=1,upper=Nobjects> underdef_objs[NunderdefObjs];

  vector[Nexperiments] experiment_shift;

  int<lower=1,upper=Nexperiments> observation2experiment[Nobservations];
  int<lower=1,upper=Niactions> observation2iaction[Nobservations];
  int<lower=0,upper=1> observation_reliable[Nobservations];
  int<lower=1,upper=Nmsprotocols> experiment2msproto[Nmsprotocols > 0 ? Nexperiments : 0];

  // map from labelXreplicateXobject to observed/missed data
  int<lower=0> Nquanted;        // total number of quantified subobjectsXexperiments
  int<lower=1,upper=Nobservations>  quant2observation[Nquanted];
  int<lower=1,upper=Nsubobjects> quant2suo[Nsubobjects > 0 ? Nquanted : 0];
  int<lower=0> Nmissed;         // total number of missed subobjectsXexperiments
  int<lower=1,upper=Nobservations> miss2observation[Nmissed];
  int<lower=1,upper=Nsubobjects> miss2suo[Nsubobjects > 0 ? Nmissed : 0];
  vector<lower=0>[Nquanted] qData; // quanted data
  vector<lower=0, upper=1>[Nmissed] missing_sigmoid_scale; // sigmoid scales for indiv. observations (<1 for higher uncertainty)

  // linear model specification
  int<lower=0> Neffects;        // number of effects (that define conditions)
  int<lower=0> NbatchEffects;   // number of batch effects (that define assay experimental variation, but not biology)
  int<lower=0> NsubBatchEffects;// number of subobject batch effects (that define assay experimental variation, but not biology)
  int<lower=0> NobjEffects;
  int<lower=1,upper=Neffects> obj_effect2effect[NobjEffects];
  int<lower=0,upper=1> effect_is_positive[Neffects];
  vector[Neffects] effect_mean;

  int<lower=0> NobjBatchEffects;
  int<lower=1,upper=NbatchEffects> obj_batch_effect2batch_effect[NobjBatchEffects];
  int<lower=0,upper=1> batch_effect_is_positive[NbatchEffects];

  int<lower=0> NsuoBatchEffects;
  int<lower=1,upper=NsubBatchEffects> suo_subbatch_effect2subbatch_effect[NsuoBatchEffects];
  int<lower=0,upper=1> subbatch_effect_is_positive[NsubBatchEffects];

  // iactXobjeff (interaction X object_effect) sparse matrix
  int<lower=0> iactXobjeff_Nw;
  vector[iactXobjeff_Nw] iactXobjeff_w;
  int<lower=0, upper=iactXobjeff_Nw+1> iactXobjeff_u[Niactions+1];
  int<lower=0, upper=NobjEffects> iactXobjeff_v[iactXobjeff_Nw];

  int<lower=0> obsXobjeff_Nw;
  vector[obsXobjeff_Nw] obsXobjeff_w;
  int<lower=0, upper=obsXobjeff_Nw+1> obsXobjeff_u[Nobservations+1];
  int<lower=0, upper=NobjEffects> obsXobjeff_v[obsXobjeff_Nw];

  // obsXobj_batcheff (observation X batch_effect) sparse matrix
  int<lower=0> obsXobjbatcheff_Nw;
  vector[obsXobjbatcheff_Nw] obsXobjbatcheff_w;
  int<lower=0, upper=obsXobjbatcheff_Nw+1> obsXobjbatcheff_u[Nobservations + 1];
  int<lower=0, upper=NobjBatchEffects> obsXobjbatcheff_v[obsXobjbatcheff_Nw];

  // obsXobj_suobatcheff (observation X sub-object batch_effect) sparse matrix
  int<lower=0> suoxobsXsuobatcheff_Nw;
  vector[suoxobsXsuobatcheff_Nw] suoxobsXsuobatcheff_w;
  int<lower=0, upper=suoxobsXsuobatcheff_Nw+1> suoxobsXsuobatcheff_u[Nsubobjects*Nobservations + 1];
  int<lower=0, upper=NsuoBatchEffects> suoxobsXsuobatcheff_v[suoxobsXsuobatcheff_Nw];

  // global model constants
  real<lower=0.0> hsprior_lambda_a_offset; // prevent lambda_a/t being too close to zero, because that negatively affects MCMC convergence
  real<lower=0.0> hsprior_lambda_t_offset;
  real global_labu_shift;   // shift to be applied to all XXX_labu variables to get the real log intensity
  vector<lower=0>[Neffects] effect_tau;
  vector<lower=0>[Neffects] effect_df;
  vector<lower=0>[Neffects] effect_df2;
  real<lower=0> effect_slab_df;
  real<lower=0> effect_slab_scale;

  real obj_labu_min; // minimal average abundance of an object
  real<lower=0> obj_labu_min_scale; // scale that defines the softness of lower abundance limit
  real<lower=0> obj_base_labu_sigma; // sigma of average abundance distribution
  real<lower=0> iact_repl_shift_tau;
  real<lower=0> iact_repl_shift_df;
  real<lower=0> batch_effect_sigma;
  real<lower=0> suo_subbatch_effect_tau;
  real<lower=0> suo_subbatch_effect_df;
  real<lower=0> suo_subbatch_effect_c;
  real<upper=0> underdef_obj_shift;

  // instrument calibrated parameters (FIXME: msproto-dependent)
  real<lower=0> zDetectionFactor;
  real zDetectionIntercept;

  real<lower=0> sigmaScaleHi;
  real<lower=0> sigmaScaleLo;
  real sigmaOffset;
  real sigmaBend;
  real sigmaSmooth;

  real zShift;
  real zScale;
}

transformed data {
  real mzShift = zShift - global_labu_shift; // zShift for the missing observation intensity
  vector[Nquanted] zScore = (log(qData) - zShift) * zScale;
  vector[Nquanted] qLogStd; // log(sd(qData))-global_labu_shift
  vector<lower=0>[Nquanted] qDataNorm; // qData/sd(qData)
  int<lower=0,upper=Nquanted> NreliableQuants = sum(observation_reliable[quant2observation]);
  int<lower=1,upper=Nquanted> reliable_quants[NreliableQuants];

  int<lower=1,upper=Niactions> quant2iaction[Nquanted] = observation2iaction[quant2observation];
  int<lower=1,upper=Nexperiments> quant2experiment[Nquanted] = observation2experiment[quant2observation];
  int<lower=1,upper=Niactions> miss2iaction[Nmissed] = observation2iaction[miss2observation];
  int<lower=1,upper=Nexperiments> miss2experiment[Nmissed] = observation2experiment[miss2observation];

  int<lower=0,upper=NobjEffects> NobjEffectsPos = sum(effect_is_positive[obj_effect2effect]);
  int<lower=0,upper=NobjEffects> NobjEffectsOther = NobjEffects - NobjEffectsPos;
  int<lower=1,upper=NobjEffects> obj_effect_reshuffle[NobjEffects];
  vector<lower=0>[NobjEffects] obj_effect_tau = effect_tau[obj_effect2effect];
  vector[NobjEffects] obj_effect_mean = effect_mean[obj_effect2effect];
  vector<lower=0>[NobjEffects] obj_effect_df = effect_df[obj_effect2effect];
  vector<lower=0>[NobjEffects] obj_effect_df2 = effect_df2[obj_effect2effect];

  int<lower=0,upper=NobjBatchEffects> NobjBatchEffectsPos = sum(batch_effect_is_positive[obj_batch_effect2batch_effect]);
  int<lower=0,upper=NobjBatchEffects> NobjBatchEffectsOther = NobjBatchEffects - NobjBatchEffectsPos;
  int<lower=1,upper=NobjBatchEffects> obj_batch_effect_reshuffle[NobjBatchEffects];

  int<lower=0,upper=NsuoBatchEffects> NsuoBatchEffectsPos = sum(subbatch_effect_is_positive[suo_subbatch_effect2subbatch_effect]);
  int<lower=0,upper=NsuoBatchEffects> NsuoBatchEffectsOther = NsuoBatchEffects - NsuoBatchEffectsPos;
  int<lower=1,upper=NsuoBatchEffects> suo_subbatch_effect_reshuffle[NsuoBatchEffects];
  int<lower=1,upper=Nobservations*Nsubobjects> quant2suoxobs[NsuoBatchEffects > 0 ? Nquanted : 0];
  int<lower=1,upper=Nobservations*Nsubobjects> miss2suoxobs[NsuoBatchEffects > 0 ? Nmissed : 0];
  real<lower=0> suo_subbatch_effect_c2 = square(suo_subbatch_effect_c);

  vector[Niactions] iactXobjbase_w = rep_vector(1.0, Niactions);
  int<lower=0> iactXobjbase_u[Niactions + 1] = one_to(Niactions + 1);

  vector[Nobservations] obsXiact_w = rep_vector(1.0, Nobservations);
  int<lower=0> obsXiact_u[Nobservations + 1] = one_to(Nobservations + 1);

  vector[Nobservations] obsXobjbase_w = rep_vector(1.0, Nobservations);
  int<lower=1> obsXobjbase_u[Nobservations + 1] = one_to(Nobservations + 1);
  int<lower=1, upper=Nobjects> obs2obj[Nobservations] = iaction2obj[observation2iaction];

  int<lower=0> NrealIactions = ndistinct(observation2iaction, Niactions);
  int<lower=0> Nobservations0 = Nobservations - NrealIactions; // number of observations degrees of freedom ()
  int<lower=0> obsXobs0_Nw = contr_poly_Nw(NrealIactions, observation2iaction);
  vector[obsXobs0_Nw] obsXobs_shift0_w;
  int<lower=1, upper=obsXobs0_Nw + 1> obsXobs_shift0_u[Nobservations + 1];
  int<lower=1, upper=Nobservations0> obsXobs_shift0_v[obsXobs0_Nw];

  int<lower=0> suoXsuo0_Nw = contr_treatment_Nw(Nobjects, suo2obj);
  vector[suoXsuo0_Nw] suoXsuo_shift0_w;
  int<lower=0, upper=suoXsuo0_Nw + 1> suoXsuo_shift0_u[Nsubobjects + 1];
  int<lower=0, upper=Nsubobjects - Nobjects> suoXsuo_shift0_v[suoXsuo0_Nw];

  matrix[Nobjects + NobjEffects, Niactions] iaction2objeffx_op;

  // prepare reshuffling of positive/other effects
  obj_effect_reshuffle = objeffects_reshuffle(obj_effect2effect, effect_is_positive);
  obj_batch_effect_reshuffle = objeffects_reshuffle(obj_batch_effect2batch_effect, batch_effect_is_positive);
  suo_subbatch_effect_reshuffle = objeffects_reshuffle(suo_subbatch_effect2subbatch_effect, subbatch_effect_is_positive);

  // process the intensity data to optimize likelihood calculation
  {
    for (i in 1:Nquanted) {
      qLogStd[i] = intensity_log_std(zScore[i], sigmaScaleHi, sigmaScaleLo, sigmaOffset, sigmaBend, sigmaSmooth);
      qDataNorm[i] = exp(log(qData[i]) - qLogStd[i]);
      qLogStd[i] -= global_labu_shift; // obs_labu is modeled without obj_base
    }
  }

  // collect quantifications of reliable object observations
  {
    int j = 1;
    for (i in 1:Nquanted) {
      if (observation_reliable[quant2observation[i]]) {
        reliable_quants[j] = i;
        j += 1;
      }
    }
  }

  // prepare obsXobs_shift0
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
            int nW = 0;
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
              real w = iact_obsXobs0[iaction2nobs_2ndpass[iact_ix], j];
              //if (w != 0) { // don't skip non-structural zeros for contr_poly (otherwise there's mismatch with contr_poly_Nw)
                nW += 1;
                obsXobs_shift0_v[obsXobs_shift0_offset + nW] = iaction2obs_shift0_offset[iact_ix] + j;
                obsXobs_shift0_w[obsXobs_shift0_offset + nW] = w;
              //}
            }
            obsXobs_shift0_u[i+1] = obsXobs_shift0_u[i] + nW;
            obsXobs_shift0_offset += nW;
        }
    }
  }
  //print("obsXobs_shift0=", csr_to_dense_matrix(Nobservations, Nobservations0, obsXobs_shift0_w, obsXobs_shift0_v, obsXobs_shift0_u));

  {
    matrix[Niactions, Nobjects+NobjEffects] objeffx2iaction_op;
    // = append_col(
    //  csr_to_dense_matrix(Niactions, Nobjects, iactXobjbase_w, iaction2obj, iactXobjbase_u),
    //  csr_to_dense_matrix(Niactions, NobjEffects, iactXobjeff_w, iactXobjeff_v, iactXobjeff_u));
    // wierd way of convertion
    for (i in 1:Nobjects) {
      vector[Nobjects] obj = rep_vector(0.0, Nobjects);
      vector[Niactions] iact;
      obj[i] = 1.0;
      iact = csr_matrix_times_vector(Niactions, Nobjects, iactXobjbase_w, iaction2obj, iactXobjbase_u, obj);
      for (j in 1:Niactions) {
        objeffx2iaction_op[j, i] = iact[j];
      }
    }
    for (i in 1:NobjEffects) {
      vector[NobjEffects] objeff = rep_vector(0.0, NobjEffects);
      vector[Niactions] iact;
      objeff[i] = 1.0;
      iact = csr_matrix_times_vector(Niactions, NobjEffects, iactXobjeff_w, iactXobjeff_v, iactXobjeff_u, objeff);
      for (j in 1:Niactions) {
        objeffx2iaction_op[j, i + Nobjects] = iact[j];
      }
    }
    // OLS operator
    iaction2objeffx_op = crossprod(objeffx2iaction_op)\(objeffx2iaction_op');
    //print("iaction2objeffx_op=", iaction2objeffx_op);
  }

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
            int nW = 0;
            matrix[obj_nsuo, obj_nsuo-1] obj_suoXsuo0 = contr_treatment(obj_nsuo);

            if (obj2nsuo_2ndpass[obj_ix] == 0) {
                // reserve (nsuo-1) suo_shift0 variables
                obj2suo_shift0_offset[obj_ix] = suo_shift0_offset;
                suo_shift0_offset += obj_nsuo-1;
                //print("suo_shift0_offset=", suo_shift0_offset);
            }
            obj2nsuo_2ndpass[obj_ix] += 1;
            // add 2npass-th row of obj_suoXsuo0 to the suoXsuo0
            for (j in 1:cols(obj_suoXsuo0)) {
              real w = obj_suoXsuo0[obj2nsuo_2ndpass[obj_ix], j];
              if (w != 0) {
                nW += 1;
                suoXsuo_shift0_v[suoXsuo_shift0_offset + nW] = obj2suo_shift0_offset[obj_ix] + j;
                suoXsuo_shift0_w[suoXsuo_shift0_offset + nW] = w;
              }
            }
            suoXsuo_shift0_u[i+1] = suoXsuo_shift0_u[i] + nW;
            suoXsuo_shift0_offset += nW;
        }
    }

    if (NsuoBatchEffects > 0) {
      // all references to the 1st protocol are redirected to index 1 (this shift is always 0)
      for (i in 1:Nquanted) {
        quant2suoxobs[i] = (quant2observation[i]-1)*Nsubobjects + quant2suo[i];
      }
      for (i in 1:Nmissed) {
        miss2suoxobs[i] = (miss2observation[i]-1)*Nsubobjects + miss2suo[i];
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

  real<lower=1.0> suo_shift_sigma;
  vector[Nsubobjects > 0 ? Nsubobjects-Nobjects : 0] suo_shift0_unscaled; // subobject shift within object

  //real<lower=0.0> obj_effect_tau;
  real<lower=0.0> effect_slab_c_t;
  vector<lower=0.0>[NobjEffects] obj_effect_lambda_t;
  vector<lower=0.0>[NobjEffects] obj_effect_lambda_a;
  vector<lower=0.0>[NobjEffects] obj_effect_eta_t;
  vector<lower=0.0>[NobjEffects] obj_effect_eta_a;
  vector<lower=0.0>[NobjEffectsPos] obj_effect_unscaled_pos;
  vector[NobjEffectsOther] obj_effect_unscaled_other;

  //real<lower=0> obj_repl_effect_sigma;
  //vector<lower=0>[Nobjects*Nexperiments] repl_shift_lambda;
  vector<lower=0>[Nobservations0 > 0 ? Niactions : 0] iact_repl_shift_lambda_t;
  vector<lower=0>[Nobservations0 > 0 ? Niactions : 0] iact_repl_shift_lambda_a;
  vector[Nobservations0] obs_shift0;

  //real<lower=0> obj_batch_effect_sigma;
  vector<lower=0.0>[NobjBatchEffectsPos] obj_batch_effect_unscaled_pos;
  vector[NobjBatchEffectsOther] obj_batch_effect_unscaled_other;

  vector<lower=0>[NsuoBatchEffects] suo_subbatch_effect_lambda_t;
  vector<lower=0>[NsuoBatchEffects] suo_subbatch_effect_lambda_a;
  vector<lower=0.0>[NsuoBatchEffectsPos] suo_subbatch_effect_unscaled_pos;
  vector[NsuoBatchEffectsOther] suo_subbatch_effect_unscaled_other;
}

transformed parameters {
  vector[Nobjects] obj_base_labu;
  vector[NobjEffects] obj_effect;
  real<lower=0> effect_slab_c;
  vector<lower=0>[NobjEffects] obj_effect_sigma; // AKA lambda_tilde*tau in rstanarm
  vector[NobjBatchEffects] obj_batch_effect;
  //vector<lower=0>[NobjBatchEffects] obj_batch_effect_sigma;
  vector[NsuoBatchEffects] suo_subbatch_effect;
  vector<lower=0>[NsuoBatchEffects] suo_subbatch_effect_sigma;

  vector<lower=0>[Nobservations0 > 0 ? Niactions : 0] iact_repl_shift_sigma;

  vector[Nobservations] obs_labu; // iaction_labu + objXexp_repl_shift * obj_repl_shift_sigma
  vector[Nobservations0 > 0 ? Nobservations : 0] obs_repl_shift; // replicate shifts for all potential observations (including missing)
  vector[NobjBatchEffects > 0 ? Nobservations : 0] obs_batch_shift;
  vector[NsuoBatchEffects > 0 ? Nobservations*Nsubobjects : 0] suoxobs_subbatch_shift;

  vector[Nsubobjects] suo_shift_unscaled; // subcomponent shift within object
  vector[Nsubobjects] suo_shift; // subcomponent shift within object

  // correct baseline abundances of underdefined objects
  obj_base_labu = obj_base_labu0 * obj_base_labu_sigma;
  for (i in 1:NunderdefObjs) {
    obj_base_labu[underdef_objs[i]] = obj_base_labu[underdef_objs[i]] + underdef_obj_shift;
  }

  // calculate effects lambdas and scale effects
  {
      vector[NobjEffects] obj_effect_sigma_pre; // AKA lambda_eta2 in rstanarm
      real effect_slab_c2 = square(effect_slab_scale) * effect_slab_c_t; // AKA c2 in rstanarm
      effect_slab_c = effect_slab_scale * sqrt(effect_slab_c_t);
      obj_effect_sigma_pre = square(obj_effect_lambda_a .* obj_effect_eta_a) .* obj_effect_lambda_t .* obj_effect_eta_t;
      obj_effect_sigma = sqrt(effect_slab_c2 * obj_effect_sigma_pre ./ (effect_slab_c2 + square(obj_effect_tau) .* obj_effect_sigma_pre)) .* obj_effect_tau;
  }
  obj_effect = obj_effect_mean + append_row(obj_effect_unscaled_pos, obj_effect_unscaled_other)[obj_effect_reshuffle] .* obj_effect_sigma;

  // calculate observations log abundance
  obs_labu = csr_matrix_times_vector(Nobservations, Nobjects, obsXobjbase_w, obs2obj, obsXobjbase_u, obj_base_labu) +
             csr_matrix_times_vector(Nobservations, NobjEffects, obsXobjeff_w, obsXobjeff_v, obsXobjeff_u, obj_effect);

  // calculate obs_shift and obs_labu
  if (Nobservations0 > 0) {
    iact_repl_shift_sigma = iact_repl_shift_lambda_a .* sqrt(iact_repl_shift_lambda_t) * iact_repl_shift_tau;
    obs_repl_shift = csr_matrix_times_vector(Nobservations, Nobservations0, obsXobs_shift0_w, obsXobs_shift0_v, obsXobs_shift0_u, obs_shift0) .*
        csr_matrix_times_vector(Nobservations, Niactions, obsXiact_w, observation2iaction, obsXiact_u, iact_repl_shift_sigma);
    obs_labu += obs_repl_shift;
  }
  // calculate obs_batch_shift (doesn't make sense to add to obs_labu)
  if (NbatchEffects > 0) {
    obj_batch_effect = append_row(obj_batch_effect_unscaled_pos, obj_batch_effect_unscaled_other)[obj_batch_effect_reshuffle] * batch_effect_sigma;
    obs_batch_shift = csr_matrix_times_vector(Nobservations, NobjBatchEffects, obsXobjbatcheff_w, obsXobjbatcheff_v, obsXobjbatcheff_u, obj_batch_effect);
  }
  // calculate suo_labu_shift
  if (Nsubobjects > 1) {
    suo_shift_unscaled = csr_matrix_times_vector(Nsubobjects, Nsubobjects - Nobjects, suoXsuo_shift0_w, suoXsuo_shift0_v, suoXsuo_shift0_u, suo_shift0_unscaled);
  } else if (Nsubobjects == 1) {
    suo_shift_unscaled = rep_vector(0.0, Nsubobjects);
  }
  suo_shift = suo_shift_unscaled * suo_shift_sigma;

  // calculate suoXobs_subbatch_shift (doesn't make sense to add to obs_labu)
  if (NsuoBatchEffects > 0) {
    vector[NsuoBatchEffects] suo_subbatch_effect_sigma_pre; // AKA lambda_eta2 in rstanarm
    suo_subbatch_effect_sigma_pre = square(suo_subbatch_effect_lambda_a) .* suo_subbatch_effect_lambda_t;
    suo_subbatch_effect_sigma = sqrt(suo_subbatch_effect_c2 * suo_subbatch_effect_sigma_pre ./ (suo_subbatch_effect_c2 + square(suo_subbatch_effect_tau) * suo_subbatch_effect_sigma_pre)) * suo_subbatch_effect_tau;
    suo_subbatch_effect = append_row(suo_subbatch_effect_unscaled_pos, suo_subbatch_effect_unscaled_other)[suo_subbatch_effect_reshuffle] .* suo_subbatch_effect_sigma;
    suoxobs_subbatch_shift = csr_matrix_times_vector(Nobservations*Nsubobjects, NsuoBatchEffects, suoxobsXsuobatcheff_w, suoxobsXsuobatcheff_v, suoxobsXsuobatcheff_u,
                                                     suo_subbatch_effect);
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

    //underdef_obj_shift ~ normal(0.0, 10.0);

    //repl_shift_lambda ~ student_t(2, 0.0, repl_shift_tau);
    //obj_repl_effect ~ normal(0.0, obj_repl_effect_lambda);
    if (Nobservations0 > 0) {
      vector[Nobservations] obs_repl_shift_unscaled; // unscaled observations shifts

      iact_repl_shift_lambda_t - hsprior_lambda_t_offset ~ inv_gamma(0.5 * iact_repl_shift_df, 0.5 * iact_repl_shift_df);
      iact_repl_shift_lambda_a - hsprior_lambda_a_offset ~ std_normal();

      //print("iact_repl_shift_sigma=", iact_repl_shift_sigma);
      //print("obsXiact=", csr_to_dense_matrix(Nobservations, Niactions,
      //          obsXiact_w, observation2iaction, obsXiact_u));
      obs_repl_shift_unscaled = csr_matrix_times_vector(Nobservations, Nobservations0, obsXobs_shift0_w, obsXobs_shift0_v, obsXobs_shift0_u, obs_shift0);
      obs_repl_shift_unscaled ~ std_normal();
    }
    //to_vector(repl_shift) ~ normal(0.0, repl_shift_lambda);

    //obj_batch_effect_lambda ~ student_t(2, 0.0, obj_batch_effect_tau);
    if (NbatchEffects > 0) {
      obj_batch_effect_unscaled_pos ~ std_normal();
      obj_batch_effect_unscaled_other ~ std_normal();
    }
    if (Nsubobjects > 0) {
      suo_shift_sigma ~ cauchy(0, 1);
      suo_shift_unscaled ~ std_normal();
      if (NsuoBatchEffects > 0) {
        suo_subbatch_effect_lambda_t - hsprior_lambda_t_offset ~ inv_gamma(0.5 * suo_subbatch_effect_df, 0.5 * suo_subbatch_effect_df);
        suo_subbatch_effect_lambda_a - hsprior_lambda_a_offset ~ std_normal();
        //obj_batch_effect ~ normal(0.0, obj_batch_effect_lambda);
        suo_subbatch_effect_unscaled_pos ~ std_normal();
        suo_subbatch_effect_unscaled_other ~ std_normal();
      }
    }

    // soft lower limit of protein abundance for each observation
    1 ~ bernoulli_logit((obs_labu - obj_labu_min) * obj_labu_min_scale);

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
            q_labu += suo_shift[quant2suo];
            m_labu += suo_shift[miss2suo];

            if (NsuoBatchEffects > 0) {
                q_labu += suoxobs_subbatch_shift[quant2suoxobs];
                m_labu += suoxobs_subbatch_shift[miss2suoxobs];
            }
        }
        if (NbatchEffects > 0) {
          q_labu += obs_batch_shift[quant2observation];
          m_labu += obs_batch_shift[miss2observation];
        }

        // model quantitations and missing data
        logcompressv(exp(q_labu - qLogStd) - qDataNorm, 0.25) ~ double_exponential(0.0, 1);
        // soft-lower-limit for subobject identification of reliable observations
        // for non-reliable observations (false identifications) we rely on double exponentual to handle outliers
        1 ~ bernoulli_logit(q_labu[reliable_quants] * (zScale * zDetectionFactor) + (-mzShift * zScale * zDetectionFactor + zDetectionIntercept));
        0 ~ bernoulli_logit(missing_sigmoid_scale .* (m_labu * (zScale * zDetectionFactor) + (-mzShift * zScale * zDetectionFactor + zDetectionIntercept)));
    }
}

generated quantities {
    vector[Nobjects] obj_base_labu_replCI;
    vector[NobjEffects] obj_effect_replCI;
    vector[Niactions] iaction_labu;
    vector[Niactions] iaction_labu_replCI;
    vector[Nsubobjects] suo_llh;

    // calculate interactions log abundance
    iaction_labu = csr_matrix_times_vector(Niactions, Nobjects, iactXobjbase_w, iaction2obj, iactXobjbase_u, obj_base_labu) +
                   csr_matrix_times_vector(Niactions, NobjEffects, iactXobjeff_w, iactXobjeff_v, iactXobjeff_u, obj_effect);
    for (i in 1:Niactions) {
      iaction_labu_replCI[i] = normal_rng(iaction_labu[i], iact_repl_shift_sigma[i]);
    }
    {
      vector[Nobjects + NobjEffects] obj_effx = iaction2objeffx_op * iaction_labu_replCI;
      for (i in 1:Nobjects) {
        obj_base_labu_replCI[i] = obj_effx[i];
      }
      for (i in 1:NobjEffects) {
        obj_effect_replCI[i] = obj_effx[i+Nobjects];
      }
    }

    // per-subobject loglikelihood (the code copied from "model" section)
    if (Nsubobjects > 0) {
        vector[Nquanted] q_labu;
        vector[Nmissed] m_labu;

        // prepare predicted abundances
        q_labu = obs_labu[quant2observation] + experiment_shift[quant2experiment] + suo_shift[quant2suo];
        m_labu = obs_labu[miss2observation] + experiment_shift[miss2experiment] + suo_shift[miss2suo];

        if (NsuoBatchEffects > 0) {
            q_labu += suoxobs_subbatch_shift[quant2suoxobs];
            m_labu += suoxobs_subbatch_shift[miss2suoxobs];
        }
        if (NbatchEffects > 0) {
          q_labu += obs_batch_shift[quant2observation];
          m_labu += obs_batch_shift[miss2observation];
        }

        // calculate log-likelihood per subobject
        suo_llh = rep_vector(0.0, Nsubobjects);
        for (i in 1:Nquanted) {
          suo_llh[quant2suo[i]] += double_exponential_lpdf(qDataNorm[i] | exp(q_labu[i] - qLogStd[i]), 1) +
              bernoulli_logit_lpmf(1 | q_labu[i] * (zScale * zDetectionFactor) + (-mzShift * zScale * zDetectionFactor + zDetectionIntercept));
        }
        for (i in 1:Nmissed) {
          suo_llh[miss2suo[i]] += bernoulli_logit_lpmf(0 | missing_sigmoid_scale[i] * m_labu[i] * (zScale * zDetectionFactor) + (-mzShift * zScale * zDetectionFactor + zDetectionIntercept));
        }
    }
}
