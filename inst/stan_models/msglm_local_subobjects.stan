functions {
#include include/array_utils.stan
#include include/contrast_matrices.stan
#include include/msnoise.stan
}

data {
  int<lower=1> Nconditions;     // number of experimental conditions
  int<lower=0> Nobjects;        // number of objects (proteins/peptides/sites etc)

  int<lower=0> Nsubobjects;     // number of objects subcomponents (peptides of proteins etc), 0 if not supported
  int<lower=1,upper=Nobjects> subobj2obj[Nsubobjects];

  int<lower=0> Niactions;       // number of interactions (observed objectXcondition pairs)
  int<lower=1,upper=Nobjects> iaction2obj[Niactions];

  int<lower=0> Nmsprotocols;    // number of MS protocols used
  int<lower=1> Nmschannels;     // number of mschannels
  vector[Nmschannels] mschannel_shift;

  int<lower=0> Nobservations;   // number of observations of interactions (objectXmschannel pairs for all iactions and mschannels of its condition)
  int<lower=1,upper=Nmschannels> observation2mschannel[Nobservations];
  int<lower=1,upper=Niactions> observation2iaction[Nobservations];
  int<lower=1,upper=Nmsprotocols> mschannel2msproto[Nmsprotocols > 0 ? Nmschannels : 0]; // TODO support by the model

  int<lower=0> Nsubobservations;// number of subobject observations (observation X subobject)
  int<lower=1,upper=Nobservations> subobs2obs[Nsubobservations];
  int<lower=1,upper=Nsubobjects> subobs2subobj[Nsubobservations];

  // map from labelXreplicateXobject to observed/missed data
  int<lower=0> Nquanted;        // total number of quantified subobjectsXmschannels
  int<lower=1,upper=Nsubobservations>  quant2subobs[Nquanted];
  int<lower=0,upper=1> quant_isreliable[Nquanted];
  int<lower=0> Nmissed;         // total number of missed subobjectsXmschannels
  int<lower=1,upper=Nsubobservations> miss2subobs[Nmissed];
  vector<lower=0>[Nquanted] qData; // quanted data
  vector<lower=0, upper=1>[Nmissed] missing_sigmoid_scale; // sigmoid scales for indiv. observations (<1 for higher uncertainty)

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

  int<lower=0> NquantBatchEffects;// number of quantification batch effects (that define subobject-level quantification variation, but not biology)
  int<lower=0> NsubobjBatchEffects;
  int<lower=1,upper=NquantBatchEffects> subobj_batch_effect2quant_batch_effect[NsubobjBatchEffects];
  int<lower=0,upper=1> quant_batch_effect_is_positive[NquantBatchEffects];

  // iactXobjeff (interaction X object_effect) sparse matrix
  int<lower=0> iactXobjeff_Nw;
  vector[iactXobjeff_Nw] iactXobjeff_w;
  int<lower=1, upper=iactXobjeff_Nw+1> iactXobjeff_u[Niactions+1];
  int<lower=1, upper=NobjEffects> iactXobjeff_v[iactXobjeff_Nw];

  int<lower=0> obsXobjeff_Nw;
  vector[obsXobjeff_Nw] obsXobjeff_w;
  int<lower=1, upper=obsXobjeff_Nw+1> obsXobjeff_u[Nobservations+1];
  int<lower=1, upper=NobjEffects> obsXobjeff_v[obsXobjeff_Nw];

  // obsXobj_batcheff (observation X object_batch_effect) sparse matrix
  int<lower=0> obsXobjbatcheff_Nw;
  vector[obsXobjbatcheff_Nw] obsXobjbatcheff_w;
  int<lower=1, upper=obsXobjbatcheff_Nw+1> obsXobjbatcheff_u[NbatchEffects > 0 ? Nobservations + 1 : 1];
  int<lower=1, upper=NobjBatchEffects> obsXobjbatcheff_v[obsXobjbatcheff_Nw];

  // subobsXsubobjbatcheff (sub-object observation X subobject_batch_effect) sparse matrix
  int<lower=0> subobsXsubobjbatcheff_Nw;
  vector[subobsXsubobjbatcheff_Nw] subobsXsubobjbatcheff_w;
  int<lower=1, upper=subobsXsubobjbatcheff_Nw+1> subobsXsubobjbatcheff_u[NquantBatchEffects > 0 ? Nsubobservations + 1 : 1];
  int<lower=1, upper=NsubobjBatchEffects> subobsXsubobjbatcheff_v[subobsXsubobjbatcheff_Nw];

  // global model constants
  real obj_labu_shift;   // shift to be applied to all XXX_labu variables to get the real log intensity
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
  real<lower=0> iact_repl_shift_tau;
  real<lower=0> iact_repl_shift_df;
  real<lower=0> batch_effect_sigma;
  real<lower=0> quant_batch_effect_tau;
  real<lower=0> quant_batch_effect_df;
  real<lower=0> quant_batch_effect_c;

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
  real mzShift = zShift - obj_labu_shift; // zShift for the missing observation intensity
  vector[Nquanted] zScore = (log2(qData) - zShift) * zScale;
  vector[Nquanted] qLog2Std; // log2(sd(qData))-obj_labu_shift
  vector<lower=0>[Nquanted] qDataNorm; // qData/sd(qData)
  int<lower=0,upper=Nquanted> NreliableQuants = sum(quant_isreliable);
  int<lower=1,upper=Nquanted> reliable_quants[NreliableQuants];

  int<lower=1,upper=Nsubobjects> quant2subobj[Nquanted] = subobs2subobj[quant2subobs];
  int<lower=1,upper=Nobservations> quant2obs[Nquanted] = subobs2obs[quant2subobs];
  int<lower=1,upper=Niactions> quant2iaction[Nquanted] = observation2iaction[quant2obs];
  int<lower=1,upper=Nmschannels> quant2mschannel[Nquanted] = observation2mschannel[quant2obs];

  int<lower=1,upper=Nsubobjects> miss2subobj[Nmissed] = subobs2subobj[miss2subobs];
  int<lower=1,upper=Nobservations> miss2obs[Nmissed] = subobs2obs[miss2subobs];
  int<lower=1,upper=Niactions> miss2iaction[Nmissed] = observation2iaction[miss2obs];
  int<lower=1,upper=Nmschannels> miss2mschannel[Nmissed] = observation2mschannel[miss2obs];

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

  int<lower=0,upper=NsubobjBatchEffects> NsubobjBatchEffectsPos = sum(quant_batch_effect_is_positive[subobj_batch_effect2quant_batch_effect]);
  int<lower=0,upper=NsubobjBatchEffects> NsubobjBatchEffectsOther = NsubobjBatchEffects - NsubobjBatchEffectsPos;
  int<lower=1,upper=NsubobjBatchEffects> subobj_batch_effect_reshuffle[NsubobjBatchEffects] =
      nonzeros_first(quant_batch_effect_is_positive[subobj_batch_effect2quant_batch_effect]);
  real<lower=0> subobj_batch_effect_c2 = square(quant_batch_effect_c);

  vector[Niactions] iactXobjbase_w = rep_vector(1.0, Niactions);
  int<lower=0> iactXobjbase_u[Niactions + 1] = one_to(Niactions + 1);

  vector[Nobservations] obsXiact_w = rep_vector(1.0, Nobservations);
  int<lower=0> obsXiact_u[Nobservations + 1] = one_to(Nobservations + 1);

  vector[Nobservations] obsXobjbase_w = rep_vector(1.0, Nobservations);
  int<lower=1> obsXobjbase_u[Nobservations + 1] = one_to(Nobservations + 1);
  int<lower=1, upper=Nobjects> obs2obj[Nobservations] = iaction2obj[observation2iaction];

  int<lower=0> NrealIactions = ndistinct(observation2iaction, Niactions);
  int<lower=0> Nobservations0 = Nobservations - NrealIactions; // number of observations degrees of freedom ()
  int<lower=0> obs_shiftXobs_shift0_Nw = contr_poly_Nw(Niactions, observation2iaction);
  vector[obs_shiftXobs_shift0_Nw] obs_shiftXobs_shift0_w;
  int<lower=1, upper=obs_shiftXobs_shift0_Nw + 1> obs_shiftXobs_shift0_u[Nobservations + 1];
  int<lower=1, upper=Nobservations0> obs_shiftXobs_shift0_v[obs_shiftXobs_shift0_Nw];

  int<lower=0> subobj_shiftXsubobj_shift0_Nw = contr_treatment_Nw(Nobjects, subobj2obj);
  vector[subobj_shiftXsubobj_shift0_Nw] subobj_shiftXsubobj_shift0_w;
  int<lower=0, upper=subobj_shiftXsubobj_shift0_Nw + 1> subobj_shiftXsubobj_shift0_u[Nsubobjects + 1];
  int<lower=0, upper=Nsubobjects - Nobjects> subobj_shiftXsubobj_shift0_v[subobj_shiftXsubobj_shift0_Nw];

  matrix[Nobjects + NobjEffects, Niactions] iaction2objeffx_op;
  //matrix[Niactions, Nobjects] iactXobjbase = csr_to_dense_matrix(Niactions, Nobjects, iactXobjbase_w, iaction2obj, iactXobjbase_u);
  //matrix[Niactions, NobjEffects] iactXobjeff = csr_to_dense_matrix(Niactions, NobjEffects, iactXobjeff_w, iactXobjeff_v, iactXobjeff_u);

  // process the intensity data to optimize likelihood calculation
  {
    for (i in 1:Nquanted) {
      qLog2Std[i] = intensity_log2_std(zScore[i], sigmaScaleHi, sigmaScaleLo, sigmaOffset, sigmaBend, sigmaSmooth);
      qDataNorm[i] = exp2(log2(qData[i]) - qLog2Std[i]);
      qLog2Std[i] -= obj_labu_shift; // obs_labu is modeled without obj_base
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

  // prepare obs_shiftXobs_shift0
  {
    int iaction2nobs[Niactions];

    int iaction2nobs_2ndpass[Niactions];
    int iaction2obs_shift0_offset[Niactions];

    int obs_shiftXobs_shift0_offset;
    int obs_shift0_offset;

    iaction2nobs = rep_array(0, Niactions);
    for (i in 1:Nobservations) {
      int iact_ix;
      iact_ix = observation2iaction[i];
      iaction2nobs[iact_ix] += 1;
    }
    //print("iaction2nobs=", iaction2nobs);
    obs_shift0_offset = 0;
    obs_shiftXobs_shift0_offset = 0;
    obs_shiftXobs_shift0_u[1] = 1;
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
                obs_shiftXobs_shift0_v[obs_shiftXobs_shift0_offset + nW] = iaction2obs_shift0_offset[iact_ix] + j;
                obs_shiftXobs_shift0_w[obs_shiftXobs_shift0_offset + nW] = w;
              //}
            }
            obs_shiftXobs_shift0_u[i+1] = obs_shiftXobs_shift0_u[i] + nW;
            obs_shiftXobs_shift0_offset += nW;
        }
    }
  }
  //print("obs_shiftXobs_shift0=", csr_to_dense_matrix(Nobservations, Nobservations0, obs_shiftXobs_shift0_w, obs_shiftXobs_shift0_v, obs_shiftXobs_shift0_u));

  {
    matrix[Niactions, Nobjects+NobjEffects] objeffx2iaction_op;
    // = append_col(iactXobjbase, iactXobjeff);
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

  subobj_shiftXsubobj_shift0_u = rep_array(0, Nsubobjects+1);
  if (Nsubobjects > 0) {
    int obj2nsuo[Nobjects];

    int obj2nsubobj_2ndpass[Nobjects];
    int obj2subobj_shift0_offset[Nobjects];

    int subobj_shiftXsubobj_shift0_offset;
    int subobj_shift0_offset;

    obj2nsuo = rep_array(0, Nobjects);
    for (i in 1:Nsubobjects) {
      int obj_ix;
      obj_ix = subobj2obj[i];
      obj2nsuo[obj_ix] += 1;
    }
    //print("object2nsuo=", object2nsuo);
    subobj_shift0_offset = 0;
    subobj_shiftXsubobj_shift0_offset = 0;
    subobj_shiftXsubobj_shift0_u[1] = 1;
    obj2nsubobj_2ndpass = rep_array(0, Nobjects);
    obj2subobj_shift0_offset = rep_array(0, Nobjects);
    for (i in 1:Nsubobjects) {
        int obj_ix;
        int obj_nsuo;
        obj_ix = subobj2obj[i];
        //print("obj_ix[", i, "]=", obj_ix);
        obj_nsuo = obj2nsuo[obj_ix];
        if (obj_nsuo > 1) {
            int nW = 0;
            matrix[obj_nsuo, obj_nsuo-1] obj_suoXsuo0 = contr_treatment(obj_nsuo);

            if (obj2nsubobj_2ndpass[obj_ix] == 0) {
                // reserve (nsuo-1) subobj_shift0 variables
                obj2subobj_shift0_offset[obj_ix] = subobj_shift0_offset;
                subobj_shift0_offset += obj_nsuo-1;
                //print("subobj_shift0_offset=", subobj_shift0_offset);
            }
            obj2nsubobj_2ndpass[obj_ix] += 1;
            // add 2npass-th row of obj_suoXsuo0 to the suoXsuo0
            for (j in 1:cols(obj_suoXsuo0)) {
              real w = obj_suoXsuo0[obj2nsubobj_2ndpass[obj_ix], j];
              if (w != 0) {
                nW += 1;
                subobj_shiftXsubobj_shift0_v[subobj_shiftXsubobj_shift0_offset + nW] = obj2subobj_shift0_offset[obj_ix] + j;
                subobj_shiftXsubobj_shift0_w[subobj_shiftXsubobj_shift0_offset + nW] = w;
              }
            }
            subobj_shiftXsubobj_shift0_u[i+1] = subobj_shiftXsubobj_shift0_u[i] + nW;
            subobj_shiftXsubobj_shift0_offset += nW;
        }
    }
  }
}

parameters {
  //real obj_base;
  //real<lower=0, upper=5.0> obj_shift_sigma;
  //vector<lower=0>[Nconditions] condition_repl_effect_sigma;

  vector[Nobjects] obj_base_labu0; // baseline object abundance

  real<lower=1.0> subobj_shift_sigma;
  vector[Nsubobjects > 0 ? Nsubobjects-Nobjects : 0] subobj_shift0_unscaled; // subobject shift within object

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
  vector<lower=0>[Nobservations0 > 0 ? Niactions : 0] iact_repl_shift_lambda_t;
  vector<lower=0>[Nobservations0 > 0 ? Niactions : 0] iact_repl_shift_lambda_a;
  vector[Nobservations0] obs_shift0;

  //real<lower=0> obj_batch_effect_sigma;
  vector<lower=0.0>[NobjBatchEffectsPos] obj_batch_effect_unscaled_pos;
  vector[NobjBatchEffectsOther] obj_batch_effect_unscaled_other;

  vector<lower=0>[NsubobjBatchEffects] subobj_batch_effect_lambda_t;
  vector<lower=0>[NsubobjBatchEffects] subobj_batch_effect_lambda_a;
  vector<lower=0.0>[NsubobjBatchEffectsPos] subobj_batch_effect_unscaled_pos;
  vector[NsubobjBatchEffectsOther] subobj_batch_effect_unscaled_other;
}

transformed parameters {
  vector[Nobjects] obj_base_labu = obj_base_labu0 * obj_base_labu_sigma;
  vector[NobjEffects] obj_effect;
  real<lower=0> effect_slab_c;
  vector<lower=0>[NobjEffects] obj_effect_sigma; // AKA lambda_tilde*tau in rstanarm
  vector[NobjBatchEffects] obj_batch_effect;
  //vector<lower=0>[NobjBatchEffects] obj_batch_effect_sigma;
  vector[NsubobjBatchEffects] subobj_batch_effect;
  vector<lower=0>[NsubobjBatchEffects] subobj_batch_effect_sigma;

  vector<lower=0>[Nobservations0 > 0 ? Niactions : 0] iact_repl_shift_sigma;

  vector[Nobservations] obs_labu; // iaction_labu + iact_repl_shift * obj_repl_shift_sigma
  vector[Nobservations0 > 0 ? Nobservations : 0] obs_repl_shift; // replicate shifts for all potential observations (including missing)
  vector[NobjBatchEffects > 0 ? Nobservations : 0] obs_batch_shift;
  vector[NsubobjBatchEffects > 0 ? Nsubobservations : 0] subobs_batch_shift;

  vector[Nsubobjects] subobj_shift_unscaled; // subcomponent shift within object
  vector[Nsubobjects] subobj_shift; // subcomponent shift within object

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
    obs_repl_shift = csr_matrix_times_vector(Nobservations, Nobservations0, obs_shiftXobs_shift0_w, obs_shiftXobs_shift0_v, obs_shiftXobs_shift0_u, obs_shift0) .*
        csr_matrix_times_vector(Nobservations, Niactions, obsXiact_w, observation2iaction, obsXiact_u, iact_repl_shift_sigma);
    obs_labu += obs_repl_shift;
  }
  // calculate obs_batch_shift (doesn't make sense to add to obs_labu)
  if (NbatchEffects > 0) {
    obj_batch_effect = append_row(obj_batch_effect_unscaled_pos, obj_batch_effect_unscaled_other)[obj_batch_effect_reshuffle] * batch_effect_sigma;
    obs_batch_shift = csr_matrix_times_vector(Nobservations, NobjBatchEffects, obsXobjbatcheff_w, obsXobjbatcheff_v, obsXobjbatcheff_u, obj_batch_effect);
  }
  // calculate subobj_shift
  if (Nsubobjects > 1) {
    subobj_shift_unscaled = csr_matrix_times_vector(Nsubobjects, Nsubobjects - Nobjects,
                                                    subobj_shiftXsubobj_shift0_w, subobj_shiftXsubobj_shift0_v,
                                                    subobj_shiftXsubobj_shift0_u, subobj_shift0_unscaled);
  } else if (Nsubobjects == 1) {
    subobj_shift_unscaled = rep_vector(0.0, Nsubobjects);
  }
  subobj_shift = subobj_shift_unscaled * subobj_shift_sigma;

  // calculate suoXobs_subbatch_shift (doesn't make sense to add to obs_labu)
  if (NsubobjBatchEffects > 0) {
    vector[NsubobjBatchEffects] subobj_batch_effect_sigma_pre; // AKA lambda_eta2 in rstanarm
    subobj_batch_effect_sigma_pre = square(subobj_batch_effect_lambda_a) .* subobj_batch_effect_lambda_t;
    subobj_batch_effect_sigma = sqrt(subobj_batch_effect_c2 * subobj_batch_effect_sigma_pre ./ (subobj_batch_effect_c2 + square(quant_batch_effect_tau) * subobj_batch_effect_sigma_pre)) * quant_batch_effect_tau;
    subobj_batch_effect = append_row(subobj_batch_effect_unscaled_pos, subobj_batch_effect_unscaled_other)[subobj_batch_effect_reshuffle] .* subobj_batch_effect_sigma;
    subobs_batch_shift = csr_matrix_times_vector(Nsubobservations, NsubobjBatchEffects, subobsXsubobjbatcheff_w, subobsXsubobjbatcheff_v, subobsXsubobjbatcheff_u,
                                                 subobj_batch_effect);
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
    if (Nobservations0 > 0) {
      vector[Nobservations] obs_repl_shift_unscaled; // unscaled observations shifts

      iact_repl_shift_lambda_t - hsprior_lambda_t_offset ~ inv_gamma(0.5 * iact_repl_shift_df, 0.5 * iact_repl_shift_df);
      iact_repl_shift_lambda_a - hsprior_lambda_a_offset ~ std_normal();

      //print("iact_repl_shift_sigma=", iact_repl_shift_sigma);
      //print("obsXiact=", csr_to_dense_matrix(Nobservations, Niactions,
      //          obsXiact_w, observation2iaction, obsXiact_u));
      obs_repl_shift_unscaled = csr_matrix_times_vector(Nobservations, Nobservations0, obs_shiftXobs_shift0_w, obs_shiftXobs_shift0_v, obs_shiftXobs_shift0_u, obs_shift0);
      obs_repl_shift_unscaled ~ std_normal();
    }
    //to_vector(repl_shift) ~ normal(0.0, repl_shift_lambda);

    //obj_batch_effect_lambda ~ student_t(2, 0.0, obj_batch_effect_tau);
    if (NbatchEffects > 0) {
      obj_batch_effect_unscaled_pos ~ std_normal();
      obj_batch_effect_unscaled_other ~ std_normal();
    }
    if (Nsubobjects > 0) {
      subobj_shift_sigma ~ cauchy(0, 1);
      subobj_shift_unscaled ~ std_normal();
      if (NsubobjBatchEffects > 0) {
        subobj_batch_effect_lambda_t - hsprior_lambda_t_offset ~ inv_gamma(0.5 * quant_batch_effect_df, 0.5 * quant_batch_effect_df);
        subobj_batch_effect_lambda_a - hsprior_lambda_a_offset ~ std_normal();
        //obj_batch_effect ~ normal(0.0, obj_batch_effect_lambda);
        subobj_batch_effect_unscaled_pos ~ std_normal();
        subobj_batch_effect_unscaled_other ~ std_normal();
      }
    }

    // soft lower limit of protein abundance for each observation
    1 ~ bernoulli_logit((obs_labu - obj_labu_min) * obj_labu_min_scale);

    // calculate the likelihood
    {
        vector[Nquanted] q_labu;
        vector[Nmissed] m_labu;

        q_labu = obs_labu[quant2obs] + mschannel_shift[quant2mschannel];
        m_labu = obs_labu[miss2obs] + mschannel_shift[miss2mschannel];
        if (Nsubobjects > 0) {
            // adjust by subcomponent shift
            q_labu += subobj_shift[quant2subobj];
            m_labu += subobj_shift[miss2subobj];

            if (NsubobjBatchEffects > 0) {
                q_labu += subobs_batch_shift[quant2subobs];
                m_labu += subobs_batch_shift[miss2subobs];
            }
        }
        if (NbatchEffects > 0) {
          q_labu += obs_batch_shift[quant2obs];
          m_labu += obs_batch_shift[miss2obs];
        }

        // model quantitations and missing data
        logcompressv(exp2(q_labu - qLog2Std) - qDataNorm, 0.25) ~ double_exponential(0.0, 1);
        // soft-lower-limit for subobject intensities of reliable quantifications
        1 ~ bernoulli_logit(q_labu[reliable_quants] * (zScale * zDetectionFactor) + (-mzShift * zScale * zDetectionFactor + zDetectionIntercept));
        0 ~ bernoulli_logit(missing_sigmoid_scale .* (m_labu * (zScale * zDetectionFactor) + (-mzShift * zScale * zDetectionFactor + zDetectionIntercept)));
    }
}

generated quantities {
    // calculate interactions log abundance
    vector[Niactions] iaction_labu =
        csr_matrix_times_vector(Niactions, Nobjects, iactXobjbase_w, iaction2obj, iactXobjbase_u, obj_base_labu) +
        csr_matrix_times_vector(Niactions, NobjEffects, iactXobjeff_w, iactXobjeff_v, iactXobjeff_u, obj_effect);
        //obj_base_labu[iaction2obj] + iactXobjeff * obj_effect;
    vector[Niactions] iaction_labu_replCI = to_vector(normal_rng(iaction_labu, iact_repl_shift_sigma));
    vector[Nobjects] obj_base_labu_replCI;
    vector[NobjEffects] obj_effect_replCI;
    vector[Nsubobjects] subobj_llh;

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
        q_labu = obs_labu[quant2obs] + mschannel_shift[quant2mschannel] + subobj_shift[quant2subobj];
        m_labu = obs_labu[miss2obs] + mschannel_shift[miss2mschannel] + subobj_shift[miss2subobj];

        if (NsubobjBatchEffects > 0) {
            q_labu += subobs_batch_shift[quant2subobs];
            m_labu += subobs_batch_shift[miss2subobs];
        }
        if (NbatchEffects > 0) {
          q_labu += obs_batch_shift[quant2obs];
          m_labu += obs_batch_shift[miss2obs];
        }

        // calculate log-likelihood per subobject
        subobj_llh = rep_vector(0.0, Nsubobjects);
        for (i in 1:Nquanted) {
          subobj_llh[quant2subobj[i]] += double_exponential_lpdf(logcompress(exp2(q_labu[i] - qLog2Std[i]) - qDataNorm[i], 0.25) | 0, 1) +
              bernoulli_logit_lpmf(1 | q_labu[i] * (zScale * zDetectionFactor) + (-mzShift * zScale * zDetectionFactor + zDetectionIntercept));
        }
        for (i in 1:Nmissed) {
          subobj_llh[miss2subobj[i]] += bernoulli_logit_lpmf(0 | missing_sigmoid_scale[i] * m_labu[i] * (zScale * zDetectionFactor) + (-mzShift * zScale * zDetectionFactor + zDetectionIntercept));
        }
    }
}
