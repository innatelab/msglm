functions {
    real intensity_log_std(real z, real scaleHi, real scaleLo, real offset, real bend, real smooth) {
        return 0.5*(scaleHi+scaleLo)*(z-bend) + 0.5*(scaleHi-scaleLo)*sqrt((z-bend)*(z-bend)+smooth) + offset;
    }
}

data {
  int<lower=1> Nmschannels;    // number of MS channels
  int<lower=1> Nobjects;       // number of objects
  int<lower=1> Nconditions;
  vector[Nmschannels]  mschannel_shift; // fixed mschannel shifts
  int<lower=1,upper=Nconditions> mschannel2condition[Nmschannels];

  matrix<lower=0>[Nobjects, Nmschannels] qData;

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
    matrix[Nobjects, Nmschannels] zScore;
    matrix[Nobjects, Nmschannels] qLogStd;
    matrix<lower=0>[Nobjects, Nmschannels] qDataScaled;
    matrix<lower=0>[Nobjects, Nmschannels] meanDenomScaled;
    int<lower=0> NobsPerObject[Nobjects];
    matrix[Nconditions, Nconditions-1] condition_shift_transform;

    condition_shift_transform = rep_matrix(0.0, Nconditions, Nconditions-1);
    for (i in 1:Nconditions-1) {
        condition_shift_transform[i, i] = 1.0;
        condition_shift_transform[Nconditions, i] = -1.0;
    }

    zScore = (log(qData) - zShift) * zScale;

    NobsPerObject = rep_array(0, Nobjects);
    for (i in 1:Nmschannels) {
        for (j in 1:Nobjects) {
            if (qData[j,i] > 0.0) {
                NobsPerObject[j] = NobsPerObject[j] + 1;
            }
        }
    }

    // process the intensity data to optimize likelihood calculation
    for (i in 1:Nmschannels) {
        for (j in 1:Nobjects) {
            qLogStd[j, i] = intensity_log_std(zScore[j, i], sigmaScaleHi, sigmaScaleLo, sigmaOffset, sigmaBend, sigmaSmooth);
        }
    }
    for (i in 1:Nmschannels) {
        for (j in 1:Nobjects) {
            if (qData[j, i] > 0.0) {
                real qScale;
                qScale = exp(-qLogStd[j, i]);
                meanDenomScaled[j, i] = qScale/NobsPerObject[j];
                qDataScaled[j, i] = qScale*qData[j, i];
            } else {
                meanDenomScaled[j, i] = 0.0;
                qDataScaled[j, i] = 0.0;
            }
        }
    }
}

parameters {
    real<lower=0> data_sigma_a;
    real<lower=0> data_sigma_t;
    real<lower=0> condition_sigma_a;
    real<lower=0> condition_sigma_t;
    vector[Nconditions-1] condition_shift0_unscaled;
}

transformed parameters {
    real<lower=0> data_sigma;
    real<lower=0> condition_sigma;
    vector[Nconditions] condition_shift_unscaled;
    vector[Nconditions] condition_shift;

    data_sigma = data_sigma_a ./ sqrt(data_sigma_t);
    condition_sigma = condition_sigma_a ./ sqrt(condition_sigma_t);
    condition_shift_unscaled = condition_shift_transform * condition_shift0_unscaled;
    condition_shift = condition_shift_unscaled * condition_sigma;
}

model {
    vector[Nmschannels] total_mschannel_shift;
    matrix[Nmschannels, Nmschannels] sum_mschannels;

    data_sigma_t ~ gamma(1.0, 1.0); // 1.0 = 2/2
    data_sigma_a ~ normal(0.0, 1.0); // 1.0 = 2/2
    //data_sigma ~ student_t(2, 0.0, 1.0);
    condition_sigma_t ~ gamma(1.0, 1.0); // 1.0 = 2/2
    condition_sigma_a ~ normal(0.0, 1.0); // 1.0 = 2/2
    //condition_sigma ~ student_t(2, 0.0, 1.0);
    condition_shift_unscaled ~ normal(0.0, 1.0);

    total_mschannel_shift = mschannel_shift + condition_shift[mschannel2condition];
    sum_mschannels = exp(rep_matrix(total_mschannel_shift', Nmschannels) - rep_matrix(total_mschannel_shift, Nmschannels));
    //print("avg_exps=", average_mschannels);
    to_vector(qDataScaled) ~ double_exponential(to_vector((qData * sum_mschannels) .* meanDenomScaled), data_sigma);
}
