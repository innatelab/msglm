functions {
    real intensity_log_std(real z, real scaleHi, real scaleLo, real offset, real bend, real smooth) {
        return 0.5*(scaleHi+scaleLo)*(z-bend) + 0.5*(scaleHi-scaleLo)*sqrt((z-bend)*(z-bend)+smooth) + offset;
    }
}

data {
  int<lower=1> Nobjects;       // number of objects

  int<lower=1> Nmschannels;    // number of MS channels
  int<lower=1> Nshifts;        // number of mschannel_shifts to estimate
  int<lower=1,upper=Nshifts> mschannel2shift[Nmschannels];

  vector[Nmschannels]  mschannel_preshift; // fixed mschannel pre-shifts
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
    matrix[Nshifts, Nshifts-1] shift_transform;

    shift_transform = rep_matrix(0.0, Nshifts, Nshifts-1);
    for (i in 1:Nshifts-1) {
        shift_transform[i, i] = 1.0;
        shift_transform[Nshifts, i] = -1.0;
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
    real<lower=0> shift_sigma_a;
    real<lower=0> shift_sigma_t;
    vector[Nshifts-1] shift0_unscaled;
}

transformed parameters {
    real<lower=0> data_sigma;
    real<lower=0> shift_sigma;
    vector[Nshifts] shift_unscaled;
    vector[Nshifts] shift;

    data_sigma = data_sigma_a ./ sqrt(data_sigma_t);
    shift_sigma = shift_sigma_a ./ sqrt(shift_sigma_t);
    shift_unscaled = shift_transform * shift0_unscaled;
    shift = shift_unscaled * shift_sigma;
}

model {
    vector[Nmschannels] total_mschan_shift;
    matrix[Nmschannels, Nmschannels] sum_mschans;

    data_sigma_t ~ gamma(1.0, 1.0); // 1.0 = 2/2
    data_sigma_a ~ normal(0.0, 1.0); // 1.0 = 2/2
    //data_sigma ~ student_t(2, 0.0, 1.0);
    shift_sigma_t ~ gamma(1.0, 1.0); // 1.0 = 2/2
    shift_sigma_a ~ normal(0.0, 1.0); // 1.0 = 2/2
    //shift_sigma ~ student_t(2, 0.0, 1.0);
    shift_unscaled ~ normal(0.0, 1.0);

    total_mschan_shift = mschannel_preshift + shift[mschannel2shift];
    // operator to sum channel intensities in each mschannel and copy it to all mschannels
    sum_mschans = exp(rep_matrix(total_mschan_shift', Nmschannels) - rep_matrix(total_mschan_shift, Nmschannels));
    //print("avg_exps=", average_mschans);
    to_vector(qDataScaled) ~ double_exponential(to_vector((qData * sum_mschans) .* meanDenomScaled), data_sigma);
}
