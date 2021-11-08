//functions {
    real intensity_log2_std(real z, data real scaleHi, data real scaleLo,
                            data real offs, data real bend, data real smooth) {
        return 0.5*(scaleHi+scaleLo)*(z-bend) + 0.5*(scaleHi-scaleLo)*sqrt((z-bend)*(z-bend)+smooth) + offs;
    }

    // compresses x: pdf(Normal,x) >= outlierProb*pdf(Cauchy,x) -> cauchy_compress(x) = x
    //               otherwise x is transformed so that pdf(Normal, compress(x)) = outlierProb*pdf(Cauchy, x)
    // k modulates the smoothness of log_sum_exp(),
    // k=4 make log_mix(..., x) match pdf(std_normal, cauchy_compress(x))
    real cauchy_compress(real x, data real a, data real k) {
      return fabs(x) - log1p_exp(k * (fabs(x) - sqrt(2 * log1p(square(x)) + a))) / k;
    }
    vector cauchy_compressv(vector x, data real a, data real k) {
      return fabs(x) - inv(k) * log1p_exp(k * (fabs(x) - sqrt(2 * log1p(square(x)) + a)));
    }
    real cauchy_compress_a(data real outlierProb) {
      return log(pi()/2) - 2 * log(outlierProb);
    }
//}