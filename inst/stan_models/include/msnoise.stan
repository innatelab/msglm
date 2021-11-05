//functions {
    real intensity_log2_std(real z, data real scaleHi, data real scaleLo,
                            data real offs, data real bend, data real smooth) {
        return 0.5*(scaleHi+scaleLo)*(z-bend) + 0.5*(scaleHi-scaleLo)*sqrt((z-bend)*(z-bend)+smooth) + offs;
    }

    // compresses x: x~0 -> logcompress(x)~x,
    // logcompress(x0, s, k) = x0, where k = logcompress_k(x0, s)
    // abs(x)>>0 -> logcompress(x)~sign(x)*pow(log(abs(x)), a)
    real logcompress(real x, data real s, data real a, data real k) {
      real t = fabs(s*x);
      return x * (1 + k*pow(log1p(t), a)) / (1 + t);
    }

    vector logcompressv(vector x, data real s, data real a, data real k) {
      vector[size(x)] t = fabs(s*x);
      return x .* (1 + k*pow(log1p(t), a)) ./ (1 + t);
    }

    real logcompress_k(data real x, data real s, data real a) {
      real t = fabs(s*x);
      return t/pow(log1p(t), a);
    }
//}