//functions {
    real intensity_log2_std(real z, real scaleHi, real scaleLo, real offs, real bend, real smooth) {
        return 0.5*(scaleHi+scaleLo)*(z-bend) + 0.5*(scaleHi-scaleLo)*sqrt((z-bend)*(z-bend)+smooth) + offs;
    }

    // compresses x: x~0 -> logcompress(x)~x, abs(x)>>0 -> logcompress(x)~sign(x)*log2(abs(x))
    real logcompress(real x, data real s) {
      return x * (1 + log1p(fabs(s*x))) / (1 + fabs(s*x));
    }

    vector logcompressv(vector x, data real s) {
      return x .* (1 + log1p(fabs(s*x))) ./ (1 + fabs(s*x));
    }
//}