# convert signal to zscore
instrument.signal2zscore <- function(signal, instr_calib) {
  if_else(!is.na(signal) & is.finite(signal),
          (log(signal) - instr_calib$zShift)*instr_calib$zScale,
          NA_real_ )
}

# inverse of sigma
instrument.signal_precision <- function(signal, instr_calib) {
  if_else(!is.na(signal) * is.finite(signal), instrument.zscore_precision(instrument.signal2zscore(signal, instr_calib), instr_calib))
  # throw(BoundsError("$res precision for signal $signal (zscore = $(signal2zscore(params, signal)))"))
}

instrument.zscore_logsd <- function(z, instr_calib) {
  zd = z - instr_calib$sigmaBend
  return (0.5 * (instr_calib$sigmaScaleHi + instr_calib$sigmaScaleLo) * zd +
          0.5 * (instr_calib$sigmaScaleHi - instr_calib$sigmaScaleLo) * sqrt(zd*zd + instr_calib$sigmaSmooth) +
          instr_calib$sigmaOffset)
}

instrument.zscore_sd <- function(z, instr_calib) {
  exp(instrument.zscore_logsd(z, instr_calib))
}
instrument.zscore_precision <- function(z, instr_calib) {
  exp(-instrument.zscore_logsd(z, instr_calib))
}

instrument.signal_sd <- function(signal, instr_calib) {
  if_else(!is.na(signal) & is.finite(signal),
          instrument.zscore_sd(instrument.signal2zscore(signal, instr_calib), instr_calib),
          NA_real_)
}

instrument.signal_likelihood_log <- function(signal, expected, signalPrecision) {
  abs(expected-signal) * signalPrecision # -Distributions.logtwo + log(signalPrecision) #this part is relatively expensive, but is not dependent on expected, so ignore
}

instrument.signal_likelihood_log <- function(signal, expected, instr_calib) {
  instrument.signal_likelihood_log(signal, expected, instrument.signal_precision(signal, instr_calib))
}

# FIXME use Rmath provided one
log1pexp <- function(x) { log1p(exp(x)) }

instrument.detection_likelihood_log <- function(is_detected, expected_log, instr_calib) {
  z = expected_log * instr_calib$signalLogDetectionFactor + instr_calib$signalLogDetectionIntercept
  return (sapply(z, function(x) if_else(is_detected, -log1pexp(-x), #+params.logDetectionMax # invlogit(z)*detMax
                    -log1pexp(x)))) #logsumexp( -Distributions.log1pexp(z)+params.logDetectionMax, params.log1mDetectionMax ) ) # invlogit(-z)*detMax+(1-detMax)
}

