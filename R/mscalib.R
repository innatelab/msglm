# convert signal to zscore
instrument.signal2zscore <- function(signal, instr_calib) {
  logbase <- mscalib_logintensityBase(instr_calib, silent = TRUE)
  if_else(!is.na(signal) & is.finite(signal),
          (log(signal, base=logbase) - instr_calib$zShift)*instr_calib$zScale,
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
  logbase <- mscalib_logintensityBase(instr_calib, silent = TRUE)
  logbase^(instrument.zscore_logsd(z, instr_calib))
}
instrument.zscore_precision <- function(z, instr_calib) {
  logbase <- mscalib_logintensityBase(instr_calib, silent = TRUE)
  logbase^(-instrument.zscore_logsd(z, instr_calib))
}

instrument.signal_logsd <- function(signal, instr_calib) {
  if_else(!is.na(signal) & is.finite(signal),
          instrument.zscore_logsd(instrument.signal2zscore(signal, instr_calib), instr_calib), NA_real_)
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
  return (purrr::map2_dbl(is_detected, z,
              ~ if_else(.x, -log1pexp(-.y), #+params.logDetectionMax # invlogit(z)*detMax
                        -log1pexp(.y)))) #logsumexp( -Distributions.log1pexp(z)+params.logDetectionMax, params.log1mDetectionMax ) ) # invlogit(-z)*detMax+(1-detMax)
}

#' Get the base of the logarithm that is used
#' for log-tranforming the intensities for the given
#' MS noise model.
#'
#' @return logintensityBase property of mscalib
#' @export
mscalib_logintensityBase <- function(mscalib, silent=FALSE) {
  if (rlang::has_name(mscalib, "logintensityBase")) {
    return(mscalib$logintensityBase)
  } else {
    if (!silent) warning("mscalib$logintensityBase not specified, defaulting to e=", exp(1))
    return(exp(1))
  }
}

#' Convert MS noise calibration model for the log_a-transformed intensities
#' to the one for log_b-transformed intensities (b is `new_base`).
#'
#' @param mscalib noise model for natural log intensities
#' @param new_base new base for the log-transformed intensities
#'
#' @return updated mscalib noise model
#' @export
mscalib_convert_logintensityBase <- function(mscalib, new_base, verbose=FALSE) {
  old_base <- mscalib_logintensityBase(mscalib, silent=!verbose)
  if (old_base == new_base) {
    if (verbose) message("Same logintensityBase=", old_base, " no conversion")
    mscalib$logintensityBase <- new_base # make sure now logintensityBase is explicitly set
    return(mscalib)
  }
  k = log(old_base)/log(new_base)
  if (verbose) message("Conversion scaling coefficient k=", k)
  mscalib_new <- list(
    logintensityBase = new_base,
    zShift = mscalib$zShift * k,
    zScale = mscalib$zScale / k,
    zDetectionFactor = mscalib$zDetectionFactor,
    zDetectionIntercept = mscalib$zDetectionIntercept,
    detectionMax = mscalib$detectionMax,
    sigmaBend = mscalib$sigmaBend,
    sigmaScaleHi = mscalib$sigmaScaleHi * k,
    sigmaScaleLo = mscalib$sigmaScaleLo * k,
    sigmaSmooth = mscalib$sigmaSmooth * k^2,
    sigmaOffset = mscalib$sigmaOffset * k
  )
  mscalib_new$signalLogDetectionFactor <- mscalib_new$zScale * mscalib$zDetectionFactor
  # signalLog2DetectionIntercept actually should stay the same
  mscalib_new$signalLogDetectionIntercept <- mscalib_new$zDetectionIntercept -
      mscalib_new$zShift * mscalib_new$zScale * mscalib_new$zDetectionFactor
  return (mscalib_new)
}
