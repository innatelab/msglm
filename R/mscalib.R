#' Read JSON file with MS noise model calibration data
#'
#' @param filename MS calibration JSON filename
#'
#' @return mscalib object
#' @export
read_mscalib_json <- function(filename) {
  mscalib_json <- read_json(filename)
  if (!rlang::has_name(mscalib_json, "mscalib")) stop("Key 'mscalib' not found in ", filename)
  return(structure(mscalib_json$mscalib, class="mscalib"))
}

# convert signal to zscore
signal2zscore <- function(mscalib, signal) {
  logbase <- logintensityBase(mscalib, silent = TRUE)
  if_else(!is.na(signal) & is.finite(signal),
          (log(signal, base=logbase) - mscalib$zShift)*mscalib$zScale,
          NA_real_ )
}

zscore_logsd <- function(mscalib, z) {
  zd = z - mscalib$sigmaBend
  return (0.5 * (mscalib$sigmaScaleHi + mscalib$sigmaScaleLo) * zd +
          0.5 * (mscalib$sigmaScaleHi - mscalib$sigmaScaleLo) * sqrt(zd*zd + mscalib$sigmaSmooth) +
          mscalib$sigmaOffset)
}

zscore_sd <- function(mscalib, z) {
  logbase <- logintensityBase(mscalib, silent = TRUE)
  logbase^(zscore_logsd(mscalib, z))
}
zscore_precision <- function(mscalib, z) {
  logbase <- logintensityBase(mscalib, silent = TRUE)
  logbase^(-zscore_logsd(mscalib, z))
}

signal_logsd <- function(mscalib, signal) {
  if_else(!is.na(signal) & is.finite(signal),
          zscore_logsd(mscalib, signal2zscore(mscalib, signal)), NA_real_)
}
signal_sd <- function(mscalib, signal) {
  if_else(!is.na(signal) & is.finite(signal),
          zscore_sd(mscalib, signal2zscore(mscalib, signal)),
          NA_real_)
}
# inverse of sigma
signal_precision <- function(mscalib, signal) {
  if_else(!is.na(signal) & is.finite(signal),
          zscore_precision(mscalib, signal2zscore(mscalib, signal)),
          NA_real_)
  # throw(BoundsError("$res precision for signal $signal (zscore = $(signal2zscore(mscalib, signal)))"))
}

signal_loglikelihood <- function(model, signal, expected) UseMethod("signal_loglikelihood")

signal_loglikelihood.default <- function(signalPrecision, signal, expected) {
  abs(expected-signal) * signalPrecision # -Distributions.logtwo + log(signalPrecision) #this part is relatively expensive, but is not dependent on expected, so ignore
}

signal_loglikelihood.mscalib <- function(mscalib, signal, expected) {
  signal_loglikelihood.default(signal_precision(mscalib, signal), signal, expected)
}

# FIXME use Rmath provided one
log1pexp <- function(x) { log1p(exp(x)) }

detection_loglikelihood <- function(mscalib, is_detected, expected_log) {
  z = expected_log * mscalib$signalLogDetectionFactor + mscalib$signalLogDetectionIntercept
  return (purrr::map2_dbl(is_detected, z,
              ~ if_else(.x, -log1pexp(-.y), #+mscalib.logDetectionMax # invlogit(z)*detMax
                        -log1pexp(.y)))) #logsumexp( -Distributions.log1pexp(z)+mscalib.logDetectionMax, params.log1mDetectionMax ) ) # invlogit(-z)*detMax+(1-detMax)
}

#' Get the base of the logarithm that is used
#' for log-tranforming the intensities for the given
#' MS noise model.
#'
#' @return logintensityBase property of mscalib
#' @export
logintensityBase <- function(mscalib, silent=FALSE) {
  if (rlang::has_name(mscalib, "logintensityBase")) {
    return(mscalib$logintensityBase)
  } else {
    if (!silent) warning(checkmate::vname(mscalib), "$logintensityBase not specified, defaulting to e=", exp(1))
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
convert_logintensityBase <- function(mscalib, new_base, verbose=FALSE) {
  old_base <- logintensityBase(mscalib, silent=!verbose)
  if (old_base == new_base) {
    if (verbose) message("Same logintensityBase=", old_base, " no conversion")
    mscalib$logintensityBase <- new_base # make sure now logintensityBase is explicitly set
    return(mscalib)
  }
  k = log(old_base)/log(new_base)
  if (verbose) message("Conversion scaling coefficient k=", k)
  mscalib_new <- structure(list(
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
  ), class="mscalib")
  mscalib_new$signalLogDetectionFactor <- mscalib_new$zScale * mscalib$zDetectionFactor
  # signalLog2DetectionIntercept actually should stay the same
  mscalib_new$signalLogDetectionIntercept <- mscalib_new$zDetectionIntercept -
      mscalib_new$zShift * mscalib_new$zScale * mscalib_new$zDetectionFactor
  return (mscalib_new)
}
