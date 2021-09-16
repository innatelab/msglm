context("mscalib")

require(jsonlite)

test_that("logintensityBase conversion is correct", {
    quant_calib_ln <- jsonlite::read_json(test_path("test_mscalib.json"))$mscalib
    expect_warning(msglm::mscalib_logintensityBase(quant_calib_ln))
    expect_equal(msglm::mscalib_logintensityBase(quant_calib_ln, silent=TRUE), exp(1))

    quant_calib_log2 <- msglm::mscalib_convert_logintensityBase(quant_calib_ln, 2)
    expect_equal(msglm::mscalib_logintensityBase(quant_calib_log2), 2)

    intensities <- seq(4, 15, by=0.1)
    ln_logsds <- msglm:::instrument.signal_logsd(intensities, quant_calib_ln)
    log2_logsds <- msglm:::instrument.signal_logsd(intensities, quant_calib_log2)
    expect_equal(ln_logsds/log(2), log2_logsds)

    ln_sds <- msglm:::instrument.signal_sd(intensities, quant_calib_ln)
    log2_sds <- msglm:::instrument.signal_sd(intensities, quant_calib_log2)
    expect_equal(ln_sds, log2_sds)
    expect_equal(ln_sds, exp(ln_logsds))
    expect_equal(log2_sds, 2^(log2_logsds))

    ln_det_llh <- msglm:::instrument.detection_likelihood_log(rep_len(TRUE, length(intensities)), log(intensities), quant_calib_ln)
    log2_det_llh <- msglm:::instrument.detection_likelihood_log(rep_len(TRUE, length(intensities)), log2(intensities), quant_calib_log2)
    expect_equal(ln_det_llh, log2_det_llh)

    ln_miss_llh <- msglm:::instrument.detection_likelihood_log(rep_len(FALSE, length(intensities)), log(intensities), quant_calib_ln)
    log2_miss_llh <- msglm:::instrument.detection_likelihood_log(rep_len(FALSE, length(intensities)), log2(intensities), quant_calib_log2)
    expect_equal(ln_miss_llh, log2_miss_llh)
})