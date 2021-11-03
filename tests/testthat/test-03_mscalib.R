context("mscalib")

require(jsonlite)

test_that("logintensityBase conversion is correct", {
    mscalib_ln <- msglm::read_mscalib_json(test_path("test_mscalib.json"))
    expect_warning(msglm::logintensityBase(mscalib_ln))
    expect_equal(msglm::logintensityBase(mscalib_ln, silent=TRUE), exp(1))

    mscalib_log2 <- msglm::convert_logintensityBase(mscalib_ln, 2)
    expect_equal(msglm::logintensityBase(mscalib_log2), 2)

    intensities <- seq(4, 15, by=0.1)
    ln_logsds <- msglm:::signal_logsd(mscalib_ln, intensities)
    log2_logsds <- msglm:::signal_logsd(mscalib_log2, intensities)
    expect_equal(ln_logsds/log(2), log2_logsds)

    ln_sds <- msglm:::signal_sd(mscalib_ln, intensities)
    log2_sds <- msglm:::signal_sd(mscalib_log2, intensities)
    expect_equal(ln_sds, log2_sds)
    expect_equal(ln_sds, exp(ln_logsds))
    expect_equal(log2_sds, 2^(log2_logsds))

    ln_signal_llh <- msglm:::signal_loglikelihood(mscalib_ln, intensities, intensities)
    log2_signal_llh <- msglm:::signal_loglikelihood(mscalib_log2, intensities, intensities)
    expect_equal(ln_signal_llh, log2_signal_llh)

    ln_det_llh <- msglm:::detection_loglikelihood(mscalib_ln, rep_len(TRUE, length(intensities)), log(intensities))
    log2_det_llh <- msglm:::detection_loglikelihood(mscalib_log2, rep_len(TRUE, length(intensities)), log2(intensities))
    expect_equal(ln_det_llh, log2_det_llh)

    ln_miss_llh <- msglm:::detection_loglikelihood(mscalib_ln, rep_len(FALSE, length(intensities)), log(intensities))
    log2_miss_llh <- msglm:::detection_loglikelihood(mscalib_log2, rep_len(FALSE, length(intensities)), log2(intensities))
    expect_equal(ln_miss_llh, log2_miss_llh)
})