context("contrasts")

test_that("basic vars_pvalue() is correct", {
  # compose frame with parameters of normally distributed vars
  vars_ref <- tibble(
    varspec = c("c0", "r1", "l1", "r2", "l2", "r3", "l3"),
    mean = c(0, 1, -1, 3, -3, 3, -3),
    sd = c(1, 1, 1, 1, 1, 0.5, 0.5)
  )
  # generate draws.array object by randomly sampling from vars_ref
  vars_draws_raw <- do.call(cbind,
                            dplyr::group_map(dplyr::rowwise(vars_ref),
                                             ~rnorm(100000, mean=.x$mean, sd=.x$sd)) %>%
                              rlang::set_names(vars_ref$varspec))
  vars_draws <- array(vars_draws_raw,
                      dim = c(nrow(vars_draws_raw), 1L, ncol(vars_draws_raw)),
                      dimnames = list(draw=NULL, chain=NULL, variable=colnames(vars_draws_raw))) %>%
    posterior::as_draws()

  for (offs in c(-1, 0, 1)) {
    test_that(paste0("vars_pvalues() with offset=", offs, " are calculated correctly"), {
      vars_ref <- dplyr::rowwise(vars_ref) %>% dplyr::group_modify(~dplyr::mutate(.x,
           prob_nonneg = pnorm(offs, mean, sd, lower.tail = FALSE),
           prob_nonpos = pnorm(offs, mean, sd, lower.tail = TRUE))) %>%
        dplyr::mutate(p_value = 2*pmin(prob_nonneg, prob_nonpos))

    vars_stats <- posterior::summarise_draws(vars_draws)
    varspecs_df <- dplyr::select(vars_stats, varspec=variable) %>%
      dplyr::mutate(var_index = dplyr::row_number(),
                    index_varspec = dplyr::row_number())
    if (offs != 0) {
      varspecs_df$prior_mean = offs
    }
    vars_pvals.draws <- msglm::vars_pvalues(vars_draws, vars_stats, varspecs_df, method="draws")
    vars_pvals.normal <- msglm::vars_pvalues(vars_draws, vars_stats, varspecs_df, method="normal")

    expect_equal(vars_pvals.draws$varspec, vars_ref$varspec)
    expect_equal(vars_pvals.draws$prob_nonneg, vars_ref$prob_nonneg, tolerance = 0.01)
    expect_equal(vars_pvals.draws$prob_nonpos, vars_ref$prob_nonpos, tolerance = 0.01)
    expect_equal(vars_pvals.draws$p_value, vars_ref$p_value, tolerance = 0.01)

    expect_equal(vars_pvals.normal$varspec, vars_ref$varspec)
    expect_equal(vars_pvals.normal$prob_nonneg, vars_ref$prob_nonneg, tolerance = 0.01)
    expect_equal(vars_pvals.normal$prob_nonpos, vars_ref$prob_nonpos, tolerance = 0.01)
    expect_equal(vars_pvals.normal$p_value, vars_ref$p_value, tolerance = 0.01)
  })}
})
