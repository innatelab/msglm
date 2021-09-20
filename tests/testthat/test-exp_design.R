context("exp_design")

test_that("msglm_model()", {
  simplest_model <- msglm_model(constant_matrix(1, list(condition = "mock", effect = c())),
                                conditions = data.frame(condition = "mock"),
                                effects = data.frame(effect = character()))
  checkmate::expect_class(simplest_model, "msglm_model")
})


test_that("set_batch_effects() correctly checks the input parameters", {
  expect_error(set_batch_effects(list()), "Must inherit from class 'msglm_model'")

  simplest_model <- msglm_model(constant_matrix(1, list(condition = "mock", effect = c())),
                                conditions = data.frame(condition = "mock"),
                                effects = data.frame(effect = character()))
  simplest_model <- set_batch_effects(simplest_model, constant_matrix(1, list(msrun = "mock_1", batch_effect = "contamination")))
  checkmate::expect_class(simplest_model, "msglm_model")
})
