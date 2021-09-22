context("exp_design")

test_that("msglm_model()", {
  simplest_model <- msglm_model(constant_matrix(1, list(condition = "mock", effect = c())),
                                conditions = data.frame(condition = "mock"),
                                effects = data.frame(effect = character()))
  checkmate::expect_class(simplest_model, "msglm_model")
  expect_equal(simplest_model$modelobject, simplest_model$quantobject)

  twoeffects_model <- msglm_model(constant_matrix(0, list(condition = "mock", effect = c("a", "b"))),
                                                    conditions = data.frame(condition = "mock"),
                                                  effects = data.frame(effect = c("a", "b")))
  checkmate::expect_class(twoeffects_model, "msglm_model")

  expect_error(msglm_model(constant_matrix(1, list(condition = "mock", effect = c("a", "a"))),
                            conditions = data.frame(condition = "mock"),
                            effects = data.frame(effect = c("a", "a"))),  "Must have unique colnames, but element 2 is duplicated")

  twoconditions_model <- msglm_model(constant_matrix(1, list(condition = c("mock", "b"), effect = "b")),
                                  conditions = data.frame(condition = c("mock", "b")),
                                  effects = data.frame(effect = "b"))
  checkmate::expect_class(twoconditions_model, "msglm_model")

  expect_equal(twoconditions_model$effects$prior_mean, 0.0)
  expect_equal(c(twoconditions_model$effects$prior_tau, twoconditions_model$effects$prior_df1, twoconditions_model$effects$prior_df2), rep(1.0, 3))

})


test_that("set_batch_effects() correctly checks the input parameters", {
  expect_error(set_batch_effects(list()), "Must inherit from class 'msglm_model'")

  simplest_model <- msglm_model(constant_matrix(1, list(condition = "mock", effect = c())),
                                conditions = data.frame(condition = "mock"),
                                effects = data.frame(effect = character()))
  simplest_model <- set_batch_effects(simplest_model, constant_matrix(1, list(msrun = "mock_1", batch_effect = "contamination")))
  checkmate::expect_class(simplest_model, "msglm_model")
})
