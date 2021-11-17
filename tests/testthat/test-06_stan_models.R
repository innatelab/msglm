context("stan models")

for (model in msglm:::msglm_model_names()) {
    test_that(paste0("stan model '", model, "' compiles"), {
        expect_class(msglm:::msglm_stan_model(model), "CmdStanModel")
    })
}
