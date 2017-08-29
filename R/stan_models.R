require(rstan)

stan_models_path <- file.path(base_scripts_path, "msglm/inst/stan_models") # FIXME
msglm.stan_model <- stan_model(file.path(stan_models_path, "msglm.stan"), "msglm", save_dso = TRUE, auto_write = TRUE)
msglm_local.stan_model <- stan_model(file.path(stan_models_path, "msglm_local.stan"), "msglm_local", save_dso = TRUE, auto_write = TRUE)
msglm_normalize.stan_model <- stan_model(file.path(stan_models_path , "msglm_normalize.stan"), "msglm_normalize", save_dso = TRUE, auto_write = TRUE)
