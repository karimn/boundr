#!/usr/bin/Rscript

"Usage:
  test_boundr single
  test_boundr multi <cores> <runs> [--append]
" -> opt_desc

script_options <- if (interactive()) {
  docopt::docopt(opt_desc, "multi 12 2 --append")
  # docopt::docopt(opt_desc, "single")
} else {
  docopt::docopt(opt_desc)
}

# Setup -------------------------------------------------------------------

library(magrittr)
library(tidyverse)
library(rlang)
library(rstan)
library(bayesplot)

library(econometr)
library(boundr)

script_options %<>%
  modify_at(c("cores", "runs"), as.integer)

options(mc.cores = max(1, parallel::detectCores()))
rstan_options(auto_write = TRUE)

# source("util.R")
# source("sim_bounded_util.R")

# Models ------------------------------------------------------------------

discrete_variables <- list2(
  define_response(
    "b",

    "program branch" = ~ 1,
    "control branch" = ~ 0,
  ),

  define_response(
    "g",

    "treatment sector" = ~ 1,
    "control sector" = ~ 0,
  ),

  define_response(
    "z",

    "village assigned treatment" = ~ 1,
    "village assigned control" = ~ 0,
  ),

  define_response(
    "m",
    input = c("b", "g", "z"),

    "never" = ~ 0,
    "program complier" = ~ b,
    "program defier" = ~ 1 - b,
    "wedge complier" = ~ g,
    "wedge defier" = ~ 1 - g,
    "treatment complier" = ~ z,
    "treatment defier" = ~ 1 - z,
    "always" = ~ 1,
  ),
)

test_model <- define_structural_causal_model(
  !!!discrete_variables,

  define_discretized_response_group(
    "y",
    cutpoints = c(-100, -20, 20, 100),
    # cutpoints = c(-100, 20, 100),

    input = c("b", "g", "m"),

    "never below" = ~ 0,
    "program complier" = ~ b,
    "program defier" = ~ 1 - b,
    "wedge complier" = ~ g,
    "wedge defier" = ~ 1 - g,
    "migration complier" = ~ m,
    "migration defier" = ~ 1 - m,
    "always below" = ~ 1,

    pruning_data = tribble(
      ~ hi,                  ~ "always below",  ~ "program complier", ~ "wedge complier", ~ "migration complier", ~ "program defier", ~ "wedge defier", ~ "migration defier", ~ "never below",

      "always below",        TRUE,              TRUE,                 TRUE,              TRUE,                   TRUE,               TRUE,             TRUE,                TRUE,
      "program complier",    FALSE,             TRUE,                 FALSE,             FALSE,                  FALSE,              FALSE,            FALSE,               TRUE,
      "wedge complier",       FALSE,             FALSE,                TRUE,              FALSE,                  FALSE,              FALSE,            FALSE,               TRUE,
      "migration complier",  FALSE,             FALSE,                FALSE,             TRUE,                   FALSE,              FALSE,            FALSE,               TRUE,
      "program defier",      FALSE,             FALSE,                FALSE,             FALSE,                  TRUE,               FALSE,            FALSE,               TRUE,
      "wedge defier",         FALSE,             FALSE,                FALSE,             FALSE,                  FALSE,              TRUE,             FALSE,               TRUE,
      "migration defier",    FALSE,             FALSE,                FALSE,             FALSE,                  FALSE,              FALSE,            TRUE,                TRUE,
      "never below",         FALSE,             FALSE,                FALSE,             FALSE,                  FALSE,              FALSE,            FALSE,               TRUE,
    ) %>%
      pivot_longer(-hi, names_to = "low", values_to = "allow") %>%
      filter(allow) %>%
      select(-allow)
  ),

  exogenous_prob = tribble(
    ~ b, ~ g, ~ z, ~ ex_prob,
    0,   0,   0,   0.4,
    1,   0,   0,   0.2,
    1,   1,   0,   0.2,
    1,   1,   1,   0.2
  ),
)

test_model2 <- define_structural_causal_model(
  !!!discrete_variables,

  define_discretized_response_group(
    "y",
    cutpoints = c(-100, -20, -10, 10, 20, 100),
    # cutpoints = c(-100, 20, 100),

    input = c("b", "g", "m"),

    "never below" = ~ 0,
    "program complier" = ~ b,
    "program defier" = ~ 1 - b,
    "wedge complier" = ~ g,
    "wedge defier" = ~ 1 - g,
    "migration complier" = ~ m,
    "migration defier" = ~ 1 - m,
    "always below" = ~ 1,

    pruning_data = tribble(
      ~ hi,                  ~ "always below",  ~ "program complier", ~ "wedge complier", ~ "migration complier", ~ "program defier", ~ "wedge defier", ~ "migration defier", ~ "never below",

      "always below",        TRUE,              TRUE,                 TRUE,              TRUE,                   TRUE,               TRUE,             TRUE,                TRUE,
      "program complier",    FALSE,             TRUE,                 FALSE,             FALSE,                  FALSE,              FALSE,            FALSE,               TRUE,
      "wedge complier",       FALSE,             FALSE,                TRUE,              FALSE,                  FALSE,              FALSE,            FALSE,               TRUE,
      "migration complier",  FALSE,             FALSE,                FALSE,             TRUE,                   FALSE,              FALSE,            FALSE,               TRUE,
      "program defier",      FALSE,             FALSE,                FALSE,             FALSE,                  TRUE,               FALSE,            FALSE,               TRUE,
      "wedge defier",         FALSE,             FALSE,                FALSE,             FALSE,                  FALSE,              TRUE,             FALSE,               TRUE,
      "migration defier",    FALSE,             FALSE,                FALSE,             FALSE,                  FALSE,              FALSE,            TRUE,                TRUE,
      "never below",         FALSE,             FALSE,                FALSE,             FALSE,                  FALSE,              FALSE,            FALSE,               TRUE,
    ) %>%
      pivot_longer(-hi, names_to = "low", values_to = "allow") %>%
      filter(allow) %>%
      select(-allow)
  ),

  exogenous_prob = tribble(
    ~ b, ~ g, ~ z, ~ ex_prob,
    0,   0,   0,   0.4,
    1,   0,   0,   0.2,
    1,   1,   0,   0.2,
    1,   1,   1,   0.2
  ),
)

test_model3 <- define_structural_causal_model(
  !!!discrete_variables,

  define_response(
    "y",

    input = c("b", "g", "m"),

    "never below" = ~ 0,
    "program complier" = ~ b,
    "program defier" = ~ 1 - b,
    "wedge complier" = ~ g,
    "wedge defier" = ~ 1 - g,
    "migration complier" = ~ m,
    "migration defier" = ~ 1 - m,
    "always below" = ~ 1,
  ),

  exogenous_prob = tribble(
    ~ b, ~ g, ~ z, ~ ex_prob,
    0,   0,   0,   0.4,
    1,   0,   0,   0.2,
    1,   1,   0,   0.2,
    1,   1,   1,   0.2
  ),
)

test_model4 <- define_structural_causal_model(
  !!!discrete_variables,

  define_discretized_response_group(
    "y",
    cutpoints = c(-100, -20, 100),

    input = c("b", "g", "m"),

    "never below" = ~ 0,
    "program complier" = ~ b,
    "program defier" = ~ 1 - b,
    "wedge complier" = ~ g,
    "wedge defier" = ~ 1 - g,
    "migration complier" = ~ m,
    "migration defier" = ~ 1 - m,
    "always below" = ~ 1,

    pruning_data = tribble(
      ~ hi,                  ~ "always below",  ~ "program complier", ~ "wedge complier", ~ "migration complier", ~ "program defier", ~ "wedge defier", ~ "migration defier", ~ "never below",

      "always below",        TRUE,              TRUE,                 TRUE,              TRUE,                   TRUE,               TRUE,             TRUE,                TRUE,
      "program complier",    FALSE,             TRUE,                 FALSE,             FALSE,                  FALSE,              FALSE,            FALSE,               TRUE,
      "wedge complier",       FALSE,             FALSE,                TRUE,              FALSE,                  FALSE,              FALSE,            FALSE,               TRUE,
      "migration complier",  FALSE,             FALSE,                FALSE,             TRUE,                   FALSE,              FALSE,            FALSE,               TRUE,
      "program defier",      FALSE,             FALSE,                FALSE,             FALSE,                  TRUE,               FALSE,            FALSE,               TRUE,
      "wedge defier",         FALSE,             FALSE,                FALSE,             FALSE,                  FALSE,              TRUE,             FALSE,               TRUE,
      "migration defier",    FALSE,             FALSE,                FALSE,             FALSE,                  FALSE,              FALSE,            TRUE,                TRUE,
      "never below",         FALSE,             FALSE,                FALSE,             FALSE,                  FALSE,              FALSE,            FALSE,               TRUE,
    ) %>%
      pivot_longer(-hi, names_to = "low", values_to = "allow") %>%
      filter(allow) %>%
      select(-allow)
  ),

  exogenous_prob = tribble(
    ~ b, ~ g, ~ z, ~ ex_prob,
    0,   0,   0,   0.4,
    1,   0,   0,   0.2,
    1,   1,   0,   0.2,
    1,   1,   1,   0.2
  ),
)

# Estimands ---------------------------------------------------------------

with_discretized_estimands <- list2(
  build_diff_estimand(
    build_atom_estimand("m", b = 1, g = 1, z = 1),
    build_atom_estimand("m", b = 0, g = 0, z = 0)
  ),

  build_discretized_diff_estimand(
    build_discretized_atom_estimand("y", b = 0, g = 0, z = 0, m = 1),
    build_discretized_atom_estimand("y", b = 0, g = 0, z = 0, m = 0)
  ),

  build_discretized_diff_estimand(
    build_discretized_atom_estimand("y", b = 1, g = 1, z = 1),
    build_discretized_atom_estimand("y", b = 0, g = 0, z = 0)
  ),

  build_discretized_diff_estimand(
    build_discretized_atom_estimand("y", b = 0, g = 0, z = 0, m = 1, cond = m == 1 & b == 0 & g == 0 & z == 0),
    build_discretized_atom_estimand("y", b = 0, g = 0, z = 0, m = 0, cond = m == 1 & b == 0 & g == 0 & z == 0)
  ),
)


test_estimands <- build_estimand_collection(
  model = test_model,
  utility = c(0, 1, 1.5),

  !!!with_discretized_estimands
)

test_estimands2 <- build_estimand_collection(
  model = test_model2,
  utility = c(0, 1, 1.5, 1.75, 1.8),

  !!!with_discretized_estimands
)

test_estimands3 <- build_estimand_collection(
  model = test_model3,

  build_diff_estimand(
    build_atom_estimand("m", b = 1, g = 1, z = 1),
    build_atom_estimand("m", b = 0, g = 0, z = 0)
  ),

  build_diff_estimand(
    build_atom_estimand("y", b = 0, g = 0, z = 0, m = 1),
    build_atom_estimand("y", b = 0, g = 0, z = 0, m = 0)
  ),

  build_diff_estimand(
    build_atom_estimand("y", b = 1, g = 1, z = 1),
    build_atom_estimand("y", b = 0, g = 0, z = 0)
  ),

  build_diff_estimand(
    build_atom_estimand("y", b = 0, g = 0, z = 0, m = 1, cond = m == 1 & b == 0 & g == 0 & z == 0),
    build_atom_estimand("y", b = 0, g = 0, z = 0, m = 0, cond = m == 1 & b == 0 & g == 0 & z == 0)
  ),
)

test_estimands4 <- build_estimand_collection(
  model = test_model4,
  utility = c(0, 1),

  !!!with_discretized_estimands
)

# Single Run --------------------------------------------------------------

if (script_options$single) {
  entity_data <- create_prior_predicted_simulation(test_model, sample_size = 4000, chains = 4, iter = 1000,
                                                     discrete_beta_hyper_sd = 2, discretized_beta_hyper_sd = 2, tau_level_sigma = 1,
                                                     num_entities = 3) %>%
    unnest(entity_data) %>%
    select(entity_index, sim) %>%
    deframe()

  known_results <- entity_data %>%
    map_dfr(boundr:::get_known_estimands, test_estimands, .id = "entity_index") %>%
    group_by_at(vars(-entity_index, -prob)) %>%
    summarize(prob = mean(prob)) %>%
    ungroup() %>%
    select(estimand_name, cutpoint, prob)

  test_sim_data <- entity_data %>%
    map_dfr(create_simulation_analysis_data, .id = "entity_index") %>%
    # mutate(y = if_else(y_1 == 0, 30, -30))
    # mutate(y = if_else(y_2 == 0, 30, if_else(y_1 == 0, sample(c(-1, 1), n(), replace = TRUE), -30)))
    mutate(y = if_else(y_2 == 0, 30, if_else(y_1 == 0, runif(n(), -19, 19), -30)))

  # test_model %>%
  #   get_linear_programming_bounds(test_sim_data, "y_1", b = 1, g = 1, z = 1, m = 0)

  test_sampler <- create_sampler(
    test_model4,
    model_levels = "entity_index",
    analysis_data = test_sim_data,
    estimands = test_estimands4,
    # y = y < -20,
    y = y,

    discrete_beta_hyper_sd = 2,
    discretized_beta_hyper_sd = 2,
    tau_level_sigma = 1,
    calculate_marginal_prob = TRUE
  )

  test_prior_fit <- test_sampler %>%
    sampling(
      chains = 4,
      iter = 1000,
      # control = lst(adapt_delta = 0.99, max_treedepth = 12),
      pars = c("iter_estimand", "marginal_p_r"),
      run_type = "prior-predict",
    )

  test_fit <- test_sampler %>%
    sampling(
      chains = 4,
      iter = 1000,
      # control = lst(adapt_delta = 0.99, max_treedepth = 12),
      pars = c("iter_estimand", "marginal_p_r"),
    )

  test_prior_results <- test_prior_fit %>%
    get_estimation_results(no_sim_diag = FALSE) %T>%
    print(n = 1000)

  test_results <- test_fit %>%
    get_estimation_results(no_sim_diag = FALSE) %T>%
    print(n = 1000)

  known_results %>%
    inner_join(test_results, by = c("cutpoint", "estimand_name")) %>%
    mutate_at(vars(starts_with("per_")), ~ . - prob) %>%
    mutate(coverage = if_else(per_0.1 > 0 | per_0.9 < 0, "outside", "inside") %>% factor()) %>%
    select(estimand_name, coverage) %>%
    print(n = 1000)

  test_prior_marginal_prob <- test_prior_fit %>% get_marginal_latent_type_prob()
}

# Multiple Runs -----------------------------------------------------------

if (script_options$multi) {
  num_runs <- script_options$runs

  test_parallel_map <- function(.x, .f, ..., cores = script_options$cores) {
    pbmcapply::pbmclapply(.x, as_mapper(.f), ..., ignore.interactive = TRUE, mc.silent = TRUE, mc.cores = cores)
  }

  test_run_data <- create_prior_predicted_simulation(test_model, sample_size = 4000, chains = 4, iter = 1000,
                                                     discrete_beta_hyper_sd = 2, discretized_beta_hyper_sd = 2, tau_level_sigma = 1,
                                                     num_entities = 3, num_sim = num_runs) %>%
    deframe() %>%
    test_parallel_map(cores = script_options$cores %/% 4,
    # map(
      function(entity_data) tryCatch({
        browser()
        entity_data %<>% deframe()

        known_results <- entity_data %>%
          map_dfr(boundr:::get_known_estimands, test_estimands, .id = "entity_index") %>%
          group_by_at(vars(-entity_index, -prob)) %>%
          summarize(prob = mean(prob)) %>%
          ungroup()

        entity_data %>%
          map_dfr(create_simulation_analysis_data, .id = "entity_index") %>%
          # mutate(y = if_else(y_2 == 0, 30, if_else(y_1 == 0, 0, -30))) %>%
          mutate(y = if_else(y_2 == 0, 30, if_else(y_1 == 0, runif(n(), -19, 19), -30))) %>%
          create_sampler(
            test_model,
            model_levels = "entity_index",
            analysis_data = .,
            estimands = test_estimands,
            y = y,

            discrete_beta_hyper_sd = 2,
            discretized_beta_hyper_sd = 2,
            tau_level_sigma = 1,
            calculate_marginal_prob = FALSE
          ) %>%
          sampling(
            pars = "iter_estimand",
            chains = 4, iter = 1000
          ) %>%
          get_estimation_results() %>%
          select(estimand_name, cutpoint, starts_with("per_")) %>% # rhat, starts_with("ess"),
          inner_join(
            select(known_results, estimand_name, cutpoint, prob),
            by = c("estimand_name", "cutpoint")
          ) %>%
          mutate_at(vars(starts_with("per_")), ~ . - prob) %>%
          mutate(coverage = if_else(per_0.1 > 0 | per_0.9 < 0, "outside", "inside") %>% factor())
      }, error = function(err) browser())
    ) %>%
    compact() %>%
    bind_rows(.id = "iter_id")

  test_run_data_file <- file.path("temp-data", "test_run3.rds")

  if (script_options$append && file.exists(test_run_data_file)) {
    test_run_data %<>%
      bind_rows(read_rds(test_run_data_file))
  }

  write_rds(test_run_data, test_run_data_file)

  test_run_data %>%
    group_by_at(vars(estimand_name, any_of("cutpoint"))) %>%
    summarize(coverage = mean(fct_match(coverage, "inside"))) %>%
    ungroup() %>%
    arrange_at(vars(estimand_name, any_of("cutpoint"))) %>%
    print(n = 1000)
}

# Manual Run --------------------------------------------------------------


# Diagnostics -------------------------------------------------------------

# color_scheme_set("darkgray")
# test_posterior <- as.array(test_fit)
# test_np <- nuts_params(test_fit)
# mcmc_pairs(test_posterior, np = test_np, regex_pars =  c("level_beta_sigma", "obs_beta\\[[2-4]"))
