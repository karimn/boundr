#!/usr/bin/Rscript

"Usage:
  test_boundr single [--different-priors --num-entities=<num-entities> --true-hyper-sd=<sd>]
  test_boundr multi <cores> <runs> [--num-entities=<num-entities> --true-hyper-sd=<sd> --constrained-prior --density-plots --output=<output-name>] [--append | [--different-priors --alt-hyper-sd=<sd>]]

Options:
  --output=<output-name>  Output name to use in file names [default: test_run]
  --num-entities=<num-entities>  Number of entities in model [default: 3]
  --true-hyper-sd=<sd>  True SD hyperparameter for prior [default: 1]
  --alt-hyper-sd=<sd>  Alternative SD hyperparameter to use for fit [default: 2.5]
" -> opt_desc

script_options <- if (interactive()) {
  docopt::docopt(opt_desc, "multi 12 3 --density-plots --num-entities=1 --output=test.rds --different-priors")
  # docopt::docopt(opt_desc, "single --different-priors --num-entities=1")
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
  modify_at(c("cores", "runs", "num-entities"), as.integer) %>%
  modify_at(c("true-hyper-sd", "alt-hyper-sd"), as.numeric)

options(mc.cores = max(1, parallel::detectCores()))
rstan_options(auto_write = TRUE)

true_discretized_beta_hyper_sd <- if (script_options$`constrained-prior`) {
  lst(default = script_options$`true-hyper-sd`,
      "migration complier" = ~ mutate(., sd = if_else(fct_match(r_m, c("always", "program defier", "wedge defier")),
                                                      0.5,
                                                      script_options$`true-hyper-sd`)))
} else script_options$`true-hyper-sd`

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
    # "treatment defier" = ~ 1 - z,
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
  entity_data <- create_prior_predicted_simulation(test_model4, sample_size = 4000, chains = 4, iter = 1000,
                                                     discrete_beta_hyper_sd = script_options$`true-hyper-sd`, discretized_beta_hyper_sd = true_discretized_beta_hyper_sd, tau_level_sigma = 1,
                                                     num_entities = script_options$`num-entities`) %>%
    unnest(entity_data) %>%
    select(entity_index, sim) %>%
    deframe()

  known_results <- entity_data %>%
    map_dfr(boundr:::get_known_estimands, test_estimands4, .id = "entity_index") %>%
    group_by_at(vars(-entity_index, -prob)) %>%
    summarize(prob = mean(prob)) %>%
    ungroup() %>%
    select(estimand_name, cutpoint, prob)

  test_sim_data <- entity_data %>%
    map_dfr(create_simulation_analysis_data, .id = "entity_index") %>%
    mutate(y = if_else(y_1 == 0, 30, -30))
    # mutate(y = if_else(y_2 == 0, 30, if_else(y_1 == 0, sample(c(-1, 1), n(), replace = TRUE), -30)))
    # mutate(y = if_else(y_2 == 0, 30, if_else(y_1 == 0, runif(n(), -19, 19), -30)))

  # test_model %>%
  #   get_linear_programming_bounds(test_sim_data, "y_1", b = 1, g = 1, z = 1, m = 0)

  test_sampler <- create_sampler(
    test_model4,
    model_levels = "entity_index",
    analysis_data = test_sim_data,
    estimands = test_estimands4,
    # y = y < -20,
    y = y,

    discrete_beta_hyper_sd = if (script_options$`different-priors`) script_options$`alt-hyper-sd` else script_options$`true-hyper-sd`,
    discretized_beta_hyper_sd = if (script_options$`different-priors`) script_options$`alt-hyper-sd` else true_discretized_beta_hyper_sd,
    # discretized_beta_hyper_sd = lst(default = 2,
    #                                 "migration complier" = ~ mutate(., sd = if_else(fct_match(r_m, c("always", "program defier", "wedge defier")),
    #                                                                                 0.1,
    #                                                                                 2))),
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

  test_prior_results <- test_prior_fit %>%
    get_estimation_results(no_sim_diag = FALSE, quants = seq(0, 1, 0.1)) %T>%
    print(n = 1000)

  test_prior_marginal_prob <- test_prior_fit %>% get_marginal_latent_type_prob()

  test_fit <- test_sampler %>%
    sampling(
      chains = 4,
      iter = 1000,
      # control = lst(adapt_delta = 0.99, max_treedepth = 12),
      pars = c("iter_estimand", "marginal_p_r"),
    )

  test_results <- test_fit %>%
    get_estimation_results(no_sim_diag = FALSE, quants = seq(0, 1, 0.1)) %T>%
    print(n = 1000)

  known_results %>%
    inner_join(test_results, by = c("cutpoint", "estimand_name")) %>%
    mutate_at(vars(starts_with("per_")), ~ . - prob) %>%
    mutate(coverage = if_else(per_0.1 > 0 | per_0.9 < 0, "outside", "inside") %>% factor()) %>%
    select(estimand_name, coverage) %>%
    print(n = 1000)

  bind_rows(prior = test_prior_results, posterior = test_results, .id = "fit_type") %>%
    filter(estimand_name == "Pr[Y^{y}_{b=0,g=0,z=0,m=0} < c | M = 1, B = 0, G = 0, Z = 0]") %>%
    select(fit_type, starts_with("per_")) %>%
    pivot_longer(names_to = "quant", cols = -fit_type, names_prefix = "per_", names_ptypes = list("quant" = numeric())) %>%
    ggplot() +
    geom_col(aes(quant, value, group = fit_type, fill = fit_type), position = position_dodge()) +
    geom_vline(xintercept = known_results %>% filter(estimand_name == "Pr[Y^{y}_{b=0,g=0,z=0,m=0} < c | M = 1, B = 0, G = 0, Z = 0]") %>% pull(prob)) +
    scale_fill_discrete("") +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(legend.position = "top")
}

# Multiple Runs -----------------------------------------------------------

if (script_options$multi) {
  num_runs <- script_options$runs

  test_parallel_map <- function(.x, .f, ..., cores = script_options$cores) {
    pbmcapply::pbmclapply(.x, as_mapper(.f), ..., ignore.interactive = TRUE, mc.silent = TRUE, mc.cores = cores)
  }

  test_sim_data <- create_prior_predicted_simulation(test_model4, sample_size = 4000, chains = 4, iter = 1000,
                                                     discrete_beta_hyper_sd = script_options$`true-hyper-sd`, discretized_beta_hyper_sd = true_discretized_beta_hyper_sd, tau_level_sigma = 1,
                                                     num_entities = script_options$`num-entities`, num_sim = num_runs) %>%
    deframe()

  test_run_data <- test_sim_data %>%
    test_parallel_map(cores = script_options$cores %/% 4,
    # map(
      function(entity_data, discrete_beta_hyper_sd, discretized_beta_hyper_sd, save_iter_data) tryCatch({
        entity_data %<>% deframe()

        known_results <- entity_data %>%
          map_dfr(boundr:::get_known_estimands, test_estimands4, .id = "entity_index") %>%
          group_by_at(vars(-entity_index, -prob)) %>%
          summarize(prob = mean(prob)) %>%
          ungroup()

        sampler <- entity_data %>%
          map_dfr(create_simulation_analysis_data, .id = "entity_index") %>%
          # mutate(y = if_else(y_2 == 0, 30, if_else(y_1 == 0, 0, -30))) %>%
          # mutate(y = if_else(y_2 == 0, 30, if_else(y_1 == 0, runif(n(), -19, 19), -30))) %>%
          mutate(y = if_else(y_1 == 0, 30, -30)) %>%
          create_sampler(
            test_model4,
            model_levels = "entity_index",
            analysis_data = .,
            estimands = test_estimands4,
            y = y,

            discrete_beta_hyper_sd = discrete_beta_hyper_sd,
            discretized_beta_hyper_sd = discretized_beta_hyper_sd,
            tau_level_sigma = 1,
            calculate_marginal_prob = FALSE
          )

        sampler %>%
          sampling(
            pars = "iter_estimand",
            chains = 4, iter = 1000
          ) %>%
          get_estimation_results(quants = seq(0, 1, 0.1)) %>%
          select(estimand_name, cutpoint, starts_with("per_"), if (save_iter_data) "iter_data") %>%
          inner_join(
            select(known_results, estimand_name, cutpoint, prob),
            by = c("estimand_name", "cutpoint")
          ) %>%
          # mutate_at(vars(starts_with("per_")), ~ . - prob) %>%
          mutate(coverage = if_else((per_0.1 - prob) > 0 | (per_0.9 - prob) < 0, "outside", "inside") %>% factor())
      }, error = function(err) browser()),
    discrete_beta_hyper_sd = if (script_options$`different-priors`) script_options$`alt-hyper-sd` else script_options$`true-hyper-sd`,
    discretized_beta_hyper_sd = if (script_options$`different-priors`) script_options$`alt-hyper-sd` else true_discretized_beta_hyper_sd,
    save_iter_data = script_options$`density-plots`
    ) %>%
    compact() %>%
    bind_rows(.id = "iter_id")

  dummy_data <- test_sim_data[[1]] %>%
    select(entity_index, sim) %>%
    deframe() %>%
    map_dfr(create_simulation_analysis_data, .id = "entity_index") %>%
    # mutate(y = if_else(y_2 == 0, 30, if_else(y_1 == 0, runif(n(), -19, 19), -30))) %>%
    mutate(y = if_else(y_1 == 0, 30, -30)) # This data isn't really used in prior prediction

  true_prior_sampler <- create_sampler(
    test_model4,
    model_levels = "entity_index",
    analysis_data = dummy_data,
    estimands = test_estimands4,
    y = y,

    discrete_beta_hyper_sd = script_options$`true-hyper-sd`,
    discretized_beta_hyper_sd = true_discretized_beta_hyper_sd,
    tau_level_sigma = 1,
    calculate_marginal_prob = FALSE
  )

  true_prior_fit <- true_prior_sampler %>%
    sampling(
      pars = "iter_estimand",
      chains = 4, iter = 1000,
      run_type = "prior-predict",
    )

  true_prior_results <- true_prior_fit %>% get_estimation_results(quants = seq(0, 1, 0.1))
  prior_results <- NULL

  if (script_options$`different-priors`) {
    prior_sampler <- create_sampler(
      test_model4,
      model_levels = "entity_index",
      analysis_data = dummy_data,
      estimands = test_estimands4,
      y = y,

      discrete_beta_hyper_sd = script_options$`alt-hyper-sd`,
      discretized_beta_hyper_sd = script_options$`alt-hyper-sd`,
      tau_level_sigma = 1,
      calculate_marginal_prob = FALSE
    )

    prior_fit <- prior_sampler %>%
      sampling(
        pars = "iter_estimand",
        chains = 4, iter = 1000,
        run_type = "prior-predict",
      )

    prior_results <- prior_fit %>% get_estimation_results(quants = seq(0, 1, 0.1))

    test_run_data %>%
      filter(estimand_name %in% c("Pr[Y^{y}_{b=0,g=0,z=0,m=0} < c | M = 1, B = 0, G = 0, Z = 0]",
                                  "Pr[Y^{y}_{b=0,g=0,z=0,m=1} < c | M = 1, B = 0, G = 0, Z = 0]")) %>%
      pack(post = starts_with("per_")) %>%
      left_join(
        prior_results %>%
          pack(wrong_prior = starts_with("per_")),
        by = c("estimand_name", "cutpoint")
      ) %>%
      left_join(
        true_prior_results %>%
          pack(true_prior = starts_with("per_")),
        by = c("estimand_name", "cutpoint")
      ) %>%
      unpack(c(true_prior, wrong_prior, post), names_sep = "_") %>%
      mutate_at(vars(contains("prior_per_")), ~ . - prob) %>%
      mutate(wrong_prior_coverage = if_else(wrong_prior_per_0.1 > 0 | wrong_prior_per_0.9 < 0, "outside", "inside") %>% factor(),
             true_prior_coverage = if_else(true_prior_per_0.1 > 0 | true_prior_per_0.9 < 0, "outside", "inside") %>% factor()) %>%
      group_by_at(vars(estimand_name, any_of("cutpoint"))) %>%
      summarize_at(vars(ends_with("coverage")), ~ mean(fct_match(., "inside"))) %>%
      ungroup() %>%
      arrange_at(vars(estimand_name, any_of("cutpoint"))) %>%
      print(n = 1000, width = 160)
  } else {
    test_run_data_file <- file.path("temp-data", str_c(script_options$output, ".rds"))

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

  if (script_options$`density-plots`) {
    density_plots <- bind_rows(prior = if (script_options$`different-priors`) {
                prior_results %>%
                  filter(estimand_name == "Pr[Y^{y}_{b=0,g=0,z=0,m=0} < c | M = 1, B = 0, G = 0, Z = 0]") %>%
                  slice(rep(1, num_runs)) %>%
                  mutate(iter_id = seq(num_runs))
              },
              posterior = test_run_data %>%
                filter(estimand_name == "Pr[Y^{y}_{b=0,g=0,z=0,m=0} < c | M = 1, B = 0, G = 0, Z = 0]") %>%
                mutate(iter_id = as.integer(iter_id)),
              "true prior" = true_prior_results %>%
                filter(estimand_name == "Pr[Y^{y}_{b=0,g=0,z=0,m=0} < c | M = 1, B = 0, G = 0, Z = 0]") %>%
                slice(rep(1, num_runs)) %>%
                mutate(iter_id = seq(num_runs)),
              .id = "fit_type") %>%
      select(iter_id, prob, fit_type, iter_data) %>% {
      # select(iter_id, prob, fit_type, starts_with("per_")) %>% {
        # pivot_longer(., names_to = "quant", cols = -c(iter_id, prob, fit_type), names_prefix = "per_", names_ptypes = list("quant" = numeric())) %>%
        mutate(., iter_data = map_if(iter_data, ~ !is_null(.x), select, -iter_id)) %>%
          unnest(iter_data) %>%
          ggplot() +
          geom_density(aes(iter_estimand, group = fit_type, color = fit_type)) +
          # geom_col(aes(quant, value, group = fit_type, fill = fit_type), position = position_dodge(), width = 0.05) +
          # geom_col(aes(value, quant, group = fit_type, fill = fit_type), position = position_dodge(), width = 0.05) +
          # geom_line(aes(quant, value, group = fit_type, color = fit_type)) +
          # geom_line(aes(value, quant, group = fit_type, color = fit_type)) +
          # geom_hline(aes(yintercept = prob), data = .) +
          geom_vline(aes(xintercept = prob), data = . %>% distinct(iter_id, prob)) +
          scale_fill_discrete("", aesthetics = c("color", "fill")) +
          # scale_x_continuous(breaks = seq(0, 1, 0.2)) +
          labs(x = "", y = "") +
          facet_wrap(vars(iter_id), scales = "free") +
          theme_minimal() +
          theme(legend.position = "top")
      }

    ggsave(file.path("temp-img", str_c(script_options$output, ".png")), density_plots)

    if (!interactive() && require(tcltk)) {
      x11()
      plot(density_plots)
      capture <- tk_messageBox(message = "Hit spacebar to close plots.")
    } else {
      plot(density_plots)
    }
  }
}

# Diagnostics -------------------------------------------------------------

# color_scheme_set("darkgray")
# test_posterior <- as.array(test_fit)
# test_np <- nuts_params(test_fit)
# mcmc_pairs(test_posterior, np = test_np, regex_pars =  c("level_beta_sigma", "obs_beta\\[[2-4]"))
