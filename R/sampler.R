#' @include util.R

#' S4 class for MCMC sampling
#'
#' @slot stan_data list.
#' @slot analysis_data data.frame.
#' @slot model_levels character.
#' @slot estimand_levels character.
#' @slot between_entity_diff_levels character.
#' @slot cv_level factor.
#' @slot unique_entity_ids data.frame.
#' @slot endogenous_latent_type_variables data.frame.
#' @slot estimands EstimandCollection.
#'
#' @export
setClass("Sampler",
         contains = "stanmodel",
         slots = c(stan_data = "list",
                   analysis_data = "data.frame",
                   model_levels = "character",
                   estimand_levels = "character",
                   between_entity_diff_levels = "character",
                   cv_level = "factor",
                   unique_entity_ids = "data.frame",
                   endogenous_latent_type_variables = "data.frame",
                   estimands = "EstimandCollection"))

#' S4 class for sampling results
#'
#' @slot sampler S4 \code{Sampler} that was used to produce this results object.
#'
#' @export
setClass("SamplingResults",
         contains = "stanfit",
         slots = c("sampler" = "Sampler"))

# Execute MCMC Sampling
#'
#' @importMethodsFrom rstan sampling
#' @export
setMethod("sampling", "Sampler", function(object, ..., run_type = c("fit", "prior-predict"), save_background_joint_prob = FALSE) {
  run_type <- arg_match(run_type) %>%
    factor(levels = c("prior-predict", "fit"))

  args <- list2(...)

  if ("data" %in% names(args)) {
     stop("Sample data cannot be specified. Data is prepared in the sampler constructor.")
  }

  pars <- c("iter_estimand", "iter_level_entity_estimand", "iter_level_entity_estimand_sd", "log_lik", "iter_between_level_entity_diff_estimand", "marginal_p_r")

  if (save_background_joint_prob) {
    pars %<>%  c("r_prob")
  }

  initializer <- function(chain_id) {
    num_discrete_types <- object@stan_data$num_discrete_r_types
    num_discretized_types <- object@stan_data$num_discretized_r_types

    list2(
      toplevel_discrete_p_r = c(MCMCpack::rdirichlet(1, rep(1, num_discrete_types))),
      toplevel_discretized_p_r = if (num_discretized_types > 0) t(MCMCpack::rdirichlet(num_discrete_types, rep(1, num_discretized_types))) else array(0, dim = c(0, num_discrete_types)),
    )
  }

  args <- lst(
    data = object@stan_data %>%
      list_modify(run_type = as.integer(run_type)),
    include = TRUE,
    init = initializer,
    pars = pars,
  ) %>%
    list_modify(!!!args) # Allow some arguments to be overridden, such as "pars"

  # callNextMethod(object, ..., data = object@stan_data)
  results <- exec(rstan::sampling, as(object, "stanmodel"), !!!args) %>%
    as("SamplingResults")

  results@sampler <- object

  return(results)
})

#' Title
#'
#' @importMethodsFrom rstan vb
#' @export
setMethod("vb", "Sampler", function(object, ..., run_type = c("fit", "prior-predict"), save_background_joint_prob = FALSE) {
  run_type <- arg_match(run_type) %>%
    factor(levels = c("prior-predict", "fit"))

  args <- list2(...)

  if ("data" %in% names(args)) {
     stop("Sample data cannot be specified. Data is prepared in the sampler constructor.")
  }

  pars <- c("iter_estimand", "iter_level_entity_estimand", "log_lik")

  if (save_background_joint_prob) {
    pars %<>%  c("r_prob")
  }

  args <- lst(
    data = object@stan_data %>%
      list_modify(run_type = as.integer(run_type)),
    include = TRUE,
    pars = pars,
  ) %>%
    list_modify(!!!args) # Allow some arguments to be overridden, such as "pars"

  # callNextMethod(object, ..., data = object@stan_data)
  results <- exec(rstan::vb, as(object, "stanmodel"), !!!args) %>%
    as("SamplingResults")

  results@sampler <- object

  return(results)
})

setGeneric("get_sampler", function(r) {
  standardGeneric("get_sampler")
})

setMethod("get_sampler", "SamplingResults", function(r) r@sampler)

setGeneric("get_estimation_results", function(r, no_levels = FALSE, no_sim_diag = TRUE, level_hist = FALSE) {
  standardGeneric("get_estimation_results")
})

#' Prepare estimation results
#'
#' @param r \code{SamplingResults}
#' @param no_levels Estimate only at the top-level
#' @param no_sim_diag Do not calculate simulation diagnostics (Rhat and ESS diagnostics)
#' @param level_hist
#'
#' @return Nested \code{tibble} with estimation results
#' @export
setMethod("get_estimation_results", "SamplingResults", function(r, no_levels = FALSE, no_sim_diag = TRUE, level_hist = FALSE) {
  estimands <- r@sampler@estimands

  if (is_null(estimands)) {
    return(NULL)
  } else {
    results <- if (no_levels || any(is.na(r@sampler@model_levels))) {
      estimands%>%
        extract_from_fit(r, no_sim_diag = no_sim_diag)
    } else {
      between_entity_diff_info <- if (!is_null(r@sampler@between_entity_diff_levels)) {
        r@sampler@analysis_data %>%
          select(all_of(r@sampler@between_entity_diff_levels)) %>%
          map(base::levels) %>%
          imap_dfr(~ tibble(level = ..2, left = ..1[-1], right = ..1[1], estimand_id = list(seq(..3))), num_estimands(estimands, "atom")) %>% {
            if (is_empty(.)) {
              return(NULL)
            } else {
              mutate(., diff_index = seq(n())) %>%
                unnest(estimand_id)
            }
          }
      }

      estimands %>%
        extract_from_fit(r, levels = r@sampler@estimand_levels, unique_ids = r@sampler@unique_entity_ids, between_entity_diff_info = between_entity_diff_info, no_sim_diag = no_sim_diag)

    }

    if ("level_estimands" %in% names(results) && level_hist) {
      results %<>%
        mutate(
          level_hist = map(level_estimands, unnest, iter_data) %>%
            map(group_by, level, iter_id) %>%
            map2(est_type,
                 ~ group_modify(
                   all_of(.x),
                   function(.data, key) {
                     breaks <- if (fct_match(.y, "atom")) {
                       seq(0, 1, 0.1)
                     } else if (fct_match(.y, "diff")) {
                       seq(-1, 1, 0.1)
                     } else 10

                     hist(.data$iter_estimand,
                          breaks = breaks,
                          include.lowest = TRUE,
                          right = TRUE,
                          plot = FALSE) %>%
                       magrittr::extract(c("breaks", "counts", "density")) %>%
                       modify_at("breaks", ~ .[-1]) %>%
                       as_tibble()
                   })
                 ) %>%
            map(ungroup) %>%
            map(group_by, level, iter_id) %>%
            map(mutate, freq = counts / sum(counts)) %>%
            map(ungroup) %>%
            map(nest, break_data = c(iter_id, counts:freq)) %>%
            map(mutate, quant = map(break_data, quantilize_est, freq, quant_probs = c(0.1, 0.25, 0.5, 0.75, 0.9))) %>%
            map(select, -break_data) %>%
            map(unnest, quant)
        )

    }

    results %>%
      select(estimand_name, any_of(c("cutpoint", "level_estimands", "between_entity_estimands", "rhat")), starts_with("ess"), starts_with("per_"), mean, iter_data) %>%
      select_if(~ !all(is.na(.x)))
  }
})

setGeneric("get_level_estimand_sd", function(est) {
  standardGeneric("get_level_estimand_sd")
})

setMethod("get_level_estimand_sd", "SamplingResults", function(est) {
  tryCatch(
    as.array(est, par = "iter_level_entity_estimand_sd"),
    error = function(err) {
      stop("Failed to find iter_level_entity_estimand_sd parameter.")

      return(NULL)
    }) %>%
    plyr::adply(3, diagnose, no_sim_diag = TRUE) %>%
    tidyr::extract(parameters, c("estimand_id", "level_index"), "(\\d+),(\\d+)", convert = TRUE) %>%
    mutate(iter_data = map(iter_data, ~ tibble(iter_estimand = c(.), iter_id = seq(NROW(.) * NCOL(.))))) %>%
    summ_iter_data() %>%
    full_join(est@sampler@estimands@estimands, ., by = c("estimand_id")) %>%
    left_join(tibble(level = est@sampler@estimand_levels) %>% mutate(level_index = seq(n())), by = "level_index")
})

setGeneric("get_marginal_latent_type_prob", function(r, no_sim_diag = TRUE) {
  standardGeneric("get_marginal_latent_type_prob")
})

#' Extract latent type marginal probabilities
#'
#' @param r \code{SamplingResults}
#' @param no_sim_diag Do not generate simulation diagnostics (Rhat and ESS diagnostics)
#'
#' @return \code{tibble} with marginal probabilities
#' @export
setMethod("get_marginal_latent_type_prob", "SamplingResults", function(r, no_sim_diag = TRUE) {
  r %>%
    as.array(par = "marginal_p_r") %>%
    plyr::adply(3, diagnose, no_sim_diag = no_sim_diag) %>%
    tidyr::extract(parameters, "marginal_latent_type_index", "(\\d+)", convert = TRUE) %>%
    mutate(iter_data = map(iter_data, ~ tibble(iter_p_r = c(.), iter_id = seq(NROW(.) * NCOL(.))))) %>%
    full_join(r@sampler@endogenous_latent_type_variables, ., by = "marginal_latent_type_index") %>%
    mutate(estimand_quantiles = map(iter_data, quantilize_est, iter_p_r, wide = TRUE, quant_probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)),
           mean = map_dbl(iter_data, ~ mean(.$iter_p_r))) %>%
    unnest(estimand_quantiles)
})

# Hint for using importMethodsFrom: add @import to boundr-package.R first then write method. Otherwise NAMESPACE won't be updated with
# the needed generic functions.

#' Approximate leave-one-out cross-validation
#'
#' @param x S4 \code{SamplingResults} object.
#' @param ... Ignored
#' @param save_psis Save intermediate \code{psis} results.
#' @param cores Number of cores to use for parallelization.
#'
#' @return \code{loo} object
#'
#' @export
setMethod("loo", "SamplingResults", function(x, ..., save_psis = FALSE, cores = getOption("mc.cores", 1)) {
  ll <- loo::extract_log_lik(x, parameter_name = "log_lik", merge_chains = FALSE)
  r_eff <- loo::relative_eff(exp(ll), cores = cores)
  loo::loo.array(ll, r_eff = r_eff, cores = cores, save_psis = save_psis)
})

setGeneric("get_latent_type_prob", function(r, ...) standardGeneric("get_latent_type_prob"))

#' Extract sampled joint latent type probabilities
#'
#' @param r S4 \code{SamplingResults} object.
#' @param ... Ignored
#' @param no_sim_diag Do not generate sampling diagnostics [default: TRUE].
#'
#' @return Nested tibble
#' @export
setMethod("get_latent_type_prob", "SamplingResults", function(r, ..., no_sim_diag = TRUE) {
  r %>%
    as.array(par = "r_prob") %>%
    plyr::adply(3, diagnose, no_sim_diag = no_sim_diag) %>%
    tidyr::extract(parameters, "latent_type_index", "(\\d+)", convert = TRUE) %>%
    mutate(
      iter_data = map(iter_data, ~ tibble(iter_r_prob = c(.), iter_id = seq(NROW(.) * NCOL(.)))),
      latent_type_index = rep(seq(r@sampler@stan_data$num_r_types), r@sampler@stan_data$num_unique_entities),
      unique_entity_id = rep(seq(r@sampler@stan_data$num_unique_entities), each = r@sampler@stan_data$num_r_types)
    ) %>%
    left_join(mutate(r@sampler@unique_entity_ids, unique_entity_id = seq(n())), by = "unique_entity_id") %>%
    as_tibble()
})
