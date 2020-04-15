#' S4 class for structural models
#'
#' @slot types_data data.frame.
#' @slot responses list.
#' @slot discretized_responses list.
#' @slot exogenous_variables character.
#' @slot exogenous_prob data.frame.
#' @slot candidate_groups data.frame.
#' @slot endogenous_latent_type_variables data.frame.
#' @slot nester function.
#'
#' @export
setClass(
  "StructuralCausalModel",
  slots = c(types_data = "data.frame",
            responses = "list",
            discretized_responses = "list",
            exogenous_variables = "character",
            exogenous_prob = "data.frame",
            candidate_groups = "data.frame",
            endogenous_latent_type_variables = "data.frame",
            nester = "function"),
  contains = "BaseModel"
)

#' S4 class used to generate synthetic data for testing purposes
#'
#' @export
setClass(
  "StructuralCausalModelSimulation",
  contains = "StructuralCausalModel"
)

setGeneric("get_discretized_pruning_data", function(r) standardGeneric("get_discretized_pruning_data"))

setMethod("get_discretized_pruning_data", "StructuralCausalModel", function(r) {
  if (length(r@discretized_responses) > 0) {
    first(r@discretized_responses)@pruning_data
  } else NULL
})

setGeneric("num_responses", function(r) {
  standardGeneric("num_responses")
})

setMethod("num_responses", "StructuralCausalModel", function(r) {
  length(r@responses)
})

setMethod("get_responses", "StructuralCausalModel", function(r) {
  map(r@responses, get_responses) %>%
    flatten() %>%
    purrr::compact()
})

setGeneric("create_simulation_analysis_data", function(r) standardGeneric("create_simulation_analysis_data"))

#' Use a simulation run's data to use as analysis data
#'
#' @param r S4 \code{StructuralCausalModelSimulation} object
#'
#' @return A \code{tibble} data set.
#' @export
setMethod("create_simulation_analysis_data", "StructuralCausalModelSimulation", function(r) {
  discretized_response_names <- names(r@discretized_responses)
  discrete_response_names <- names(purrr::discard(r@responses, ~ is(., "DiscretizedResponse")))

  r@types_data %>%
    unnest(outcomes) %>%
    select(all_of(discrete_response_names), matches(str_interp("${discretized_response_names}_\\d")), num_obs) %>%
    map(~ rep(.x, .y), num_obs = .$num_obs) %>%
    as_tibble()
})

setMethod("get_discretized_response_info", "StructuralCausalModel", function(r, name) {
  r@discretized_responses[[name]] %>% get_discretized_response_info()
})

setMethod("get_discretized_cutpoints", "StructuralCausalModel", function(r, name = 1) {
  r@discretized_responses[[name]] %>% get_discretized_cutpoints()
})

#' Get names of discretized variables
#'
#' Continuous variables are discretized over a set of cutpoints. This method returns the names of the discrete variables that correspond to these cutpoints.
#'
#' @param r S4 object for structural causal model.
#'
#' @return vector of names.
#' @export
setMethod("get_discretized_variable_names", "StructuralCausalModel", function(r, name) {
  r@discretized_responses[[name]] %>% get_discretized_variable_names()
})

setMethod("discretize_continuous_variables", "StructuralCausalModel", function(r, ...) {
  continuous_var <- list2(...) %>%
    magrittr::extract(names(r@discretized_responses)) # Ignore arguments not corresponding to discretized variables

  imap(continuous_var, function(var_val, var_name) {
    r@discretized_responses[[var_name]] %>%
      discretize_continuous_variables(var_val)
  }) %>%
    flatten()
})

setMethod("set_obs_outcomes", "StructuralCausalModel", function (r, ..., cond = TRUE) {
  intervention <- list2(...)

  was_nested <- FALSE

  if ("outcomes" %in% names(r@types_data)) {
    r@types_data %<>% unnest(outcomes)
    was_nested <- TRUE
  }

  # Prune out based on conditional (abduction)
  r@types_data %<>%
    filter(!!cond)

  # Apply intervention
  if (!is_empty(intervention)) {
    r@types_data %<>% mutate(!!!intervention)
  }

  leaf_responses <- r@responses

  leaf_responses %>%
    purrr::discard(~ get_output_variable_name(.x) %in% names(intervention)) %>%
    purrr::iwalk(function(response, r_type_name) {
      input_arg <- get_input_variable_names(response)

      output_val <- if (!is_null(input_arg)) {
        response_input <- select(r@types_data, all_of(input_arg)) %>%
          mutate_if(is.factor, as.character) %>%
          as.list()

        set_obs_outcomes(response,
                         curr_r_type = r@types_data %>% pull(str_c("r_", r_type_name)),
                         !!!response_input)
      } else {
        set_obs_outcomes(response, curr_r_type = r@types_data %>% pull(str_c("r_", r_type_name)))
      }

      r@types_data <<- r@types_data %>%
        mutate(!!sym(get_output_variable_name(response)) := output_val)
    })

  if (was_nested) {
    r@types_data %<>%
      r@nester()
  }

  return(r)
})

setGeneric("get_ex_prob", function(r) standardGeneric("get_ex_prob"))

setMethod("get_ex_prob", "StructuralCausalModel", function(r) r@exogenous_prob$ex_prob)

setGeneric("num_endogenous_types", function(r) standardGeneric("num_endogenous_types"))

setMethod("num_endogenous_types", "StructuralCausalModel", function(r) nrow(r@types_data))

setGeneric("num_discrete_types", function(r) standardGeneric("num_discrete_types"))

setMethod("num_discrete_types", "StructuralCausalModel", function(r) max(r@types_data$discrete_r_type_id))

setGeneric("num_discretized_types", function(r) standardGeneric("num_discretized_types"))

setMethod("num_discretized_types", "StructuralCausalModel", function(r) {
  if (!is_empty(r@discretized_responses)) {
    return(length(r@discretized_responses[[1]]@child_responses[[1]]@finite_states))
  } else return(0)
})

setGeneric("get_endogenous_responses", function(model) standardGeneric("get_endogenous_responses"))

setMethod("get_endogenous_responses", "StructuralCausalModel", function(model) model@responses %>% magrittr::extract(setdiff(names(.), model@exogenous_variables)))

setMethod("get_candidates", "StructuralCausalModel", function(r, analysis_data = NULL) {
  if (missing(analysis_data)) {
    analysis_data <- if ("outcomes" %in% names(r@types_data)) {
      unnest(r@types_data, outcomes)
    } else {
      r@types_data
    }
  }

  get_endogenous_responses(r) %>%
    map(get_candidates, analysis_data) %>%
    pmap(~ list2(...) %>%
           imap(~ rename(.x, !!str_c("r_", .y) := type)) %>%
           reduce(function(model, to_filter_on) semi_join(model, to_filter_on, by = names(to_filter_on)), .init = r@types_data) %>%
           pull(latent_type_index))
})

setGeneric("get_prob_indices", function(r, outcome, ..., as_data = FALSE) standardGeneric("get_prob_indices"))

#' Get latent type indices for given intervention
#'
#' @param r S4 \code{StructuralCausalModel} object
#' @param outcome Outcome to predict probability equals to 1.
#' @param ... Intervention statement
#' @param cond Conditional statement (for abduction) [default = TRUE]
#'
#' @return Vector of indices
#' @export
setMethod("get_prob_indices", "StructuralCausalModel", function(r, outcome, ..., as_data = FALSE) {
  r %<>% set_obs_outcomes(...) #, cond = cond)

  outcome_exper <- str_c(outcome, " == 1") %>%
      parse_expr() %>%
      as_quosure(env = global_env())

  masked_data <- r@types_data %>%
    unnest(outcomes) %>%
    mutate(mask = !!outcome_exper)

  if (as_data) {
    return(masked_data)
  } else {
    masked_data %>%
      pull(mask) %>%
      which()
  }
})

setGeneric("num_exogenous_comb", function(r) standardGeneric("num_exogenous_comb"))

setMethod("num_exogenous_comb", "StructuralCausalModel", function(r) nrow(r@exogenous_prob))

setGeneric("create_sampler", function(r,
                                      estimands = NULL, rep_estimands = NULL,
                                      model_levels = NA_character_, cv_level = NA_character_, estimand_levels = NULL, between_entity_diff_levels = NULL,
                                      analysis_data,
                                      ...,
                                      tau_level_sigma = 1, discrete_beta_hyper_sd = 1, discretized_beta_hyper_sd = 1,
                                      calculate_marginal_prob = FALSE,
                                      use_random_binpoint = TRUE,
                                      num_sim_unique_entities = 0,
                                      alternative_model_file = NULL) {
  standardGeneric("create_sampler")
})

create_sampler_creator <- function() {
  function(r,
           estimands = NULL, rep_estimands = NULL,
           model_levels = NA_character_, cv_level = NA_character_, estimand_levels = NULL, between_entity_diff_levels = NULL,
           analysis_data,
           ...,
           tau_level_sigma = 1, discrete_beta_hyper_sd = 1, discretized_beta_hyper_sd = 1,
           calculate_marginal_prob = FALSE,
           use_random_binpoint = TRUE,
           num_sim_unique_entities = 0,
           alternative_model_file = NULL) {
    new_sampler <- if (is_null(alternative_model_file)) {
      stanmodels$bounded %>% as("Sampler")
    } else {
      stan_model(alternative_model_file) %>% as("Sampler")
    }

    new_sampler@structural_model <- r
    new_sampler@endogenous_latent_type_variables <- r@endogenous_latent_type_variables

    new_sampler@cv_level <- factor(arg_match(cv_level, values = c(model_levels, NA_character_)), levels = model_levels)
    new_sampler@model_levels <- model_levels

    level_ids <- model_levels %>% purrr::set_names(seq_along(.), .)

    stopifnot(all(estimand_levels %in% purrr::discard(model_levels, is.na)))
    stopifnot(all(between_entity_diff_levels %in% purrr::discard(model_levels, is.na)))

    discretized_response_names <- names(r@discretized_responses)
    discrete_response_names <- names(purrr::discard(r@responses, ~ is(., "DiscretizedResponse")))

    discrete_rename_list <- enquos(...) %>%
      magrittr::extract(intersect(names(.), discrete_response_names))

    discretized_data_list <- enquos(...) %>%
      magrittr::extract(intersect(names(.), discretized_response_names))

    new_sampler@analysis_data <- if (missing(analysis_data)) {
      analysis_data <- create_simulation_analysis_data(r)
    } else if (!is_null(analysis_data)) {
      discretized_data_list <- transmute(analysis_data, !!!discretized_data_list)
      discretized_data_list <- discretize_continuous_variables(r, !!!discretized_data_list)

      analysis_data <- analysis_data %>%
        mutate(!!!discrete_rename_list, !!!discretized_data_list) %>%
        select(-one_of(c("candidate_group_id", "experiment_assignment_id")))
    } else {
      tibble()
    }

    if (nrow(new_sampler@analysis_data) > 0) {
      new_sampler@analysis_data %<>%
        mutate_at(vars(all_of(purrr::discard(model_levels, is.na))), ~ if (!is.factor(.x)) factor(.x) else .x) %>%
        mutate(unique_entity_id = group_indices(., !!!syms(purrr::discard(model_levels, is.na)))) %>%
        left_join(select(r@candidate_groups, -candidate_group), by = setdiff(names(r@candidate_groups), c("candidate_group", "candidate_group_id"))) %>%
        left_join(select(r@exogenous_prob, -ex_prob), by = r@exogenous_variables)

      new_sampler@unique_entity_ids <- new_sampler@analysis_data %>%
        select(unique_entity_id, purrr::discard(model_levels, is.na)) %>%
        arrange(unique_entity_id) %>%
        select(-unique_entity_id) %>%
        distinct() %>%
        mutate_if(~ !is.factor(.), forcats::as_factor) %>%
        mutate_all(fct_drop)
    } else {
      new_sampler@unique_entity_ids <- tibble()
    }

    new_sampler@estimand_levels <- if (is_null(estimand_levels)) character(0) else estimand_levels
    new_sampler@between_entity_diff_levels <- if (is_null(between_entity_diff_levels)) character(0) else between_entity_diff_levels

    stan_est_info <- if (!is_null(estimands)) {
      new_sampler@estimands <- estimands

      # if (!is_null(between_entity_diff_levels)) {
      #   new_sampler@estimators %<>%
      #     add_between_level_entity_diff_estimands(between_entity_diff_levels, new_sampler@analysis_data)
      # }

      get_stan_info(estimands) %>%
        list_modify(
          num_estimand_levels = length(new_sampler@estimand_levels),
          num_abducted_estimands = sum(.$abducted_prob_size > 0),
          abducted_estimand_ids = which(.$abducted_prob_size > 0),
          abducted_prob_size = .$abducted_prob_size %>% keep(~ . > 0),
          estimand_levels = as.array(level_ids[new_sampler@estimand_levels]),

          num_between_entity_diff_levels = length(new_sampler@between_entity_diff_levels),
          between_entity_diff_levels =  as.array(level_ids[new_sampler@between_entity_diff_levels]),
        )
    }

    new_sampler@stan_data <- {
      if (nrow(new_sampler@analysis_data) > 0) {
        new_sampler@analysis_data %>%
          select(
            obs_unique_entity_id = unique_entity_id,
            obs_candidate_group = candidate_group_id,
          )
        } else {
          tibble(
            obs_unique_entity_id = array(0, dim = 0),
            obs_candidate_group = array(0, dim = 0)
          )
        }
      } %>%
      as.list() %>%
      list_modify(!!!lst(
        num_obs = length(.$obs_unique_entity_id),
        num_r_types = num_endogenous_types(r),

        num_discrete_r_types = num_discrete_types(r),
        num_discretized_r_types = num_discretized_types(r),
        num_compatible_discretized_r_types = if (num_discretized_r_types > 0) r %>% get_discretized_pruning_data() %>% arrange(low, hi) %>% count(low) %>% pull(n) else array(0, dim = 0),
        compatible_discretized_r_types = if (num_discretized_r_types > 0) r %>% get_discretized_pruning_data() %>% arrange(low, hi) %>% mutate_all(as.integer) %>% pull(hi) else array(0, dim = 0),

        discrete_r_type_id = r@types_data$discrete_r_type_id,

        calculate_marginal_prob,

        num_candidate_groups = nrow(r@candidate_groups),
        candidate_group_size = r@candidate_groups$candidate_group %>% map_int(nrow),
        candidate_group_ids = r@candidate_groups %>%
          unnest(candidate_group) %>%
          pull(latent_type_index),

        num_unique_entities = max(1, nrow(new_sampler@unique_entity_ids), num_sim_unique_entities),

        unique_entity_candidate_groups = if (num_obs > 0) {
          new_sampler@analysis_data %>%
            count(unique_entity_id, candidate_group_id) %>%
            arrange(unique_entity_id, candidate_group_id) %>%
            pull(candidate_group_id) %>%
            as.array()
        } else {
          array(0, dim = 0)
        },


        num_unique_entity_candidate_groups = if (num_obs > 0) {
          new_sampler@analysis_data %>%
            count(unique_entity_id, candidate_group_id) %>%
            count(unique_entity_id, name = "num_candidates") %>%
            arrange(unique_entity_id) %>%
            pull(num_candidates) %>%
            as.array()
        } else {
          array(0, dim = 0)
        },

        num_unique_entity_in_candidate_groups = if (num_obs > 0) {
          new_sampler@analysis_data %>%
            count(unique_entity_id, candidate_group_id) %>%
            arrange(unique_entity_id, candidate_group_id) %>%
            pull(n)
        } else {
          array(0, dim = 0)
        },

        obs_in_unique_entity_in_candidate_groups = if (num_obs > 0) {
          new_sampler@analysis_data %>%
            transmute(unique_entity_id, candidate_group_id, obs_index = seq(n())) %>%
            arrange_at(vars(-obs_index)) %>%
            pull(obs_index)
        } else {
          array(0, dim = 0)
        },

        num_experiment_types = num_exogenous_comb(r),
        experiment_types_prob = r %>% get_ex_prob(),

        num_levels = if (num_sim_unique_entities > 0) 1 else ncol(new_sampler@unique_entity_ids),

        unique_entity_ids = if (!is_empty(new_sampler@unique_entity_ids)) {
          new_sampler@unique_entity_ids %>% mutate_all(as.integer) %>% as.matrix() %>% as.array()
        } else array(1, dim = c(num_unique_entities, num_levels)),

        num_bg_variables = n_distinct(r@endogenous_latent_type_variables$type_variable),
        num_bg_variable_types = r@endogenous_latent_type_variables %>% count(type_variable) %>% pull(n) %>% as.array(),
        num_bg_variable_type_combo_members = map_int(r@endogenous_latent_type_variables$latent_type_ids, length),
        bg_variable_type_combo_members = unlist(r@endogenous_latent_type_variables$latent_type_ids),

        cutpoints = if (is_empty(r@discretized_responses)) array(0, dim = 0) else unlist(r@discretized_responses[[1]]@cutpoints),
        num_cutpoints = length(cutpoints),

        compatible_discretized_pair_ids =
          if (num_cutpoints >= 3) {
          discretized_var_name <- names(r@discretized_responses)

          r@types_data %>% select(str_c("r_", discretized_var_name, "_1"),
                                  if (num_cutpoints >= 3) num_range(str_c(discretized_var_name, "_pair_id_"), seq(2, num_cutpoints - 2))) %>%
            mutate_all(as.integer) %>%
            as.matrix() %>%
            t()
        } else {
          array(0, dim = c(0, num_r_types))
        },

        # Default values
        num_discrete_estimands = 0,
        num_atom_estimands = 0,
        num_diff_estimands = 0,
        num_mean_diff_estimands = 0,
        num_utility_diff_estimands = 0,
        num_discretized_groups = 0,
        num_abducted_estimands = 0,
        num_between_entity_diff_levels = 0,
        num_op = 0,
        abducted_prob_size = array(1, dim = 0),
        abducted_prob_index = array(1, dim = 0),
        est_prob_size = array(1, dim = 0),
        est_prob_index = array(1, dim = 0),
        diff_estimand_atoms = array(0, dim = 0),
        mean_diff_estimand_atoms = array(0, dim = 0),
        discretized_group_ids = array(0, dim = 0),
        abducted_estimand_ids = array(0, dim = 0),
        between_entity_diff_levels = array(0, dim = 0),
        utility = array(0, dim = 0),
        utility_diff_estimand_atoms = array(0, dim = 0),
        num_discrete_utility_values = 0,

        num_responses = length(get_responses(r)),
        type_response_value = r@types_data %>%
          unnest(outcomes) %>%
          select(names(get_responses(r))) %>%
          map(matrix, nrow = length(get_ex_prob(r)), byrow = FALSE),

        # experiment_assign_entity = count(new_sampler@analysis_data, unique_entity_id, experiment_assignment_id),
        # num_experiment_entities = nrow(experiment_assign_entity),

        # num_rep_corr_estimands = length(rep_estimands),
        # rep_corr_outcomes = if (num_rep_corr_estimands > 0) rep_estimands %>% map(~ which(names(r@responses) %in% c(.x@outcome1, .x@outcome2))) %>% unlist() %>% as.array() else array(0, dim = 0),
        # rep_corr_cond = if (num_rep_corr_estimands > 0) rep_estimands %>% map_if(~ !is.na(.x@cond), ~ which(names(r@responses) == .x@cond), .else = ~ 0) %>% unlist() %>% as.array() else array(0, dim = 0),

        generate_rep = 0, # num_rep_corr_estimands > 0,

        num_estimand_levels = 0,
        estimand_levels = array(1, dim = 0),

        log_lik_level = if (!is.na(cv_level)) as.integer(new_sampler@cv_level) else 0,

        num_shards = 0,

        use_random_binpoint,

        # Priors

        discrete_beta_hyper_sd,
        discretized_beta_hyper_sd,
        tau_level_sigma = if (all(is.na(model_levels)) && num_sim_unique_entities == 0) array(0, dim = 0)
          else if (length(tau_level_sigma) == 1 && !is_named(tau_level_sigma)) as.array(rep(tau_level_sigma, num_levels))
          else as.array(unlist(tau_level_sigma[model_levels])),
      )) %>%
      list_modify(!!!stan_est_info) %>%
      map_if(is.factor, as.integer)

    return(new_sampler)
  }
}


#' Create S4 sampler object for given model and estimands
#'
#' @param estimands estimands to calculate
#' @param rep_estimands
#' @param model_levels Model level names
#' @param cv_level Cross validation level
#' @param estimand_levels Which model levels to calculate estimands for
#' @param between_entity_diff_levels Which model levels to calculate difference between estimands for
#' @param analysis_data Dataset to use for analysis
#' @param ... specify names of model variables in data
#' @param tau_level_sigma hyperparameter
#' @param discrete_beta_hyper_sd hyperparameter
#' @param discretized_beta_hyper_sd hyperparameter
#' @param calculate_marginal_prob Calculate latent type marginal probabilities
#' @param use_random_binpoint
#' @param num_sim_unique_entities
#' @param alternative_model_file
#'
#' @return \code{Sampler} S4 object
#' @export
setMethod("create_sampler", "StructuralCausalModel", create_sampler_creator())

setGeneric("get_known_estimands", function(r, estimands) {
  standardGeneric("get_known_estimands")
})

setMethod("get_known_estimands", "StructuralCausalModel", function(r, estimands) {
  estimands %>% calculate_from_known_dgp(r)
})

setGeneric("build_estimand_collection", function(model, ...) standardGeneric("build_estimand_collection"))

#' Create collection of estimands
#'
#' @param model current model
#' @param ... estimands
#' @param utility utility values for discretized cutpoints
#' @param cores Number of cores to use in preparing estimands
#'
#' @return \code{EstimandCollection} S4 object
#' @export
setMethod("build_estimand_collection", "StructuralCausalModel", function (model, ..., utility = NA_real_, cores = 1) {
  estimand_data_builder <- function(accum_data, next_est) {
    next_estimand_id <- if (is_null(accum_data)) 1 else pull(accum_data, estimand_id) %>% max() %>% add(1)
    next_estimand_group_id <- if (is_null(accum_data)) 1 else pull(accum_data, estimand_group_id) %>% c(0) %>% max(na.rm = TRUE) %>% add(1)

    next_est %>%
      set_model(model) %>%
      get_component_estimands(next_estimand_id, next_estimand_group_id) %>%
      bind_rows(accum_data)
  }

  new_est_collection <- new(
    "EstimandCollection",
    model = model,
    estimands = rlang::list2(...) %>%
      compact() %>%
      reduce(estimand_data_builder, .init = NULL) %>%
      mutate(
        est_type = case_when(
          map_lgl(est_obj, is, "AtomEstimand") ~ "atom",
          map_lgl(est_obj, is, "DiscreteDiffEstimand") ~ "diff",
          map_lgl(est_obj, is, "DiscretizedMeanEstimand") ~ "mean",
          map_lgl(est_obj, is, "DiscretizedUtilityEstimand") ~ "utility",
          map_lgl(est_obj, is, "DiscretizedMeanDiffEstimand") ~ "mean-diff",
          map_lgl(est_obj, is, "DiscretizedUtilityDiffEstimand") ~ "utility-diff",
        ) %>% factor(levels = c("atom", "diff", "mean", "utility", "mean-diff", "utility-diff"))
      ) %>% {
        if ("cutpoint" %in% names(.)) . else mutate(., cutpoint = NA_real_)
      }
  )

  new_est_collection@estimands %<>%
    arrange(est_type, estimand_id) %>%
    mutate(new_estimand_id = seq(n())) %>% {
      if (any(str_detect(names(.), "estimand_id_(left|right)"))) {
        left_join(., select(., estimand_id, new_estimand_id), by = c("estimand_id_left" = "estimand_id"), suffix = c("", "_left")) %>%
        left_join(select(., estimand_id, new_estimand_id), by = c("estimand_id_right" = "estimand_id"), suffix = c("", "_right"))
      } else .
    } %>%
    select(-starts_with("estimand_id")) %>%
    rename_at(vars(starts_with("new_estimand")), str_remove, "new_") %>%
    group_by(est_type) %>%
    mutate(within_est_type_index = seq(n())) %>%
    ungroup() %>%
    mutate(
      atom_index = if_else(fct_match(est_type, "atom"), within_est_type_index, NA_integer_),
      diff_index = if_else(fct_match(est_type, "diff"), within_est_type_index, NA_integer_),
      mean_index = if_else(fct_match(est_type, "mean"), within_est_type_index, NA_integer_),
      utility_index = if_else(fct_match(est_type, "utility"), within_est_type_index, NA_integer_),
    ) %>%
    select(-within_est_type_index) %>% {
      if (any(str_detect(names(.), "estimand_id_(left|right)"))) {
        left_join(., filter(., fct_match(est_type, "mean")) %>% select(estimand_id, mean_index), by = c("estimand_id_left" = "estimand_id"), suffix = c("", "_left")) %>%
          left_join(filter(., fct_match(est_type, "mean")) %>% select(estimand_id, mean_index), by = c("estimand_id_right" = "estimand_id"), suffix = c("", "_right")) %>%
          left_join(filter(., fct_match(est_type, "utility")) %>% select(estimand_id, utility_index), by = c("estimand_id_left" = "estimand_id"), suffix = c("", "_left")) %>%
          left_join(filter(., fct_match(est_type, "utility")) %>% select(estimand_id, utility_index), by = c("estimand_id_right" = "estimand_id"), suffix = c("", "_right"))
      } else .
    }

  if (any(!is.na(utility))) {
    stopifnot(length(utility) == (length(get_discretized_cutpoints(model)) + 1))
  } else {
    stopifnot(filter(new_est_collection@estimands, fct_match(est_type, "utility")) %>% nrow() %>% equals(0))
  }

  num_diff_estimands <- num_estimands(new_est_collection, "diff")

  discretized_group_ids <- if (any(fct_match(new_est_collection@estimands$est_type, "mean"))) new_est_collection@estimands %>%
    filter(fct_match(est_type, "mean")) %>%
    semi_join(new_est_collection@estimands, ., by = c("estimand_group_id" = "mean_estimand_group_id")) %>%
    select(estimand_group_id, estimand_id)

  new_est_collection@est_stan_info <- new_est_collection %>%
    get_stan_data_structures(cores = cores) %>%
    list_modify(!!!lst(
      num_discrete_estimands = num_estimands(new_est_collection, c("atom", "diff")),
      num_atom_estimands = num_estimands(new_est_collection, "atom"),
      num_diff_estimands,
      num_mean_diff_estimands = num_estimands(new_est_collection, "mean-diff"),
      num_utility_diff_estimands = num_estimands(new_est_collection, "utility-diff"),

      utility = if (all(!is.na(utility))) as.array(utility) else array(0, dim = 0),
      num_discrete_utility_values = length(utility),

      diff_estimand_atoms = new_est_collection@estimands %>%
        filter(fct_match(est_type, "diff")) %>%
        arrange(estimand_id) %>%
        select(matches("^estimand_id_(left|right)$")) %>% {
          if (nrow(.) > 0) as.matrix(.) %>% t() %>% c() else array(0, dim = 0)
        },

      mean_diff_estimand_atoms = new_est_collection@estimands %>%
        filter(fct_match(est_type, "mean-diff")) %>%
        arrange(estimand_id) %>%
        select(matches("mean_index_(left|right)")) %>% {
          if (nrow(.) > 0) as.matrix(.) %>% t() %>% c() else array(0, dim = 0)
        },

      utility_diff_estimand_atoms = new_est_collection@estimands %>%
        filter(fct_match(est_type, "utility-diff")) %>%
        arrange(estimand_id) %>%
        select(matches("utility_index_(left|right)")) %>% {
          if (nrow(.) > 0) as.matrix(.) %>% t() %>% c() else array(0, dim = 0)
        },

      num_discretized_groups = if (!is_empty(discretized_group_ids)) n_distinct(discretized_group_ids$estimand_group_id) else 0,
      discretized_group_ids = if (!is_empty(discretized_group_ids)) discretized_group_ids %>% pull(estimand_id) else array(0, dim = 0),
    ))

  return(new_est_collection)
})

setGeneric("get_latent_type_by_observed_outcomes", function(model, obs_data, ...) {
  standardGeneric("get_latent_type_by_observed_outcomes")
})

#' Get linear programming restriction matrix
#'
#' @param S4 \code{StructuralCausalModel} object
#'
#' @return Data set
setMethod("get_latent_type_by_observed_outcomes", "StructuralCausalModel", function(model, obs_data, ...) {
  all_variables <- names(model@responses)

  types_by_outcomes <- model@types_data %>%
    unnest(outcomes) %>%
    arrange_at(vars(all_of(all_variables))) %>%
    group_by_at(vars(all_of(all_variables))) %>%
    group_nest(.key = "latent_type_mask") %>%
    mutate(latent_type_mask = map(latent_type_mask, ~ inset(.y, .x$latent_type_index, 1), rep(0, max(model@types_data$latent_type_index))))
    # data %>%
    # do.call(rbind, .)

  obs_data %>%
    arrange_at(vars(all_of(all_variables))) %>%
    group_by_at(vars(all_of(all_variables))) %>%
    count() %>%
    group_by_at(vars(all_of(model@exogenous_variables))) %>%
    mutate(prob = n / sum(n)) %>%
    ungroup() %>%
    right_join(types_by_outcomes, by = all_variables) %>%
    mutate(prob = coalesce(prob, 0.0))
})

#' Defined the structural model's directed acyclic graph and each variable's response function
#'
#' @param ... model variables
#' @param exogenous_prob \code{data.frame} with the experiment's assignment mechanism's probabilities.
#'
#' @return A \code{StructuralCausalModel} S3 object
#' @export
define_structural_causal_model <- function(..., exogenous_prob) {
  exogenous_variables <- setdiff(names(exogenous_prob), "ex_prob")
  responses <- list2(...) %>%
    purrr::compact()

  if (!missing(exogenous_prob) && !"ex_prob" %in% names(exogenous_prob)) {
    stop("'ex_prob' undefined in exogenous probability.")
  }

  leaf_responses <- map(responses, get_responses) %>% flatten() %>% purrr::compact()
  discretized_responses <- keep(responses, ~ is(., "DiscretizedResponseGroup")) %>%
    purrr::set_names(map_chr(., get_output_variable_name))

  if (length(discretized_responses) > 1) {
    stop("More than one discretized reponse not yet supported.")
  }

  comb <- new(
    "StructuralCausalModel",
    exogenous_variables = exogenous_variables,
    exogenous_prob = exogenous_prob %>%
      mutate(experiment_assignment_id = seq(n())),
    responses = leaf_responses %>%
      purrr::set_names(map_chr(., get_output_variable_name)),
    discretized_responses = discretized_responses,
    types_data = leaf_responses %>% {
        # if (prune_discretized) purrr::discard(., ~ is(., "DiscretizedResponse")) else .
        purrr::discard(., ~ is(., "DiscretizedResponse"))
      } %>%
      purrr::set_names(map_chr(., get_output_variable_name) %>% str_c("r_", .)) %>%
      map(get_response_type_names) %>%
      map(forcats::as_factor) %>%
      do.call(expand.grid, .)
  )

  if (length(discretized_responses) > 0) {
    discretized_var_name <- first(discretized_responses) %>% get_output_variable_name()

    discretized_var_grid <- first(discretized_responses) %>%
      get_types_grid() %>%
      select(num_range(str_c("r_", discretized_var_name, "_"), seq(length(get_discretized_cutpoints(comb)))), # Arrange columns in proper order
             matches("pair_id_\\d+$")) %>%
      arrange_at(vars(num_range(str_c("r_", discretized_var_name, "_"), seq(length(get_discretized_cutpoints(comb))))))

    comb@types_data %<>%
      merge(discretized_var_grid)
  }

  comb %<>% set_obs_outcomes()

  comb@types_data %<>%
    # right_join(exogenous_prob, by = setdiff(names(exogenous_prob), c("experiment_assignment_id", "ex_prob"))) %>% # exogenous_prob is used to also exclude any exogenous variable combinations (hence right join)
    right_join(exogenous_prob, by = exogenous_variables) %>% # exogenous_prob is used to also exclude any exogenous variable combinations (hence right join)
    nest(outcomes = c(str_c("r_", exogenous_variables), map_chr(leaf_responses, get_output_variable_name), ex_prob))

  discrete_r_type_col <- comb@responses %>%
    purrr::discard(is, "DiscretizedResponse") %>%
    names() %>%
    str_c("r_", .) %>%
    intersect(names(comb@types_data))

  comb@types_data %<>%
    group_by_at(discrete_r_type_col) %>%
    group_nest(.key = "temp_id_nest") %>%
    mutate(discrete_r_type_id = seq(n())) %>%
    unnest(temp_id_nest) %>%
    mutate(latent_type_index = seq(n()))

  comb@nester <- function(data) {
    nest(data, outcomes = c(str_c("r_", exogenous_variables), map_chr(leaf_responses, get_output_variable_name), ex_prob, r_candidates, num_r_candidates))
  }

  comb@types_data %<>%
    unnest(outcomes) %>%
    mutate(
      r_candidates = get_candidates(comb),
      num_r_candidates = map_int(r_candidates, length)
    ) %>%
    comb@nester()

  comb@candidate_groups <- comb@types_data %>%
    unnest(outcomes) %>%
    select(all_of(map_chr(leaf_responses, get_output_variable_name)), latent_type_index) %>%
    nest(candidate_group = latent_type_index) %>%
    mutate(candidate_group_id = seq(n()))

  get_latent_type_ids <- function(type_variable, type) filter(comb@types_data, fct_match(!!sym(type_variable), type)) %>% pull(latent_type_index)

  r_type_cols <- comb@responses %>%
    map_chr(~ .@output) %>%
    setdiff(comb@exogenous_variables) %>%
    str_c("r_", .)

  comb@endogenous_latent_type_variables <- comb@types_data %>%
    select(all_of(r_type_cols)) %>%
    map(fct_unique) %>%
    enframe(name = "type_variable", value = "type") %>%
    unnest(type) %>%
    mutate(latent_type_ids = map2(type_variable, type, get_latent_type_ids),
           marginal_latent_type_index = seq(n()))

  return(comb)
}

create_scm_simulation <- function(scm, sample_size, prob, concentration_alpha = 1) {
  new_sim <- scm

  exogenous_variables <- scm@exogenous_variables

  new_sim@nester <- function(data) {
    nest(data, outcomes = c(str_c("r_", exogenous_variables), map_chr(scm@responses, get_output_variable_name), ex_prob, num_obs, r_candidates, num_r_candidates))
  }

  if (!missing(prob)) {
    new_sim@types_data %<>%
      mutate(prob = prob)
  } else if (!has_name(new_sim@types_data, "prob")) {
    new_sim@types_data %<>%
      mutate(prob = c(MCMCpack::rdirichlet(1, rep(concentration_alpha, n()))))
  } else {
    new_sim@types_data %<>%
      mutate(outcomes = map(outcomes, select, -num_obs))
  }

  new_sim@types_data %<>%
    mutate(num_obs = c(rmultinom(1, sample_size, prob))) %>%
    unnest(outcomes) %>%
    mutate(num_obs = round(num_obs * ex_prob)) %>%
    # mutate(num_obs = c(rmultinom(1, sample_size, prob * ex_prob))) %>%
    new_sim@nester()

  return(as(new_sim, "StructuralCausalModelSimulation"))
}

setGeneric("create_prior_predicted_simulation", function(scm, ...) {
  standardGeneric("create_prior_predicted_simulation")
})

#' Title
#'
#' @param scm
#' @param sample_size
#' @param chains
#' @param iter
#' @param ...
#' @param num_sim
#' @param num_entities
#'
#' @return
#' @export
#'
#' @examples
setMethod("create_prior_predicted_simulation", "StructuralCausalModel", function(scm, sample_size, chains, iter, ..., num_sim = 1, num_entities = 1) {
  prior_predict_sampler <- scm %>% create_sampler(analysis_data = NULL, ..., num_sim_unique_entities = num_entities)

  prior_predict_fit <- prior_predict_sampler %>% sampling(run_type = "prior-predict", pars = "r_log_prob", chains = chains, iter = iter)

  prior_predict_fit %>%
    as.data.frame(pars = "r_log_prob") %>%
    sample_n(num_sim, replace = FALSE) %>%
    mutate(iter_id = seq(n())) %>%
    pivot_longer(cols = -iter_id, names_to = "latent_type_index", values_to = "iter_prob") %>%
    tidyr::extract(latent_type_index, "latent_type_index", "(\\d+)", convert = TRUE) %>%
    mutate(
      iter_prob = exp(iter_prob),
      entity_index = ((latent_type_index - 1) %/% prior_predict_sampler@stan_data$num_r_types) + 1,
      latent_type_index = ((latent_type_index - 1) %% prior_predict_sampler@stan_data$num_r_types) + 1
    ) %>%
    nest(prob_data = c(latent_type_index, iter_prob)) %>%
    mutate(
      sim = map(prob_data, ~ create_scm_simulation(scm, sample_size %/% num_entities, .x$iter_prob)),
    ) %>%
    select(-prob_data) %>%
    group_nest(iter_id, .key = "entity_data")
    # mutate(
    #   analysis_data = map(entity_data, ~ map2_dfr(.x$sim, .x$entity_index, ~ mutate(create_simulation_analysis_data(.x), entity_index = .y)))
    # ) %>%
    # select(-iter_id)

    # group_by(iter_id) %>%
    # group_map(~ map2_dfr(.x$sim, .x$entity_index, ~ mutate(create_simulation_analysis_data(.x), entity_index = .y))) %>%
    # map(mutate, entity_index = factor(entity_index)) %>% {
    #   if (length(.) > 1) . else first(.)
    # }
})

setGeneric("get_linear_programming_bounds", function(scm, obs_data, outcome, ...) {
  standardGeneric("get_linear_programming_bounds")
})

#' Get linear programming bounds
#'
#' @param scm S4 \code{StructuralCausalModel} object.
#' @param obs_data Data used to calculate conditional probabilities.
#' @param outcome Estimation outcome,
#' @param ... Intervention.
#'
#' @return Vector with minimum and maximum bounds.
#'
#' @export
setMethod("get_linear_programming_bounds", "StructuralCausalModel", function(scm, obs_data, outcome, ...) {
  # all_variables <- names(scm@responses)

  types_by_outcomes <- get_latent_type_by_observed_outcomes(scm, obs_data)

  latent_type_restrict <- types_by_outcomes %>%
    pull(latent_type_mask) %>%
    do.call(rbind, .)

  # types_by_outcomes <- obs_data %>%
  #   arrange_at(vars(all_of(all_variables))) %>%
  #   group_by_at(vars(all_of(all_variables))) %>%
  #   count() %>%
  #   group_by_at(vars(all_of(scm@exogenous_variables))) %>%
  #   mutate(prob = n / sum(n)) %>%
  #   ungroup() %>%
  #   right_join(types_by_outcomes, by = all_variables) %>%
  #   mutate(prob = coalesce(prob, 0.0))

  intervention_latent_types <- get_prob_indices(scm, outcome, ..., as_data = TRUE)
  objective_fun <- intervention_latent_types %>%
    transmute(latent_type_index, mask = mask * ex_prob) %>%
    group_by(latent_type_index) %>%
    summarize(mask = sum(mask)) %>%
    pull(mask)

  lst(min = lpSolve::lp("min", objective_fun, rbind(latent_type_restrict, 1), "=", c(types_by_outcomes$prob, 1)),
      max = lpSolve::lp("max", objective_fun, rbind(latent_type_restrict, 1), "=", c(types_by_outcomes$prob, 1)))
})
