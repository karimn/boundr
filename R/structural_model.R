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
            endogenous_joint_discrete_latent_type_variables = "data.frame",
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
                                      tau_level_sigma = 1, discretized_beta_hyper_mean = 0, discrete_beta_hyper_sd = 1, discretized_beta_hyper_sd = 1,
                                      calculate_marginal_prob = FALSE,
                                      use_random_binpoint = FALSE,
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
           tau_level_sigma = 1, discretized_beta_hyper_mean = 0, discrete_beta_hyper_sd = 1, discretized_beta_hyper_sd = 1,
           calculate_marginal_prob = FALSE,
           use_random_binpoint = FALSE,
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
        select(!any_of(c("candidate_group_id", "experiment_assignment_id")))
    } else {
      tibble()
    }

    if (nrow(new_sampler@analysis_data) > 0) {
      new_sampler@analysis_data %<>%
        mutate(
          across(all_of(purrr::discard(model_levels, is.na)), ~ if (!is.factor(.x)) factor(.x) else .x),
          unique_entity_id = group_by(., across(all_of(purrr::discard(model_levels, is.na)))) %>% group_indices()
        ) %>%
        left_join(select(r@candidate_groups, -candidate_group), by = setdiff(names(r@candidate_groups), c("candidate_group", "candidate_group_id"))) %>%
        left_join(select(r@exogenous_prob, -ex_prob), by = r@exogenous_variables)

      if (any(is.na(new_sampler@analysis_data$candidate_group_id))) {
        stop("Some observations do not match any candidate group. Check your model's assumptions.")
      }

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

    num_discrete_r_types <- num_discrete_types(r)
    num_discretized_r_types <- num_discretized_types(r)

    discretized_beta_hyper_mean %<>% {
      if (is.numeric(.)) {
        # if (any(dim(.) != c(num_discretized_r_types, num_discrete_r_types))) {
        if (NROW(.) != num_discretized_r_types || NROW(.) != num_discrete_r_types) {
          matrix(., num_discretized_r_types, num_discrete_r_types)
        } else .
      } else if (is.list(.)) {
        default_mean <- .$default

        if (is_null(default_mean)) stop("Missing default mean value for discretized prior specification.")

        if (length(.) > 1) {
          discrete_type_var <- str_c("r_", discrete_response_names)
          discretized_type_var <- str_c("r_", discretized_response_names, "_1")

          type_combos <-r@types_data %>%
            select(any_of(c(discrete_type_var))) %>%
            distinct()

          list_modify(., default = NULL) %>%
            map_if(.,
                   ~ rlang::is_function(.x) | rlang::is_formula(.x),
                   ~ rlang::as_function(.x)(type_combos),
                   .else = ~ type_combos %>% mutate(mean = .x)) %>%
            tibble::enframe(value = "discrete_types") %>%
            mutate(name = factor(name, levels = r@types_data %>% pull(discretized_type_var) %>% levels())) %>%
            complete(name) %>%
            rename(!!str_remove(discretized_type_var, "_1$") := "name") %>%
            mutate(discrete_types = map_if(discrete_types, is_null, ~ rep(default_mean, nrow(type_combos)), .else = ~ arrange(.x) %>% pull(mean))) %>%
            pull(discrete_types) %>%
            do.call(rbind, .)
        } else {
          matrix(default_mean, num_discretized_r_types, num_discrete_r_types)
        }
      }
    }

    discretized_beta_hyper_sd %<>% {
      if (is.numeric(.)) {
        if (NROW(.) != num_discretized_r_types || NROW(.) != num_discrete_r_types) {
          matrix(., num_discretized_r_types, num_discrete_r_types)
        } else .
      } else if (is.list(.)) {
        default_sd <- .$default

        if (is_null(default_sd)) stop("Missing default SD value for discretized prior specification.")

        if (length(.) > 1) {
          discrete_type_var <- str_c("r_", discrete_response_names)
          discretized_type_var <- str_c("r_", discretized_response_names, "_1")

          type_combos <-r@types_data %>%
            select(any_of(c(discrete_type_var))) %>%
            distinct()

          list_modify(., default = NULL) %>%
            map_if(.,
                   ~ rlang::is_function(.x) | rlang::is_formula(.x),
                   ~ rlang::as_function(.x)(type_combos),
                   .else = ~ type_combos %>% mutate(sd = .x)) %>%
            tibble::enframe(value = "discrete_types") %>%
            mutate(name = factor(name, levels = r@types_data %>% pull(discretized_type_var) %>% levels())) %>%
            complete(name) %>%
            rename(!!str_remove(discretized_type_var, "_1$") := "name") %>%
            mutate(discrete_types = map_if(discrete_types, is_null, ~ rep(default_sd, nrow(type_combos)), .else = ~ arrange(.x) %>% pull(sd))) %>%
            pull(discrete_types) %>%
            do.call(rbind, .)
        } else {
          matrix(default_sd, num_discretized_r_types, num_discrete_r_types)
        }
      }
    }

    discretized_beta_hyper_sd <- if (is.numeric(discretized_beta_hyper_sd)) {
      matrix(discretized_beta_hyper_sd, num_discretized_r_types, num_discrete_r_types)
    } else if (is.list(discretized_beta_hyper_sd)) {
      stop("Not yet supported")
    }

    new_sampler@stan_data <- if (nrow(new_sampler@analysis_data) > 0) {
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

    new_sampler@stan_data %<>%
      as.list() %>%
      list_modify(!!!lst(
        num_obs = length(.$obs_unique_entity_id),
        num_r_types = num_endogenous_types(r),

        num_discrete_r_types,
        num_discretized_r_types,
        num_compatible_discretized_r_types = if (num_discretized_r_types > 0) {
          r %>%
            get_discretized_pruning_data() %>% {
              if (!is_empty(.)) {
                arrange(., low, hi) %>%
                  count(low) %>%
                  pull(n)
              } else {
                rep(num_discretized_r_types, num_discretized_r_types)
              }
            }
        } else array(0, dim = 0),

        compatible_discretized_r_types = if (num_discretized_r_types > 0) {
          r %>%
            get_discretized_pruning_data() %>% {
              if (!is_empty(.)) {
                arrange(., low, hi) %>%
                  mutate_all(as.integer) %>%
                  pull(hi)
              } else {
                rep(seq(num_discretized_r_types), num_discretized_r_types)
              }
            }
        } else array(0, dim = 0),

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
            distinct(unique_entity_id, candidate_group_id) %>%
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

        num_discrete_bg_variables = r@endogenous_latent_type_variables %>%
          filter(!discretized) %$%
          n_distinct(type_variable),
        num_discrete_bg_variable_types = r@endogenous_latent_type_variables %>%
          filter(!discretized) %>%
          count(type_variable) %>%
          pull(n) %>%
          as.array(),

        # Which joint probabilities include each of the marginal types
        num_discrete_bg_variable_type_combo_members = r@endogenous_latent_type_variables %>%
          filter(!discretized) %>%
          pull(latent_type_ids) %>%
          map_int(length),
        discrete_bg_variable_type_combo_members = r@endogenous_latent_type_variables %>%
          filter(!discretized) %>%
          pull(latent_type_ids) %>%
          unlist(),

        num_discretized_bg_variables = r@endogenous_latent_type_variables %>%
          filter(discretized) %$%
          n_distinct(type_variable),

        num_discretized_bg_variable_type_combo_members = r@endogenous_latent_type_variables %>%
          filter(discretized) %>%
          unnest(latent_type_ids) %>%
          pull(latent_type_ids) %>% # there's a second one in there
          map_int(nrow),

        discretized_bg_variable_type_combo_members = if (num_discretized_bg_variables > 0) {
          r@endogenous_latent_type_variables %>%
            filter(discretized) %>%
            unnest(latent_type_ids) %>%
            unnest(latent_type_ids) %>%
            pull(latent_type_index)
        } else array(0, dim = 0),

        num_joint_discrete_combo_members = r@endogenous_joint_discrete_latent_type_variables$latent_type_ids %>%
          first() %>%
          length(),
        joint_discrete_combo_members = r@endogenous_joint_discrete_latent_type_variables$latent_type_ids %>%
          unlist(),

        cutpoints = if (is_empty(r@discretized_responses)) array(0, dim = 0) else unlist(r@discretized_responses[[1]]@cutpoints),
        num_cutpoints = length(cutpoints),

        compatible_discretized_pair_ids = if (num_cutpoints >= 3) {
          discretized_var_name <- names(r@discretized_responses)

          r@types_data %>%
            select(
              str_c("r_", discretized_var_name, "_1"),
              if (num_cutpoints >= 3) num_range(str_c(discretized_var_name, "_pair_id_"), seq(2, num_cutpoints - 2))
            ) %>%
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

        generate_rep = 0,

        num_estimand_levels = 0,
        estimand_levels = array(1, dim = 0),

        log_lik_level = if (!is.na(cv_level)) as.integer(new_sampler@cv_level) else 0,

        num_shards = 0,

        use_random_binpoint,

        # Priors

        discretized_beta_hyper_mean,

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

setGeneric("get_known_marginal_latent_type_prob", function(model, ...) {
  standardGeneric("get_known_marginal_latent_type_prob")
})

setMethod("get_known_marginal_latent_type_prob", "StructuralCausalModel", function(model, ...) {
  discrete_marginal_prob <- model@endogenous_joint_discrete_latent_type_variables %>%
    unnest(latent_type_ids) %>%
    left_join(select(model@types_data, latent_type_index, prob), by = c("latent_type_ids" = "latent_type_index")) %>%
    group_by(discrete_r_type_id) %>%
    summarize(marginal_prob = sum(prob)) %>%
    ungroup()

  discrete_type_variables <- model@endogenous_latent_type_variables %>%
    filter(!discretized) %>%
    pull(type_variable) %>%
    unique()

  model@endogenous_latent_type_variables %>%
    mutate(latent_type_ids = map_if(latent_type_ids, discretized, unnest, latent_type_ids, .else = ~ tibble(latent_type_index = .))) %>%
    select(type_variable, type, discretized, latent_type_ids) %>%
    unnest(latent_type_ids) %>%
    left_join(select(model@types_data, latent_type_index, prob), by = "latent_type_index") %>%
    group_by(type_variable, type, discretized, discrete_r_type_id) %>%
    summarize(marginal_prob = sum(prob)) %>%
    ungroup() %>%
    left_join(discrete_marginal_prob, by = "discrete_r_type_id", suffix = c("", "_discrete")) %>%
    mutate(
      marginal_prob_discrete = coalesce(marginal_prob_discrete, 1),
      marginal_prob = marginal_prob / marginal_prob_discrete
    ) %>%
    select(-marginal_prob_discrete) %>%
    left_join(distinct_at(model@types_data, vars(discrete_r_type_id, all_of(discrete_type_variables))), by = "discrete_r_type_id")
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
setMethod("get_latent_type_by_observed_outcomes", "StructuralCausalModel", function(model, obs_data, cond = TRUE, ...) {
  all_variables <- names(model@responses)

  types_by_outcomes <- model@types_data %>%
    unnest(outcomes) %>%
    arrange_at(vars(all_of(all_variables))) %>%
    group_by_at(vars(all_of(all_variables))) %>%
    group_nest(.key = "latent_type_mask") %>%
    mutate(latent_type_mask = map(latent_type_mask,
                                  ~ inset(.y, .x$latent_type_index, 1),
                                  rep(0, max(model@types_data$latent_type_index))))

  cond_types <- model@types_data %>%
    unnest(outcomes) %>%
    filter({{ cond }}) %>%
    select(latent_type_index)

  cond_prob <- if (!missing(obs_data)) {
    if (!is_logical(cond) || !cond) {
      stop("Conditional statements not supported using data.")
    }

    obs_data %>%
      arrange_at(vars(all_of(all_variables))) %>%
      group_by_at(vars(all_of(all_variables))) %>%
      count() %>%
      group_by_at(vars(all_of(model@exogenous_variables))) %>%
      mutate(prob = n / sum(n)) %>%
      ungroup()
  } else if (is(model, "StructuralCausalModelSimulation")) {
    model@types_data %>%
      semi_join(cond_types, by = "latent_type_index") %>%
      unnest(outcomes) %>%
      mutate(prob = prob * ex_prob) %>%
      select(all_of(all_variables), prob) %>%
      group_by_at(vars(all_of(model@exogenous_variables))) %>%
      mutate(prob = prob / sum(prob)) %>%
      ungroup()
  } else {
    stop("Must provide data if not using a simulation model.")
  }

  cond_prob %>%
    group_by_at(vars(all_of(all_variables))) %>%
    summarize(prob = sum(prob)) %>%
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
  } else if (length(discretized_responses) > 0) {
    first(discretized_responses)@pruning_data
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

  get_latent_type_ids <- function(type_variable, type, discretized = FALSE) {
    # TODO Use the discretized parameter to calculate conditional on discrete types
    filter(comb@types_data, fct_match(!!sym(type_variable), type)) %>% pull(latent_type_index)
  }

  r_type_cols <- comb@responses %>%
    map_df(~ tibble(type_variable = .x@output, discretized = is(.x, "DiscretizedResponse"))) %>%
    filter(!type_variable %in% comb@exogenous_variables) %>%
    mutate(type_variable = str_c("r_", type_variable))

  comb@endogenous_latent_type_variables <- comb@types_data %>%
    select(all_of(r_type_cols$type_variable)) %>%
    map(fct_unique) %>%
    enframe(name = "type_variable", value = "type") %>%
    unnest(type) %>%
    left_join(r_type_cols, by = "type_variable") %>%
    mutate(latent_type_ids = pmap(lst(type_variable, type, discretized), get_latent_type_ids) %>%
             map_if(discretized,
                    ~ filter(comb@types_data, latent_type_index %in% .) %>%
                      select(discrete_r_type_id, latent_type_index) %>%
                      nest(latent_type_ids = latent_type_index))) %>%
    group_by(discretized) %>%
    mutate(marginal_latent_type_index = seq(n())) %>%
    ungroup()

  comb@endogenous_joint_discrete_latent_type_variables <- comb@types_data %>%
    select(latent_type_index, discrete_r_type_id, all_of(filter(r_type_cols, !discretized) %>% pull(type_variable))) %>%
    nest(latent_type_ids = latent_type_index) %>%
    mutate(latent_type_ids = map(latent_type_ids, pull, latent_type_index))

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
setMethod("create_prior_predicted_simulation", "StructuralCausalModel", function(scm, sample_size, chains, iter, ..., num_sim = 1, num_entities = 1, same_dist = FALSE) {
  prior_predict_sampler <- scm %>% create_sampler(analysis_data = NULL, ..., num_sim_unique_entities = num_entities)

  prior_predict_fit <- prior_predict_sampler %>% sampling(run_type = "prior-predict", pars = "r_log_prob", chains = chains, iter = iter)

  prior_predict_fit %>%
    as.data.frame(pars = "r_log_prob") %>% {
      if (!same_dist) {
        sample_n(., num_sim, replace = FALSE)
      } else {
        sample_n(., 1, replace = FALSE) %>%
          slice(rep(1, num_sim))
      }
    } %>%
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
})

setGeneric("get_linear_programming_objective", function(scm, outcome, ..., cond = TRUE, as_data = FALSE) {
  standardGeneric("get_linear_programming_objective")
})

setMethod("get_linear_programming_objective", "StructuralCausalModel", function(scm, outcome, ..., cond = TRUE, as_data = FALSE) {
  endog_type_variables <- scm@endogenous_latent_type_variables %$% unique(type_variable)

  cond_types <- scm@types_data %>%
    unnest(outcomes) %>%
    filter({{ cond }}) %>%
    select(latent_type_index)

  get_prob_indices(scm, outcome, ..., as_data = TRUE) %>%
    mutate(latent_type_index, mask = mask * ex_prob) %>%
    group_by_at(vars(latent_type_index, all_of(endog_type_variables))) %>%
    summarize(mask = sum(mask)) %>%
    ungroup() %>% {
      if (as_data) {
        filter(., mask > 0) %>%
          semi_join(cond_types, by = "latent_type_index")
      } else {
        pull(., mask)
      }
    }
})

setGeneric("get_linear_programming_bounds", function(scm, outcome, obs_data, ..., cond = TRUE) {
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
setMethod("get_linear_programming_bounds", "StructuralCausalModel", function(scm, outcome, obs_data, ..., cond = TRUE) {
  types_by_outcomes <- get_latent_type_by_observed_outcomes(scm, obs_data, {{ cond }})

  latent_type_restrict <- types_by_outcomes %>%
    pull(latent_type_mask) %>%
    do.call(rbind, .)

  objective_fun <- scm %>% get_linear_programming_objective(outcome, ...)

  lst(min = lpSolve::lp("min", objective_fun, rbind(latent_type_restrict, 1), "=", c(types_by_outcomes$prob, 1)),
      max = lpSolve::lp("max", objective_fun, rbind(latent_type_restrict, 1), "=", c(types_by_outcomes$prob, 1)))
})
