summ_iter_data <- . %>% 
  mutate(estimand_quantiles = map_if(iter_data, ~ !is_null(.), quantilize_est, iter_estimand, wide = TRUE, quant_probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))) %>% 
  unnest(estimand_quantiles) %>% 
  mutate(mean = map_dbl(iter_data, ~ mean(.$iter_estimand)))  

# Structural Model Classes ------------------------------------------------

setClass(
  "StructuralCausalModel",
  slots = c(types_data = "data.frame", 
            responses = "list", 
            discretized_responses = "list", 
            exogenous_variables = "character", 
            exogenous_prob = "data.frame", 
            candidate_groups = "data.frame",
            endogenous_latent_type_variables = "data.frame",
            nester = "function")
)

setClass(
  "StructuralCausalModelSimulation",
  contains = "StructuralCausalModel"
)

# Estimation Classes ------------------------------------------------------

setClass(
  "EstimationBase",
  slots = lst(model = "StructuralCausalModel")
)

setClass(
  "ReplicationEstimand",
  contains = "EstimationBase",
  slots = c("name" = "character")
)

setClass(
  "ReplicationCorrelationEstimand",
  contains = "ReplicationEstimand", 
  slots = c("outcome1" = "character", "outcome2" = "character", cond = "ANY", condition_string = "character")
)

setClass(
  "Estimand",
  contains = "EstimationBase",
  slots = list(name = "character", outcome_group = "character", condition_string = "character")
)

setClass(
  "DiffEstimand",
  contains = "Estimand",
  slots = list(left = "Estimand", right = "Estimand")
)

setClass(
  "DiscreteEstimand",
  contains = "Estimand"
)

# setClass(
#   "DiscreteCorrelation",
#   contains = "DiscreteEstimand"
# )

setClass(
  "AtomEstimand",
  contains = "DiscreteEstimand",
  slots = list(intervention = "list", outcome = "ANY", cond = "ANY")
)

setClass(
  "DiscreteDiffEstimand",
  contains = c("DiscreteEstimand", "DiffEstimand")
)

setClass(
  "DiscretizedEstimand",
  contains = "DiscreteEstimand",
  slots = c("direction" = "factor", "cutpoint" = "numeric")
)

setClass(
  "DiscretizedDiffEstimand",
  contains = c("DiscretizedEstimand", "DiscreteDiffEstimand")
)

setClass(
  "DiscretizedAtomEstimand", 
  contains = c("DiscretizedEstimand", "AtomEstimand"),
)

setClass(
  "DiscretizedAtomEstimatorCollection",
  contains = c("DiscretizedEstimand", "AtomEstimand")
)

setClass(
  "DiscretizedDiffEstimatorCollection",
  contains = c("DiscretizedEstimand", "DiscreteDiffEstimand")
)

setClass(
  "DiscretizedMeanEstimator",
  contains = "Estimand",
  slots = c("group" = "DiscretizedAtomEstimatorCollection")
)

setClass(
  "DiscretizedUtilityEstimand",
  contains = "Estimand",
  slots = c("group" = "DiscretizedAtomEstimatorCollection")
)

setClass(
  "DiscretizedMeanDiffEstimator",
  contains = "DiffEstimand"
)

setClass(
  "DiscretizedUtilityDiffEstimand",
  contains = "DiffEstimand"
)

setClass(
  "EstimandCollection",
  contains = "EstimationBase",
  slots = list(estimands = "data.frame", 
               est_stan_info = "list")
)

EstimandResults <- setClass(
  "EstimandResults",
  contains = "tbl_df"
)

setMethod("initialize", "DiscretizedMeanDiffEstimator", function(.Object, ...) {
  .Object <- callNextMethod(.Object, ...)
  
  left_right <- list2(...)
 
  .Object@name <- left_right %>% 
    map_chr(~ .@name) %>% 
    str_c(collapse = " - ") 
  
  .Object@outcome_group <- left_right[[1]]@outcome_group
 
  return(.Object) 
})

setMethod("initialize", "DiscretizedUtilityDiffEstimand", function(.Object, ...) {
  .Object <- callNextMethod(.Object, ...)
  
  left_right <- list2(...)
 
  .Object@name <- left_right %>% 
    map_chr(~ .@name) %>% 
    str_c(collapse = " - ") 
  
  .Object@outcome_group <- left_right[[1]]@outcome_group
 
  return(.Object) 
})

# Sampler -----------------------------------------------------------------

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

setClass("SamplingResults",
         contains = "stanfit",
         slots = c("sampler" = "Sampler"))

setMethod("sampling", "Sampler", function(object, ..., run_type = c("fit", "prior-predict"), save_background_joint_prob = FALSE) {
  run_type <- arg_match(run_type) %>% 
    factor(levels = c("prior-predict", "fit"))
    
  args <- list2(...)
  
  if ("data" %in% names(args)) {
     stop("Sample data cannot be specified. Data is prepared in the sampler constructor.")
  }
  
  # pars <- c("iter_estimand", "iter_level_entity_estimand", "iter_level_entity_estimand_sd", "log_lik", "iter_between_level_entity_diff_estimand", "rep_corr_estimand", "marginal_p_r")
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
  results <- exec(sampling, as(object, "stanmodel"), !!!args) %>% 
    as("SamplingResults")
  
  results@sampler <- object
  
  return(results)
})

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
  results <- exec(vb, as(object, "stanmodel"), !!!args) %>% 
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

setMethod("get_marginal_latent_type_prob", "SamplingResults", function(r, no_sim_diag = TRUE) {
  r %>% 
    as.array(par = "marginal_p_r") %>% 
    plyr::adply(3, diagnose, no_sim_diag = no_sim_diag) %>% 
    tidyr::extract(parameters, "marginal_latent_type_index", "(\\d+)", convert = TRUE) %>% 
    mutate(iter_data = map(iter_data, ~ tibble(iter_p_r= c(.), iter_id = seq(NROW(.) * NCOL(.))))) %>% 
    full_join(r@sampler@endogenous_latent_type_variables, ., by = "marginal_latent_type_index") %>% 
    mutate(estimand_quantiles = map(iter_data, quantilize_est, iter_p_r, wide = TRUE, quant_probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)),
           mean = map_dbl(iter_data, ~ mean(.$iter_p_r))) %>% 
    unnest(estimand_quantiles) 
})

# Response and Response Type Combinations ---------------------------------

setClass(
  "BaseResponse",
  slots = c(output = "character")
)

setClass(
  "Response", 
  contains =  "BaseResponse",
  slots = lst(input = "character", finite_states = "list")
)

setClass("BaseDiscretizedResponse", 
         contains = "BaseResponse") 

setClass(
  "DiscretizedResponse", 
  contains =  c("BaseDiscretizedResponse", "Response"),
  slots = c("cutpoint" = "numeric")
)

setClass("DiscretizedResponseGroup", 
         contains = "BaseDiscretizedResponse", 
         slots = c(child_responses = "list", cutpoints = "list", direction = "factor", pruning_data = "data.frame"))

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

setGeneric("get_responses", function(r) {
  standardGeneric("get_responses")
})

setMethod("get_responses", "StructuralCausalModel", function(r) {
  map(r@responses, get_responses) %>% 
    flatten() %>% 
    compact()
})

setGeneric("create_simulation_analysis_data", function(r) standardGeneric("create_simulation_analysis_data"))

setMethod("create_simulation_analysis_data", "StructuralCausalModelSimulation", function(r) {
  discretized_response_names <- names(r@discretized_responses)
  discrete_response_names <- names(purrr::discard(r@responses, ~ is(., "DiscretizedResponse")))
  
  r@types_data %>% 
    unnest(outcomes) %>% 
    select(all_of(discrete_response_names), matches(str_interp("${discretized_response_names}_\\d")), num_obs) %>%
    map(~ rep(.x, .y), num_obs = .$num_obs) %>%
    as_tibble()
})

setMethod("get_responses", "Response", function(r) {
  return(r)
})

setMethod("get_responses", "DiscretizedResponseGroup", function(r) {
  map(r@child_responses, get_responses) %>% 
    flatten() %>% 
    compact()
})

setGeneric("get_cutpoint", function(r) {
  standardGeneric("get_cutpoint")
})

setMethod("get_cutpoint", "DiscretizedResponse", function(r) r@cutpoint)

setGeneric("get_discretized_response_info", function(r, ...) {
  standardGeneric("get_discretized_response_info")
})

setMethod("get_discretized_response_info", "DiscretizedResponseGroup", function(r, ...) {
  lst(
    cutpoint = r@child_responses %>% map_dbl(get_cutpoint),
    outcome = r@child_responses %>% map_chr(get_output_variable_name),
    direction = r@direction
  )
})

setMethod("get_discretized_response_info", "StructuralCausalModel", function(r, name) {
  r@discretized_responses[[name]] %>% get_discretized_response_info()
})

setGeneric("get_discretized_cutpoints", function(r, ...) {
  standardGeneric("get_discretized_cutpoints")
})

setMethod("get_discretized_cutpoints", "DiscretizedResponseGroup", function(r, ...) {
  r@child_responses %>% { set_names(map_dbl(., get_cutpoint), map_chr(., get_output_variable_name)) }
})

setMethod("get_discretized_cutpoints", "StructuralCausalModel", function(r, name = 1) {
  r@discretized_responses[[name]] %>% get_discretized_cutpoints()
})

setGeneric("get_discretized_variable_names", function(r, ...) {
  standardGeneric("get_discretized_variable_names")
})

setMethod("get_discretized_variable_names", "DiscretizedResponseGroup", function(r, ...) {
  r@child_responses %>% map_chr(get_output_variable_name)
})

setMethod("get_discretized_variable_names", "StructuralCausalModel", function(r, name) {
  r@discretized_responses[[name]] %>% get_discretized_variable_names()
})

setGeneric("discretize_continuous_variables", function(r, ...) {
  standardGeneric("discretize_continuous_variables")
})

setMethod("discretize_continuous_variables", "DiscretizedResponseGroup", function(r, var_val) {
  map_dbl(r@child_responses, get_cutpoint) %>% 
    map(~ if (fct_match(r@direction, "<")) var_val < .x else if (fct_match(r@direction, ">")) var_val > .x) %>% 
    set_names(map_chr(r@child_responses, ~ .@output))
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

# setGeneric("prune_types", function(r, types) {
#   standardGeneric("prune_types")
# })
# 
# setMethod("prune_types", "DiscretizedResponseGroup", function(r, types) {
#   col_names <- str_c("r_", map_chr(r@child_responses, get_output_variable_name))
#   num_child_response <- length(r@child_responses)
#   
#   reduce2(col_names[-num_child_response], col_names[-1], function(curr_pruned, left_col, right_col) { 
#     curr_pruned %>% 
#       semi_join(r@pruning_data, by = unlist(lst(!!sym(left_col) := "low", !!sym(right_col) := "hi")))
#   },
#   .init = types)
# })


setGeneric("get_types_grid", function(r) {
  standardGeneric("get_types_grid")
})

setMethod("get_types_grid", "DiscretizedResponseGroup", function(r) {
  num_cutpoints <- length(r@child_responses)
  discretized_types <- first(r@child_responses)@finite_states %>% names() 
  
  if (num_cutpoints > 1) {
    types_grid <- reduce(
      seq(num_cutpoints - 1, ), 
      function(cumul, curr_i) { 
        inner_join(rename_all(r@pruning_data, ~ str_c(., curr_i)), 
                   cumul, 
                   set_names(str_c("hi", curr_i + 1), str_c("low", curr_i))) 
      },
      .init = rename_all(r@pruning_data, ~ str_c(., num_cutpoints))) %>% 
      select(-hi1) %>% 
      distinct() %>%
      set_names(map_chr(r@child_responses, get_output_variable_name) %>% str_c("r_", .) %>% rev()) %>% 
      mutate_all(factor, levels = discretized_types)
    
    pair_ids <- types_grid %>%
      rev() %>% { 
      pmap(list(hi = rev(.[-1]), low = rev(.)[-1]), 
           function(hi, low, pruning_data) tibble(hi, low) %>% left_join(pruning_data, by = c("hi", "low")) %>% pull(pair_id), 
           pruning_data = r@pruning_data %>% arrange(low, hi) %>% mutate(pair_id = seq(n()))) 
    } %>% 
      set_names(str_remove(names(.), "^r_") %>% str_replace("(\\d+)$", "pair_id_\\1"))
    
    types_grid %>%
      mutate(!!!pair_ids)
  } else {
    r@pruning_data %>% 
      select(hi) %>% 
      distinct() %>% 
      set_names(map_chr(r@child_responses, get_output_variable_name) %>% str_c("r_", .)) %>% 
      mutate_all(factor, levels = discretized_types)
  }
})

setGeneric("get_output_variable_name", function(response) standardGeneric("get_output_variable_name"))

setMethod("get_output_variable_name", "BaseResponse", function(response) response@output)

setGeneric("get_input_variable_names", function(response) standardGeneric("get_input_variable_names"))

setMethod("get_input_variable_names", "Response", function(response) if (length(response@input) > 1 || !is.na(response@input[1L])) response@input)

setGeneric(
  "get_response_type_names", function(response) standardGeneric("get_response_type_names")
)

setMethod("get_response_type_names", "Response", function(response) names(response@finite_states))

setGeneric(
  "set_obs_outcomes", function(r, ...) standardGeneric("set_obs_outcomes")
)

setMethod("set_obs_outcomes", "Response", function (r, curr_r_type, ...) {
  case_when(!!!imap(r@finite_states, function(fun, name) !!fct_match(curr_r_type, name) ~ !!exec(fun, !!!list2(...), r = curr_r_type))) %>% 
    as.integer()
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
    iwalk(function(response, r_type_name) { 
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

setGeneric("get_candidates", function(r, analysis_data) standardGeneric("get_candidates"))

setMethod("get_candidates", "Response", function(r, analysis_data) {
  map(r@finite_states, function(fun) exec(fun, !!!select(analysis_data, r@input), r = rep(NA, nrow(analysis_data)))) %>%  # For each type/class in the current r, produce a response given input values 
    map_if(~ length(.) == 1, ~ rep(., nrow(analysis_data))) %>% 
    bind_cols() %>% 
    mutate(row_index = seq(n())) %>% 
    pivot_longer(cols = -row_index, names_to = "type", values_to = "output") %>% 
    mutate(type = factor(type, levels = names(r@finite_states))) %>% 
    nest_join(analysis_data %>% 
                select(output = r@output) %>% 
                mutate(row_index = seq(n())), 
              ., 
              by = c("row_index", "output"), 
              name = "candidates") %>% 
    pull(candidates)
})

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

setGeneric("get_prob_indices", function(r, outcome, ..., cond = TRUE) standardGeneric("get_prob_indices"))

setMethod("get_prob_indices", "StructuralCausalModel", function(r, outcome, ..., cond = TRUE) {
  r %<>% set_obs_outcomes(..., cond = cond)
  
  r@types_data %>%
    unnest(outcomes) %>%
    mutate(mask = !!outcome) %>%
    pull(mask) %>%
    which() 
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
  compiled_model <- stan_model(file.path("boundr", "bounded.stan"))
  
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
      compiled_model %>% as("Sampler")
    } else {
      stan_model(alternative_model_file) %>% as("Sampler")
    }
    
    new_sampler@endogenous_latent_type_variables <- r@endogenous_latent_type_variables
    
    new_sampler@cv_level <- factor(arg_match(cv_level, values = c(model_levels, NA_character_)), levels = model_levels)
    new_sampler@model_levels <- model_levels
    
    level_ids <- model_levels %>% set_names(seq_along(.), .)
    
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
        mutate_if(~ !is.factor(.), as_factor) %>% 
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

setMethod("create_sampler", "StructuralCausalModel", create_sampler_creator())

setGeneric("get_known_estimands", function(r, estimands) {
  standardGeneric("get_known_estimands")
})

setMethod("get_known_estimands", "StructuralCausalModel", function(r, estimands) {
  estimands %>% calculate_from_known_dgp(r) 
})


define_response <- function(output, input = NA_character_, ...) {
  input_arg <- if (length(input) > 1 || !is.na(input[1L])) input %>% set_names(rep(NA, length(.)), .)
  
  input_arg %<>% c("r" = NA) 
  
  new("Response",
      output = output,
      input = input,
      finite_states = list2(...) %>% 
        compact() %>% 
        map(~ new_function(input_arg, f_rhs(.x))))
}

define_discretized_response <- function(output, cutpoint, input = NA_character_, ...) {
  input_arg <- if (length(input) > 1 || !is.na(input[1L])) input %>% set_names(rep(NA, length(.)), .)
  
  input_arg %<>% c("r" = NA) 
  
  new("DiscretizedResponse",
      output = output,
      input = input,
      finite_states = list2(...) %>% 
        compact() %>% 
        map(~ new_function(input_arg, f_rhs(.x))),
      cutpoint = cutpoint)
}

define_discretized_response_group <- function(output, cutpoints, direction = c("<", ">"), input = NA_character_, ..., pruning_data = NULL) {
  stopifnot(!is.unsorted(cutpoints) || length(cutpoints) > 2)
  
  child_responses <- map2(sprintf("%s_%d", output, seq(length(cutpoints) - 2)), cutpoints[-c(1, length(cutpoints))], 
                          define_discretized_response, input = input, ...)
  
  finite_state_names <- names(child_responses[[1]]@finite_states)
  
  new_disc <- new(
    "DiscretizedResponseGroup",
    output = output,
    child_responses = child_responses,
    pruning_data = pruning_data %>% 
      filter_all(~ . %in% finite_state_names) %>% 
      mutate_all(factor, levels = finite_state_names),
    cutpoints = as.list(cutpoints),
    direction = arg_match(direction) %>% factor(levels = c("<", ">"))
  )
 
  new_disc@cutpoints <- new_disc@cutpoints %>% 
    set_names(c(NA, map_chr(new_disc@child_responses, get_output_variable_name), NA))
  
  return(new_disc) 
}

define_structural_causal_model <- function(..., exogenous_prob, exogenous_variables = setdiff(names(exogenous_prob), "ex_prob")) {
  responses <- list2(...) %>% 
    compact()
  
  if (!missing(exogenous_prob) && !"ex_prob" %in% names(exogenous_prob)) {
    stop("'ex_prob' undefined in exogenous probability.")
  }
  
  leaf_responses <- map(responses, get_responses) %>% flatten() %>% compact()
  discretized_responses <- keep(responses, ~ is(., "DiscretizedResponseGroup")) %>% 
    set_names(map_chr(., get_output_variable_name))
  
  if (length(discretized_responses) > 1) {
    stop("More than one discretized reponse not yet supported.")
  }
  
  comb <- new(
    "StructuralCausalModel",
    exogenous_variables = exogenous_variables,
    exogenous_prob = exogenous_prob %>% 
      mutate(experiment_assignment_id = seq(n())),
    responses = leaf_responses %>% 
      set_names(map_chr(., get_output_variable_name)),
    discretized_responses = discretized_responses,
    types_data = leaf_responses %>% {
        # if (prune_discretized) purrr::discard(., ~ is(., "DiscretizedResponse")) else . 
        purrr::discard(., ~ is(., "DiscretizedResponse")) 
      } %>% 
      set_names(map_chr(., get_output_variable_name) %>% str_c("r_", .)) %>%
      map(get_response_type_names) %>% 
      map(as_factor) %>% 
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


create_prior_predicted_simulation <- function(scm, sample_size, chains, iter, ..., num_sim = 1, num_entities = 1) {
  prior_predict_sampler <- scm %>% create_sampler(analysis_data = NULL, ..., num_sim_unique_entities = num_entities)
  
  prior_predict_fit <- prior_predict_sampler %>% sampling(run_type = "prior-predict", pars = "r_prob", chains = chains, iter = iter)
  
  prior_predict_fit %>% 
    as.data.frame(pars = "r_prob") %>% 
    sample_n(num_sim, replace = FALSE) %>% 
    mutate(iter_id = seq(n())) %>% 
    pivot_longer(cols = -iter_id, names_to = "latent_type_index", values_to = "iter_prob") %>% 
    tidyr::extract(latent_type_index, "latent_type_index", "(\\d+)", convert = TRUE) %>%
    mutate(
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
}

# Estimation ---------------------------------------------------------------

diagnose <- function(cell, no_sim_diag) { 
  tibble(iter_data = list(cell)) %>% 
    {
      if (!no_sim_diag) {
        mutate(., 
               ess_bulk = ess_bulk(cell), 
               ess_tail = ess_tail(cell),
               rhat = Rhat(cell)) 
      } else .
    }
}

setGeneric("calculate_from_known_dgp",
           signature = "est",
           function(est, joint_dist, ...) standardGeneric("calculate_from_known_dgp"))

setMethod("calculate_from_known_dgp", "AtomEstimand", function(est, joint_dist, prob_var = prob) {
  joint_dist@types_data %<>%  
    unnest(outcomes) %>%
    mutate(full_joint_prob = !!est@cond * ex_prob * {{ prob_var }}, 
           full_joint_prob = full_joint_prob / sum(full_joint_prob))
  
  joint_dist %<>% 
    set_obs_outcomes(!!!est@intervention) 
  
  joint_dist@types_data %>% 
    filter(!!est@outcome) %$% 
    sum(full_joint_prob)
})

setMethod("calculate_from_known_dgp", "DiscreteDiffEstimand", function(est, joint_dist, prob_var = prob) {
  return(calculate_from_known_dgp(est@left, joint_dist, {{ prob_var }}) - calculate_from_known_dgp(est@right, joint_dist, {{ prob_var }}))
})


setMethod("calculate_from_known_dgp", "EstimandCollection", function(est, joint_dist, as_df = TRUE) {
  est_calculator <- function(joint_dist, prob_var, calculated_var, set_model = FALSE) {
    if (set_model) {
      joint_dist %<>% 
        right_join(model, by = "latent_type_index")
    }
    
    est@estimands %>% 
      filter(fct_match(est_type, c("atom", "diff"))) %>% 
      mutate(
        {{ calculated_var }} := map_dbl(est_obj, calculate_from_known_dgp, joint_dist, prob_var = {{ prob_var }})
      ) %>% 
      select(-est_obj)
  }
  
  calculated_est <- est_calculator(joint_dist, prob_var = prob, calculated_var = prob)

  if (!as_df) {
    calculated_est <- calculated_est$r_prob %>% set_names(calculated_est$name)
  }
  
  return(calculated_est)
})

setGeneric("get_discrete_estimand_info", function(est) {
  standardGeneric("get_discrete_estimand_info")
})

setMethod("get_discrete_estimand_info", "DiscreteEstimand", function(est) {
  tibble(estimand_name = est@name, outcome_group = est@outcome_group)
})

setMethod("get_discrete_estimand_info", "DiscretizedDiffEstimand", function(est) {
  tibble(estimand_name = est@name, outcome_group = est@outcome_group, cutpoint = est@cutpoint)
})

setMethod("get_discrete_estimand_info", "DiscretizedAtomEstimand", function(est) {
  tibble(estimand_name = est@name, outcome_group = est@outcome_group, cutpoint = est@cutpoint)
})

setGeneric("extract_from_fit", 
           signature = c("est"), 
           function(est, fit, levels, unique_ids, between_entity_diff_info, no_sim_diag = FALSE) standardGeneric("extract_from_fit"))

setMethod("extract_from_fit", "EstimandCollection", function(est, fit, levels, unique_ids, between_entity_diff_info, no_sim_diag = FALSE) {
  discrete_est_info <- est@estimands 
  
  discrete_estimation <- fit %>% 
    as.array(par = "iter_estimand") %>% 
    plyr::adply(3, diagnose, no_sim_diag) %>% 
    tidyr::extract(parameters, "estimand_id", "(\\d+)", convert = TRUE) %>% 
    mutate(iter_data = map(iter_data, ~ tibble(iter_estimand = c(.), iter_id = seq(NROW(.) * NCOL(.))))) %>% 
    full_join(discrete_est_info, ., by = c("estimand_id")) 
  
  stopifnot(!any(map_lgl(discrete_estimation$iter_data, is_empty)))

  if (!missing(levels) && !is_empty(levels) && !missing(unique_ids)) { 
    discrete_level_estimation <- tryCatch(
      as.array(fit, par = "iter_level_entity_estimand"), 
      error = function(err) {
        stop("Failed to find iter_level_entity_estimand parameter.")
        
        return(NULL)
      })
    
    long_entity_ids <- map_dfr(
      levels, 
      ~ unique_ids %>% 
        select(all_of(.x)) %>% 
        distinct() %>% 
        rename("entity_name" = .x) %>% 
        mutate_all(lst(entity_index = as.integer)) %>% 
        mutate(
          entity_name = as.character(entity_name),
          level = .x),
      .id = "level_index") %>% 
      arrange(level_index, entity_index) %>% 
      mutate(long_entity_index = seq(n()))
    
    if (!is_null(discrete_level_estimation)) {
      discrete_estimation <- discrete_level_estimation %>% 
        plyr::adply(3, diagnose, no_sim_diag) %>% 
        tidyr::extract(parameters, c("estimand_id", "long_entity_index"), "(\\d+),(\\d+)", convert = TRUE) %>% 
        mutate(iter_data = map(iter_data, ~ tibble(iter_estimand = c(.), iter_id = seq(NROW(.) * NCOL(.))))) %>% 
        summ_iter_data() %>%
        left_join(long_entity_ids, by = c("long_entity_index")) %>% 
        group_nest(estimand_id, .key = "level_estimands") %>%
        left_join(discrete_estimation, ., by = c("estimand_id"))
    }
    
    if (!is_null(between_entity_diff_info)) {
      between_entity_diff_estimation <- tryCatch(
        as.array(fit, par = "iter_between_level_entity_diff_estimand"), 
        error = function(err) {
          stop("Failed to find iter_between_level_entity_diff_estimand parameter.")
          
          return(NULL)
        })
      
      discrete_estimation <- between_entity_diff_estimation %>% 
        plyr::adply(3, diagnose, no_sim_diag) %>% 
        tidyr::extract(parameters, c("estimand_id", "diff_index"), "(\\d+),(\\d+)", convert = TRUE) %>% 
        mutate(iter_data = map(iter_data, ~ tibble(iter_estimand = c(.), iter_id = seq(NROW(.) * NCOL(.))))) %>% 
        summ_iter_data() %>% 
        left_join(between_entity_diff_info, by = c("diff_index", "estimand_id")) %>% 
        group_nest(estimand_id, .key = "between_entity_estimands") %>% 
        left_join(discrete_estimation, ., by = "estimand_id")
    }
  }
  
  discrete_estimation %>% 
    summ_iter_data() %>%
    new_tibble(nrow = nrow(.), class = "EstimandResults")
})

setGeneric("latex_tablular", function(est) standardGeneric("latex_tablular"))

setMethod("latex_tablular", "EstimandResults", function(est) {
  quants <- str_subset(names(est), "^per_0\\.\\d+$") %>% 
    str_extract("0\\.\\d+") %>% 
    as.numeric()
  num_quants <- length(quants)
  
  results_latex <- est %>% 
    # mutate_at(vars(starts_with("per_0."), rhat, ess_bulk, ess_tail), ~ sprintf("%.2f", .)) %>% 
    # unite("latex", starts_with("per_0."), rhat, ess_bulk, ess_tail, sep = " & ") %>% 
    mutate_at(vars(starts_with("per_0.")), ~ sprintf("%.2f", .)) %>% 
    unite("latex", starts_with("per_0."), sep = " & ") %>% 
    mutate(
      row_color = if_else((row_number() %% 2) == 0, "\\rowcolor{gray!20}", ""),
      latex = str_c(row_color, "$", estimand_name, "$ & ", latex, "\\\\")
    ) %>% 
    pull(latex) %>% 
    str_c(collapse = "\n")
  
    #                  & \\multicolumn{${num_quants}}{c}{Percentiles} & \\multicolumn{3}{c}{Simulation Diagnostics} \\\\
    #                    \\cmidrule(l){2-${num_quants + 1}} \\cmidrule(l){${num_quants + 2}-${num_quants + 4}}
    # Estimand         & ${str_c(quants * 100, collapse = '\\\\% &')}\\% & $\\widehat{R}$ & Bulk ESS & Tail ESS \\\\
 
  str_c( 
    str_interp(
      
    "\\begin{tabular}{l*{${num_quants}}{c}}
     \\toprule
                     & \\multicolumn{${num_quants}}{c}{Percentiles} \\\\
                       \\cmidrule(l){2-${num_quants + 1}} 
    Estimand         & ${str_c(quants * 100, collapse = '\\\\% &')}\\% \\\\
    \\midrule"
  ),
  
  results_latex,
  
  "\\bottomrule
   \\end{tabular}"
   
 )

  # cat(str_interp("& \\multicolumn{${length(quants)}}{c}{Percentiles} \\\\\n"))
  # cat(str_interp("\\cmidrule(l){2-${ength(quants)}}\n"))
  # cat(str_c("Treatment Effect & ", str_c(quants * 100, collapse = "\\% &"), " \\% \\\\\n"))
  # cat("\\midrule\n")
})

setMethod("plot", c("EstimandResults"), function(x, 
                                                 y = "estimand_id", 
                                                 estimands = NULL,
                                                 levels = NULL, 
                                                 prior_results = NULL,
                                                 zero_line = TRUE,
                                                 wrap_after_friendly_name = FALSE) {
  plot_type <- "linerange" 
  
  estimand_friendly_renamer <- identity
  
  if (!is_null(prior_results)) {
    x %<>% bind_rows(fit = ., prior = prior_results, .id = "run_type") 
  }
  
  if (!is_null(estimands)) {
    x %<>% filter(estimand_name %in% estimands)
    
    if (is_named(estimands)) {
      estimand_friendly_renamer <- function(est_names) {
        rename_sep <- if (wrap_after_friendly_name) "\n" else ", "
        
        rename_list <- estimands %>% 
          set_names(str_c(names(.), ., sep = rename_sep)) 
        
        exec(fct_recode, est_names, !!!rename_list)
      }
    }
  } 
  
  ordered_estimand_labels <- if (is_null(estimands)) unique(x$estimand_name) else estimands
 
  if (is_null(levels)) {
    ordered_estimand_labels %<>% rev()
  } 
  
  x %<>% 
    mutate(estimand_id = factor(estimand_id, 
                                levels = estimand_id, 
                                labels = estimand_name) %>% 
             exec(fct_relevel, ., !!!ordered_estimand_labels) %>%
             estimand_friendly_renamer() %>% 
             fct_relabel(. %>% str_replace_all(c("\\[" = "\\\\[", "\\]" = "\\\\]"))))
  
  if (!is_null(levels)) {
    x %<>%
      select(-one_of("iter_data", "ess_bulk", "ess_tail", "rhat"), -starts_with("per_0."), -mean) %>% # Not using population level iter_data 
      unnest(level_estimands) %>% 
      filter(fct_match(level, levels))
  }
  
  plot_obj <- if (plot_type == "density") {
    if (!is_null(levels)) stop("Density plot for levels not yet supported.")
    
    x %>% 
      unnest(iter_data) %>%
      ggplot(aes(x = iter_estimand, y = estimand_id)) +
      ggridges::geom_density_ridges(quantile_lines = TRUE, quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9)) +
      labs(x = "") +
      scale_y_discrete("", labels = latex2exp::TeX) +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 14), strip.text = element_blank()) +
      NULL
  } else {
    plot_position <- ggstance::position_dodgev(height = if (!is_null(prior_results)) 0.75 else 0)
    
    x %>% 
      ggplot(aes_string(y = y, 
                 group = if (!is_null(prior_results)) "run_type" else if (!is_null(levels)) "entity_name" else NA,
                 color = if (!is_null(prior_results)) "run_type" else NULL)) +
      ggstance::geom_crossbarh(aes(x = per_0.5, xmin = per_0.25, xmax = per_0.75), width = 0.25, position = plot_position) +
      ggstance::geom_crossbarh(aes(x = per_0.5, xmin = per_0.1, xmax = per_0.9), width = 0.4, position = plot_position) +
      ggstance::geom_linerangeh(aes(xmin = per_0.05, xmax = per_0.95), fatten = 3, position = plot_position) +
      scale_y_discrete("", labels = latex2exp::TeX) +
      labs(
        x = "",
        caption = "Line range: 90% credible interval. Outer box: 80% credible interval. Inner box: 50% credible interval. 
                   Thick vertical line: median.") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 12), strip.text = element_text(size = 12)) + 
      NULL
  }
  
  if (!is_null(prior_results)) {
    plot_obj <- plot_obj + 
      scale_color_manual("", values = c("fit" = "black", "prior"= "grey"), labels = c("fit" = "Fit", "prior"= "Prior Prediction")) +
      theme(legend.position = "bottom")
  }
  
  if (zero_line) {
    plot_obj <- plot_obj + geom_vline(xintercept = 0, linetype = "dotted") 
  }
  
  return(plot_obj)
})

setGeneric("plot_cdf", function(data, ..., just_data = FALSE, per = seq(0.0, 0.8, 0.2)) {
  standardGeneric("plot_cdf")
})

setMethod("plot_cdf", "EstimandResults", function(data, ..., just_data = FALSE, per = seq(0.0, 0.8, 0.2)) {
  estimands <- list2(...)
  facet <- FALSE
  
  prep_est_data <- function(curr_estimand_name, data) {
    data %>% 
      filter(estimand_name %in% curr_estimand_name) %>% 
      arrange(cutpoint) %>% 
      bind_rows(
        group_by(., estimand_name) %>% 
          filter(estimand_id == max(estimand_id)) %>% 
          mutate(cutpoint = Inf) %>% 
          ungroup()
      ) %>% 
      mutate(percentile = rep(per, each = 2))  
  }
  
  if (is.list(estimands)) {
    facet <- TRUE
  } else {
    estimands %<>% list() 
  }
  
  data %<>% 
    map(estimands, prep_est_data, data = .) %>% 
    bind_rows(.id = "est_group") %>% 
    tidyr::extract(estimand_name, into = c("outcome", "intervention_name", "intervention_value", "condition"), 
                   "Y\\^\\{(\\w+)\\}_\\{(\\w)=(\\d+)\\}(?:.+\\|\\s([^\\]]+))?", remove = FALSE) %>% 
    mutate(
      condition = if_else(!is.na(condition) > 0, str_c("|", condition), ""),
      intervention = str_glue("E\\[Y^{{{outcome}}}_{{{intervention_name}={intervention_value}}} < c{condition}\\]"),
      est_group = str_glue("Pr\\[Y^{{{outcome}}}_{{{intervention_name}}} < c{condition}\\]") 
  )
  
  if (just_data) {
     return(data)
  }
  
  plot_obj <- ggplot(data, aes(percentile, per_0.5)) +
    geom_step(aes(color = intervention_value)) +
    geom_area(aes(ymax = per_0.5, fill = intervention_value), alpha = 0.125, position = position_identity(),
              data = . %>%
                group_by(est_group, intervention_value) %>% 
                bind_rows(mutate(., cutpoint = dplyr::lead(cutpoint), percentile = dplyr::lead(percentile))) %>% 
                arrange(estimand_id) %>% 
                ungroup()
    ) +
    scale_y_continuous("Cumulative Probability") +
    scale_x_continuous(latex2exp::TeX("Percentiles of observed $Y$"), breaks = data$percentile, labels = . %>% multiply_by(100)) +
    scale_color_brewer("Intervention", palette = "Dark2", labels = latex2exp::TeX) +
    scale_color_brewer("Intervention", palette = "Dark2", labels = latex2exp::TeX, aesthetics = "fill") +
    theme(legend.position = "right", panel.grid.minor.x = element_blank()) + 
    NULL
  
  if (facet) {
    plot_obj <- plot_obj +
      facet_wrap(vars(est_group), ncol = 3, labeller = as_labeller(latex2exp::TeX, default = label_parsed)) 
  }
  
  return(plot_obj)
})

setGeneric("plot_cdf_diff", function(data, ..., just_data = FALSE, per = seq(0.0, 0.8, 0.2)) {
  standardGeneric("plot_cdf_diff")
})

# setMethod("plot_cdf_diff", "EstimandResults", function(data, ..., just_data = FALSE, per = seq(0.0, 0.8, 0.2)) {
setMethod("plot_cdf_diff", "EstimandResults", function(data, ..., just_data = FALSE) {
  estimands <- list2(...)
  
  data <- map(estimands, function(curr_estimand_name, data) {
    data %>% 
      filter(estimand_name == curr_estimand_name) %>% 
      arrange(year, cutpoint) %>% 
      group_by(year) %>% 
      bind_rows(
        filter(., estimand_id == max(estimand_id)) %>% 
          mutate(estimand_id = max(estimand_id) + 1, cutpoint = Inf)
      ) %>% 
      ungroup()
      # mutate(percentile = per) 
  },
  data = data) %>% 
    set_names(estimands) %>% 
    bind_rows(.id = "est_group") %>% 
    mutate(est_group = str_replace_all(est_group, c("\\(?E" = "Pr", "([\\[\\]])" = "\\\\\\1", "\\)$" = "")))
  
  if (just_data) return(data)
  
  plot_obj <- data %>% { 
    ggplot(., aes(percentile, per_0.5)) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_step(direction = "hv") +
      geom_ribbon(aes(percentile, ymin = per_0.25, ymax = per_0.75), alpha = 0.25,
                  data = . %>% 
                    group_by(est_group) %>% 
                    bind_rows(mutate(., cutpoint = dplyr::lead(cutpoint), percentile = dplyr::lead(percentile))) %>% 
                    arrange(estimand_id) %>% 
                    ungroup()
      ) +
      geom_ribbon(aes(percentile, ymin = per_0.1, ymax = per_0.9), alpha = 0.25,
                  data = . %>% 
                    group_by(est_group) %>% 
                    bind_rows(mutate(., cutpoint = dplyr::lead(cutpoint), percentile = dplyr::lead(percentile))) %>% 
                    arrange(estimand_id) %>% 
                    ungroup()
      ) +
      geom_ribbon(aes(percentile, ymin = per_0.05, ymax = per_0.95), alpha = 0.25,
                  data = . %>% 
                    group_by(est_group) %>% 
                    bind_rows(mutate(., cutpoint = dplyr::lead(cutpoint), percentile = dplyr::lead(percentile))) %>% 
                    arrange(estimand_id) %>% 
                    ungroup()
      ) +
      scale_y_continuous("Cumulative Probability Difference") +
      scale_x_continuous(latex2exp::TeX("Percentiles of observed $Y$"), breaks = .$percentile, labels = . %>% multiply_by(100)) +
      labs(caption = "The black line shows the posterior median while the grey ribbons represent the 50, 80, and 90% credible intervals.") +
      theme(legend.position = "top", panel.grid.minor.x = element_blank()) + 
      NULL
  }
  
  if (length(estimands) > 1) {
    plot_obj <- plot_obj +
      facet_wrap(vars(est_group), ncol = 2, labeller = as_labeller(latex2exp::TeX, default = label_parsed), scales = "free") 
  }
  
  return(plot_obj)
})

setGeneric("obs_outcomes_setter<-", 
           function(est, value) standardGeneric("obs_outcomes_setter<-"),
           signature = "est")

setMethod("obs_outcomes_setter<-", "EstimationBase", function(est, value) {
  est@obs_outcomes_setter <- value
  
  return(est)
})

setMethod("obs_outcomes_setter<-", "DiscreteDiffEstimand", function(est, value) {
  est <- callNextMethod()
  
  obs_outcomes_setter(est@left) <- value
  obs_outcomes_setter(est@right) <- value
  
  return(est)
})

setGeneric("get_range", function(est) {
  standardGeneric("get_range")
})

setMethod("get_range", "Estimand", function(est) c(0, 1))

setMethod("get_range", "DiscretizedEstimand", function(est) est@model@discretized_types[[1]]@cutpoints %>% { c(first(.), last(.)) } )

setGeneric("set_model", 
           function(est, model) standardGeneric("set_model"),
           signature = "est")

setMethod("set_model", "EstimationBase", function(est, model) {
  est@model <- model
  
  return(est)
})

setMethod("set_model", "DiffEstimand", function(est, model) {
  est <- callNextMethod()
  
  est@left %<>% set_model(model)
  est@right %<>% set_model(model)
  
  return(est)
})

setMethod("set_model", "DiscretizedAtomEstimatorCollection", function(est, model) {
  est <- callNextMethod()
  
  discretized_type <- model@discretized_responses[[est@outcome_group]]
  est@cutpoint <- get_discretized_cutpoints(discretized_type)
  est@direction <- discretized_type@direction
  est@outcome = names(est@cutpoint)
  
  intervention_string <- est@intervention %>% 
    imap_chr(~ str_interp("${.y}=${.x}")) %>% 
    str_c(collapse = ",")
  
  est@name = str_interp("Pr[Y^{${est@outcome_group}}_{${intervention_string}} ${est@direction} c${est@condition_string}]")  
  
  return(est)
})

setGeneric("get_component_estimands",
           function(est, next_estimand_id, next_estimand_group_id) standardGeneric("get_component_estimands"))

setMethod("get_component_estimands", "Estimand", function(est, next_estimand_id, next_estimand_group_id) {
  est %>% 
    get_discrete_estimand_info() %>% 
    mutate(
      est_obj = list(est),
      estimand_id = next_estimand_id,
      estimand_group_id = NA_integer_
    )
}) 

setMethod("get_component_estimands", "DiscreteDiffEstimand", function(est, next_estimand_id) {
  left_est_data <- est@left %>% get_component_estimands(next_estimand_id)
  next_estimand_id <- left_est_data %>% pull(estimand_id) %>% max() %>% add(1)
  
  right_est_data <- est@right %>% get_component_estimands(next_estimand_id)
  next_estimand_id <- right_est_data %>% pull(estimand_id) %>% max() %>% add(1)
  
  est %>% 
    get_discrete_estimand_info() %>% 
    mutate(
      est_obj = list(est),
      estimand_id = next_estimand_id,
      estimand_group_id = NA_integer_,
      estimand_id_left = left_est_data$estimand_id, 
      estimand_id_right = right_est_data$estimand_id, 
    ) %>% 
    bind_rows(left_est_data, right_est_data)
}) 

setMethod("get_component_estimands", "DiscretizedAtomEstimatorCollection", function(est, next_estimand_id, next_estimand_group_id) {
  discretized_atoms <- est@cutpoint %>% 
    imap(function(cutpoint, outcome) {
      new(
        "DiscretizedAtomEstimand",
        name = est@name,
        outcome_group = est@outcome_group,
        cutpoint = cutpoint,
        direction = est@direction,
        intervention = est@intervention,
        outcome = str_c(outcome, " == 1") %>%
          parse_expr() %>%
          as_quosure(env = global_env()),
        cond = est@cond
      ) 
    }) %>% 
    map(set_model, est@model) %>% 
    map_df(~ mutate(get_discrete_estimand_info(.), est_obj = list(.))) %>% 
    mutate(
      estimand_id = seq(next_estimand_id, next_estimand_id + n() - 1),
      estimand_group_id = next_estimand_group_id
    )
  
  next_estimand_id <- discretized_atoms %>% 
    pull(estimand_id) %>% 
    max() %>% 
    add(1)
  
  atom_mean_est <- tibble(
    estimand_name = str_remove(est@name, "\\s*[<>]\\s*c\\s*") %>% str_replace("Pr\\[", "E["),
    outcome_group = est@outcome_group,
    estimand_id = next_estimand_id,
    mean_estimand_group_id = next_estimand_group_id,
    est_obj = list(new(
      "DiscretizedMeanEstimator",
      name = estimand_name,
      outcome_group = outcome_group,
      group = est
    ))
  )
  
  next_estimand_id <- atom_mean_est %>% 
    pull(estimand_id) %>% 
    max() %>% 
    add(1)
  
  utility_est <- tibble(
    estimand_name = str_remove(est@name, "\\s*[<>]\\s*c\\s*") %>% str_replace("Pr\\[", "EU["),
    outcome_group = est@outcome_group,
    estimand_id = next_estimand_id,
    mean_estimand_group_id = next_estimand_group_id,
    est_obj = list(new(
      "DiscretizedUtilityEstimand",
      name = estimand_name,
      outcome_group = outcome_group,
      group = est
    ))
  )
  
  bind_rows(discretized_atoms, atom_mean_est, utility_est)
})

setMethod("get_component_estimands", "DiscretizedDiffEstimatorCollection", function(est, next_estimand_id, next_estimand_group_id) {
  left_est_data <- est@left %>% get_component_estimands(next_estimand_id, next_estimand_group_id)
  next_estimand_id <- left_est_data %>% pull(estimand_id) %>% max() %>% add(1)
  next_estimand_group_id <- left_est_data %>% pull(estimand_group_id) %>% max(na.rm = TRUE) %>% add(1)
  
  right_est_data <- est@right %>% get_component_estimands(next_estimand_id, next_estimand_group_id)
  next_estimand_id <- right_est_data %>% pull(estimand_id) %>% max() %>% add(1)
  next_estimand_group_id <- right_est_data %>% pull(estimand_group_id) %>% max(na.rm = TRUE) %>% add(1)
  
  mean_diff_data <- list(left_est_data, right_est_data) %>% 
    map(filter, map_lgl(est_obj, is, "DiscretizedMeanEstimator")) %>% 
    map_df(select, estimand_id, est_obj) %>% {
      left_right_list <- set_names(.$est_obj, c("left", "right"))
      mean_diff_est_obj <- exec(new, "DiscretizedMeanDiffEstimator", !!!left_right_list)
    
      transmute(., estimand_id, name = c("left", "right")) %>% 
        pivot_wider(values_from = estimand_id, names_prefix = "estimand_id_") %>% 
        mutate(est_obj = list(mean_diff_est_obj),
               estimand_id = next_estimand_id,
               outcome_group = mean_diff_est_obj@outcome_group,
               estimand_name = mean_diff_est_obj@name)  
    }
  
  next_estimand_id <- mean_diff_data %>% pull(estimand_id) %>% max() %>% add(1)
  
  utility_diff_data <- list(left_est_data, right_est_data) %>% 
    map(filter, map_lgl(est_obj, is, "DiscretizedUtilityEstimand")) %>% 
    map_df(select, estimand_id, est_obj) %>% {
      left_right_list <- set_names(.$est_obj, c("left", "right"))
      utility_diff_est_obj <- exec(new, "DiscretizedUtilityDiffEstimand", !!!left_right_list)
    
      transmute(., estimand_id, name = c("left", "right")) %>% 
        pivot_wider(values_from = estimand_id, names_prefix = "estimand_id_") %>% 
        mutate(est_obj = list(utility_diff_est_obj),
               estimand_id = next_estimand_id,
               outcome_group = utility_diff_est_obj@outcome_group,
               estimand_name = utility_diff_est_obj@name)  
    }
  
  next_estimand_id <- utility_diff_data %>% pull(estimand_id) %>% max() %>% add(1)
 
  inner_join(
    select(left_est_data, est_obj, cutpoint, estimand_id) %>% filter(!is.na(cutpoint)), 
    select(right_est_data, est_obj, cutpoint, estimand_id) %>% filter(!is.na(cutpoint)),
    by = c("cutpoint"), suffix = c("_left", "_right")
  ) %>% 
    mutate(
      estimand_id = seq(next_estimand_id, next_estimand_id + n() - 1),
      estimand_group_id = next_estimand_group_id,
      est_obj = pmap(lst(est_obj_left, est_obj_right, cutpoint), function(est_obj_left, est_obj_right, cutpoint) {
        new(
          "DiscretizedDiffEstimand",
          name = str_c(est_obj_left@name, " - ", est_obj_right@name),
          outcome_group = est_obj_left@outcome_group,
          cutpoint = cutpoint,
          left = est_obj_left,
          right = est_obj_right
        )
      }),
      
      outcome_group = map_chr(est_obj, ~ .@outcome_group),
      estimand_name = map_chr(est_obj, ~ .@name),
    ) %>% 
    select(-matches("(est_obj)_(left|right)$")) %>% 
    bind_rows(left_est_data, right_est_data, mean_diff_data, utility_diff_data)
})

setGeneric("num_estimands",
           function(est, est_class = "ANY") standardGeneric("num_estimands"))

setMethod("num_estimands", "EstimandCollection", function(est, est_class = NA) {
  if (all(is.na(est_class))) {
    nrow(est@estimands)
  } else {
    filter(est@estimands, fct_match(est_type, est_class)) %>% 
      nrow()
  }
})

setGeneric("get_stan_data_structures",
           function(est, ...) standardGeneric("get_stan_data_structures"),
           signature = "est")

setMethod("get_stan_data_structures", "AtomEstimand", function(est) {
  abducted_mask <- est@model@types_data %>%
    unnest(outcomes) %>%
    mutate(abducted_mask = !!est@cond) %>%
    pull(abducted_mask) 
  
  est_prob_index <- est@model %>%
    get_prob_indices(est@outcome, !!!est@intervention)
  
  lst(
    abducted_prob_size = if (all(abducted_mask)) 0 else sum(abducted_mask),
    abducted_prob_index = if (abducted_prob_size > 0) which(abducted_mask), # Row major index
    
    est_prob_index = if (abducted_prob_size > 0) intersect(est_prob_index, abducted_prob_index) else est_prob_index,
    est_prob_size = as.array(length(est_prob_index)),
  )
})

setMethod("get_stan_data_structures", "Estimand", function(est) {
  return(NULL)
})

setMethod("get_stan_data_structures", "EstimandCollection", function(est, cores = 1) {
  map_fun <- if (cores > 1) partial(pbmcapply::pbmclapply, ignore.interactive = TRUE, mc.silent = TRUE, mc.cores = cores) else map
  
  stan_data_structures <- est@estimands$est_obj %>% 
    map_fun(get_stan_data_structures) %>% 
    reduce(function(accum, to_add) list_merge(accum, !!!to_add)) 
  
  return(stan_data_structures)
})

setGeneric("get_stan_info", function(est) standardGeneric("get_stan_info"))

setMethod("get_stan_info", "EstimandCollection", function(est) est@est_stan_info)

setGeneric("add_between_level_entity_diff_estimands", function(est, levels, analysis_data) {
  standardGeneric("add_between_level_entity_diff_estimands")
})

setMethod("add_between_level_entity_diff_estimands", "EstimandCollection", function(est, levels, analysis_data) {
  # Do nothing
})

build_atom_estimand <- function(outcome, ..., cond, cond_desc) {
  intervention <- list2(...)
  intervention_string <- if (!is_empty(intervention)) { 
    intervention %>% 
      imap_chr(~ str_interp("${.y}=${.x}")) %>% 
      str_c(collapse = ",") %>% 
      str_c("_{", ., "}")
  } else ""
  
  condition_string <- if (!missing(cond_desc)) {
      str_c(" | ", cond_desc)
    } else if (!missing(cond)) { 
      as_label(enquo(cond)) %>% 
        str_split("\\s*&\\s*", simplify = TRUE) %>%
        c() %>% 
        map_chr(str_to_upper) %>%
        map(str_replace_all, fixed("=="), "=") %>%
        str_c(., collapse = ", ") %>% 
        str_c(" | ", .)  
    } else ""
  
  new(
    "AtomEstimand", 
    name = str_interp("E[${str_to_upper(outcome)}${intervention_string}${condition_string}]"), 
    outcome_group = outcome,
    intervention = intervention, 
    outcome = str_c(outcome, " == 1") %>% 
      parse_expr() %>% 
      as_quosure(env = global_env()),
    cond = if (!missing(cond)) enquo(cond) else TRUE,
    condition_string = condition_string)
}


# build_discrete_correlation_estimand <- function(outcome1, outcome2, cond) {
#   
# }


build_replication_correlation_estimand <- function(outcome1, outcome2, cond = NA_character_) {
  condition_string <- if (!is.na(cond)) { 
    str_c(" | ", str_to_upper(cond))
    # as_label(enquo(cond)) %>% 
    #   str_split("\\s*&\\s*", simplify = TRUE) %>%
    #   c() %>% 
    #   map_chr(str_to_upper) %>%
    #   # map(str_replace_all, fixed("=="), "=") %>%
    #   # str_c(., collapse = ", ") %>% 
    #   str_c(" | ", .)
  } else ""
  
  new(
    "ReplicationCorrelationEstimand",
    
    name = str_interp("SampleCor[${str_to_upper(outcome1)}, ${str_to_upper(outcome2)}${condition_string}]"), 
    outcome1 = outcome1,
    outcome2 = outcome2,
    # cond = if (!missing(cond)) enquo(cond) else TRUE,
    cond = cond,
    condition_string = condition_string
  )
}

build_diff_estimand <- function(left, right) {
  diff_outcome_group <- union(left@outcome_group, right@outcome_group)
  
  stopifnot(length(diff_outcome_group) == 1)
  
  new("DiscreteDiffEstimand", name = str_interp("${left@name[1]} - ${right@name[1]}"), outcome_group = diff_outcome_group, left = left, right = right)
}

build_discretized_atom_estimand <- function(outcome_group, ..., cond, cond_desc) {
  new(
    "DiscretizedAtomEstimatorCollection",
    outcome_group = outcome_group,
    intervention = list(...), 
    cond = if (!missing(cond)) enquo(cond) else TRUE,
    condition_string = if (!missing(cond_desc)) {
      str_c(" | ", cond_desc)
    } else if (!missing(cond)) { 
      as_label(enquo(cond)) %>% 
        str_split("\\s*&\\s*", simplify = TRUE) %>%
        c() %>% 
        map_chr(str_to_upper) %>%
        map(str_replace_all, fixed("=="), "=") %>%
        str_c(., collapse = ", ") %>% 
        str_c(" | ", .)  
    } else ""
  )
}

build_discretized_diff_estimand <- function(left, right) {
  new(
    "DiscretizedDiffEstimatorCollection",
    name = str_interp("${left@name[1]} - ${right@name[1]}"),
    left = left,
    right = right,
    outcome_group = NA_character_)
}

  
build_estimand_collection <- function (model, ..., utility = NA_real_, cores = 1) {
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
          map_lgl(est_obj, is, "DiscretizedMeanEstimator") ~ "mean",
          map_lgl(est_obj, is, "DiscretizedUtilityEstimand") ~ "utility",
          map_lgl(est_obj, is, "DiscretizedMeanDiffEstimator") ~ "mean-diff",
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
}

build_discretized_estimands <- function(types, outcome, fun) {
  types %>% 
    get_discretized_response_info(outcome) %>% {
      pmap(.[c("cutpoint", "outcome")], fun, outcome_group = outcome, direction = .$direction) 
    } %>% 
    flatten()
}
