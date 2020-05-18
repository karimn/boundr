#' @include structural_model.R
NULL

#' S4 base class for all responses
#'
#' @slot output Variable name.
#'
#' @export
setClass(
  "BaseResponse",
  slots = c(output = "character"),
  contains = "BaseModel"
)

#' S4 class for a discrete variable's responses
#'
#' @slot input Variable name.
#' @slot finite_states list of latent types and their corresponding structural functions.
#'
#' @export
setClass(
  "Response",
  contains =  "BaseResponse",
  slots = c(input = "character", finite_states = "list")
)

#' S4 base class for all discretized variable responses.
#'
#' @export
setClass("BaseDiscretizedResponse",
         contains = "BaseResponse")

setClass(
  "DiscretizedResponse",
  contains =  c("BaseDiscretizedResponse", "Response"),
  slots = c("cutpoint" = "numeric")
)

#' Create a collection of variables and their responses for a discretized continuous variable
#'
#' @slot child_responses Individual variable for each cutpoint.
#' @slot cutpoints list of cutpoints to discretize and the boundary values.
#' @slot direction Calculate proportion below or above cutpoints.
#' @slot pruning_data \code{data.frame} to exclude some response type combinations.
#'
#' @export
setClass("DiscretizedResponseGroup",
         contains = "BaseDiscretizedResponse",
         slots = c(child_responses = "list", cutpoints = "list", direction = "factor", pruning_data = "data.frame"))

setMethod("get_responses", "Response", function(r) {
  return(list(r))
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

setMethod("get_discretized_response_info", "DiscretizedResponseGroup", function(r, ...) {
  lst(
    cutpoint = r@child_responses %>% map_dbl(get_cutpoint),
    outcome = r@child_responses %>% map_chr(get_output_variable_name),
    direction = r@direction
  )
})

setMethod("get_discretized_cutpoints", "DiscretizedResponseGroup", function(r, ...) {
  r@child_responses %>% { purrr::set_names(map_dbl(., get_cutpoint), map_chr(., get_output_variable_name)) }
})

#' Get names of discretized variables
#'
#' Continuous variables are discretized over a set of cutpoints. This method returns the names of the discrete variables that correspond to these cutpoints.
#'
#' @param r S4 object for discretized response group.
#'
#' @return vector of names.
#' @export
setMethod("get_discretized_variable_names", "DiscretizedResponseGroup", function(r, ...) {
  r@child_responses %>% map_chr(get_output_variable_name)
})

setMethod("discretize_continuous_variables", "DiscretizedResponseGroup", function(r, var_val) {
  map_dbl(r@child_responses, get_cutpoint) %>%
    map(~ if (fct_match(r@direction, "<")) var_val < .x else if (fct_match(r@direction, ">")) var_val > .x) %>%
    purrr::set_names(map_chr(r@child_responses, ~ .@output))
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

  if (!is_empty(r@pruning_data)) {
    if (num_cutpoints > 1) {
      types_grid <- reduce(
        seq(num_cutpoints - 1, ),
        function(cumul, curr_i) {
          inner_join(rename_all(r@pruning_data, ~ str_c(., curr_i)),
                     cumul,
                     purrr::set_names(str_c("hi", curr_i + 1), str_c("low", curr_i)))
        },
        .init = rename_all(r@pruning_data, ~ str_c(., num_cutpoints))) %>%
        select(-hi1) %>%
        distinct() %>%
        purrr::set_names(map_chr(r@child_responses, get_output_variable_name) %>% str_c("r_", .) %>% rev()) %>%
        mutate_all(factor, levels = discretized_types)

      pair_ids <- types_grid %>%
        rev() %>% {
        pmap(list(hi = rev(.[-1]), low = rev(.)[-1]),
             function(hi, low, pruning_data) tibble(hi, low) %>% left_join(pruning_data, by = c("hi", "low")) %>% pull(pair_id),
             pruning_data = r@pruning_data %>% arrange(low, hi) %>% mutate(pair_id = seq(n())))
      } %>%
        purrr::set_names(str_remove(names(.), "^r_") %>% str_replace("(\\d+)$", "pair_id_\\1"))

      types_grid %>%
        mutate(!!!pair_ids)
    } else {
      r@pruning_data %>%
        select(hi) %>%
        distinct() %>%
        purrr::set_names(map_chr(r@child_responses, get_output_variable_name) %>% str_c("r_", .)) %>%
        mutate_all(factor, levels = discretized_types)
    }
  } else {
    types_grid <- r@child_responses %>%
      purrr::set_names(map_chr(., get_output_variable_name) %>% str_c("r_", .)) %>%
      map(get_response_type_names) %>%
      map(factor, levels = discretized_types) %>%
      do.call(expand.grid, .)

    pair_ids <- types_grid %>% {
        list(hi = rev(.[-1]), low = rev(.)[-1])
      } %>%
      pmap(function(hi, low, pruning_data) tibble(hi, low) %>% mutate(pair_id = group_indices(., hi, low)) %>% pull(pair_id)) %>%
      purrr::set_names(str_remove(names(.), "^r_") %>% str_replace("(\\d+)$", "pair_id_\\1"))

    types_grid %>%
      mutate(!!!pair_ids)
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

setMethod("set_obs_outcomes", "Response", function (r, curr_r_type, ...) {
  case_when(!!!imap(r@finite_states, function(fun, name) !!fct_match(curr_r_type, name) ~ !!exec(fun, !!!list2(...), r = curr_r_type))) %>%
    as.integer()
})

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

#' Define an observable variable and its response function
#'
#' @param output Name of variable
#' @param input Antecedent variable names
#' @param ... finite responses to input
#'
#' @return A \code{Response} S4 object
#' @export
define_response <- function(output, input = NA_character_, ...) {
  input_arg <- if (length(input) > 1 || !is.na(input[1L])) input %>% purrr::set_names(rep(NA, length(.)), .)

  input_arg %<>% c("r" = NA)

  new("Response",
      output = output,
      input = input,
      finite_states = list2(...) %>%
        compact() %>%
        map(~ rlang::new_function(input_arg, rlang::f_rhs(.x))))
}

define_discretized_response <- function(output, cutpoint, input = NA_character_, ...) {
  input_arg <- if (length(input) > 1 || !is.na(input[1L])) input %>% purrr::set_names(rep(NA, length(.)), .)

  input_arg %<>% c("r" = NA)

  new("DiscretizedResponse",
      output = output,
      input = input,
      finite_states = list2(...) %>%
        compact() %>%
        map(~ rlang::new_function(input_arg, rlang::f_rhs(.x))),
      cutpoint = cutpoint)
}

#' Define a discretized variable for a continuous variable
#'
#' @param output Variable name.
#' @param cutpoints list of cutpoints to discretize over and the boundary values.
#' @param direction Whether discretization calculates proportion above or below cutpoints.
#' @param input Antecedent variable names.
#' @param ... finite responses to input
#' @slot pruning_data \code{data.frame} to exclude some response type combinations.
#'
#' @return A \code{DiscretizedResponseGroup} S4 object
#' @export
define_discretized_response_group <- function(output, cutpoints, direction = c("<", ">"), input = NA_character_, ..., pruning_data = NULL) {
  stopifnot(!is.unsorted(cutpoints) || length(cutpoints) > 2)

  child_responses <- map2(sprintf("%s_%d", output, seq(length(cutpoints) - 2)), cutpoints[-c(1, length(cutpoints))],
                          define_discretized_response, input = input, ...)

  finite_state_names <- names(child_responses[[1]]@finite_states)

  new_disc <- new(
    "DiscretizedResponseGroup",
    output = output,
    child_responses = child_responses,
    pruning_data = if (!is_null(pruning_data)) {
      pruning_data %>%
        filter_all(~ . %in% finite_state_names) %>%
        mutate_all(factor, levels = finite_state_names)
    } else tibble(),
    cutpoints = as.list(cutpoints),
    direction = arg_match(direction) %>% factor(levels = c("<", ">"))
  )

  new_disc@cutpoints <- new_disc@cutpoints %>%
    purrr::set_names(c(NA, map_chr(new_disc@child_responses, get_output_variable_name), NA))

  return(new_disc)
}

