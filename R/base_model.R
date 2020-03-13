#' Base S4 class for model objects
#'
#' @export
setClass("BaseModel")

setGeneric("get_responses", function(r) {
  standardGeneric("get_responses")
})

setGeneric("get_discretized_response_info", function(r, ...) {
  standardGeneric("get_discretized_response_info")
})

setGeneric("get_discretized_cutpoints", function(r, ...) {
  standardGeneric("get_discretized_cutpoints")
})

setGeneric("get_discretized_variable_names", function(r, ...) {
  standardGeneric("get_discretized_variable_names")
})

setGeneric("discretize_continuous_variables", function(r, ...) {
  standardGeneric("discretize_continuous_variables")
})

setGeneric("set_obs_outcomes", function(r, ...) standardGeneric("set_obs_outcomes"))

setGeneric("get_candidates", function(r, analysis_data) standardGeneric("get_candidates"))
