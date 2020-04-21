#' @include util.R structural_model.R
NULL

setClass(
  "EstimationBase",
  slots = c(model = "StructuralCausalModel")
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

setClass(
  "AtomEstimand",
  contains = "DiscreteEstimand",
  slots = list(intervention = "list", outcome = "character", cond = "ANY")
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
  "DiscretizedAtomEstimandCollection",
  contains = c("DiscretizedEstimand", "AtomEstimand")
)

setClass(
  "DiscretizedDiffEstimandCollection",
  contains = c("DiscretizedEstimand", "DiscreteDiffEstimand")
)

setClass(
  "DiscretizedMeanEstimand",
  contains = "Estimand",
  slots = c("group" = "DiscretizedAtomEstimandCollection")
)

setClass(
  "DiscretizedUtilityEstimand",
  contains = "Estimand",
  slots = c("group" = "DiscretizedAtomEstimandCollection")
)

setClass(
  "DiscretizedMeanDiffEstimand",
  contains = "DiffEstimand"
)

setClass(
  "DiscretizedUtilityDiffEstimand",
  contains = "DiffEstimand"
)

#' Collection of S4 estimands
#'
#' @slot estimands \code{data.frame} for internal storage of estimands.
#' @slot est_stan_info list used in Stan model.
#'
#' @export
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

setMethod("initialize", "DiscretizedMeanDiffEstimand", function(.Object, ...) {
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

  outcome_exper <- str_c(est@outcome, " == 1") %>%
      parse_expr() %>%
      as_quosure(env = global_env())

  joint_dist@types_data %>%
    filter(!!outcome_exper) %$%
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
    calculated_est <- calculated_est$r_prob %>% purrr::set_names(calculated_est$name)
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
           function(est, fit, levels, unique_ids, between_entity_diff_info, no_sim_diag, ...) standardGeneric("extract_from_fit"))

setMethod("extract_from_fit", "EstimandCollection", function(est, fit, levels, unique_ids, between_entity_diff_info, no_sim_diag, quants, ...) {
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
        summ_iter_data(quants = quants) %>%
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
        summ_iter_data(quants = quants) %>%
        left_join(between_entity_diff_info, by = c("diff_index", "estimand_id")) %>%
        group_nest(estimand_id, .key = "between_entity_estimands") %>%
        left_join(discrete_estimation, ., by = "estimand_id")
    }
  }

  discrete_estimation %>%
    summ_iter_data(quants = quants) %>%
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
          purrr::set_names(str_c(names(.), ., sep = rename_sep))

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
    purrr::set_names(estimands) %>%
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

setMethod("set_model", "DiscretizedAtomEstimandCollection", function(est, model) {
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

setMethod("get_component_estimands", "DiscretizedAtomEstimandCollection", function(est, next_estimand_id, next_estimand_group_id) {
  discretized_atoms <- est@cutpoint %>%
    imap(function(cutpoint, outcome) {
      new(
        "DiscretizedAtomEstimand",
        name = est@name,
        outcome_group = est@outcome_group,
        cutpoint = cutpoint,
        direction = est@direction,
        intervention = est@intervention,
        outcome = outcome,
        # str_c(outcome, " == 1") %>%
        #   parse_expr() %>%
        #   as_quosure(env = global_env()),
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
      "DiscretizedMeanEstimand",
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

setMethod("get_component_estimands", "DiscretizedDiffEstimandCollection", function(est, next_estimand_id, next_estimand_group_id) {
  left_est_data <- est@left %>% get_component_estimands(next_estimand_id, next_estimand_group_id)
  next_estimand_id <- left_est_data %>% pull(estimand_id) %>% max() %>% add(1)
  next_estimand_group_id <- left_est_data %>% pull(estimand_group_id) %>% max(na.rm = TRUE) %>% add(1)

  right_est_data <- est@right %>% get_component_estimands(next_estimand_id, next_estimand_group_id)
  next_estimand_id <- right_est_data %>% pull(estimand_id) %>% max() %>% add(1)
  next_estimand_group_id <- right_est_data %>% pull(estimand_group_id) %>% max(na.rm = TRUE) %>% add(1)

  mean_diff_data <- list(left_est_data, right_est_data) %>%
    map(filter, map_lgl(est_obj, is, "DiscretizedMeanEstimand")) %>%
    map_df(select, estimand_id, est_obj) %>% {
      left_right_list <- purrr::set_names(.$est_obj, c("left", "right"))
      mean_diff_est_obj <- exec(new, "DiscretizedMeanDiffEstimand", !!!left_right_list)

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
      left_right_list <- purrr::set_names(.$est_obj, c("left", "right"))
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

#' Create S4 instance of estimand for a discrete variable
#'
#' @param outcome Variable name
#' @param ... Intervention
#' @param cond Conditional
#' @param cond_desc Description of conditional
#'
#' @return \code{AtomEstimand} S4 object
#' @export
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
    outcome = outcome,
    # outcome = str_c(outcome, " == 1") %>%
    #   parse_expr() %>%
    #   as_quosure(env = global_env()),
    cond = if (!missing(cond)) enquo(cond) else TRUE,
    condition_string = condition_string)
}

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

#' Create an S4 instance of difference estimand for discrete variables
#'
#' @param left Left estimand to difference from
#' @param right  Right estimand to difference out
#'
#' @return \code{DiscreteDiffEstimand} S4 object
#' @export
build_diff_estimand <- function(left, right) {
  diff_outcome_group <- union(left@outcome_group, right@outcome_group)

  stopifnot(length(diff_outcome_group) == 1)

  new("DiscreteDiffEstimand", name = str_interp("${left@name[1]} - ${right@name[1]}"), outcome_group = diff_outcome_group, left = left, right = right)
}

#' Create an S4 instance for an estimand of a discretized variable
#'
#' @param outcome_group Name of discretized variable group .
#' @param ... Intervention.
#' @param cond Conditional.
#' @param cond_desc String description of conditional.
#'
#' @return \code{DiscretizedAtomEstimandCollection} S4 object
#' @export
build_discretized_atom_estimand <- function(outcome_group, ..., cond, cond_desc) {
  new(
    "DiscretizedAtomEstimandCollection",
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

#' Create an S4 instance for a different estimand for a discretized variable
#'
#' @param left Left estimand, to difference from
#' @param right  Right estimand, to difference out
#'
#' @return \code{DiscretizedDiffEstimandCollection} S4 object
#' @export
build_discretized_diff_estimand <- function(left, right) {
  new(
    "DiscretizedDiffEstimandCollection",
    name = str_interp("${left@name[1]} - ${right@name[1]}"),
    left = left,
    right = right,
    outcome_group = NA_character_)
}
