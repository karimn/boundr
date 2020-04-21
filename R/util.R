quantiles_pivot_wider <- function(quantiles_data) {
  pivot_wider(quantiles_data, names_from = per, values_from = est, names_prefix = "per_")
}

quantilize_est <- function(iter, var, wide = TRUE, quant_probs = c(0.05, 0.1, 0.5, 0.9, 0.95)) {
  quant_data <- iter %>%
    pull({{ var }}) %>%
    quantile(probs = quant_probs, names = FALSE) %>%
    enframe(name = NULL, value = "est") %>%
    mutate(per = quant_probs)

  if (wide) {
    quant_data %<>%
      quantiles_pivot_wider()
  }

  return(quant_data)
}

summ_iter_data <- function(.data, quants = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)) {
  .data %>%
    mutate(estimand_quantiles = map_if(iter_data, ~ !is_null(.), quantilize_est, iter_estimand, wide = TRUE, quant_probs = quants)) %>%
    unnest(estimand_quantiles) %>%
    mutate(mean = map_dbl(iter_data, ~ mean(.$iter_estimand)))
}

diagnose <- function(cell, no_sim_diag) {
  tibble(iter_data = list(cell)) %>% {
    if (!no_sim_diag) {
      mutate(.,
             ess_bulk = ess_bulk(cell),
             ess_tail = ess_tail(cell),
             rhat = Rhat(cell))
    } else .
  }
}
