summ_iter_data <- function(.data) {
  .data %>%
    mutate(estimand_quantiles = map_if(iter_data, ~ !is_null(.), quantilize_est, iter_estimand, wide = TRUE, quant_probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))) %>%
    unnest(estimand_quantiles) %>%
    mutate(mean = map_dbl(iter_data, ~ mean(.$iter_estimand)))
}
