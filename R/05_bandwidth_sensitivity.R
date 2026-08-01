# =============================================================================
# Daniell bandwidth sensitivity analysis (Scenario 3)
# =============================================================================

simulate_bandwidth_cell <- function(m,
                                    n = 200L,
                                    T = 12L,
                                    k = 30L,
                                    nsim = NSIM_POWER,
                                    nperm = NPERM,
                                    alpha = ALPHA,
                                    seed = BASE_SEED) {
  set.seed(seed)
  p_values <- numeric(nsim)

  for (simulation in seq_len(nsim)) {
    DFT_list <- generate_scenario_1_3_DFT(n, T, k, scenario = 3L)
    p_values[simulation] <- permutation_test_weighted(
      DFT_list,
      m = m,
      nperm = nperm,
      seed = seed + simulation
    )$p_value
  }

  power <- mean(p_values <= alpha)
  data.frame(
    Scenario = "Scenario_3",
    T = T,
    n = n,
    k = k,
    m = m,
    default_m = article_bandwidth(n, k),
    nsim = nsim,
    nperm = nperm,
    rejections = sum(p_values <= alpha),
    power = power,
    Monte_Carlo_SE = sqrt(power * (1 - power) / nsim),
    seed = seed
  )
}

run_bandwidth_sensitivity <- function() {
  m_values <- c(0L, 5:14, 34L)
  output <- lapply(seq_along(m_values), function(i) {
    simulate_bandwidth_cell(
      m = m_values[i],
      seed = BASE_SEED + 3000000L + i * 10000L
    )
  })
  output <- do.call(rbind, output)
  write.csv(output,
            file.path(GENERATED_RESULTS_DIR, "bandwidth_sensitivity_scenario_3.csv"),
            row.names = FALSE)
  output
}
