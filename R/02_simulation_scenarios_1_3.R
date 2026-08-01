# =============================================================================
# Simulation Scenarios 1-3 (manuscript Tables 2-7)
# =============================================================================

simulate_stationary_AR1 <- function(n, phi, stationary_variance = 1) {
  innovation_sd <- sqrt(stationary_variance * (1 - phi^2))
  as.numeric(arima.sim(model = list(ar = phi), n = n, sd = innovation_sd))
}

generate_PC_series <- function(n, T, a, phi) {
  time_index <- seq_len(n)
  periodic_amplitude <- a * (1 + cos(2 * pi * time_index / T))
  periodic_amplitude * simulate_stationary_AR1(n, phi, stationary_variance = 1)
}

scenario_parameters <- function(scenario) {
  switch(
    as.character(as.integer(scenario)),
    "1" = list(a = c(2, 4, 6), phi = c(0.8, 0.8, 0.8)),
    "2" = list(a = c(2, 2, 2), phi = c(0.2, 0.5, 0.8)),
    "3" = list(a = c(2, 4, 6), phi = c(0.2, 0.5, 0.8)),
    stop("scenario must be 1, 2, or 3")
  )
}

generate_scenario_1_3_DFT <- function(n,
                                      T,
                                      k,
                                      scenario,
                                      replicates_per_group = 3L,
                                      perturbation_variance = 0.001) {
  parameters <- scenario_parameters(scenario)
  lambda_j <- base_frequency_grid(T, k)
  output <- list()
  index <- 1L

  for (g in seq_along(parameters$a)) {
    for (r in seq_len(replicates_per_group)) {
      a_gr <- parameters$a[g] + rnorm(1L, 0, sqrt(perturbation_variance))
      phi_gr <- parameters$phi[g] + rnorm(1L, 0, sqrt(perturbation_variance))
      phi_gr <- min(max(phi_gr, -0.99), 0.99)

      x <- generate_PC_series(n, T, a_gr, phi_gr)
      x <- x - mean(x)
      output[[index]] <- compute_DFT_matrix(x, T, lambda_j)
      index <- index + 1L
    }
  }

  names(output) <- paste("Series", seq_along(output))
  output
}

simulate_power_cell_1_3 <- function(n,
                                    T,
                                    k,
                                    scenario,
                                    nsim = NSIM_POWER,
                                    nperm = NPERM,
                                    alpha = ALPHA,
                                    seed = BASE_SEED) {
  set.seed(seed)
  m <- article_bandwidth(n, k)
  p_values <- numeric(nsim)

  for (simulation in seq_len(nsim)) {
    DFT_list <- generate_scenario_1_3_DFT(n, T, k, scenario)
    p_values[simulation] <- permutation_test_weighted(
      DFT_list,
      m = m,
      nperm = nperm,
      seed = seed + simulation
    )$p_value
  }

  power <- mean(p_values <= alpha)
  data.frame(
    Scenario = paste0("Scenario_", scenario),
    T = T,
    n = n,
    k = k,
    m = m,
    nsim = nsim,
    nperm = nperm,
    rejections = sum(p_values <= alpha),
    power = power,
    Monte_Carlo_SE = sqrt(power * (1 - power) / nsim),
    seed = seed
  )
}

run_scenarios_1_3 <- function() {
  design <- expand.grid(
    scenario = 1:3,
    T = T_VALUES,
    n = N_VALUES,
    k = K_VALUES,
    KEEP.OUT.ATTRS = FALSE
  )

  run_one <- function(i) {
    row <- design[i, ]
    simulate_power_cell_1_3(
      n = row$n,
      T = row$T,
      k = row$k,
      scenario = row$scenario,
      seed = BASE_SEED + i * 10000L
    )
  }

  if (.Platform$OS.type != "windows" && N_CORES > 1L) {
    output <- parallel::mclapply(seq_len(nrow(design)), run_one, mc.cores = N_CORES)
  } else {
    output <- lapply(seq_len(nrow(design)), run_one)
  }

  output <- do.call(rbind, output)
  write.csv(output,
            file.path(GENERATED_RESULTS_DIR, "simulation_power_scenarios_1_3.csv"),
            row.names = FALSE)
  output
}
