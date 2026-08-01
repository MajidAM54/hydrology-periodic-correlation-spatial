# =============================================================================
# Simulation Scenario 4: square-wave PAR(1) (manuscript Tables 8-9)
# =============================================================================

generate_PAR1_square <- function(n,
                                 T,
                                 phi_low,
                                 phi_high,
                                 innovation_variance = 0.001,
                                 innovation = c("Gaussian", "t3")) {
  innovation <- match.arg(innovation)

  epsilon <- if (innovation == "Gaussian") {
    rnorm(n, 0, sqrt(innovation_variance))
  } else {
    # Var(t_3) = 3. This scaling matches the Gaussian innovation variance.
    rt(n, df = 3) * sqrt(innovation_variance / 3)
  }

  x <- numeric(n)
  phi_average <- (phi_low + phi_high) / 2
  x[1L] <- epsilon[1L] / sqrt(max(1 - phi_average^2, 1e-8))

  for (tt in 2:n) {
    phase <- ((tt - 1L) %% T) + 1L
    phi_t <- if (phase <= floor(T / 2)) phi_low else phi_high
    x[tt] <- phi_t * x[tt - 1L] + epsilon[tt]
  }
  x
}

generate_scenario_4_DFT <- function(n,
                                    T,
                                    k,
                                    innovation = "Gaussian",
                                    replicates_per_group = 3L,
                                    coefficient_perturbation_variance = 0.0001,
                                    innovation_variance = 0.001,
                                    coefficient_pairs = rbind(
                                      c(0.10, 0.30),
                                      c(0.40, 0.60),
                                      c(0.70, 0.90)
                                    )) {
  lambda_j <- base_frequency_grid(T, k)
  output <- list()
  index <- 1L

  for (g in seq_len(nrow(coefficient_pairs))) {
    for (r in seq_len(replicates_per_group)) {
      phi_low <- coefficient_pairs[g, 1L] +
        rnorm(1L, 0, sqrt(coefficient_perturbation_variance))
      phi_high <- coefficient_pairs[g, 2L] +
        rnorm(1L, 0, sqrt(coefficient_perturbation_variance))
      phi_low <- min(max(phi_low, -0.99), 0.99)
      phi_high <- min(max(phi_high, -0.99), 0.99)

      x <- generate_PAR1_square(
        n = n,
        T = T,
        phi_low = phi_low,
        phi_high = phi_high,
        innovation_variance = innovation_variance,
        innovation = innovation
      )
      x <- x - mean(x)
      output[[index]] <- compute_DFT_matrix(x, T, lambda_j)
      index <- index + 1L
    }
  }

  names(output) <- paste("Series", seq_along(output))
  output
}

simulate_power_cell_4 <- function(n,
                                  T,
                                  k,
                                  innovation = "Gaussian",
                                  nsim = NSIM_POWER,
                                  nperm = NPERM,
                                  alpha = ALPHA,
                                  seed = BASE_SEED) {
  set.seed(seed)
  m <- article_bandwidth(n, k)
  p_values <- numeric(nsim)

  for (simulation in seq_len(nsim)) {
    DFT_list <- generate_scenario_4_DFT(n, T, k, innovation = innovation)
    p_values[simulation] <- permutation_test_weighted(
      DFT_list,
      m = m,
      nperm = nperm,
      seed = seed + simulation
    )$p_value
  }

  power <- mean(p_values <= alpha)
  data.frame(
    Scenario = "Scenario_4",
    Innovation = innovation,
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

run_scenario_4 <- function() {
  design <- expand.grid(T = T_VALUES, n = N_VALUES, k = K_VALUES,
                        KEEP.OUT.ATTRS = FALSE)

  run_one <- function(i) {
    row <- design[i, ]
    simulate_power_cell_4(
      n = row$n,
      T = row$T,
      k = row$k,
      innovation = "Gaussian",
      seed = BASE_SEED + 500000L + i * 10000L
    )
  }

  if (.Platform$OS.type != "windows" && N_CORES > 1L) {
    output <- parallel::mclapply(seq_len(nrow(design)), run_one, mc.cores = N_CORES)
  } else {
    output <- lapply(seq_len(nrow(design)), run_one)
  }
  output <- do.call(rbind, output)
  write.csv(output,
            file.path(GENERATED_RESULTS_DIR, "simulation_power_scenario_4.csv"),
            row.names = FALSE)

  sensitivity_design <- expand.grid(
    T = c(12L, 3L),
    Innovation = c("Gaussian", "t3"),
    KEEP.OUT.ATTRS = FALSE
  )
  sensitivity <- lapply(seq_len(nrow(sensitivity_design)), function(i) {
    simulate_power_cell_4(
      n = 200L,
      T = sensitivity_design$T[i],
      k = 30L,
      innovation = sensitivity_design$Innovation[i],
      seed = BASE_SEED + 900000L + i * 10000L
    )
  })
  sensitivity <- do.call(rbind, sensitivity)
  write.csv(sensitivity,
            file.path(GENERATED_RESULTS_DIR, "scenario_4_nongaussian_sensitivity.csv"),
            row.names = FALSE)

  list(power = output, nongaussian = sensitivity)
}
