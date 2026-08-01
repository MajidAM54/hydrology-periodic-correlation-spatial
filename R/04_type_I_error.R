# =============================================================================
# Empirical Type I error under no spatial difference
# =============================================================================

null_parameters_1_3 <- function(scenario) {
  switch(
    as.character(as.integer(scenario)),
    "1" = list(a = 4, phi = 0.8),
    "2" = list(a = 2, phi = 0.5),
    "3" = list(a = 4, phi = 0.5),
    stop("scenario must be 1, 2, or 3")
  )
}

generate_null_DFT_1_3 <- function(n, T, k, scenario, G = 9L) {
  pars <- null_parameters_1_3(scenario)
  lambda_j <- base_frequency_grid(T, k)
  output <- lapply(seq_len(G), function(i) {
    x <- generate_PC_series(n, T, pars$a, pars$phi)
    x <- x - mean(x)
    compute_DFT_matrix(x, T, lambda_j)
  })
  names(output) <- paste("Series", seq_len(G))
  output
}

generate_null_DFT_4 <- function(n, T, k, G = 9L) {
  lambda_j <- base_frequency_grid(T, k)
  output <- lapply(seq_len(G), function(i) {
    x <- generate_PAR1_square(
      n = n,
      T = T,
      phi_low = 0.40,
      phi_high = 0.60,
      innovation_variance = 0.001,
      innovation = "Gaussian"
    )
    x <- x - mean(x)
    compute_DFT_matrix(x, T, lambda_j)
  })
  names(output) <- paste("Series", seq_len(G))
  output
}

simulate_type1_cell <- function(n,
                                T,
                                k,
                                scenario,
                                nsim = NSIM_TYPE1,
                                nperm = NPERM,
                                alpha = ALPHA,
                                seed = BASE_SEED) {
  set.seed(seed)
  m <- article_bandwidth(n, k)
  p_values <- numeric(nsim)

  for (simulation in seq_len(nsim)) {
    DFT_list <- if (scenario <= 3L) {
      generate_null_DFT_1_3(n, T, k, scenario)
    } else {
      generate_null_DFT_4(n, T, k)
    }
    p_values[simulation] <- permutation_test_weighted(
      DFT_list,
      m = m,
      nperm = nperm,
      seed = seed + simulation
    )$p_value
  }

  size <- mean(p_values <= alpha)
  data.frame(
    Scenario = paste0("Scenario_", scenario),
    T = T,
    n = n,
    k = k,
    m = m,
    nsim = nsim,
    nperm = nperm,
    rejections = sum(p_values <= alpha),
    empirical_type_I_error = size,
    Monte_Carlo_SE = sqrt(size * (1 - size) / nsim),
    seed = seed
  )
}

run_type_I_error <- function() {
  design <- expand.grid(
    scenario = 1:4,
    T = T_VALUES,
    n = N_VALUES,
    k = K_VALUES,
    KEEP.OUT.ATTRS = FALSE
  )

  run_one <- function(i) {
    row <- design[i, ]
    simulate_type1_cell(
      n = row$n,
      T = row$T,
      k = row$k,
      scenario = row$scenario,
      seed = BASE_SEED + 2000000L + i * 10000L
    )
  }

  if (.Platform$OS.type != "windows" && N_CORES > 1L) {
    output <- parallel::mclapply(seq_len(nrow(design)), run_one, mc.cores = N_CORES)
  } else {
    output <- lapply(seq_len(nrow(design)), run_one)
  }

  output <- do.call(rbind, output)
  write.csv(output,
            file.path(GENERATED_RESULTS_DIR, "empirical_type_I_error.csv"),
            row.names = FALSE)
  output
}
