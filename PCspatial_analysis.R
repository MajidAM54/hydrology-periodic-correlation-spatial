

library(MASS)
library(dplyr)
library(tidyr)


compute_DFT_matrix <- function(x, T, lambda_j) {
  n <- length(x)
  k <- length(lambda_j)
  T_int <- as.integer(T)
  D <- matrix(0, nrow = 2 * T_int, ncol = k)
  for (j in 1:k) {
    lam <- lambda_j[j]
    freqs <- lam + 2 * pi * (0:(T_int - 1)) / T_int
    d <- numeric(T_int)
    for (idx in 1:T_int) {
      d[idx] <- sum(x * exp(-1i * freqs[idx] * (1:n))) / sqrt(n)
    }
    D[, j] <- c(Re(d), Im(d))
  }
  return(D)
}

daniell_smooth_matrix <- function(D_mat, m) {
  if (m == 0) return(D_mat)
  p <- nrow(D_mat)
  k <- ncol(D_mat)
  smoothed <- matrix(0, p, k)
  for (j in 1:k) {
    idx_start <- max(1, j - m)
    idx_end   <- min(k, j + m)
    smoothed[, j] <- rowMeans(D_mat[, idx_start:idx_end, drop = FALSE])
  }
  return(smoothed)
}

compute_weighted_Wilks_Lambda <- function(station_DFT_list, m) {
  G <- length(station_DFT_list)
  k <- ncol(station_DFT_list[[1]])
  p <- nrow(station_DFT_list[[1]])
  
  smoothed_list <- lapply(station_DFT_list, daniell_smooth_matrix, m = m)
  
  
  grand_sum <- rep(0, p)
  for (g in 1:G) grand_sum <- grand_sum + rowSums(smoothed_list[[g]])
  grand_mean <- grand_sum / (G * k)
  
 
  station_means <- matrix(0, G, p)
  freq_means <- matrix(0, k, p)
  for (g in 1:G) station_means[g, ] <- rowMeans(smoothed_list[[g]])
  for (j in 1:k) {
    temp <- sapply(smoothed_list, function(mat) mat[, j])
    freq_means[j, ] <- rowMeans(temp)
  }
  
  # H and W matrices (equations 14‑15)
  H <- matrix(0, p, p)
  W <- matrix(0, p, p)
  for (g in 1:G) {
    for (j in 1:k) {
      resid_H <- smoothed_list[[g]][, j] - freq_means[j, ]
      H <- H + tcrossprod(resid_H)
      resid_W <- smoothed_list[[g]][, j] - station_means[g, ] - freq_means[j, ] + grand_mean
      W <- W + tcrossprod(resid_W)
    }
  }
  
  
  log_det_W <- determinant(W, logarithm = TRUE)$modulus
  log_det_WH <- determinant(W + H, logarithm = TRUE)$modulus
  log_Wilks <- log_det_W - log_det_WH
  Wilks <- exp(log_Wilks)
  return(list(Wilks_Lambda = Wilks, log_Wilks = log_Wilks, H = H, W = W))
}

permutation_test_weighted <- function(station_DFT_list, m, nperm = 1999) {
  G <- length(station_DFT_list)
  k <- ncol(station_DFT_list[[1]])
  p <- nrow(station_DFT_list[[1]])
  
  obs_log <- compute_weighted_Wilks_Lambda(station_DFT_list, m)$log_Wilks
  null_log <- numeric(nperm)
  
  for (b in 1:nperm) {
    perm_list <- vector("list", G)
    for (j in 1:k) {
      col_vals <- lapply(station_DFT_list, function(mat) mat[, j])
      perm_idx <- sample(G)
      perm_cols <- col_vals[perm_idx]
      for (g in 1:G) {
        if (j == 1) perm_list[[g]] <- matrix(0, p, k)
        perm_list[[g]][, j] <- perm_cols[[g]]
      }
    }
    null_log[b] <- compute_weighted_Wilks_Lambda(perm_list, m)$log_Wilks
  }
  p_val <- mean(null_log <= obs_log)
  return(list(obs_stat = exp(obs_log), p_value = p_val, null_distribution = exp(null_log)))
}

pairwise_test_weighted <- function(DFT1, DFT2, m, nperm = 1999) {
  station_list <- list(DFT1, DFT2)
  res <- permutation_test_weighted(station_list, m = m, nperm = nperm)
  return(res$p_value)
}

generate_PC_series <- function(n, T, a, phi, target_sd = 1) {
  sigma_epsilon <- target_sd * sqrt(1 - phi^2)
  t <- 1:n
  g_t <- a * (1 + cos(2 * pi * t / T))
  Y <- arima.sim(model = list(ar = phi), n = n, sd = sigma_epsilon)
  X <- g_t * Y
  return(X)
}

generate_PAR1_square <- function(n, T, phi_low, phi_high, sigma2 = 0.001) {
  X <- numeric(n)
  Z <- rnorm(n, 0, sqrt(sigma2))
  X[1] <- Z[1] / sqrt(1 - ((phi_low + phi_high)/2)^2)
  for (t in 2:n) {
    phase <- ((t - 1) %% T) + 1
    phi_t <- ifelse(phase <= T/2, phi_low, phi_high)
    X[t] <- phi_t * X[t-1] + Z[t]
  }
  return(X)
}

simulate_power_scenario3 <- function(n, T, k, nsim = 100, nperm = 199,
                                     group_a = c(2,4,6),
                                     group_phi = c(0.2,0.5,0.8),
                                     n_series_per_group = c(3,3,3)) {
  
  m <- min(floor(sqrt(n)), floor((k-1)/2))
  lambda_j <- (2 * pi / T) * (0:(k-1)) / k
  G <- length(group_a)
  reject <- 0
  
  for (sim in 1:nsim) {
    station_list <- list()
    idx <- 1
    for (g in 1:G) {
      for (rep in 1:n_series_per_group[g]) {
        a_pert <- group_a[g] + rnorm(1, 0, 0.001)
        phi_pert <- group_phi[g] + rnorm(1, 0, 0.001)
        phi_pert <- min(max(phi_pert, -0.99), 0.99)
        x <- generate_PC_series(n, T, a_pert, phi_pert, target_sd = 1)
        x <- x - mean(x)
        D <- compute_DFT_matrix(x, T, lambda_j)
        station_list[[idx]] <- D
        idx <- idx + 1
      }
    }
    p_val <- permutation_test_weighted(station_list, m = m, nperm = nperm)$p_value
    if (p_val < 0.05) reject <- reject + 1
  }
  return(reject / nsim)
}
