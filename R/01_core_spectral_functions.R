# =============================================================================
# Core functions for the DFT-based multivariate complete-block test
# =============================================================================

article_bandwidth <- function(n, k) {
  min(floor(sqrt(n)), floor((k - 1L) / 2L))
}

base_frequency_grid <- function(T, k) {
  (2 * pi / T) * (0:(k - 1L)) / k
}

compute_DFT_matrix <- function(x, T, lambda_j) {
  x <- as.numeric(x)
  n <- length(x)
  T <- as.integer(T)
  k <- length(lambda_j)
  time_index <- seq_len(n)

  D <- matrix(0, nrow = 2L * T, ncol = k)
  for (j in seq_len(k)) {
    frequencies <- lambda_j[j] + 2 * pi * (0:(T - 1L)) / T
    d <- complex(length = T)
    for (ell in seq_len(T)) {
      d[ell] <- sum(x * exp(-1i * frequencies[ell] * time_index)) / sqrt(n)
    }
    D[, j] <- c(Re(d), Im(d))
  }
  D
}

# The Daniell half-width m is implemented as a symmetric moving average.
# Near the frequency-grid boundaries, the window is truncated and the
# remaining equal weights are renormalized to sum to one. Because equations
# (14)-(15) sum the smoothed periodogram-type matrices over all k centers,
# this function returns the aggregate weight received by each frequency block.
daniell_aggregate_weights <- function(k, m) {
  k <- as.integer(k)
  m <- as.integer(m)
  if (k < 1L) stop("k must be positive.")
  if (m < 0L) stop("m must be non-negative.")

  aggregate <- numeric(k)
  for (center in seq_len(k)) {
    index <- seq.int(max(1L, center - m), min(k, center + m))
    aggregate[index] <- aggregate[index] + 1 / length(index)
  }
  aggregate
}

validate_DFT_list <- function(station_DFT_list) {
  if (length(station_DFT_list) < 2L) stop("At least two series are required.")
  first_dim <- dim(station_DFT_list[[1L]])
  if (length(first_dim) != 2L) stop("Each DFT object must be a matrix.")
  valid <- vapply(station_DFT_list, function(x) identical(dim(x), first_dim), logical(1L))
  if (!all(valid)) stop("All DFT matrices must have identical dimensions.")
  invisible(first_dim)
}

# Equations (6)-(7): block-specific periodogram-type quadratic forms.
compute_HW_blocks <- function(station_DFT_list) {
  dimensions <- validate_DFT_list(station_DFT_list)
  p <- dimensions[1L]
  k <- dimensions[2L]
  G <- length(station_DFT_list)

  array_D <- array(unlist(station_DFT_list), dim = c(p, k, G))
  array_D <- aperm(array_D, c(3L, 1L, 2L))  # G x p x k

  station_means <- matrix(0, nrow = G, ncol = p)
  for (g in seq_len(G)) station_means[g, ] <- rowMeans(array_D[g, , , drop = FALSE][1, , ])

  frequency_means <- matrix(0, nrow = p, ncol = k)
  for (j in seq_len(k)) frequency_means[, j] <- colMeans(array_D[, , j, drop = FALSE][, , 1])

  grand_mean <- colMeans(station_means)
  H_blocks <- array(0, dim = c(p, p, k))
  W_blocks <- array(0, dim = c(p, p, k))

  for (j in seq_len(k)) {
    for (g in seq_len(G)) {
      d_gj <- array_D[g, , j]
      delta_H <- d_gj - frequency_means[, j]
      residual <- d_gj - station_means[g, ] - frequency_means[, j] + grand_mean
      H_blocks[, , j] <- H_blocks[, , j] + tcrossprod(delta_H)
      W_blocks[, , j] <- W_blocks[, , j] + tcrossprod(residual)
    }
  }

  list(H_blocks = H_blocks, W_blocks = W_blocks)
}

# Equations (14)-(16): apply Daniell weighting to the quadratic-form blocks,
# then calculate Wilks' Lambda.
compute_weighted_Wilks <- function(station_DFT_list,
                                   m,
                                   ridge_relative = RIDGE_RELATIVE) {
  dimensions <- validate_DFT_list(station_DFT_list)
  p <- dimensions[1L]
  k <- dimensions[2L]

  blocks <- compute_HW_blocks(station_DFT_list)
  frequency_weights <- daniell_aggregate_weights(k, m)

  H_weighted <- matrix(0, nrow = p, ncol = p)
  W_weighted <- matrix(0, nrow = p, ncol = p)
  for (j in seq_len(k)) {
    H_weighted <- H_weighted + frequency_weights[j] * blocks$H_blocks[, , j]
    W_weighted <- W_weighted + frequency_weights[j] * blocks$W_blocks[, , j]
  }

  H_weighted <- (H_weighted + t(H_weighted)) / 2
  W_weighted <- (W_weighted + t(W_weighted)) / 2
  total <- H_weighted + W_weighted

  reference_scale <- mean(diag(total))
  if (!is.finite(reference_scale) || reference_scale <= 0) reference_scale <- 1
  ridge <- diag(ridge_relative * reference_scale, nrow = p)

  log_det_W <- as.numeric(determinant(W_weighted + ridge, logarithm = TRUE)$modulus)
  log_det_total <- as.numeric(determinant(total + ridge, logarithm = TRUE)$modulus)
  log_lambda <- log_det_W - log_det_total

  list(
    Wilks_Lambda = exp(log_lambda),
    log_Wilks = log_lambda,
    H_weighted = H_weighted,
    W_weighted = W_weighted,
    Daniell_aggregate_weights = frequency_weights
  )
}

permutation_test_weighted <- function(station_DFT_list,
                                      m,
                                      nperm = NPERM,
                                      finite_sample_correction = FINITE_SAMPLE_CORRECTION,
                                      seed = NULL,
                                      keep_null = FALSE) {
  dimensions <- validate_DFT_list(station_DFT_list)
  p <- dimensions[1L]
  k <- dimensions[2L]
  G <- length(station_DFT_list)
  if (!is.null(seed)) set.seed(seed)

  observed <- compute_weighted_Wilks(station_DFT_list, m = m)
  null_log <- numeric(nperm)

  for (b in seq_len(nperm)) {
    permuted <- lapply(seq_len(G), function(g) matrix(0, nrow = p, ncol = k))
    for (j in seq_len(k)) {
      permutation <- sample.int(G)
      for (g in seq_len(G)) {
        permuted[[g]][, j] <- station_DFT_list[[permutation[g]]][, j]
      }
    }
    null_log[b] <- compute_weighted_Wilks(permuted, m = m)$log_Wilks
  }

  # Smaller Wilks' Lambda values provide stronger evidence against H0.
  count <- sum(null_log <= observed$log_Wilks + 1e-12)
  p_value <- if (finite_sample_correction) {
    (1 + count) / (nperm + 1)
  } else {
    count / nperm
  }

  result <- list(
    Wilks_Lambda = observed$Wilks_Lambda,
    log_Wilks = observed$log_Wilks,
    p_value = p_value,
    count_no_larger = count,
    nperm = nperm
  )
  if (keep_null) result$null_log_Wilks <- null_log
  result
}

pairwise_permutation_holm <- function(station_DFT_list,
                                      m,
                                      nperm = NPERM,
                                      alpha = ALPHA,
                                      seed = BASE_SEED) {
  station_names <- names(station_DFT_list)
  if (is.null(station_names) || any(!nzchar(station_names))) {
    station_names <- paste("Station", seq_along(station_DFT_list))
  }

  pairs <- combn(seq_along(station_DFT_list), 2L)
  output <- vector("list", ncol(pairs))
  for (r in seq_len(ncol(pairs))) {
    i <- pairs[1L, r]
    j <- pairs[2L, r]
    test <- permutation_test_weighted(
      station_DFT_list[c(i, j)],
      m = m,
      nperm = nperm,
      seed = seed + r
    )
    output[[r]] <- data.frame(
      station1 = station_names[i],
      station2 = station_names[j],
      comparison = paste(station_names[i], "vs", station_names[j]),
      Wilks_Lambda = test$Wilks_Lambda,
      raw_p = test$p_value,
      stringsAsFactors = FALSE
    )
  }

  output <- do.call(rbind, output)
  output$holm_p <- p.adjust(output$raw_p, method = "holm")
  output$significant_0_05 <- output$holm_p <= alpha
  output
}

make_DFT_list_from_dataframe <- function(data, T, k, center = TRUE) {
  data <- as.data.frame(data, check.names = FALSE)
  data[] <- lapply(data, as.numeric)
  if (anyNA(data)) stop("Input data contain missing or non-numeric values.")

  n <- nrow(data)
  lambda_j <- base_frequency_grid(T, k)
  DFT_list <- lapply(names(data), function(name) {
    x <- data[[name]]
    if (center) x <- x - mean(x)
    compute_DFT_matrix(x, T = T, lambda_j = lambda_j)
  })
  names(DFT_list) <- names(data)

  list(
    station_DFT_list = DFT_list,
    n = n,
    T = T,
    k = k,
    m = article_bandwidth(n, k),
    lambda_j = lambda_j
  )
}

make_zone_average_DFT <- function(station_DFT_list, zones) {
  output <- lapply(zones, function(members) {
    matrices <- station_DFT_list[members]
    if (length(matrices) != length(members) || any(vapply(matrices, is.null, logical(1L)))) {
      stop("A zone contains an unknown station label.")
    }
    Reduce("+", matrices) / length(matrices)
  })
  names(output) <- names(zones)
  output
}

read_station_csv <- function(path) {
  data <- read.csv(path, check.names = FALSE)
  data <- data[, !grepl("^Month_Index$", names(data)), drop = FALSE]
  data[] <- lapply(data, as.numeric)
  if (anyNA(data)) stop("Missing/non-numeric values found in ", path)
  data
}
