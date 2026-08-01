# =============================================================================
# Conventional comparison: Friedman + Conover post-hoc + Holm correction
# =============================================================================

extract_conover_pairs <- function(test_object) {
  p_matrix <- test_object$p.value
  index <- which(!is.na(p_matrix), arr.ind = TRUE)
  output <- lapply(seq_len(nrow(index)), function(i) {
    station1 <- rownames(p_matrix)[index[i, 1L]]
    station2 <- colnames(p_matrix)[index[i, 2L]]
    data.frame(
      station1 = station1,
      station2 = station2,
      comparison = paste(station1, "vs", station2),
      raw_p = p_matrix[index[i, 1L], index[i, 2L]],
      stringsAsFactors = FALSE
    )
  })
  output <- do.call(rbind, output)
  output$holm_p <- p.adjust(output$raw_p, method = "holm")
  output$significant_0_05 <- output$holm_p <= ALPHA
  output
}

cohen_kappa_binary <- function(x, y) {
  x <- as.logical(x)
  y <- as.logical(y)
  if (length(x) != length(y)) stop("Decision vectors must have equal lengths.")
  observed <- mean(x == y)
  px <- mean(x)
  py <- mean(y)
  expected <- px * py + (1 - px) * (1 - py)
  if (abs(1 - expected) < .Machine$double.eps) return(NA_real_)
  (observed - expected) / (1 - expected)
}

pair_key <- function(a, b) paste(sort(c(a, b)), collapse = "__")

align_pairwise_decisions <- function(conventional, spectral) {
  conventional$key <- mapply(pair_key, conventional$station1, conventional$station2)
  spectral$key <- mapply(pair_key, spectral$station1, spectral$station2)
  merge(
    conventional[, c("key", "significant_0_05")],
    spectral[, c("key", "significant_0_05")],
    by = "key",
    suffixes = c("_conventional", "_spectral")
  )
}

run_conventional_one <- function(data, dataset_name, spectral_pairwise_file) {
  if (!requireNamespace("PMCMRplus", quietly = TRUE)) {
    stop("Package 'PMCMRplus' is required. Run R/00_install_packages.R first.")
  }

  data <- as.data.frame(data, check.names = FALSE)
  data[] <- lapply(data, as.numeric)
  if (anyNA(data)) stop("Missing/non-numeric values in ", dataset_name)

  friedman <- friedman.test(as.matrix(data))
  conover <- PMCMRplus::frdAllPairsConoverTest(
    as.matrix(data),
    p.adjust.method = "none"
  )
  pairwise <- extract_conover_pairs(conover)

  write.csv(pairwise,
            file.path(GENERATED_RESULTS_DIR,
                      paste0(dataset_name, "_Friedman_Conover_Holm_pairwise.csv")),
            row.names = FALSE)

  spectral <- read.csv(spectral_pairwise_file, stringsAsFactors = FALSE)
  aligned <- align_pairwise_decisions(pairwise, spectral)
  agreement_vector <- aligned$significant_0_05_conventional ==
    aligned$significant_0_05_spectral

  summary <- data.frame(
    Dataset = dataset_name,
    Number_of_blocks = nrow(data),
    Number_of_stations = ncol(data),
    Friedman_statistic = unname(friedman$statistic),
    Friedman_df = unname(friedman$parameter),
    Friedman_p_value = friedman$p.value,
    Significant_pairwise_after_Holm = sum(pairwise$significant_0_05),
    Total_pairwise = nrow(pairwise),
    Agreement_count = sum(agreement_vector),
    Agreement_total = length(agreement_vector),
    Cohen_kappa = cohen_kappa_binary(
      aligned$significant_0_05_conventional,
      aligned$significant_0_05_spectral
    )
  )

  write.csv(aligned,
            file.path(GENERATED_RESULTS_DIR,
                      paste0(dataset_name, "_conventional_vs_spectral_decisions.csv")),
            row.names = FALSE)
  list(summary = summary, pairwise = pairwise, aligned = aligned)
}

run_conventional_comparisons <- function() {
  precipitation <- read_station_csv(file.path(DATA_DIR, "precipitation_monthly_180.csv"))
  runoff <- read_station_csv(file.path(DATA_DIR, "runoff_monthly_180.csv"))

  precipitation_result <- run_conventional_one(
    precipitation,
    "precipitation",
    file.path(GENERATED_RESULTS_DIR, "precipitation_spectral_pairwise_Holm.csv")
  )
  runoff_result <- run_conventional_one(
    runoff,
    "runoff",
    file.path(GENERATED_RESULTS_DIR, "runoff_spectral_pairwise_Holm.csv")
  )

  summary <- rbind(precipitation_result$summary, runoff_result$summary)
  write.csv(summary,
            file.path(GENERATED_RESULTS_DIR, "conventional_comparison_summary.csv"),
            row.names = FALSE)
  summary
}
