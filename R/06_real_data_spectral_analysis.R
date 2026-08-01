# =============================================================================
# Real-data spectral analysis: global, pairwise, and zone comparisons
# =============================================================================

moving_average <- function(x, width = 24L) {
  as.numeric(stats::filter(x, rep(1 / width, width), sides = 1))
}

run_one_real_dataset <- function(data,
                                 dataset_name,
                                 T = REAL_T,
                                 k = REAL_K,
                                 nperm = REAL_NPERM,
                                 seed = BASE_SEED) {
  prepared <- make_DFT_list_from_dataframe(
    data,
    T = T,
    k = k,
    center = CENTER_REAL_SERIES
  )

  global <- permutation_test_weighted(
    prepared$station_DFT_list,
    m = prepared$m,
    nperm = nperm,
    seed = seed
  )

  pairwise <- pairwise_permutation_holm(
    prepared$station_DFT_list,
    m = prepared$m,
    nperm = nperm,
    seed = seed + 10000L
  )

  global_table <- data.frame(
    Dataset = dataset_name,
    n = prepared$n,
    T = T,
    k = k,
    m = prepared$m,
    Number_of_stations = ncol(data),
    Wilks_Lambda = global$Wilks_Lambda,
    p_value = global$p_value,
    nperm = nperm
  )

  write.csv(global_table,
            file.path(GENERATED_RESULTS_DIR, paste0(dataset_name, "_spectral_global.csv")),
            row.names = FALSE)
  write.csv(pairwise,
            file.path(GENERATED_RESULTS_DIR, paste0(dataset_name, "_spectral_pairwise_Holm.csv")),
            row.names = FALSE)

  list(prepared = prepared, global = global_table, pairwise = pairwise)
}

run_precipitation_zone_analysis <- function(precipitation_prepared,
                                            nperm = REAL_NPERM,
                                            seed = BASE_SEED) {
  zones <- list(
    West = c("Station 1", "Station 2", "Station 3"),
    Central = c("Station 4", "Station 5", "Station 6"),
    East = c("Station 7", "Station 8", "Station 9")
  )

  zone_DFT <- make_zone_average_DFT(precipitation_prepared$station_DFT_list, zones)
  global <- permutation_test_weighted(
    zone_DFT,
    m = precipitation_prepared$m,
    nperm = nperm,
    seed = seed
  )
  pairwise <- pairwise_permutation_holm(
    zone_DFT,
    m = precipitation_prepared$m,
    nperm = nperm,
    seed = seed + 10000L
  )

  global_table <- data.frame(
    Dataset = "precipitation_zones",
    n = precipitation_prepared$n,
    T = precipitation_prepared$T,
    k = precipitation_prepared$k,
    m = precipitation_prepared$m,
    Number_of_zones = length(zones),
    Wilks_Lambda = global$Wilks_Lambda,
    p_value = global$p_value,
    nperm = nperm
  )

  write.csv(global_table,
            file.path(GENERATED_RESULTS_DIR, "precipitation_zone_global.csv"),
            row.names = FALSE)
  write.csv(pairwise,
            file.path(GENERATED_RESULTS_DIR, "precipitation_zone_pairwise_Holm.csv"),
            row.names = FALSE)

  list(global = global_table, pairwise = pairwise)
}

write_real_data_descriptive_tables <- function(precipitation, runoff) {
  zones <- list(West = 1:3, Central = 4:6, East = 7:9)
  precip_rows <- lapply(names(zones), function(region) {
    regional_series <- rowMeans(precipitation[, zones[[region]], drop = FALSE])
    data.frame(
      Region = region,
      Mean_monthly_precipitation_mm = mean(regional_series),
      SD_of_regional_monthly_series = sd(regional_series)
    )
  })
  precip_zones <- do.call(rbind, precip_rows)

  upstream <- rowMeans(runoff[, 1:3, drop = FALSE])
  runoff_groups <- data.frame(
    Group = c("Stations 1-3", "Station 4"),
    Mean_monthly_runoff_m3_s = c(mean(upstream), mean(runoff[, 4])),
    SD_of_monthly_series = c(sd(upstream), sd(runoff[, 4]))
  )

  write.csv(precip_zones,
            file.path(GENERATED_RESULTS_DIR, "precipitation_regional_descriptive.csv"),
            row.names = FALSE)
  write.csv(runoff_groups,
            file.path(GENERATED_RESULTS_DIR, "runoff_group_descriptive.csv"),
            row.names = FALSE)

  list(precipitation = precip_zones, runoff = runoff_groups)
}

make_real_data_figures <- function(precipitation, runoff) {
  precip_zone_series <- data.frame(
    West = rowMeans(precipitation[, 1:3, drop = FALSE]),
    Central = rowMeans(precipitation[, 4:6, drop = FALSE]),
    East = rowMeans(precipitation[, 7:9, drop = FALSE])
  )
  precip_smoothed <- as.data.frame(lapply(precip_zone_series, moving_average, width = 24L))

  png(file.path(GENERATED_FIGURES_DIR, "precipitation_regions_24month_moving_average.png"),
      width = 1800, height = 1200, res = 200)
  matplot(seq_len(nrow(precip_smoothed)), precip_smoothed,
          type = "l", lty = 1:3, lwd = 2,
          xlab = "Time (Month)", ylab = "Monthly precipitation (mm)")
  legend("topleft", legend = names(precip_smoothed), lty = 1:3, lwd = 2, bty = "n")
  dev.off()

  runoff_group_series <- data.frame(
    Stations_1_3 = rowMeans(runoff[, 1:3, drop = FALSE]),
    Station_4 = runoff[, 4]
  )
  png(file.path(GENERATED_FIGURES_DIR, "runoff_upstream_vs_outlet.png"),
      width = 1800, height = 1200, res = 200)
  matplot(seq_len(nrow(runoff_group_series)), runoff_group_series,
          type = "l", lty = 1:2, lwd = 2,
          xlab = "Time (Month)", ylab = expression("Monthly runoff (m"^3*"/s)"))
  legend("topleft", legend = c("Average of Stations 1-3", "Station 4"),
         lty = 1:2, lwd = 2, bty = "n")
  dev.off()
}

run_real_data_spectral_analysis <- function() {
  precipitation <- read_station_csv(file.path(DATA_DIR, "precipitation_monthly_180.csv"))
  runoff <- read_station_csv(file.path(DATA_DIR, "runoff_monthly_180.csv"))

  precipitation_result <- run_one_real_dataset(
    precipitation,
    dataset_name = "precipitation",
    seed = BASE_SEED + 4000000L
  )
  runoff_result <- run_one_real_dataset(
    runoff,
    dataset_name = "runoff",
    seed = BASE_SEED + 4100000L
  )
  zone_result <- run_precipitation_zone_analysis(
    precipitation_result$prepared,
    seed = BASE_SEED + 4200000L
  )

  descriptive <- write_real_data_descriptive_tables(precipitation, runoff)
  make_real_data_figures(precipitation, runoff)

  list(
    precipitation = precipitation_result,
    runoff = runoff_result,
    precipitation_zones = zone_result,
    descriptive = descriptive
  )
}
