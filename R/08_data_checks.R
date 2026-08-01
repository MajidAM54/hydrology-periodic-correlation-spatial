# =============================================================================
# Input-data integrity checks
# =============================================================================

check_one_dataset <- function(path, expected_rows, expected_series) {
  raw <- read.csv(path, check.names = FALSE)
  month_index_ok <- "Month_Index" %in% names(raw) &&
    identical(as.integer(raw$Month_Index), seq_len(expected_rows))
  data <- raw[, setdiff(names(raw), "Month_Index"), drop = FALSE]

  data.frame(
    File = path,
    Rows = nrow(data),
    Expected_rows = expected_rows,
    Series = ncol(data),
    Expected_series = expected_series,
    Missing_values = sum(is.na(data)),
    Month_index_sequential = month_index_ok,
    Row_check = nrow(data) == expected_rows,
    Series_check = ncol(data) == expected_series
  )
}

run_data_checks <- function() {
  checks <- rbind(
    check_one_dataset(file.path(DATA_DIR, "precipitation_monthly_180.csv"), 180L, 9L),
    check_one_dataset(file.path(DATA_DIR, "runoff_monthly_180.csv"), 180L, 4L)
  )
  write.csv(checks,
            file.path(GENERATED_RESULTS_DIR, "input_data_integrity_checks.csv"),
            row.names = FALSE)
  if (!all(checks$Row_check & checks$Series_check &
           checks$Missing_values == 0 & checks$Month_index_sequential)) {
    stop("One or more input-data checks failed. See input_data_integrity_checks.csv")
  }
  checks
}
