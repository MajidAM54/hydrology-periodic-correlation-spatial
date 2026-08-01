# =============================================================================
# Master reproducibility script
# Usage:
#   Rscript run_all.R quick   # smoke test
#   Rscript run_all.R real    # complete real-data and conventional analyses
#   Rscript run_all.R full    # complete manuscript simulation design + real data
# =============================================================================

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) >= 1L) Sys.setenv(REPRO_MODE = tolower(arguments[1L]))

source(file.path("R", "00_config.R"))
source(file.path("R", "01_core_spectral_functions.R"))
source(file.path("R", "02_simulation_scenarios_1_3.R"))
source(file.path("R", "03_simulation_scenario_4.R"))
source(file.path("R", "04_type_I_error.R"))
source(file.path("R", "05_bandwidth_sensitivity.R"))
source(file.path("R", "06_real_data_spectral_analysis.R"))
source(file.path("R", "07_conventional_friedman_conover.R"))
source(file.path("R", "08_data_checks.R"))

message("Reproducibility mode: ", REPRO_MODE)
message("Checking input data...")
data_checks <- run_data_checks()

if (REPRO_MODE != "real") {
  message("1/5: Running Scenarios 1-3...")
  scenario_1_3_results <- run_scenarios_1_3()

  message("2/5: Running Scenario 4 and t3 sensitivity...")
  scenario_4_results <- run_scenario_4()

  message("3/5: Running empirical Type I error study...")
  type_I_results <- run_type_I_error()

  message("4/5: Running bandwidth sensitivity study...")
  bandwidth_results <- run_bandwidth_sensitivity()
}

message("5/5: Running real-data spectral and conventional analyses...")
real_data_results <- run_real_data_spectral_analysis()
conventional_results <- run_conventional_comparisons()
validation_results <- run_validation()

message("Completed. Outputs are in results/generated and figures/generated.")
