
REPRO_MODE <- tolower(Sys.getenv("REPRO_MODE", unset = "quick"))
REPRO_MODE <- match.arg(REPRO_MODE, c("quick", "real", "full"))

BASE_SEED <- 20260801L
ALPHA <- 0.05

T_VALUES <- c(12L, 3L)
N_VALUES <- c(100L, 150L, 200L, 500L)
K_VALUES <- c(10L, 15L, 20L, 30L, 50L)

if (REPRO_MODE == "full") {
  NSIM_POWER <- 1000L
  NSIM_TYPE1 <- 1000L
  NPERM <- 2000L
  N_CORES <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
} else if (REPRO_MODE == "real") {
  NSIM_POWER <- 0L
  NSIM_TYPE1 <- 0L
  NPERM <- 2000L
  N_CORES <- 1L
} else {
  NSIM_POWER <- 10L
  NSIM_TYPE1 <- 10L
  NPERM <- 99L
  N_CORES <- 1L
}


FINITE_SAMPLE_CORRECTION <- FALSE

REAL_T <- 12L
REAL_K <- 50L
REAL_NPERM <- if (REPRO_MODE == "quick") 199L else 2000L
CENTER_REAL_SERIES <- TRUE

# Small numerical stabilization for determinant calculations. This does not
# change the model or permutation scheme; it avoids failures from nearly
# singular matrices in finite samples.
RIDGE_RELATIVE <- 1e-10

DATA_DIR <- "data"
GENERATED_RESULTS_DIR <- file.path("results", "generated")
GENERATED_FIGURES_DIR <- file.path("figures", "generated")

dir.create(GENERATED_RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(GENERATED_FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
