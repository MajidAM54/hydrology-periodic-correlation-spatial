# Install the only non-base package used by the conventional comparison.
required <- c("PMCMRplus")
missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing) > 0L) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}
message("Required packages are available.")
