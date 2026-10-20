# Package initialization hooks
# This file runs when the package is loaded

.onLoad <- function(libname, pkgname) {
  # --- STRICTLY pinned (do NOT auto-install) ---
  pinned_pkgs <- list(
    xgboost = "1.7.11.1",
    SHAPforxgboost = "0.1.3"
  )

  # Enforce versions (fail fast, no auto-install)
  for (pkg in names(pinned_pkgs)) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      packageStartupMessage(sprintf(
        "ERROR: Required package '%s' not installed.\nPlease install version %s",
        pkg, pinned_pkgs[[pkg]]
      ))
      warning(sprintf(
        "Package '%s' version %s is required but not installed.",
        pkg, pinned_pkgs[[pkg]]
      ))
    }

    installed <- as.character(utils::packageVersion(pkg))
    if (!identical(installed, pinned_pkgs[[pkg]])) {
      warning(sprintf(
        "Package version mismatch for '%s': required %s but found %s.\nPlease install the exact version.",
        pkg, pinned_pkgs[[pkg]], installed
      ))
    }
  }

  packageStartupMessage("Version checks passed: xgboost 1.7.11.1, SHAPforxgboost 0.1.3")
}
