# diffval (development version)

# diffval 1.2.0

* Added function `internal_assignment()`
* Improved documentation, linking to Monteiro-Henriques (2025,
  <https://doi.org/10.3897/VCS.140466>).

# diffval 1.1.0

## Major changes

* Function `bigdata_tdv()` added, allowing Total Differential Value calculation
  for big matrices, optionally using fork with package _parallel_.
* Two auxiliary functions added to assist `bigdata_tdv()` in `utils.R`.

## Bug fixes

* Improved sorting in `tabulation()` function.
* TDV calculation was corrected for some cases in `optim_tdv_hill_climb()` and
  in the associated auxiliary function in `utils.R`.

# diffval 1.0.0

* Updated version, for package release.

# diffval 0.0.0.9025

* Added a `NEWS.md` file to track changes to the package.
