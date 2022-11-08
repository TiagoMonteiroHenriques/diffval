#' *Taxus baccata* forests
#'
#' A binary phytosociological table containing relevés of *Taxus baccata*
#'   forests, in the northwest of the Iberian Peninsula.
#'
#' @format A matrix with 209 rows and 33 columns. Each column correspond to a
#'   phytosociological relevé and each row correspond to a taxon. Values in the
#'   matrix denote presences (1) and absences (0).
#'
#' @source Portela-Pereira E., Monteiro-Henriques T., Casas C., Forner N.,
#'   Garcia-Cabral I., Fonseca J.P. & Neto C. 2021. *Teixedos no noroeste da*
#'   *Península Ibérica*. Finisterra 56(117): 127‐150.
#'   \doi{10.18055/FINIS18102}.
#'
#' @examples
#' # Getting the Taxus baccata forests data set
#' data(taxus_bin)
#'
#' # Inspect the first rows and columns of taxus_bin
#' head(taxus_bin, n = c(5, 5))
#'
"taxus_bin"
