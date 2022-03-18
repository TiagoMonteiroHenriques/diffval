# select_taxa.R
#'
#' @title Select taxa by a 'constancy' range
#'
#' @description Trims a phytosociological table using the given 'constancy' range.
#'
#' @param m A `matrix`, i.e. a phytosociological table of 0s (absences) and 1s (presences), where rows correspond to taxa and columns correspond to relevés.
#' @param const A `vector` with two numeric values, which will be used as minimum (the first value) and maximum (the second value) 'constancy' values to trim table `m`.
#' @param min.pres 	A `numeric` integer giving, optionally, the minimum number presences of each taxon (in all relevés). If provided, it overrides the 'constancy' range.
#'
#' @details The function accepts a phytosociological table (`m`), and minimum and maximum 'constancy' values (`const`). It returns a trimmed phytosociological table, removing taxa outside the given 'constancy' range and removing relevés that become empty. Optionally, a minimum number of presences (in relevés) might be used (`min.pres`) instead of the 'constancy' range.
#'
#' @return A `list` with the following components:
#'
#'
#'\describe{
#'   \item{remo.rel}{The indices of the removed relevé(s).}
#'   \item{remo.sp}{The indices of the removed taxon(taxa).}
#'   \item{m.select}{The trimmed table.}
#'}
#'
#' @author Tiago Monteiro-Henriques. E-mail: \email{tiagomonteirohenriques@@gmail.com}.
#'
#' @export
#'
select_taxa <- function(m, const=c(0,100), min.pres=NULL) {
  stopifnot(is.matrix(m))
  mode(m) <- "integer"
  if (!identical(c(0L,1L), sort(unique(as.vector(m))))) {stop("Matrix m should contain only 0's and 1's")}
  if (min(rowSums(m))==0) {stop("At least one taxa is not present in any relev\u00e9")}
  if (min(colSums(m))==0) {stop("At least one relev\u00e9 contains no taxa")}
  mpa <- (m>0)+0
  m.col <- ncol(m)
  m.row <- nrow(m)
  if (is.null(min.pres)) { #selects based on constancy values
    inf <- const[1]/100
    sup <- const[2]/100
    sp.ind <- rowSums(mpa) > round(m.col*inf) & rowSums(mpa)<round(m.col*sup)
    m.select <- m[sp.ind,]
    rel.ind <- colSums(m.select) > 0
    m.select <- m.select[, rel.ind]
    cat(paste0(m.row-sum(sp.ind), " taxon(taxa) removed and ", m.col-sum(rel.ind), " empty relev\u00e9(s) removed\n"))
    return(list(remo.rel=which(!rel.ind), remo.sp=which(!sp.ind), m.select=m.select))
  } else { #selects using number of presences in relevés
    sp.ind <- rowSums(mpa) >= min.pres
    m.select <- m[sp.ind,]
    rel.ind <- colSums(m.select)>0
    m.select <- m.select[, rel.ind]
    cat(paste0(m.row-sum(sp.ind), " taxon(taxa) removed and ", m.col-sum(rel.ind), " empty relev\u00e9(s) removed\n"))
    return(list(remo.rel=which(!rel.ind), remo.sp=which(!sp.ind), m.select=m.select))
  }
}
