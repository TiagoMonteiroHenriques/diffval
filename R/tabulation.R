# tabulation.R
#'
#' @title Rearrange a phytosociological table, showing differential taxa on top.
#'
#' @description This function reorders a phytosociological table rows based on, firstly, the number of groups in which a taxon occurs in, and secondly,
#' on the within-group relative frequency. The columns are reordered using the given k-partition `p`.
#'
#' @param m A `matrix`, i.e. a phytosociological table of 0s (absences) and 1s (presences), where rows correspond to taxa and columns correspond to relevés.
#' @param p A `vector` of integer numbers with the partition of the relevés (i.e. a k-partition, consisting in a vector with values from 1 to k, with length equal to the number of columns of m, ascribing each relevé to one of the k groups).
#' @param taxa.names A `character vector` (with length equal to the number of rows of `m`) with the taxa names.
#' @param plot.im By default, `NULL`, returns without plotting. If `plot.im` = "normal", plots an image of the tabulated matrix. If `plot.im` = "condensed", plots an image of the tabulated matrix but presenting sets of differential taxa as solid coloured blocks.
#' @param palette A `character` with the name of the colour palette (one of \code{\link[grDevices]{hcl.pals}}`()`) to be passed to \code{\link[grDevices]{hcl.colors}}. Defaults to "Vik".
#'
#' @details The function accepts a phytosociological table (`m`), a k-partition of its columns (`p`) and the names of the taxa (corresponding to
#' the rows of `m`), returning a rearranged/reordered matrix (and plotting optionally).
#'
#' @return If `plot.im` = `NULL`, a `list` with the following components:
#'
#' \describe{
#'   \item{taxa.names}{The given `taxa.names`}
#'   \item{taxa.ord}{A `vector` with the order of the rows/taxa.}
#'   \item{tabulated}{The rearranged/reordered `m` `matrix`.}
#'   \item{condensed}{The matrix used to create the "condensed" image.}
#' }
#'
#' If `plot.im` = "normal", it returns the above list and, additionally, plots an image of the tabulated matrix.
#' If `plot.im` = "condensed", it returns the above list and, additionally, plots an image of the tabulated matrix, but presenting the sets of differential taxa as solid coloured blocks.
#'
#' @author Tiago Monteiro-Henriques. E-mail: \email{tiagomonteirohenriques@@gmail.com}.
#'
#' @export
#'
tabulation <- function (m, p, taxa.names, plot.im = NULL, palette = "Vik") {
  stopifnot(is.matrix(m))
  nr <- ncol(m) # no. of relevés
  if (!identical(length(p), nr)) {stop("Object p must be a partition of the columns of m")}
  k <- max(p)
  mode(p) <- "integer"
  if (!identical(sort(unique(p)), 1:k)) {stop("Object p is not a valid partition of the columns of m")}
  mode(m) <- "integer"
  if (!identical(c(0L,1L), sort(unique(as.vector(m))))) {stop("Matrix m should contain only 0's and 1's")}
  if (min(rowSums(m)) == 0) {stop("At least one taxa is not present in any relev\u00e9")}
  if (min(colSums(m)) == 0) {stop("At least one relev\u00e9 contains no taxa")}
  ns <- nrow(m) #no. of taxa
  tp <- table(p)
  if (length(taxa.names)!=ns) {stop("The length of taxa.names must match the number of rows of m")}
  if (!identical(length(p), nr)) {stop("Object p must be a partition of the columns of m")}
  res <- tdv(m, p, output.type = "full")
  ot <- (res$afg > 0) * (1:k)

  order.0 <- apply(ot, 2, function(x) {as.numeric(paste0(x[x>0], collapse=""))})
  #order.1 <- res$t2 #no. of groups containing the taxon
  #order.2 <- apply(ot, 2, function(x) {prod(x[x>0])}) #product of group numbers when the taxon is present
  order.3 <- 1 - (colSums(res$arf) / res$t2) #1 - the sum of all the adjusted relative frequencies for the taxon

  sort.rel <- order(p)
  #taxa.ord <- order(order.1, order.2, order.3)
  taxa.ord <- order(order.0, order.3)

  mat1 <- rbind(sort(p), 0, m[taxa.ord, sort.rel])
  colnames(mat1) <- sort.rel
  ht <- colnames(m)[sort.rel]
  rownames(mat1) <- c("group","space",as.character(taxa.names)[taxa.ord])
  mat2 <- t(res$adp*res$arf)[taxa.ord,]
  rownames(mat2) <- c(as.character(taxa.names)[taxa.ord])
  if (!is.null(plot.im)) {
    if (plot.im == "normal") {
      mat1.im <- mat1
      mat1.im[3:(ns+2),] <- mat1.im[3:(ns+2),] * matrix(sort(p)+1,ns,nr,byrow=TRUE)
      mat1.im[mat1.im==0] <- 1
      mat1.im[1,] <- mat1.im[1,]+1
      mat1.im[2,] <- 0
      graphics::image(t(mat1.im[(ns+2):1,]),col = c("black","white", grDevices::hcl.colors(k, palette)), xaxt="n", yaxt="n")
    }
    if (plot.im == "condensed") {
      mat2.im <- mat2>0
      mat2.im <- rbind((1:k)+1, 0, mat2.im)
      mat2.im[3:(ns+2),] <- mat2.im[3:(ns+2),]*matrix((1:k)+1,ns,k,byrow=TRUE)
      mat2.im[mat2.im==0] <- 1
      mat2.im[2,] <- 0
      graphics::image(t(mat2.im[(ns+2):1,]),col = c("black","white", grDevices::hcl.colors(k, palette)), xaxt="n", yaxt="n")
    }
  }
  return(list('taxa.names' = taxa.names, 'taxa.ord'=taxa.ord, header = ht, tabulated = mat1[-2,], condensed = mat2))
}
