# tabulation.R
#'
#' @title Rearrange a phytosociological table, showing differential taxa on top.
#'
#' @description This function reorders a phytosociological table rows using,
#' firstly, the increasing number of groups in which a taxon occurs, and
#' secondly, the decreasing sum of the inner frequency of presences of each
#' taxon (see \code{\link{tdv}}).
#' The columns are also reordered, simply using the increasing number of the
#' respective group membership.
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
#' If `plot.im` = "normal", it returns the above list and, additionally, plots
#' an image of the tabulated matrix.
#' If `plot.im` = "condensed", it returns the above list and, additionally,
#' plots an image of the tabulated matrix, but presenting the sets of
#' differential taxa as solid coloured blocks of equal width.
#'
#' @author Tiago Monteiro-Henriques. E-mail: \email{tiagomonteirohenriques@@gmail.com}.
#'
#' @examples
#'
#' #getting the Taxus baccata forests data set
#' data(taxus_bin)
#' #creating a group partition, as presented in the original article of the data set
#' groups <- rep(c(1,2,3), c(3,11,19))
#'
#' #removing taxa occurring in only one relevé in order to
#' #reproduce exactly the example in the original article of the data set
#' taxus_bin_wmt <- taxus_bin[rowSums(taxus_bin) > 1,]
#'
#' #sorts the phytosociological table, putting exclusive taxa in the top and
#' #plots an image of it
#' tabul <- tabulation(taxus_bin_wmt, groups, taxa.names = rownames(taxus_bin_wmt), plot.im = "normal")
#'
#' #inspect the first rows and columns of the reordered phytosociological table
#' tabul$tabulated[1:6, 1:10]
#'
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
  order.1 <- apply(ot, 2, function(x) {as.numeric(paste0(x[x>0], collapse=""))}) #this will give
  #shortlex/radix order, i.e. it sorts taxa firstly by the no. of groups containing the taxon
  #and secondly by the lexicographical order of the groups the taxa belongs to
  order.2 <- - colSums(res$ifp) # - the sum of all the relative frequencies for the taxon (no
  #need to divide by res$e, as will be used only as order to ties of order.1
  taxa.ord <- order(order.1, order.2)
  sort.rel <- order(p) #to sort relevés by the respective group numbers

  mat1 <- rbind(sort(p), 0, m[taxa.ord, sort.rel])
  colnames(mat1) <- sort.rel
  ht <- colnames(m)[sort.rel]
  rownames(mat1) <- c("group","space",as.character(taxa.names)[taxa.ord])
  mat2 <- t(res$ofda * res$ifp)[taxa.ord,]
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
      mat2.im <- mat2 > 0
      mat2.im <- rbind((1:k)+1, 0, mat2.im)
      mat2.im[3:(ns+2),] <- mat2.im[3:(ns+2),]*matrix((1:k)+1,ns,k,byrow=TRUE)
      mat2.im[mat2.im==0] <- 1
      mat2.im[2,] <- 0
      graphics::image(t(mat2.im[(ns+2):1,]),col = c("black","white", grDevices::hcl.colors(k, palette)), xaxt="n", yaxt="n")
    }
  }
  return(list('taxa.names' = taxa.names, 'taxa.ord' = taxa.ord, header = ht, tabulated = mat1[-2,], condensed = mat2))
}
