#' Rearrange a phytosociological table, showing differential taxa on top
#'
#' This function reorders a phytosociological table rows using, firstly, the
#'   increasing number of groups in which a taxon occurs, and secondly, the
#'   decreasing sum of the inner frequency of presences of each taxon
#'   (see [tdv()]). The columns are also reordered, simply using the increasing
#'   number of the respective group membership.
#'
#' @param m_bin A matrix. A phytosociological table of 0s (absences) and 1s
#'   (presences), where rows correspond to taxa and columns correspond to
#'   relevés.
#' @param p A vector of integer numbers with the partition of the relevés (i.e.,
#'   a k-partition, consisting in a vector with values from 1 to k, with length
#'   equal to the number of columns of `m_bin`, ascribing each relevé to one of
#'   the k groups).
#' @param taxa_names A character vector (with length equal to the number of rows
#'   of `m_bin`) with the taxa names.
#' @param plot_im By default, `NULL`, returns without plotting. If
#'   `plot_im = "normal"`, plots an image of the tabulated matrix. If
#'   `plot_im = "condensed"`, plots an image of the tabulated matrix but
#'   presenting sets of differential taxa as solid coloured blocks.
#' @param palette A character with the name of the colour palette (one of
#'   [grDevices::hcl.pals()] to be passed to [grDevices::hcl.colors()]. Defaults
#'   to "Vik".
#' @param greyout A logical. If `TRUE` (the default), non-differential taxa are
#'   greyed out (using the colour defined by `greyout_colour`). If `FALSE`,
#'   non-differential taxa is depicted with the respective group colours.
#' @param greyout_colour A character with the name of the colour to use for
#'   non-differential taxa. Defaults to "grey".
#'
#' @details The function accepts a phytosociological table (`m_bin`), a
#'   k-partition of its columns (`p`) and the names of the taxa (corresponding
#'   to the rows of `m_bin`), returning a rearranged/reordered matrix (and
#'   plotting optionally).
#'
#' @return If `plot_im = NULL`, a list with the following components:
#'   \describe{
#'     \item{taxa.names}{The given `taxa_names`}
#'     \item{taxa.ord}{A vector with the order of the rows/taxa.}
#'     \item{tabulated}{The rearranged/reordered `m_bin` matrix.}
#'     \item{condensed}{The matrix used to create the "condensed" image.}
#'   }
#'
#'   If `plot_im = "normal"`, it returns the above list and, additionally, plots
#'   an image of the tabulated matrix.
#'   If `plot_im = "condensed"`, it returns the above list and, additionally,
#'   plots an image of the tabulated matrix, but presenting the sets of
#'   differential taxa as solid coloured blocks of equal width.
#'
#' @author Tiago Monteiro-Henriques. E-mail: \email{tmh.dev@@icloud.com}.
#'
#' @examples
#' # Getting the Taxus baccata forests data set
#' data(taxus_bin)
#'
#' # Creating a group partition, as presented in the original article of the
#' # data set
#' groups <- rep(c(1, 2, 3), c(3, 11, 19))
#'
#' # Removing taxa occurring in only one relevé in order to
#' # reproduce exactly the example in the original article of the data set
#' taxus_bin_wmt <- taxus_bin[rowSums(taxus_bin) > 1, ]
#'
#' # Sorting the phytosociological table, putting exclusive taxa in the top and
#' # plotting an image of it
#' tabul <- tabulation(
#'   m_bin = taxus_bin_wmt,
#'   p = groups,
#'   taxa_names = rownames(taxus_bin_wmt),
#'   plot_im = "normal",
#'   palette = "Zissou 1"
#' )
#'
#' # Inspect the first rows and columns of the reordered phytosociological table
#' head(tabul$tabulated, n = c(5, 5))
#'
#' @export
tabulation <- function(m_bin,
                       p,
                       taxa_names,
                       plot_im = NULL,
                       palette = "Vik",
                       greyout = TRUE,
                       greyout_colour = "grey") {
  stopifnot(is.matrix(m_bin))
  nr <- ncol(m_bin) # No. of relevés
  if (!identical(length(p), nr)) {
    stop("Object `p` must be a partition of the columns of `m_bin`.")
  }
  k <- max(p)
  mode(p) <- "integer"
  if (!identical(sort(unique(p)), 1:k)) {
    stop("Object `p` is not a valid partition of the columns of `m_bin`.")
  }
  mode(m_bin) <- "integer"
  if (!identical(c(0L, 1L), sort(unique(as.vector(m_bin))))) {
    stop("Matrix `m_bin` should contain only 0's and 1's.")
  }
  if (min(rowSums(m_bin)) == 0) {
    stop("At least one taxa is not present in any relev\u00e9.")
  }
  if (min(colSums(m_bin)) == 0) {
    stop("At least one relev\u00e9 contains no taxa.")
  }
  ns <- nrow(m_bin) # No. of taxa
  if (length(taxa_names) != ns) {
    stop("The length of `taxa_names` must match the number of rows of `m_bin`.")
  }
  if (!identical(length(p), nr)) {
    stop("Object `p` must be a partition of the columns of `m_bin`.")
  }
  res <- tdv(m_bin, p, output_type = "full")
  ot <- (res$afg > 0) * (1:k)
  order_0 <- res$e
  order_1_mat <- t(apply(ot, 2, function(x) {
    groups <- x[x > 0]
    c(groups, rep(Inf, k - length(groups))) # fill to equal length
  })) # data frame to be sorted column by column after order_0, giving the
  # shortlex order of the rows
  order_2 <- -colSums(res$ifp) # Minus the sum of all the relative frequencies
  # for the taxon (no need to divide by res$e, as will be used only as ordering
  # factor for the ties within order_0 plus order_1_mat
  taxa_ord <- do.call(
    order, c(list(order_0), as.data.frame(order_1_mat), list(order_2))
  ) # order_0 plus order_1_mat
  # give the shortlex order, i.e. order_0 first sorts taxa by the number of
  # groups in which they occur, and order_1 then sorts ties lexicographically by
  # the groups the taxa belong to
  sort_rel <- order(p) # To sort relevés by the respective group numbers

  mat1 <- rbind(sort(p), 0, m_bin[taxa_ord, sort_rel])
  colnames(mat1) <- sort_rel
  ht <- colnames(m_bin)[sort_rel]
  rownames(mat1) <- c("group", "space", as.character(taxa_names)[taxa_ord])
  mat2 <- t(res$ofda * res$ifp)[taxa_ord, , drop = FALSE]
  rownames(mat2) <- c(as.character(taxa_names)[taxa_ord])
  if (!is.null(plot_im)) {
    if (plot_im == "normal") {
      mat1_im <- mat1
      mat1_im[3:(ns + 2), ] <- mat1_im[3:(ns + 2), ] *
        matrix(sort(p) + 1, ns, nr, byrow = TRUE)
      mat1_im[mat1_im == 0] <- 1
      mat1_im[1, ] <- mat1_im[1, ] + 1
      mat1_im[2, ] <- 0
      if (greyout) {
        if (all(!rowSums(mat2) == 0)) {
          group_colour <- c(
            "black",
            "white",
            grDevices::hcl.colors(max(k, 2), palette)[1:k]
          )
        } else {
          group_colour <- c(
            "black",
            "white",
            grDevices::hcl.colors(max(k, 2), palette)[1:k],
            greyout_colour
          )
          mat1_im[-c(1:2), ][rowSums(mat2) == 0, ][mat1_im[-c(1:2), ][rowSums(mat2) == 0, ] > 1] <- k + 2
        }
        graphics::image(t(mat1_im[(ns + 2):1, ]),
          col = group_colour,
          xaxt = "n",
          yaxt = "n"
        )
      } else {
        graphics::image(t(mat1_im[(ns + 2):1, ]),
          col = c(
            "black",
            "white",
            grDevices::hcl.colors(max(k, 2), palette)[1:k]
          ),
          xaxt = "n",
          yaxt = "n"
        )
      }
    }
    if (plot_im == "condensed") {
      mat2_im <- mat2[, sort(p)] > 0
      mat2_im <- rbind(sort(p) + 1, 0, mat2_im)
      mat2_im[3:(ns + 2), ] <- mat2_im[3:(ns + 2), ] *
        matrix(sort(p) + 1, ns, nr, byrow = TRUE)

      mat2_im[mat2_im == 0] <- 1
      mat2_im[2, ] <- 0
      if (greyout) {
        if (all(!rowSums(mat2) == 0)) {
          group_colour <- c(
            "black",
            "white",
            grDevices::hcl.colors(max(k, 2), palette)[1:k]
          )
        } else {
          group_colour <- c(
            "black",
            "white",
            grDevices::hcl.colors(max(k, 2), palette)[1:k],
            greyout_colour
          )
          mat2_im[-c(1:2), ][rowSums(mat2) == 0, ] <- k + 2
        }
        graphics::image(t(mat2_im[(ns + 2):1, ]),
          col = group_colour,
          xaxt = "n",
          yaxt = "n"
        )
      } else {
        graphics::image(t(mat2_im[(ns + 2):1, ]),
          col = c(
            "black",
            "white",
            grDevices::hcl.colors(max(k, 2), palette)[1:k]
          ),
          xaxt = "n",
          yaxt = "n"
        )
      }
    }
  }
  list(
    taxa.names = taxa_names,
    taxa.ord   = taxa_ord,
    header     = ht,
    tabulated  = mat1[-2, ],
    condensed  = mat2
  )
}
