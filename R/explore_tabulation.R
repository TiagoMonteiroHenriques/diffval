#' Interactively explore a tabulation of a phytosociological matrix
#'
#' This function plots an interactive image of a tabulation.
#'
#' @param tab A list as returned by the [tabulation()] function.
#' @param palette A character with the name of the colour palette (one of
#'   [grDevices::hcl.pals()] to be passed to [grDevices::hcl.colors()]. Defaults
#'   to "Vik".
#'
#' @details The function explore.tabulation accepts an object returned by the
#'   [tabulation()] function, plotting a condensed image of the
#'   respective tabulated matrix, permitting the user to click on the coloured
#'   blocks and receive the respective list of taxa names on the console.
#'
#' @return Returns invisibly, although it prints taxa names on the console upon
#'   the user click on the figure.
#'
#' @author Tiago Monteiro-Henriques. E-mail: \email{tmh.dev@@icloud.com}.
#'
#' @examples
#' # Getting the Taxus baccata forests data set
#' data(taxus_bin)
#' # Creating a group partition, as presented in the original article of
#' # the data set
#' groups <- rep(c(1, 2, 3), c(3, 11, 19))
#'
#' # Removing taxa occurring in only one relevé in order to
#' # reproduce exactly the example in the original article of the data set
#' taxus_bin_wmt <- taxus_bin[rowSums(taxus_bin) > 1, ]
#'
#' # Sorts the phytosociological table, putting exclusive taxa at the top and
#' # plots an image of it
#' tabul <- tabulation(
#'   m_bin = taxus_bin_wmt,
#'   p = groups,
#'   taxa_names = rownames(taxus_bin_wmt),
#'   plot_im = "normal",
#'   palette = "Zissou 1"
#' )
#'
#' # This creates an interactive plot (where you can click)
#' if (interactive()) {
#'   explore_tabulation(tabul, palette = "Zissou 1")
#' }
#'
#' @export
explore_tabulation <- function(tab, palette = "Vik") {
  mat2 <- tab$condensed
  ns <- nrow(mat2)
  k <- ncol(mat2)
  taxa_names <- tab$taxa.names
  taxa_ord <- tab$taxa.ord
  mat2_im <- mat2 > 0
  mat2_im <- rbind((1:k) + 1, 0, mat2_im)
  mat2_im[3:(ns + 2), ] <- mat2_im[3:(ns + 2), ] * matrix(
    (1:k) + 1,
    ns,
    k,
    byrow = TRUE
  )
  mat2_im[mat2_im == 0] <- 1
  mat2_im[2, ] <- 0
  graphics::image(
    t(mat2_im[(ns + 2):1, ]),
    col = c("black", "white", grDevices::hcl.colors(k, palette)),
    xaxt = "n",
    yaxt = "n"
  )
  id_x <- rep(seq(0, 1, length.out = k))
  id_y <- rev(rep(seq(0, 1, length.out = ns + 2)[1:ns], times = k))
  list_cent_label <- apply((mat2 > 0) + 0, 2, function(x) {
    ic <- 1
    i_from <- NULL
    i_to <- NULL
    sen <- "from"
    x <- c(x, 0)
    while (ic <= length(x)) { # Maybe replace with strgsplit and grep
      if (sen == "from") {
        if (x[ic] == 0) {
          ic <- ic + 1
        } else {
          i_from <- c(i_from, ic)
          sen <- "to"
          ic <- ic + 1
        }
      } else {
        if (x[ic] == 1) {
          ic <- ic + 1
        } else {
          i_to <- c(i_to, ic - 1)
          sen <- "from"
          ic <- ic + 1
        }
      }
    }
    if (is.null(cbind(i_from, i_to))) {
      return()
    }
    ind <- apply(cbind(i_from, i_to), 1, function(z) {
      seq(z[1], z[2])
    })
    if (is.matrix(ind)) {
      apply(ind, 2, function(y) {
        res_y <- mean(id_y[y]) # Find mean coordinate for plotting
        res_label <- paste(
          paste(taxa_names[taxa_ord][y], collapse = "\n"),
          "\n\n"
        )
        return(list(res.y = res_y, res.label = res_label))
      })
    } else {
      lapply(ind, function(y) {
        res_y <- mean(id_y[y]) # Find mean coordinate for plotting
        res_label <- paste(
          paste(taxa_names[taxa_ord][y], collapse = "\n"),
          "\n\n"
        )
        return(list(res.y = res_y, res.label = res_label))
      })
    }
  })
  id_x <- id_x[!sapply(list_cent_label, is.null)]
  list_cent_label <- Filter(Negate(is.null), list_cent_label)
  id_y_cent_labt <- lapply(list_cent_label, function(x) {
    as.matrix(simplify2array(x, higher = FALSE)[1, ])
  })
  id_y_cent_lab <- unlist(id_y_cent_labt)
  id_lab <- c("", unlist(sapply(list_cent_label, function(x) {
    as.matrix(simplify2array(x, higher = FALSE)[2, ])
  })))
  id_x_cent_lab <- rep.int(id_x, times = sapply(id_y_cent_labt, length))
  graphics::points(id_x_cent_lab, id_y_cent_lab, cex = 0.5, pch = 20)
  y_exit <- 0.95 # it can be improved
  graphics::points(1, y_exit, cex = 3)
  graphics::points(1, y_exit, cex = 2, col = "red")
  graphics::points(1, y_exit, cex = 1, col = "white")
  cat("Click on the plot black dots to retrieve taxon(taxa) name(s).\nClick on
      the top-right circle to exit.\n\n")
  res_click <- "start"
  while (res_click != "") {
    res_click <- id_lab[(graphics::identify(
      x = c(1, id_x_cent_lab),
      y = c(y_exit, id_y_cent_lab),
      labels = id_lab, cex = 0.5, plot = FALSE, n = 1
    ))]
    cat(res_click)
  }
  invisible()
}
