# explore_tabulation.R
#'
#' @title Interactively explore a tabulation of a phytosociologic matrix.
#'
#' @description This function plots an interactive image of a tabulation.
#'
#' @param tab A `list` as returned by the \code{\link{tabulation}} function.
#' @param palette A `character` with the name of the colour palette (one of \code{\link[grDevices]{hcl.pals}}`()`) to be passed to \code{\link[grDevices]{hcl.colors}}. Defaults to "Vik".
#'
#' @details The function explore.tabulation accepts an object returned by the \code{\link{tabulation}} function, plotting a condensed image
#' of the respective tabulated matrix, permitting the user to click on the coloured blocks and receive the respective list of taxa names on
#' the console.
#'
#' @return Returns invisibly, although it prints taxa names on the console upon the user click on the figure.
#'
#' @author Tiago Monteiro-Henriques. E-mail: \email{tiagomonteirohenriques@@gmail.com}.
#'
#' @export
#'
explore_tabulation <- function (tab, palette = "Vik") {
  mat2 <- tab$condensed
  ns <- nrow(mat2)
  k <- ncol(mat2)
  taxa.names <- tab$taxa.names
  taxa.ord <- tab$taxa.ord
  mat2.im <- mat2 > 0
  mat2.im <- rbind((1:k) + 1, 0, mat2.im)
  mat2.im[3:(ns + 2),] <- mat2.im[3:(ns + 2),] * matrix((1:k) + 1, ns, k, byrow=TRUE)
  mat2.im[mat2.im == 0] <- 1
  mat2.im[2,] <- 0
  graphics::image(t(mat2.im[(ns + 2):1,]), col = c("black","white", grDevices::hcl.colors(k, palette)), xaxt="n", yaxt="n")
  id.x <- rep(seq(0, 1, length.out = k))
  id.y <- rev(rep(seq(0, 1, length.out = ns + 2)[1:ns], times = k))
  list.cent.label <- apply((mat2 > 0) + 0, 2, function (x) {
    #x <- (mat2>0)+0
    #x <- x[,2]
    ic <- 1
    i.from <- NULL
    i.to <- NULL
    sen <- "from"
    x <- c(x,0)
    while (ic <= length(x)) { #maybe replace with strgsplit + grep!
      if (sen == "from") {
        if (x[ic] == 0) {
          ic <- ic + 1
        } else {
          i.from <- c(i.from, ic)
          sen <- "to"
          ic <- ic + 1
        }
      } else {
        if (x[ic] == 1) {
          ic <- ic + 1
        } else {
          i.to <- c(i.to, ic - 1)
          sen <- "from"
          ic <- ic + 1
        }
      }
    }
    if (is.null(cbind(i.from, i.to))) {
      return()
    }
    ind <- apply(cbind(i.from, i.to), 1, function (z) {
      seq(z[1],z[2])
    })
    if (is.matrix(ind)) {
      apply(ind, 2, function (y) {
        res.y <- mean(id.y[y]) #find mean coordinate for plotting
        res.label <- paste(paste(taxa.names[taxa.ord][y], collapse="\n"), "\n\n")
        return(list(res.y = res.y, res.label = res.label))
      })
    } else {
      lapply(ind, function (y) {
        res.y <- mean(id.y[y]) #find mean coordinate for plotting
        res.label <- paste(paste(taxa.names[taxa.ord][y], collapse="\n"), "\n\n")
        return(list(res.y = res.y, res.label = res.label))
      })
    }
  })
  id.x <- id.x[!sapply(list.cent.label, is.null)]
  list.cent.label <- Filter(Negate(is.null), list.cent.label)
  id.y.cent.labt <- lapply(list.cent.label, function (x) {
    as.matrix(simplify2array(x, higher = FALSE)[1,])
  })
  id.y.cent.lab <- unlist(id.y.cent.labt)
  id.lab <- c("",unlist(sapply(list.cent.label, function (x) {
    as.matrix(simplify2array(x,higher = FALSE)[2,])
  })))
  id.x.cent.lab <- rep.int(id.x, times = sapply(id.y.cent.labt, length))
  #; print(id.y.cent.label); print(id.label)
  graphics::points(id.x.cent.lab , id.y.cent.lab, cex = 0.5, pch = 20)
  y.exit <- 0.95 #it can be improved
  graphics::points(1, y.exit, cex = 3)
  graphics::points(1, y.exit, cex = 2, col = "red")
  graphics::points(1, y.exit, cex = 1, col = "white")
  cat("Click on the plot black dots to retrieve taxon(taxa) names.\nClick on the top-right circle to exit.\n\n")
  res.click <- "start"
  while (res.click != "") {
    res.click <- id.lab[(graphics::identify(x = c(1, id.x.cent.lab), y = c(y.exit, id.y.cent.lab), labels = id.lab, cex = 0.5, plot = FALSE, n = 1))]
    cat(res.click)
  }
  invisible()
}
