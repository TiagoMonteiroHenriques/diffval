#' The Total Differential Value of a big phytosociological data set
#'
#' Given a big phytosociological data set represented as a list, and a partition
#'   of the relevés in that list, this function calculates the respective Total
#'   Differential Value (TDV).
#'
#' @param phyto_list A list. This is a very light representation of what could
#'   be a usual phytosociological table, registering only taxa presences. Each
#'   component should uniquely represent a taxon and should contain a vector (of
#'   numeric values) with the relevé(s) id(s) where that taxon was observed.
#'   Relevé's ids are expected to be represented by consecutive integers,
#'   starting with 1. The components of the list might be named (e.g. using the
#'   taxon name) or empty (decreasing further memory burden). However, for
#'   `output_type == "normal"` taxa names are useful for output interpretation.
#' @param p A vector of integer numbers with the partition of the relevés (i.e.,
#'   a k-partition, consisting in a vector with values from 1 to k, with length
#'   equal to the number of relevés in `phyto_list`, ascribing each relevé to
#'   one of the k groups).
#' @param n_rel The number of relevés in `phyto_list`, obtained, for
#'   example, using the instruction `length(unique(unlist(phyto_list)))`.
#' @param output_type A character determining the amount of information returned
#'   by the function and also the amount of pre-validations. Possible values are
#'   "normal" (the default) and "fast".
#' @param parallel Logical. Should function [parallel::mclapply()]) be used to
#'   improve computation time by forking? Not available on Windows. Refer to
#'   that function manual for more information. Defaults to `FALSE`.
#' @param mc_cores The number of cores to be passed to [parallel::mclapply()] if
#'   `parallel = TRUE`. See [parallel::mclapply()] for more information.
#'
#' @details This function accepts a list (`phyto_list`) representing a
#'   phytosociological data set, as well as a k-partition of its relevés (`p`),
#'   returning the corresponding TDV (see [diffval::tdv()] for an explanation
#'   on TDV).
#'   Partition `p` gives the group to which each relevé is ascribed, by
#'   increasing order of relevé id.
#'   Big phytosociological tables can occupy a significant amount of computer
#'   memory, which mostly relate to the fact that the absences (usually more
#'   frequent than presences) are also recorded in memory. The use of a list,
#'   focusing only on presences, reduces significantly the amount of needed
#'   memory to store all the information that a phytosociological table contains
#'   and also the computation time of TDV, allowing computations for big data
#'   sets.
#'
#' @return If `output_type = "normal"` (the default) pre-validations are done
#'   (which can take some time) and a list is returned, with the following
#'   components (see [diffval::tdv()] for the mathematical notation):
#'
#'   \describe{
#'     \item{ifp}{A matrix with the \eqn{\frac{a}{b}}{a/b} values for each taxon
#'       in each group, for short called the 'inner frequency of presences'.}
#'     \item{ofda}{A matrix with the \eqn{\frac{c}{d}}{c/d} values for each
#'       taxon in each group, for short called the 'outer frequency of
#'       differentiating absences'.}
#'     \item{e}{A vector with the \eqn{e} values for each taxon, i.e., the
#'       number of groups containing that taxon.}
#'     \item{diffval}{A matrix with the \eqn{DiffVal} for each taxon.}
#'     \item{tdv}{A numeric with the TDV of matrix `m_bin,` given the partition
#'       `p`.}
#'   }
#'
#'   If `output_type = "fast"`, only TDV is returned and no pre-validations are
#'   done.
#'
#' @author Tiago Monteiro-Henriques. E-mail: \email{tmh.dev@@icloud.com}.
#'
#' @examples
#' # Getting the Taxus baccata forests data set
#' data(taxus_bin)
#'
#' # Creating a group partition, as the one presented in the original article of
#' # the data set
#' groups <- rep(c(1, 2, 3), c(3, 11, 19))
#'
#' # Removing taxa occurring in only one relevé, in order to reproduce exactly
#' # the example in the original article of the data set
#' taxus_bin_wmt <- taxus_bin[rowSums(taxus_bin) > 1, ]
#'
#' # Calculating TDV using tdv()
#' tdv(taxus_bin_wmt, groups)$tdv
#'
#' # Converting from the phytosociologic matrix format to the list format
#' taxus_phyto_list <- apply(taxus_bin_wmt, 1, function(x) which(as.logical(x)))
#'
#' # Getting the number of relevés in the list
#' n_rel <- length(unique(unlist(taxus_phyto_list)))
#'
#' # Calculating TDV using bigdata_tdv(), even if this is not a big matrix
#' bigdata_tdv(
#'   phyto_list = taxus_phyto_list,
#'   p = groups,
#'   n_rel = n_rel,
#'   output_type = "normal"
#' )$tdv
#'
#' @export
bigdata_tdv <- function(phyto_list,
                        p,
                        n_rel,
                        output_type = "normal",
                        parallel = FALSE,
                        mc_cores = getOption("mc.cores", 2L)) {
  if (output_type == "fast") {
    k <- max(p)
  } else {
    stopifnot("`phyto_list` must be a list." = is.list(phyto_list))
    rel_ids <- sort(unique(unlist(phyto_list)))
    if (!identical(as.integer(n_rel), length(rel_ids))) {
      stop("`nrel` is not matching the number of relev\u00e9s in `phyto_list`.")
    }
    mode(rel_ids) <- "integer"
    if (!identical(rel_ids, seq_along(rel_ids))) {
      stop("Relev\u00e9's ids should be consecutive integers starting with 1.")
    }
    if (!identical(length(p), as.integer(n_rel))) {
      stop("`p` must be a partition of the relev\u00e9s in `phyto_list`.")
    }
    k <- max(p)
    mode(p) <- "integer"
    if (!identical(sort(unique(p)), 1:k)) {
      stop("`p` is not a valid partition of the relev\u00e9s in `phyto_list`.")
    }
  }
  n_taxa <- length(phyto_list)
  k <- max(p)
  b <- tabulate(p)
  d <- n_rel - b

  if (output_type == "fast") {
    if (parallel) {
      result <- parallel::mclapply(
        phyto_list,
        dvilf,
        p = p,
        k = k,
        n_taxa = n_taxa,
        b = b,
        d = d,
        mc.cores = mc_cores
      )
      result <- unlist(result)
    } else {
      result <- sapply(
        phyto_list,
        dvilf,
        p = p,
        k = k,
        n_taxa = n_taxa,
        b = b,
        d = d
      )
    }
    return(tdv = sum(result) / n_taxa)
  }

  if (output_type == "normal") {
    if (parallel) {
      result <- parallel::mclapply(
        phyto_list,
        dv_in_list,
        p = p,
        k = k,
        n_taxa = n_taxa,
        b = b,
        d = d,
        mc.cores = mc_cores
      )
    } else {
      result <- lapply(
        phyto_list,
        dv_in_list,
        p = p,
        k = k,
        n_taxa = n_taxa,
        b = b,
        d = d
      )
    }
    # Creating an output identical to tdv()
    diffval <- sapply(result, function(x) {
      x$diffval
    })
    ifp <- t(sapply(result, function(x) {
      x$ifp
    }))
    colnames(ifp) <- 1:k
    ofda <- t(sapply(result, function(x) {
      x$ofda
    }))
    colnames(ofda) <- 1:k
    tdv <- sum(diffval) / n_taxa
    diffval <- as.matrix(diffval)
    colnames(diffval) <- "DiffVal"
    return(list(
      ifp = ifp,
      ofda = ofda,
      e = sapply(result, function(x) {
        x$e
      }), # This is returning integer, while tdv() is returning double
      diffval = diffval,
      tdv = tdv
    ))
  }
  stop('Argument `output_type` must be "fast" or "normal".')
}
