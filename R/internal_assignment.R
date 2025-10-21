#' Check the internal assignment of a given classification
#'
#' Given a phytosociological table and a partition of its columns, this function
#' checks the internal assignment of relevés to groups, based on the  presence
#' of taxa that are exclusive to each group (or group combination) as defined by
#' the partition, distinguishing relevés assigned unambiguously to their group
#' from those for which there is ambiguity.
#'
#' @param m_bin A matrix. A phytosociological table of 0s (absences) and 1s
#'   (presences), where rows correspond to taxa and columns correspond to
#'   relevés.
#' @param p A vector of integer numbers with the partition of the relevés (i.e.,
#'   a k-partition, consisting in a vector with values from 1 to k, with length
#'   equal to the number of columns of `m_bin`, ascribing each relevé to one of
#'   the k groups).
#'
#' @details The function accepts a phytosociological table (`m_bin`) and a
#'   k-partition of its columns (`p`), and assesses which relevés are assigned
#'   unambiguously to their group and which are not.
#'   The assignment of a relevé to a group is considered unambiguous when
#'   transferring it to another group would alter the pattern of differential
#'   taxa defined by `p`. Conversely, if a relevé could be moved to a different
#'   group without changing this pattern, its assignment is considered
#'   ambiguous.
#'
#' @return A list with the following components:
#'
#'   \describe{
#'     \item{rel_ambiguous_assign}{A vector containing the names of the relevés
#'     with ambiguous assignment.}
#'     \item{possible_assignments}{A data frame with all the possible
#'     assignments for the ambiguously assigned relevés.}
#'     \item{iap}{The internal assignment precision (IAP), i.e., the proportion
#'     of relevés with unambiguous assignment.}
#'     \item{iaa}{The internal assignment ambiguity (IAA), i.e., the proportion
#'     of relevés with ambiguous assignment (IAA = 1 - IAP).}
#'   }
#'
#' @author Tiago Monteiro-Henriques. E-mail: \email{tmh.dev@@icloud.com}.
#'
#' @examples
#' # Getting the Taxus baccata forests data set
#' data(taxus_bin)
#'
#' # Creating some group partitions
#' groups1 <- rep(c(1, 2, 3), c(3, 11, 19))
#' set.seed(1)
#' groups2 <- sample(rep(c(1, 2, 3), c(3, 11, 19)))
#'
#' # In this case, all relevés are unambiguously assigned to a group
#' internal_assignment(taxus_bin, groups1)
#'
#' # In this other case, some relevés could be moved to a different group, as
#' # their assignment is ambiguous
#' internal_assignment(taxus_bin, groups2)
#'
#' @export
internal_assignment <- function(m_bin, p) {
  stopifnot(is.matrix(m_bin))
  if (!identical(length(p), ncol(m_bin))) {
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
  afg <- rowsum(t(m_bin), group = p)
  nrel <- ncol(m_bin)
  k <- max(p)
  if (k == 1) {
    warning("Only one group is present!")
    return(
      list(
        releves_with_ambiguous_assignment = NA,
        possible_assignments_of_ambiguous_releves = as.data.frame(NA),
        IAA = NA,
        IAP = NA
      )
    )
  }
  key <- t(afg) > 0 # key of species presences in groups
  # possible groups to which the relevé can be assigned given the species it
  # contains
  poss_g <- sapply(1:nrel, function(x) {
    apply(
      key[names(m_bin[, x][m_bin[, x] > 0]), , drop = FALSE],
      2,
      all
    )
  })
  colnames(poss_g) <- colnames(m_bin)
  # number of groups to which each relevé can be assigned
  n_poss_g <- colSums(poss_g)
  names(n_poss_g) <- colnames(m_bin)
  amb_rel <- names(n_poss_g)[n_poss_g > 1] # relevés with ambiguous assignment
  unamb_rel <- names(n_poss_g)[n_poss_g == 1] # relevés with precise assignment
  iap <- length(unamb_rel) / nrel
  # possible assignments (for each relevé)
  poss_a <- cbind(possible_assignments = apply(poss_g, 2, function(x) {
    paste0("{", paste((1:k)[x], collapse = ", "), "}")
  }))
  # possible assignments of ambiguous relevés
  poss_a_amb <- cbind(groups = poss_a[n_poss_g > 1, ])
  list(
    releves_with_ambiguous_assignment = amb_rel,
    possible_assignments_of_ambiguous_releves = as.data.frame(poss_a_amb),
    iap = iap, # Internal assignment precision
    iaa = 1 - iap # Internal assignment ambiguity
  )
}
