#' Obtain a partition using a Greedy-type algorithm
#'
#' This function obtains a partition of the columns of a given phytosociological
#'   matrix, aiming at high values of the Total Differential Value (TDV),
#'   implementing a Greedy-type algorithm.
#'
#' @param m_bin A matrix. A phytosociological table of 0s (absences) and 1s
#'   (presences), where rows correspond to taxa and columns correspond to
#'   relevés.
#' @param k A numeric giving the number of desired groups.
#' @param verify A logical. If `TRUE` (the default) the function verifies if
#'   basic features of `m_bin` data structure are met. Otherwise if `FALSE`.
#'
#' @details Given the phytosociological table `m_bin` (rows corresponding to
#'   taxa and columns corresponding to relevés), this function uses a
#'   Greedy-type algorithm (a simplified version of the Greedy algorithm) to
#'   obtain a `k`-partition (`k`, defined by the user) of the columns of
#'   `m_bin`, aiming at high values of TDV.
#'   The algorithm operates in the following way: Firstly, `k` columns are
#'   selected randomly to work as seeds for each one of the desired `k` groups.
#'   Secondly, one of the remaining columns is selected randomly and added to
#'   the partition group which maximizes the upcoming TDV. This second step is
#'   repeated until all columns are placed in a group of the `k`-partition.
#'
#'   This function is expected to perform faster than [partition_tdv_grasp()],
#'   yet returning worse partitions in terms of TDV. For the (true) Greedy
#'   algorithm see [partition_tdv_grasp()].
#'   See [tdv()] for an explanation on the TDV of a phytosociological table.
#'
#' @return A numeric vector, which length is the same as the number of columns
#'   of `m_bin`, with numbers from 1 to `k`, representing the group to which the
#'   respective column was ascribed.
#'
#' @author Jorge Orestes Cerdeira and Tiago Monteiro-Henriques.
#'   E-mail: \email{tmh.dev@@icloud.com}.
#'
#' @examples
#' # Getting the Taxus baccata forests data set
#' data(taxus_bin)
#'
#' # Obtaining a partiton based on a Greedy-type algorithm
#' partition_tdv_grdtp(taxus_bin, 3)
#'
#' @export
partition_tdv_grdtp <- function(m_bin, k, verify = TRUE) {
  if (verify) {
    stopifnot(is.matrix(m_bin))
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
  }
  if (k <= 1) {
    stop("Given `k` size is too small.")
  }
  nr <- ncol(m_bin) # No. of relevés
  if (k > nr) {
    stop("Given `k` size is too big.")
  }
  ns <- nrow(m_bin) # No. of taxa

  # GRDTP initial partition
  par_grdtp <- rep(0, nr)
  seed <- sample(1:nr, k) # Simple random seed
  par_grdtp[seed] <- 1:k

  # Preparing mat_cur (the matrix to assist TDV calculation)
  mat_cur <- matrix(0, ns, 6 * k + 2)
  ind_a <- 1 + 0:(k - 1) * 6
  ind_b <- ind_a + 1
  ind_c <- ind_b + 1
  ind_d <- ind_c + 1
  ind_ab <- ind_d + 1
  ind_cd <- ind_ab + 1
  ind_e <- k * 6 + 1
  ind_usable <- ind_e + 1

  presences_in_groups <- rowSums(m_bin[, seed])
  absences_in_groups <- k - presences_in_groups

  # For the special case of seed (i.e., only one relevé in each of the groups)
  mat_cur[, ind_a] <- m_bin[, seed]
  mat_cur[, ind_b] <- 1
  mat_cur[, ind_c] <- absences_in_groups - !m_bin[, seed] # CF. THIS
  mat_cur[, ind_d] <- k - 1
  mat_cur[, ind_ab] <- mat_cur[, ind_a] # At this stage mat_cur[,ind_b] is
  # always 1
  mat_cur[, ind_cd] <- mat_cur[, ind_c] / mat_cur[, ind_d]
  mat_cur[, ind_e] <- presences_in_groups
  mat_cur[, ind_usable] <- as.numeric(mat_cur[, ind_e] != k)

  present_in_groups <- presences_in_groups != 0 # Present in at least one group
  # (i.e., not absent)
  usable_row <- as.logical(mat_cur[, ind_usable]) # Taxa not present in all
  # groups (as taxa that is in all groups have DiffVal = 0)
  p_a_u <- present_in_groups & usable_row # p_a_u (present and usable). Present
  # taxa in at least in one group but not in all of them.

  rel_ind <- sample((1:nr)[-seed])
  for (newcol in rel_ind) {
    if (all(!usable_row)) {
      warning("Attention, there are no usable taxa for TDV calculation!")
    }
    tdv_changes <- get_tdv_grdtp(
      m_bin = m_bin,
      k = k,
      newcol = newcol,
      mat_cur = mat_cur,
      usable_row = usable_row,
      p_a_u = p_a_u,
      ind_a = ind_a,
      ind_b = ind_b,
      ind_c = ind_c,
      ind_d = ind_d,
      ind_e = ind_e,
      ind_usable = ind_usable
    )
    g <- max.col(tdv_changes) # Defining g as the group with higher TDV
    par_grdtp[newcol] <- g # Updating par_grdtp accordingly
    # Update mat_cur, given newcol and g
    # Auxiliary indices to update "e" and "usable"
    ind_for_e <- mat_cur[usable_row, ind_a[g]] == 0 &
      m_bin[usable_row, newcol] == 1 # Previously absent in group g, but will be
    # present as is present in newcol
    ind_for_usable <- ind_for_e & mat_cur[usable_row, ind_e] == k - 1 #
    # Previously absent in group g, but will be present as is present in newcol,
    # and g was the only empty group before
    # Updating parameter "e"
    mat_cur[usable_row, ind_e][ind_for_e] <-
      mat_cur[usable_row, ind_e][ind_for_e] + 1
    # Updating "usable"
    mat_cur[usable_row, ind_usable][ind_for_usable] <- 0
    # Updating "usable_row" and "p_a_u"
    usable_row <- as.logical(mat_cur[, ind_usable]) # Taxa not present in all
    # groups (taxa that are present in all groups have DiffVal = 0)
    present_in_groups <- mat_cur[, ind_e] != 0 # Taxa not absent from all groups
    p_a_u <- present_in_groups & usable_row # Taxa present in at least one group
    # but not present in all of them

    # Outside group g
    # Updating parameter "c"
    mat_cur[usable_row, ind_c[-g]] <- mapply(
      aux_function_c,
      a_wg = mat_cur[usable_row, ind_a[g]],
      newrel = m_bin[usable_row, newcol],
      c_og = mat_cur[usable_row, ind_c[-g]],
      b_wg = mat_cur[usable_row, ind_b[g]]
    )
    # Updating parameter "d"
    mat_cur[usable_row, ind_d[-g]] <- mat_cur[usable_row, ind_d[-g]] + 1
    # Updating "c/d"
    mat_cur[p_a_u, ind_cd[-g]] <- mat_cur[p_a_u, ind_c[-g]] /
      mat_cur[p_a_u, ind_d[-g]]

    # Within group g
    # Updating parameter "a"
    mat_cur[p_a_u, ind_a[g]] <- mat_cur[p_a_u, ind_a[g]] + m_bin[p_a_u, newcol]
    # Updating parameter "b"
    mat_cur[usable_row, ind_b[g]] <- mat_cur[usable_row, ind_b[g]] + 1
    # Updating "a/b"
    mat_cur[p_a_u, ind_ab[g]] <- mat_cur[p_a_u, ind_a[g]] /
      mat_cur[p_a_u, ind_b[g]]
  }
  par_grdtp
}
