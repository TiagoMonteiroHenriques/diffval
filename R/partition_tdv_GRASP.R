#' Obtain a partition using a GRASP algorithm
#'
#' This function obtains a partition of the columns of a given phytosociological
#'   matrix, aiming at high values of the Total Differential Value (TDV) using a
#'   GRASP algorithm.
#'
#' @param m_bin A matrix. A phytosociological table of 0s (absences) and 1s
#'   (presences), where rows correspond to taxa and columns correspond to
#'   relevés.
#' @param k A numeric giving the number of desired groups.
#' @param thr A numeric giving a threshold value (from 0 to 1 ) with the
#'   probability used to compute the sample quantile, in order to get the best
#'   `m_bin` columns from which to select one to be include in the GRASP
#'   solution (in each step of the procedure).
#' @param verify A logical. If `TRUE` (the default) the function verifies if
#'   basic features of `m_bin` data structure are met. Otherwise if `FALSE`.
#'
#' @details This function uses a Greedy Randomized Adaptive Search Procedure
#'   (GRASP) to obtain a partition of `m_bin`.
#'   Given a phytosociological table (`m_bin`, with rows corresponding to taxa
#'   and columns corresponding to relevés) this function searches for a
#'   `k`-partition (`k`, defined by the user) aiming at high values of the TDV.
#'   See [tdv()] for an explanation on the TDV of a phytosociological table.
#'
#'   With `thr = 1`, the algorithm corresponds to the Greedy algorithm.
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
#' # Obtaining a partition based on the GRASP algorithm
#' partition_tdv_grasp(taxus_bin, 3)
#'
#' @export
partition_tdv_grasp <- function(m_bin, k, thr = 0.95, verify = TRUE) {
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
  nr <- ncol(m_bin) # no. of relevés
  if (k > nr) {
    stop("Given `k` size is too big.")
  }
  ns <- nrow(m_bin) # no. of taxa

  # GRASP initial partition
  par_grasp <- rep(0, nr)
  seed <- sample(1:nr, k) # Simple random seed
  par_grasp[seed] <- 1:k

  # Preparing mat_cur (the matrix to assist DiffVal and TDV calculation)
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

  ind_0 <- (1:k - 1) * 2 + 1 # k indices to store values when adding 0 values
  ind_1 <- ind_0 + 1 # k indices to store values when adding 1 values

  # while (0 %in% par_grasp & any(as.logical(mat_cur[,ind_usable]))) { # Do
  # while there is an empty column in the partition and at least one line is
  # "usable" (CF!)
  while (0 %in% par_grasp) { # Do while there is an empty column in the
    # partition. What happens if all lines are not usable? Try with a matrix of
    # only 1s?
    rel_ind <- which(par_grasp == 0)
    res_mat <- matrix(0, length(rel_ind), k)
    curr_dv <- rep(0, ns)
    if (all(!usable_row)) {
      warning("Attention, there are no usable taxa for TDV calculation!")
    }
    curr_dv[p_a_u] <- rowSums(
      mat_cur[p_a_u, ind_ab, drop = FALSE] *
        mat_cur[p_a_u, ind_cd, drop = FALSE] / mat_cur[p_a_u, ind_e]
    )
    dv_changes <- get_dv_01(
      k = k,
      mat_cur = mat_cur,
      usable_row = usable_row,
      p_a_u = p_a_u,
      ns = ns,
      ind_0 = ind_0,
      ind_1 = ind_1,
      ind_a = ind_a,
      ind_b = ind_b,
      ind_c = ind_c,
      ind_d = ind_d,
      ind_e = ind_e,
      ind_usable = ind_usable
    )
    for (rel in seq_along(rel_ind)) {
      for (g in 1:k) {
        res_mat[rel, g] <- get_tdv_newcol(m_bin,
          newcol = rel_ind[rel],
          g = g,
          p_a_u = p_a_u,
          present_in_groups = present_in_groups,
          curr_dv = curr_dv,
          dv_changes = dv_changes,
          ind_0 = ind_0,
          ind_1 = ind_1
        )
      }
    }
    res_mat[res_mat < stats::quantile(res_mat, thr)] <- 0
    res_mat_cum <- res_mat
    res_mat_cum[] <- cumsum(as.vector(res_mat))
    un_val <- stats::runif(1, 0, res_mat_cum[length(rel_ind), k])
    ind_aux <- which(res_mat_cum >= un_val, arr.ind = TRUE)[1, ]
    par_grasp[rel_ind[ind_aux[1]]] <- ind_aux[2]

    # Update mat_cur, given the selected column
    newcol <- rel_ind[ind_aux[1]]
    g <- ind_aux[2]
    # Auxiliary indices to update "e" and "usable"
    ind_for_e <- mat_cur[usable_row, ind_a[g]] == 0 &
      m_bin[usable_row, newcol] == 1 # Previously absent in group g, but will be
    # present as it is present in newcol
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
  par_grasp
}
