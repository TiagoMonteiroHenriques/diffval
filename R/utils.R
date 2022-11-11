# SOME AUXILIARY FUNCTIONS

# Four auxiliary functions to HillClimb_optim_tdv =============================

# random_neighbour_hc: this function returns randomly one of the
# stoch_neigh_size-neighbour partitions, assuring that the minimum group size
# (mgs) is respected

random_neighbour_hc <- function(p, k, mgs, stoch_neigh_size) {
  tp <- tabulate(p)
  k_int <- which(tp > mgs) # k of interest to sample
  k_sam <- tp[k_int] - mgs # k samplable
  swap <- sample(rep(k_int, k_sam), min(stoch_neigh_size, sum(k_sam))) #
  # stoch_neigh_size cannot be greater than sum(k_sam)
  niter <- table(swap)
  gn_tot <- NULL
  in_tot <- NULL
  for (gv in sort(unique(swap))) { # It assumes that tabulate function
    # also sorts!
    if (length(c(1:k)[-gv]) != 1) {
      gn_tot <- c(gn_tot, sample(c(1:k)[-gv], niter[paste(gv)], replace = TRUE))
    } else {
      gn_tot <- c(gn_tot, rep(c(1:k)[-gv], niter[paste(gv)]))
    }
    in_tot <- c(in_tot, sample(which(p == gv), niter[paste(gv)]))
  }
  p[in_tot] <- gn_tot
  p
}

# max_tdv_neighbour: this function assesses all 1-neighbouring partitions of the
# given partition and returns the(one of the) partition(s) presenting the
# greater TDV

max_tdv_neighbour <- function(mt,
                              p,
                              k,
                              ns,
                              nr,
                              mgs,
                              ofda,
                              ifp,
                              afg,
                              empty_size,
                              gct,
                              i_mul,
                              dv) {
  tp <- tabulate(p)
  mtemp <- matrix(p, nr, nr, byrow = TRUE)
  fordiag <- p
  res1 <- NULL
  res2 <- NULL
  if (min(tp) > mgs) { # All groups have more than mgs elements
    for (i in 1:(k - 1)) {
      fordiag <- fordiag + 1
      fordiag[which(fordiag == k + 1)] <- 1
      diag(mtemp) <- fordiag
      res1 <- rbind(res1, mtemp)
      res2 <- c(res2, fordiag)
    }
    res2 <- cbind(rep(p, k - 1), res2)
    colnames(res2) <- NULL
  } else { # At least one group has only mgs elements
    ind_rm <- as.numeric(sapply(which(tp == mgs), function(x) { # It returns the
      # indices of the partition corresponding to groups presenting only mgs
      # elements. as.numeric is probably not necessary.
      which(p == x)
    }))
    for (i in 1:(k - 1)) {
      fordiag <- fordiag + 1
      fordiag[which(fordiag == k + 1)] <- 1
      diag(mtemp) <- fordiag
      res1 <- rbind(res1, mtemp[-ind_rm, ])
      res2 <- c(res2, fordiag[-ind_rm])
    }
    res2 <- cbind(rep(p[-ind_rm], k - 1), res2)
    colnames(res2) <- NULL
  }
  mat_neig <- list(p.list = res1, pairs = res2) # Matrix of neighbouring
  # partitions

  mat_neig_tdv <- sapply(seq_len(nrow(mat_neig$p.list)), function(x) { # Like
    # this, it is difficult to parallelize!
    pn <- mat_neig$p.list[x, ]
    kc <- mat_neig$pairs[x, ]
    return(tdv_neig(
      mt = mt, k = k, ns = ns, nr = nr, ofda = ofda, ifp = ifp,
      afg = afg, empty_size = empty_size, gct = gct,
      i_mul = i_mul, dv = dv, pn = pn, kc = kc
    ))
  })
  list(
    tdv = tdv <- max(mat_neig_tdv),
    p   = mat_neig$p.list[which(mat_neig_tdv == tdv)[1], ]
  )
}

# tdv_neig: an auxiliary function for tdv calculation (of a 1-neighbour
# partition, pn) having as starting point a partition p (aiming at efficiency of
# the optimization functions)

# ofda, ifp, afg, empty_size, gct, i_mul and dv relate to partition p
# pn  is the neighbouring partition
# kc  gives the pair of swapping groups

tdv_neig <- function(mt,
                     k,
                     ns,
                     nr,
                     ofda,
                     ifp,
                     afg,
                     empty_size,
                     gct,
                     i_mul,
                     dv,
                     pn,
                     kc) {
  tp_n <- tabulate(pn) # Size of each group (inner) in neighbour partition
  outer_size_n <- nr - tp_n # Sum of the sizes of the outer groups in neighbour
  # partition

  # Updating ofda, ifp, afg, empty_size and gct (not changing their names)
  for (i in kc) { # (Updates afg) no. of relevés containing the taxa, within
    # each group (absolute frequency in each group), only for the two groups
    # that changed a relevé!
    afg[i, ] <- colSums(mt[which(pn == i), , drop = FALSE])
  }
  i_aff_tx <- afg[kc, ][1, ] > 0 | afg[kc, ][2, ] > 0 # Taxa affected by the
  # swap of groups (i.e., present in at least one of the swapping groups)
  empty_size[kc, ] <- (afg == 0)[kc, ] * tp_n[kc] # (Updates empty_size) no. of
  # relevés of each group, when the taxon is not present. i_aff_tx must not be
  # used here!
  gct[i_aff_tx] <- colSums(afg > 0)[i_aff_tx] # (Updates gct) no. of groups
  # containing the taxon [e]
  i_mul[i_aff_tx] <- (gct > 1)[i_aff_tx] # (Updates i_mul) indices of the taxa
  # occurring in more than one group (they must occur in at least one)
  for (g in 1:k) { # Fills matrices ofda [c/d] and ifp [a/b], only when the
    # taxon is present in the group and only for the affected taxa!
    i_tx <- afg[g, ] > 0 # Indices of the taxa present in the group g
    if (sum(i_tx & i_aff_tx) > 0) { # In the case that the group g has affected
      # taxa
      ofda[g, i_tx & !i_mul] <- 1 # ofda is 1 for the taxa occurring in one
      # group only!
      if (sum(i_mul) > 0) { # If there are taxa occurring in more than one group
        i_tx_mul <- i_tx & i_mul # Taxa of group g occurring also in other group
        # than g
        ofda[g, i_tx_mul] <- colSums(empty_size[-g, i_tx_mul, drop = FALSE]
        / outer_size_n[g]) # Size of outer empty groups divided by the sum of
        # the sizes of the outer groups [c/d]
      }
      ofda[g, i_aff_tx & !i_tx] <- 0 # Inserts a zero to the affected taxa that
      # is no more present in the group!
      ifp[g, i_aff_tx] <- afg[g, i_aff_tx] / tp_n[g] # Presences inside group g
      # divided by the group g size [a/b]
    }
  }
  dv <- colSums(ifp[, i_aff_tx,
    drop = FALSE
  ] * ofda[, i_aff_tx,
    drop = FALSE
  ]) / gct[i_aff_tx]
  sum(dv) / ns
}

# tdv_aux: an auxiliary function for direct tdv calculation (aiming at
# efficiency of the optimization functions)

tdv_aux <- function(mt, p, k, ns, nr) {
  tp <- tabulate(p) # Size of each group (inner)
  outer_size <- nr - tp # Sum of the sizes of the outer groups
  ofda <- ifp <- matrix(0, k, ns) # Matrices to store [a/b], i.e., the inner
  # frequency of presences (ifp) and [c/d], i.e., the outer frequency of
  # differenciating absences (ofda)
  afg <- rowsum(mt, group = p) # No. of relevés containing the taxon, within
  # each group (absolute frequency in each group)
  empty_size <- (afg == 0) * tp # No. of relevés of each group (group size),
  # when the taxon is not present
  gct <- colSums(afg > 0) # No. of groups containing the taxon [e]
  i_mul <- gct > 1 # Indices of the taxa occurring in more than one group (taxa
  # must occur in at least one group)
  for (g in 1:k) { # Fills matrices ofda [c/d] and ifp [a/b], only when the
    # taxon is present in the group!
    i_tx <- afg[g, ] > 0 # Indices of the taxa present in the group g
    ofda[g, i_tx & !i_mul] <- 1 # ofda is 1 for the taxa occurring in one group
    # only!
    if (sum(i_mul) > 0) { # If there are taxa occurring in more than one group
      i_tx_mul <- i_tx & i_mul # Taxa of group g occurring also in other group
      # than g
      ofda[g, i_tx_mul] <- colSums(empty_size[-g, i_tx_mul, drop = FALSE]
      / outer_size[g]) # Size of outer empty groups divided by the sum of the
      # sizes of the outer groups [c/d]
    }
    ifp[g, i_tx] <- afg[g, i_tx] / tp[g] # Presences inside group g divided by
    # the group g size [a/b]
  }
  sum(colSums(ifp * ofda) / gct) / ns # For TDV1 replace gct by gct^2
}

# Two auxiliary functions to GRASP_partition_tdv ==============================

# get_dv_01: auxiliary function for GRASP efficiency
# this function uses current calculation matrix (mat_cur) to obtain, for each
# (usable) row (and for each group g!), the result of DiffVal by introducing a
# new relevé. For each row and for each group the DiffVal is presented in two
# columns considering that new relevé brings a 0 or a 1 to that row.

# usable_row  the lines of mat_cur that are still useful for the TDV calculation
#             (i.e., the respective taxa is not in all groups), but the
#             respective taxa could still be absent from all groups
# p_a_u       (present and usable) the lines of mat_cur that are not null and
#             that are still useful for the TDV calculation (i.e., the
#             respective taxa is not in all groups)
# ind_0       k indices to store values when adding 0 values
# ind_1       k indices to store values when adding 1 values

get_dv_01 <- function(k,
                      mat_cur,
                      usable_row,
                      p_a_u, ns,
                      ind_0,
                      ind_1,
                      ind_a,
                      ind_b,
                      ind_c,
                      ind_d,
                      ind_e,
                      ind_usable) {
  res_mat_01 <- matrix(0, ns, k * 2) # To store DiffVal for each group
  # (depending if the new relevé to enter the group has a 0 or 1 in the row)
  res_mat_02 <- matrix(NA, ns, k) # To store p_a_u_new for each group
  for (g in 1:k) {

    # For the 0 column
    b_local <- mat_cur[p_a_u, ind_b, drop = FALSE]
    b_local[, g] <- b_local[, g] + 1
    ab_local <- mat_cur[p_a_u, ind_a, drop = FALSE] / b_local
    c_local <- mat_cur[p_a_u, ind_c, drop = FALSE]
    c_local[, -g] <- mapply(aux_function_c_if0,
      a_wg = mat_cur[p_a_u, ind_a[g]],
      c_og = mat_cur[p_a_u, ind_c[-g]]
    )
    d_local <- mat_cur[p_a_u, ind_d, drop = FALSE]
    d_local[, -g] <- d_local[, -g] + 1
    cd_local <- c_local / d_local

    res_mat_01[p_a_u, ind_0[g]] <- rowSums(ab_local * cd_local) /
      (mat_cur[p_a_u, ind_e]) # DiffVal

    # For the 1 column
    e_param <- mat_cur[, ind_e] # Get the "e" parameter (all rows)
    # Auxiliary indices to update "e" and "usable"
    ind_for_e <- mat_cur[usable_row, ind_a[g]] == 0 # Previously absent in group
    # g, but will be present as is present in newcol
    ind_for_usable <- ind_for_e & (e_param[usable_row] == (k - 1)) #
    # Previously absent in group g, but will be present as is present in newcol,
    # and g was the only empty group before

    # Updating parameter "e" locally
    e_param[usable_row][ind_for_e] <- (e_param[usable_row][ind_for_e] +
      1)
    # Updating "usable" locally
    usable_col <- mat_cur[, ind_usable] # Get the entire column (all rows)
    usable_col[usable_row][ind_for_usable] <- 0 # Lost rows
    # Updating "usable_row" and "p_a_u"
    usable_row_new <- as.logical(usable_col) # Taxa not present in all groups
    # (taxa that are present in all groups have DiffVal = 0)
    p_a_u_new <- usable_row_new

    a_local <- mat_cur[p_a_u_new, ind_a, drop = FALSE]
    a_local[, g] <- a_local[, g] + 1
    b_local <- mat_cur[p_a_u_new, ind_b, drop = FALSE]
    b_local[, g] <- b_local[, g] + 1
    ab_local <- a_local / b_local
    c_local <- mat_cur[p_a_u_new, ind_c, drop = FALSE]
    c_local[, -g] <- mapply(aux_function_c_if1,
      a_wg = mat_cur[p_a_u_new, ind_a[g]],
      c_og = mat_cur[p_a_u_new, ind_c[-g]],
      b_wg = mat_cur[p_a_u_new, ind_b[g]]
    )
    d_local <- mat_cur[p_a_u_new, ind_d, drop = FALSE]
    d_local[, -g] <- d_local[, -g] + 1
    cd_local <- c_local / d_local

    res_mat_01[p_a_u_new, ind_1[g]] <- rowSums(ab_local * cd_local) /
      (e_param[p_a_u_new]) # DiffVal
    res_mat_02[, g] <- p_a_u_new
  }
  list(res_mat_01, res_mat_02)
}

# get_tdv_newcol: auxiliary function for GRASP efficiency (based on get_dv_01)
# this function uses current mat_cur to calculate TDV efficiently, given a
# newcol and dv_changes, which is the output of function get_dv_01.

# newcol    the index of the column of m_bin (the relevé) to be included in
#           group g
# g         the group (of the partition) where the new relevé is to be included
# p_a_u     the lines of mat_cur that are not null and that are still useful for
#           the TDV calculation (i.e., the respective taxa is not in all groups)
# ind_0     k indices to store values when adding 0 values
# ind_1     k indices to store values when adding 1 values

get_tdv_newcol <- function(m_bin,
                           newcol,
                           g,
                           p_a_u,
                           present_in_groups,
                           curr_dv,
                           dv_changes,
                           ind_0,
                           ind_1) {
  newrel <- as.logical(m_bin[, newcol])
  # Update when 0
  curr_dv[!newrel & p_a_u] <- dv_changes[[1]][, ind_0[g]][!newrel & p_a_u]
  # Update when 1
  p_a_u_new <- dv_changes[[2]][, g]
  p_a_u_final <- p_a_u
  p_a_u_final[newrel] <- p_a_u_new[newrel]
  curr_dv[newrel & p_a_u_new] <- dv_changes[[1]][, ind_1[g]][newrel & p_a_u_new]

  present_in_groups_final <- present_in_groups | newrel

  sum(curr_dv[p_a_u_final]) / sum(present_in_groups_final) # Attention this is
  # not divided by ns (the number ot rows of m_bin) as it should if TDV is
  # needed). As ns is a constant it is not needed for the optimization.
}

# Auxiliary function for GRDTP efficiency =====================================

# get_tdv_grdtp: this function uses current calculation matrix (mat_cur) to
# obtain, for each group g, the result of TDV by introducing a new relevé in to
# that group g.

# usable_row  the lines of mat_cur that are still useful for the TDV calculation
#             (i.e., the respective taxa is not in all groups), but the
#             respective taxa could still be absent from all groups
# p_a_u       (present and usable) the lines of mat_cur that are not null and
#             that are still useful for the TDV calculation (i.e., the
#             respective taxa is not in all groups)

get_tdv_grdtp <- function(m_bin,
                          k,
                          newcol,
                          mat_cur,
                          usable_row,
                          p_a_u,
                          ind_a,
                          ind_b,
                          ind_c,
                          ind_d,
                          ind_e,
                          ind_usable) {
  res_mat <- matrix(0, 1, k) # To store TDV for each group
  where_0 <- m_bin[, newcol] == 0
  where_1 <- !where_0
  for (g in 1:k) {
    # For the 0s
    where_0_pau <- p_a_u & where_0
    a_local0 <- mat_cur[where_0_pau, ind_a, drop = FALSE]
    b_local0 <- mat_cur[where_0_pau, ind_b, drop = FALSE]
    b_local0[, g] <- b_local0[, g] + 1
    ab_local0 <- a_local0 / b_local0
    c_local0 <- mat_cur[where_0_pau, ind_c, drop = FALSE]
    c_local0[, -g] <- mapply(aux_function_c_if0,
      a_wg = a_local0[, g],
      c_og = c_local0[, -g]
    )
    d_local0 <- mat_cur[where_0_pau, ind_d, drop = FALSE]
    d_local0[, -g] <- d_local0[, -g] + 1
    cd_local0 <- c_local0 / d_local0

    # For the 1s
    e_param <- mat_cur[, ind_e] # Get the "e" parameter (all rows)
    # Auxiliary indices to update "e" and "usable"
    ind_for_e <- mat_cur[usable_row, ind_a[g]] == 0 & where_1[usable_row] #
    # Previously absent in group g, but will be present as is present in newcol
    ind_for_usable <- ind_for_e & (e_param[usable_row] == (k - 1)) #
    # Previously absent in group g, but will be present as is present in newcol,
    # and g was the only empty group before

    # Updating parameter "e" locally
    e_param[usable_row][ind_for_e] <- (e_param[usable_row][ind_for_e] + 1)
    # Updating "usable" locally
    usable_col <- mat_cur[, ind_usable] # Get the entire column (all rows)
    usable_col[usable_row][ind_for_usable] <- 0 # Lost rows
    # Updating "usable_row" and "p_a_u"
    usable_row_new <- as.logical(usable_col) # Taxa not present in all groups
    # (taxa that are present in all groups have DiffVal = 0)
    p_a_u_new <- usable_row_new # Taxa present in at least one group but not
    # present in all of them

    where_1_paunew <- p_a_u_new & where_1
    a_local1 <- mat_cur[where_1_paunew, ind_a, drop = FALSE]
    b_local1 <- mat_cur[where_1_paunew, ind_b, drop = FALSE]
    c_local1 <- mat_cur[where_1_paunew, ind_c, drop = FALSE]
    d_local1 <- mat_cur[where_1_paunew, ind_d, drop = FALSE]
    c_local1[, -g] <- mapply(aux_function_c_if1,
      a_wg = a_local1[, g],
      c_og = c_local1[, -g],
      b_wg = b_local1[, g]
    )
    # a and b can only be updated after updating c, as it uses the older values
    # of a and b
    a_local1[, g] <- a_local1[, g] + 1
    b_local1[, g] <- b_local1[, g] + 1
    ab_local1 <- a_local1 / b_local1
    d_local1[, -g] <- d_local1[, -g] + 1
    cd_local1 <- c_local1 / d_local1

    res_mat[1, g] <- sum(
      sum(rowSums(ab_local0 * cd_local0) /
        (e_param[p_a_u & where_0])),
      sum(rowSums(ab_local1 * cd_local1) /
        (e_param[p_a_u_new & where_1]))
    )
  }
  res_mat
}

# Three auxiliary functions to decrease GRASP and GRDTP computation time ======

# These are functions to recalculate parameter "c" (outside group g) given a new
# column relevé to be included in group g
# These functions could be adapted to the use of function ifelse() but
# apparently the gain in time is only for smaller datasets, as mapply()
# returned faster for big datasets.
# a_wg    the "a" parameter within the group g
# newrel  the new relevé to be included
# c_og    the "c" parameter, outside the group g #this vector could be longer,
# the other vectors are recycled by the mapply function
# b_wg    the "b" parameter, within the group g

# aux_function_c:
aux_function_c <- function(a_wg, newrel, c_og, b_wg) {
  if (a_wg == 0) {
    if (newrel == 0) {
      return(c_og + 1)
    } else {
      return(c_og - b_wg)
    }
  } else {
    return(c_og)
  }
}

# aux_function_c_if0:
aux_function_c_if0 <- function(a_wg, c_og) { # Simpler function, to use when
  # newrel is a vector of 0s
  if (a_wg == 0) {
    return(c_og + 1)
  } else {
    return(c_og)
  }
}

# aux_function_c_if1:
aux_function_c_if1 <- function(a_wg, c_og, b_wg) { # Simpler function, to use
  # when newrel is a vector of 1s
  if (a_wg == 0) {
    return(c_og - b_wg)
  } else {
    return(c_og)
  }
}

# Auxiliary function for optim_tdv_simul_anne =================================

# random_neighbour_sa: selects a random neighbour (or the same partition),
# changing just one relevé

random_neighbour_sa <- function(p, nr, k) {
  tp <- table(p)
  k_int <- which(tp > 1) # Could fail when all groups have only one element...
  change_col <- sample((1:nr)[p %in% k_int], 1) # Changing just one column
  p[change_col] <- sample(1:k, length(change_col), replace = TRUE)
  p
}

# Auxiliary functions to optim_tdv_gurobi_k_2 =================================

# These are two auxiliary functions to assist the preparation of all necessary
# objects to pass to Gurobi solver, from a binary table (e.g. a binary
# phytosociological table)
# table   a binary matrix
# t			  the size of one of the groups

# optim_tdv_gurobi_td: for the t-dependent formulation
optim_tdv_gurobi_td <- function(table, t, n, alphai) {
  m <- nrow(table) # i index

  # Number of lines for each restriction
  num_restr_1 <- 1
  num_restr_2 <- sum(table)
  num_restr_3 <- num_restr_2

  # Restrictions matrix
  mat <- matrix(0, num_restr_1 + num_restr_2 + num_restr_3, n + m + m)
  nlinha <- 0

  # Restriction 1: sum x_j = t
  nlinha <- nlinha + 1
  mat[nlinha, 1:n] <- 1

  # Restriction 2: x_j + G1_i <= 1 #for a_i_j = 1
  for (i in 1:m) {
    for (j in 1:n) {
      if (table[i, j] == 1) {
        nlinha <- nlinha + 1
        mat[nlinha, c(j, n + i)] <- c(1, 1)
      }
    }
  }

  # Restriction 3: x_j - G2_i >= 0 #for a_i_j = 1
  for (i in 1:m) {
    for (j in 1:n) {
      if (table[i, j] == 1) {
        nlinha <- nlinha + 1
        mat[nlinha, c(j, n + m + i)] <- c(1, -1)
      }
    }
  }

  # Right hand side vector
  rhs <- rep(0, num_restr_1 + num_restr_2 + num_restr_3)
  rhs[1] <- t
  rhs[1 + (1:(num_restr_2))] <- 1
  rhs[(1 + num_restr_2) + (1:(num_restr_3))] <- 0

  # Objective function vector
  obj <- c(rep(0, n), alphai / (n - t), alphai / t)

  # Directions vector/sense
  sense <- c(rep("=", 1), rep("<=", num_restr_2), rep(">=", num_restr_3))

  # Objective function aim
  modelsense <- "max"

  # Variable types
  types <- c(rep("B", n), rep("C", m), rep("C", m))

  list(
    A          = mat,
    obj        = obj,
    modelsense = modelsense,
    rhs        = rhs,
    sense      = sense,
    vtype      = types
  )
}

# optim_tdv_gurobi_ti: for the t-independent formulation
optim_tdv_gurobi_ti <- function(table, n, alphai) {
  m <- nrow(table) # i index

  # Rows of restrictions matrix
  n_restr1 <- 1 # sum x_j >= 1
  n_restr2 <- n_restr1 # sum x_j <= n-1
  n_restr3 <- sum(table) # x_{j} + G1_{i} <= 1 #for a_i_j = 1
  n_restr4 <- n_restr3 # x_{j} - G2_{i} >= 0 #for a_i_j = 1
  n_restr5 <- n * m # Y1_{i} - Z1_{ij} >= 0
  n_restr6 <- n_restr5 # x_{j} - Z1_{ij} >= 0
  # restriction 7 was eliminated
  n_restr8 <- n_restr5 # x_{j} + Y1_{i} - Z1_{ij} <= 1
  n_restr9 <- m # alpha_{i}*G1_{i} - sum_{j}(Z1_{i,j}) = 0
  n_restr10 <- n_restr5 # Y2_{i} - Z2_{ij} >= 0
  n_restr11 <- n_restr5 # x_{j} - Z2_{ij} >= 0
  # restriction 12 was eliminated
  n_restr13 <- n_restr5 # x_{j} + Y2_{i} - Z2_{ij} <= 1
  n_restr14 <- m # alpha_{i}*G2_{i} - n * Y2_{i} + sum_{j}(Z2_{i,j}) = 0

  total_rest_rows <- 2 * n_restr1 + 2 * n_restr3 + 6 * n_restr5 + 2 * m

  # Columns of  restrictions matrix
  n_var1 <- n # x
  n_var2 <- m # G1
  n_var3 <- m # G2
  n_var4 <- m # Y1
  n_var5 <- m # Y2
  n_var6 <- n * m # Z1
  n_var7 <- n * m # Z2

  col_ini_g1 <- n_var1
  col_ini_g2 <- col_ini_g1 + n_var2
  col_ini_y1 <- col_ini_g2 + n_var3
  col_ini_y2 <- col_ini_y1 + n_var4
  col_ini_z1 <- col_ini_y2 + n_var5
  col_ini_z2 <- col_ini_z1 + n_var6

  total_var_cols <- col_ini_z2 + n_var7

  # Empty restriction matrix
  mat <- matrix(0, total_rest_rows, total_var_cols)
  nlinha <- 0

  # Restriction 1: sum x_j >= 1
  nlinha <- nlinha + 1
  mat[nlinha, 1:n_var1] <- 1

  # Restriction 2: sum x_j <= n-1 #TO DO: experimentar com <= floor(n/2)
  nlinha <- nlinha + 1
  mat[nlinha, 1:n_var1] <- 1

  # Restriction 3: x_{j} + G2_{i} <= 1 #para a_i_j = 1
  for (i in 1:m) {
    for (j in 1:n) {
      if (table[i, j] == 1) {
        nlinha <- nlinha + 1
        mat[nlinha, c(j, col_ini_g2 + i)] <- c(1, 1)
      }
    }
  }

  # Restriction 4: x_{j} - G1_{i} >= 0 #para a_i_j = 1
  for (i in 1:m) {
    for (j in 1:n) {
      if (table[i, j] == 1) {
        nlinha <- nlinha + 1
        mat[nlinha, c(j, col_ini_g1 + i)] <- c(1, -1)
      }
    }
  }

  # Restriction 5: Y1_{i} - Z1_{ij} >= 0
  for (i in 1:m) {
    for (j in 1:n) {
      nlinha <- nlinha + 1
      mat[nlinha, c(
        col_ini_y1 + i,
        col_ini_z1 + j + n * (i - 1)
      )] <- c(1, -1)
    }
  }

  # Restriction 6: x_{j} - Z1_{ij} >= 0
  for (i in 1:m) {
    for (j in 1:n) {
      nlinha <- nlinha + 1
      mat[nlinha, c(
        j,
        col_ini_z1 + j + n * (i - 1)
      )] <- c(1, -1)
    }
  }

  # Restriction 8: x_{j} + Y1_{i} - Z1_{ij} <= 1
  for (i in 1:m) {
    for (j in 1:n) {
      nlinha <- nlinha + 1
      mat[nlinha, c(
        j,
        col_ini_y1 + i,
        col_ini_z1 + j + n * (i - 1)
      )] <- c(1, 1, -1)
    }
  }

  # Restriction 9: alpha_{i}*G1_{i} - sum_{j}(Z1_{i,j}) = 0
  for (i in 1:m) {
    nlinha <- nlinha + 1
    mat[nlinha, col_ini_g1 + i] <- alphai[i]
    mat[
      nlinha,
      (col_ini_z1 + 1 + n * (i - 1)):(col_ini_z1 + n + n * (i - 1))
    ] <- -1
  }

  # Restriction 10: Y2_{i} - Z2_{ij} >= 0
  for (i in 1:m) {
    for (j in 1:n) {
      nlinha <- nlinha + 1
      mat[nlinha, c(
        col_ini_y2 + i,
        col_ini_z2 + j + n * (i - 1)
      )] <- c(1, -1)
    }
  }

  # Restriction 11: x_{j} - Z2_{ij} >= 0 (old: x_{j} + Z2_{ij} <= 1)
  for (i in 1:m) {
    for (j in 1:n) {
      nlinha <- nlinha + 1
      mat[nlinha, c(
        j,
        col_ini_z2 + j + n * (i - 1)
      )] <- c(1, -1)
      # )] <- c(1, 1)
    }
  }

  # Restriction 13: x_{j} + Y2_{i} - Z2_{ij} <= 1
  for (i in 1:m) {
    for (j in 1:n) {
      nlinha <- nlinha + 1
      mat[nlinha, c(
        j,
        col_ini_y2 + i,
        col_ini_z2 + j + n * (i - 1)
      )] <- c(1, 1, -1)
    }
  }

  # Restriction 14: alpha_{i}*G2_{i} - n * Y2_{i} + sum_{j}(Z2_{i,j}) = 0
  for (i in 1:m) {
    nlinha <- nlinha + 1
    mat[nlinha, c(col_ini_g2 + i, col_ini_y2 + i)] <- c(alphai[i], -n)
    mat[
      nlinha,
      (col_ini_z2 + 1 + n * (i - 1)):(col_ini_z2 + n + n * (i - 1))
    ] <- 1
  }

  # Right hand side vector
  rhs <- rep(0, total_rest_rows)
  rhs[n_restr1] <- 1
  rhs[n_restr1 + n_restr2] <- n - 1 # Possibly try with floor(n/2)
  rhs[n_restr1 + n_restr2 + (1:n_restr3)] <- 1
  rhs[n_restr1 + n_restr2 + n_restr3 + (1:n_restr4)] <- 0
  rhs[n_restr1 + n_restr2 + n_restr3 + n_restr4 + (1:n_restr5)] <- 0
  rhs[n_restr1 + n_restr2 + n_restr3 + n_restr4 + n_restr5 + (1:n_restr6)] <- 0
  rhs[n_restr1 + n_restr2 + n_restr3 + n_restr4 + n_restr5 + n_restr6 + (1:n_restr8)] <- 1
  rhs[n_restr1 + n_restr2 + n_restr3 + n_restr4 + n_restr5 + n_restr6 + n_restr8 + (1:n_restr9)] <- 0
  rhs[n_restr1 + n_restr2 + n_restr3 + n_restr4 + n_restr5 + n_restr6 + n_restr8 + n_restr9 + (1:n_restr10)] <- 0
  rhs[n_restr1 + n_restr2 + n_restr3 + n_restr4 + n_restr5 + n_restr6 + n_restr8 + n_restr9 + n_restr10 + (1:n_restr11)] <- 0
  rhs[n_restr1 + n_restr2 + n_restr3 + n_restr4 + n_restr5 + n_restr6 + n_restr8 + n_restr9 + n_restr10 + n_restr11 + (1:n_restr13)] <- 1
  rhs[n_restr1 + n_restr2 + n_restr3 + n_restr4 + n_restr5 + n_restr6 + n_restr8 + n_restr9 + n_restr10 + n_restr11 + n_restr13 + (1:n_restr14)] <- 0

  # Directions vector/sense
  sense <- c(
    rep(">=", n_restr1),
    rep("<=", n_restr2),
    rep("<=", n_restr3),
    rep(">=", n_restr4 + n_restr5 + n_restr6),
    rep("<=", n_restr8),
    rep("=", n_restr9),
    rep(">=", n_restr10 + n_restr11),
    rep("<=", n_restr13),
    rep("=", n_restr14)
  )

  # Objective function vector
  obj <- c(
    rep(0, n_var1 + n_var2 + n_var3), rep(1, n_var4 + n_var5),
    rep(0, n_var6 + n_var7)
  )

  # Objective function aim
  modelsense <- "max"

  # Variable types
  types <- c(
    rep("B", n_var1), rep("C", n_var2 + n_var3),
    rep("C", n_var4 + n_var5), rep("C", n_var6 + n_var7)
  )

  list(
    A          = mat,
    obj        = obj,
    modelsense = modelsense,
    rhs        = rhs,
    sense      = sense,
    vtype      = types
  )
}
