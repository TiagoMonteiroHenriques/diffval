# GRDTP_partition_tdv.R
#'
#' @title Obtain a partition using a Greedy-type algorithm
#'
#' @description This function obtains a partition of the columns of a given phytosociological matrix, aiming at high values of the
#' Total Differential Value (TDV), implementing a Greedy-type algorithm.
#'
#' @param m A `matrix`, i.e. a phytosociological table of 0s (absences) and 1s (presences), where rows correspond to taxa and columns correspond to relevés.
#' @param k A `numeric` giving the number of desired groups.
#' @param verify A `logical`. If `TRUE` (the default) the function verifies if basic features of `m` data structure are met. Otherwise if `FALSE`.
#'
#' @details Given the phytosociological table `m` (rows corresponding to taxa and columns corresponding to relevés), this function
#' uses a Greedy-type algorithm (a simplified version of the Greedy algorithm) to obtain a k-partition (`k`, defined by the user) of the columns of `m`, aiming at high values of TVD.
#' The algorithm operates in the following way: Firstly, `k` columns are selected randomly to work as seeds for each one of the
#' desired `k` groups. Secondly, one of the remaining columns is selected randomly and added to the partition group which maximizes
#' the upcoming TDV. This second step is repeated until all columns are placed in a group of the k-partition.
#'
#' Being a simplified version of the Greedy algorithm, it is expected to perform faster than \code{\link{GRASP_partition_tdv}}, yet
#' returning worse partitions in terms of TDV. For the (true) Greedy algorithm see \code{\link{GRASP_partition_tdv}}.
#' See \code{\link{tdv}} for an explanation on the TDV of a phytosociological table.
#'
#' @return A `numeric vector`, which length is the same as the number of columns of m, with numbers from 1 to `k`, representing the group to which
#'the respective column was ascribed.
#'
#' @author Jorge Orestes Cerdeira and Tiago Monteiro-Henriques.
#'
#' E-mail: \email{tmh.dev@@icloud.com}.
#'
#' @examples
#'
#' #getting the Taxus baccata forests data set
#' data(taxus_bin)
#'
#' #obtaining a partiton based on a Greedy-type algorithm
#' GRDTP_partition_tdv(taxus_bin, 3)
#'
#' @export
#'
GRDTP_partition_tdv <- function(m, k, verify = TRUE) {
  if (verify) {
    stopifnot(is.matrix(m))
    mode(m) <- "integer"
    if (!identical(c(0L, 1L), sort(unique(as.vector(m))))) {stop("Matrix m should contain only 0's and 1's")}
    if (min(rowSums(m))==0) {stop("At least one taxa is not present in any relev\u00e9")}
    if (min(colSums(m))==0) {stop("At least one relev\u00e9 contains no taxa")}
  }
  if (k <= 1) {stop("Given k size is too small.")}
  nr <- ncol(m) # no. of relevés
  if (k > nr) {stop("Given k size is too big.")}
  ns <- nrow(m) # no. of taxa

  #GRDTP initial partition
  par.GRDTP <- rep(0, nr)
  seed <- sample(1:nr, k) #simple random seed
  par.GRDTP[seed] <- 1:k

  #preparing mat_cur (the matrix to assist TDV calculation)
  mat_cur <- matrix(0, ns, 6 * k + 2)
  ind_a <- 1 + 0:(k-1) * 6
  ind_b <- ind_a + 1
  ind_c <- ind_b + 1
  ind_d <- ind_c + 1
  ind_ab <- ind_d + 1
  ind_cd <- ind_ab + 1
  ind_e <- k * 6 + 1
  ind_usable <- ind_e + 1

  presences.in.groups <- rowSums(m[,seed])
  absences.in.groups <- k - presences.in.groups

  #for the special case of seed (i.e. only one relevé in each of the groups):
  mat_cur[,ind_a] <- m[,seed]
  mat_cur[,ind_b] <- 1
  mat_cur[,ind_c] <- absences.in.groups - !m[,seed] #CF. THIS
  mat_cur[,ind_d] <- k-1
  mat_cur[,ind_ab] <- mat_cur[,ind_a] #at this stage mat_cur[,ind_b] is always 1
  mat_cur[,ind_cd] <- mat_cur[,ind_c]/mat_cur[,ind_d]
  mat_cur[,ind_e] <- presences.in.groups
  mat_cur[,ind_usable] <- as.numeric(mat_cur[,ind_e] != k)

  present.in.groups <- presences.in.groups != 0 #present in at least one group (i.e. not absent)
  usable.row <- as.logical(mat_cur[,ind_usable]) #taxa not present in all groups (as taxa that is in all groups have DiffVal = 0)
  p_a_u <- present.in.groups & usable.row #(present and usable) present taxa in at least in one group but not in all of them

  rel.ind <- sample((1:nr)[-seed])
  for (newcol in rel.ind) {
    if (all(!usable.row)) {warning("Attention, there are no usable taxa for TDV calculation!")}
    TDVchanges <- get_TDV(m = m, k = k, newcol = newcol, mat_cur = mat_cur, usable.row = usable.row, p_a_u = p_a_u, ind_a =ind_a, ind_b = ind_b, ind_c = ind_c, ind_d = ind_d, ind_ab = ind_ab, ind_cd = ind_cd, ind_e = ind_e, ind_usable = ind_usable)
    g <- max.col(TDVchanges) #defining g as the group with higher tdv
    par.GRDTP[newcol] <- g #updating par.GRDTP accordingly
    #update mat_cur, given newcol and g
    #auxiliary indices to update "e" and "usable"
    ind.for.e <- mat_cur[usable.row, ind_a[g]] == 0 & m[usable.row, newcol] == 1 #previously absent in group g, but will be present as is present in newcol
    ind.for.usable <- ind.for.e & mat_cur[usable.row, ind_e] == k - 1 #previously absent in group g, but will be present as is present in newcol, and g was the only empty group before
    #updating parameter "e"
    mat_cur[usable.row, ind_e][ind.for.e] <- mat_cur[usable.row, ind_e][ind.for.e] + 1
    #updating "usable"
    mat_cur[usable.row, ind_usable][ind.for.usable] <- 0
    #updating "usable.row" and "p_a_u"
    usable.row <- as.logical(mat_cur[,ind_usable]) #taxa not present in all groups (taxa that are present in all groups have DiffVal = 0)
    present.in.groups <- mat_cur[, ind_e] != 0 #taxa not absent from all groups
    p_a_u <- present.in.groups & usable.row #taxa present in at least one group but not present in all of them
    #outside group g
    #updating parameter "c"
    mat_cur[usable.row, ind_c[-g]] <- mapply(aux_function_c, a_wg = mat_cur[usable.row, ind_a[g]], newrel = m[usable.row, newcol], c_og = mat_cur[usable.row, ind_c[-g]], b_wg = mat_cur[usable.row, ind_b[g]])
    #updating parameter "d"
    mat_cur[usable.row, ind_d[-g]] <- mat_cur[usable.row, ind_d[-g]] + 1
    #updating "c/d"
    mat_cur[p_a_u, ind_cd[-g]] <- mat_cur[p_a_u, ind_c[-g]] / mat_cur[p_a_u,ind_d[-g]]
    #within group g
    #updating parameter "a"
    mat_cur[p_a_u, ind_a[g]] <- mat_cur[p_a_u, ind_a[g]] + m[p_a_u, newcol]
    #updating parameter "b"
    mat_cur[usable.row, ind_b[g]] <- mat_cur[usable.row, ind_b[g]] + 1
    #updating "a/b"
    mat_cur[p_a_u, ind_ab[g]] <- mat_cur[p_a_u,ind_a[g]] / mat_cur[p_a_u,ind_b[g]]
  }
  #print(max(TDVchanges)/ns)
  return(par.GRDTP)
}

### get_TDV

#auxiliary function for GRDTP efficiency
#this function uses current calculation matrix (mat_cur) to obtain, for each group g, the result of TDV by introducing a new relevé in to that group g.

#usable.row   the lines of mat_cur that are still useful for the TDV calculation (i.e. the respective taxa is not in all groups), but the respective taxa could still be absent from all groups
#p_a_u        (present and usable) the lines of mat_cur that are not null and that are still useful for the TDV calculation (i.e. the respective taxa is not in all groups)

get_TDV <- function(m, k, newcol, mat_cur, usable.row, p_a_u, ind_a, ind_b, ind_c, ind_d, ind_ab, ind_cd, ind_e, ind_usable) {
  res.mat <- matrix(0, 1, k) #to store TDV for each group
  where.0 <- m[,newcol] == 0
  where.1 <- !where.0
  for (g in 1:k) {
    #for the 0s
    where.0_pau <- p_a_u & where.0
    a_local0 <- mat_cur[where.0_pau, ind_a, drop = FALSE]
    b_local0 <- mat_cur[where.0_pau, ind_b, drop = FALSE]
    b_local0[, g] <- b_local0[, g] + 1
    ab_local0 <- a_local0 / b_local0
    c_local0 <- mat_cur[where.0_pau, ind_c, drop = FALSE]
    c_local0[, -g] <- mapply(aux_function_c_if0, a_wg = a_local0[,g], c_og = c_local0[,-g])
    d_local0 <- mat_cur[where.0_pau, ind_d, drop = FALSE]
    d_local0[, -g] <- d_local0[, -g] + 1
    cd_local0 <- c_local0 / d_local0

    #for the 1s
    e_parameter <- mat_cur[, ind_e] #get the e parameter (all rows)
    #auxiliary indices to update "e" and "usable"
    ind.for.e <- mat_cur[usable.row, ind_a[g]] == 0 & where.1[usable.row] #previously absent in group g, but will be present as is present in newcol
    ind.for.usable <- ind.for.e & (e_parameter[usable.row] == (k - 1)) #previously absent in group g, but will be present as is present in newcol, and g was the only empty group before
    #updating parameter "e" locally
    e_parameter[usable.row][ind.for.e] <- (e_parameter[usable.row][ind.for.e] + 1)
    #updating "usable" locally
    usable_col <- mat_cur[, ind_usable] #get the entire column (all rows)
    usable_col[usable.row][ind.for.usable] <- 0 #lost rows
    #updating "usable.row" and "p_a_u"
    usable.row_new <- as.logical(usable_col) #taxa not present in all groups (taxa that are present in all groups have DiffVal = 0)
    #present.in.groups_new <- rep(TRUE, ns) #all will have at least presence in group g (i.e. all taxa not absent from all groups)
    p_a_u_new <- usable.row_new #as present.in.groups_new is all TRUE it is the same as usable.row_new (present and usable) taxa present in at least one group but not present in all of them

    where.1_paunew <- p_a_u_new & where.1
    a_local1 <- mat_cur[where.1_paunew, ind_a, drop = FALSE]
    b_local1 <- mat_cur[where.1_paunew, ind_b, drop = FALSE]
    c_local1 <- mat_cur[where.1_paunew, ind_c, drop = FALSE]
    d_local1 <- mat_cur[where.1_paunew, ind_d, drop = FALSE]
    c_local1[,-g] <- mapply(aux_function_c_if1, a_wg = a_local1[,g], c_og = c_local1[,-g], b_wg = b_local1[,g])
    #a and b can only be updated after updating c, as it uses the older values of a and b
    a_local1[, g] <- a_local1[, g] + 1
    b_local1[, g] <- b_local1[,g] + 1
    ab_local1 <- a_local1 / b_local1
    d_local1[,-g] <- d_local1[,-g] + 1
    cd_local1 <- c_local1 / d_local1

    res.mat[1, g] <- sum(sum(rowSums(ab_local0 * cd_local0) / (e_parameter[p_a_u & where.0])), sum(rowSums(ab_local1 * cd_local1) / (e_parameter[p_a_u_new & where.1])))
  }
  return(res.mat)
}
