# GRASP_partition_tdv.R
#'
#' @title Obtain a partition using a GRASP algorithm
#'
#' @description This function obtains a partition of the columns of a given phytosociological matrix, aiming at high values of the
#'Total Differential Value (TDV) using a GRASP algorithm.
#'
#' @param m.bin A `matrix`, i.e. a phytosociological table of 0s (absences) and 1s (presences), where rows correspond to taxa and columns correspond to relevés.
#' @param k A `numeric` giving the number of desired groups.
#' @param thr A `numeric` giving a threshold value (from 0 to 1 ) with the probability used to compute the sample quantile, in order to get the best `m.bin` columns from which to select one to be include in the GRASP solution (in each step of the procedure).
#' @param verify A `logical`. If `TRUE` (the default) the function verifies if basic features of `m.bin` data structure are met. Otherwise if `FALSE`.
#'
#' @details This function uses a Greedy Randomized Adaptive Search Procedure (GRASP) to obtain a partition of `m.bin`.
#' Given a phytosociological table (`m.bin`, with rows corresponding to taxa and columns corresponding to relevés) this function searches
#' for a k-partition (`k`, defined by the user) aiming at high values of the TDV. See \code{\link{tdv}} for an explanation on the
#' TDV of a phytosociological table.
#'
#' With `thr = 1`, the algorithm corresponds to the Greedy algorithm.
#'
#' @return A `numeric vector`, which length is the same as the number of columns of `m.bin`, with numbers from 1 to `k`, representing the group to which
#' the respective column was ascribed.
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
#' #obtaining a partiton based on the GRASP algorithm
#' GRASP_partition_tdv(taxus_bin, 3)
#'
#' @export
#'
GRASP_partition_tdv <- function(m.bin, k, thr = 0.95, verify = TRUE) {
  if (verify) {
    stopifnot(is.matrix(m.bin))
    mode(m.bin) <- "integer"
    if (!identical(c(0L,1L), sort(unique(as.vector(m.bin))))) {stop("Matrix m.bin should contain only 0's and 1's")}
    if (min(rowSums(m.bin)) == 0) {stop("At least one taxa is not present in any relev\u00e9")}
    if (min(colSums(m.bin)) == 0) {stop("At least one relev\u00e9 contains no taxa")}
  }
  if (k <= 1) {stop("Given k size is too small.")}
  nr <- ncol(m.bin) # no. of relevés
  if (k > nr) {stop("Given k size is too big.")}
  ns <- nrow(m.bin) # no. of taxa

  #GRASP initial partition
  par.GRASP <- rep(0, nr)
  seed <- sample(1:nr, k) #simple random seed
  par.GRASP[seed] <- 1:k

  #preparing mat_cur (the matrix to assist DiffVal and TDV calculation)
  mat_cur <- matrix(0, ns, 6 * k + 2)
  ind_a <- 1 + 0:(k-1) * 6
  ind_b <- ind_a + 1
  ind_c <- ind_b + 1
  ind_d <- ind_c + 1
  ind_ab <- ind_d + 1
  ind_cd <- ind_ab + 1
  ind_e <- k * 6 + 1
  ind_usable <- ind_e + 1

  presences.in.groups <- rowSums(m.bin[,seed])
  absences.in.groups <- k - presences.in.groups

  #for the special case of seed (i.e. only one relevé in each of the groups):
  mat_cur[,ind_a] <- m.bin[,seed]
  mat_cur[,ind_b] <- 1
  mat_cur[,ind_c] <- absences.in.groups - !m.bin[,seed] #CF. THIS
  mat_cur[,ind_d] <- k-1
  mat_cur[,ind_ab] <- mat_cur[,ind_a] #at this stage mat_cur[,ind_b] is always 1
  mat_cur[,ind_cd] <- mat_cur[,ind_c]/mat_cur[,ind_d]
  mat_cur[,ind_e] <- presences.in.groups
  mat_cur[,ind_usable] <- as.numeric(mat_cur[,ind_e] != k)

  present.in.groups <- presences.in.groups != 0 #present in at least one group (i.e. not absent)
  usable.row <- as.logical(mat_cur[,ind_usable]) #taxa not present in all groups (as taxa that is in all groups have DiffVal = 0)
  p_a_u <- present.in.groups & usable.row #(present and usable) present taxa in at least in one group but not in all of them

  ind_0 <- (1:k-1)*2+1 #k indices to store values when adding 0 values
  ind_1 <- ind_0+1 #k indices to store values when adding 1 values

  #while (0 %in% par.GRASP & any(as.logical(mat_cur[,ind_usable]))) { #do while there is an empty column in the partition and at least one line is "usable" (CF!)
  while (0 %in% par.GRASP) { #do while there is an empty column in the partition (what happen if all lines are not usable? Try with a matrix of only 1s?)
    rel.ind <- which(par.GRASP==0)
    res.mat <- matrix(0, length(rel.ind), k)
    currDV <- rep(0, ns)
    if (all(!usable.row)) {warning("Attention, there are no usable taxa for TDV calculation!")}
    currDV[p_a_u] <- rowSums(mat_cur[p_a_u, ind_ab, drop = FALSE] * mat_cur[p_a_u, ind_cd, drop = FALSE] / mat_cur[p_a_u, ind_e])
    DVchanges <- get_DV_01(k = k, mat_cur = mat_cur, usable.row = usable.row, p_a_u = p_a_u, ns = ns, ind_0 = ind_0, ind_1 = ind_1, ind_a =ind_a, ind_b = ind_b, ind_c = ind_c, ind_d = ind_d, ind_e = ind_e, ind_usable = ind_usable)
    for (rel in 1:length(rel.ind)) {
      for (g in 1:k) {
        res.mat[rel, g] <- get_tdv_newcol(m.bin = m.bin, newcol = rel.ind[rel], g = g, p_a_u = p_a_u, present.in.groups = present.in.groups, currDV = currDV, DVchanges = DVchanges, ind_0 = ind_0, ind_1 = ind_1)
      }
    }
    res.mat[res.mat < stats::quantile(res.mat, thr)] <- 0
    res.mat.cum <- res.mat
    res.mat.cum[] <- cumsum(as.vector(res.mat))
    #set.seed(1)
    un.val <- stats::runif(1, 0, res.mat.cum[length(rel.ind), k]) #check if '/ns' could be removed previously
    ind.aux <- which(res.mat.cum >= un.val, arr.ind=TRUE)[1,]
    par.GRASP[rel.ind[ind.aux[1]]] <- ind.aux[2]

    #update mat_cur, given the selected column
    newcol <- rel.ind[ind.aux[1]]
    g <- ind.aux[2]
    #auxiliary indices to update "e" and "usable"
    ind.for.e <- mat_cur[usable.row, ind_a[g]] == 0 & m.bin[usable.row, newcol] == 1 #previously absent in group g, but will be present as is present in newcol
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
    mat_cur[usable.row, ind_c[-g]] <- mapply(aux_function_c, a_wg = mat_cur[usable.row, ind_a[g]], newrel = m.bin[usable.row, newcol], c_og = mat_cur[usable.row, ind_c[-g]], b_wg = mat_cur[usable.row, ind_b[g]])
    #updating parameter "d"
    mat_cur[usable.row, ind_d[-g]] <- mat_cur[usable.row, ind_d[-g]] + 1
    #updating "c/d"
    mat_cur[p_a_u, ind_cd[-g]] <- mat_cur[p_a_u, ind_c[-g]] / mat_cur[p_a_u,ind_d[-g]]
    #within group g
    #updating parameter "a"
    mat_cur[p_a_u, ind_a[g]] <- mat_cur[p_a_u, ind_a[g]] + m.bin[p_a_u, newcol]
    #updating parameter "b"
    mat_cur[usable.row, ind_b[g]] <- mat_cur[usable.row, ind_b[g]] + 1
    #updating "a/b"
    mat_cur[p_a_u, ind_ab[g]] <- mat_cur[p_a_u,ind_a[g]] / mat_cur[p_a_u,ind_b[g]]
  }
  return(par.GRASP)
}
