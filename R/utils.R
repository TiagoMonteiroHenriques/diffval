#SOME AUXILIARY FUNCTIONS

#Four auxiliary functions to HillClimb_optim_tdv

#random.neighbour: this function returns randomly one of the stoch.neigh.size-neighbour partitions, assuring that the minimum group size (mgs) is respected
random.neighbour <- function (p, k, mgs, stoch.neigh.size) {
  tp <- tabulate(p)
  k.int <- which(tp > mgs) #k of interest to sample
  k.sam <- tp[k.int] - mgs #k samplable
  swap <- sample(rep(k.int, k.sam), min(stoch.neigh.size, sum(k.sam))) #neighbourhood.size cannot be greater than sum(k.sam)
  niter <- table(swap)
  gn.tot <- NULL
  in.tot <- NULL
  for (gv in sort(unique(swap))) { #it assumes that table function also sorts!
    if (length(c(1:k)[-gv]) != 1) {
      gn.tot <- c(gn.tot, sample(c(1:k)[-gv], niter[paste(gv)], replace = TRUE))
    } else {
      gn.tot <- c(gn.tot, rep(c(1:k)[-gv], niter[paste(gv)]))
    }
    in.tot <- c(in.tot, sample(which(p == gv), niter[paste(gv)]))
  }
  p[in.tot] <- gn.tot
  return(p)
}

#max.tdv.neighbour: this function assesses all 1-neighbouring partitions of the given partition and returns the(one of the) partition(s) presenting the greater TDV
max.tdv.neighbour <- function (mt, p, k, ns, nr, mgs, ofda, ifp, afg, empty.size, gct, i.mul, DV) {
  tp <- tabulate(p)
  mtemp <- matrix(p, nr, nr, byrow = TRUE)
  fordiag <- p
  res1 <- NULL
  res2 <- NULL
  if (min(tp) > mgs) { #all groups have more than mgs elements
    for (i in 1:(k - 1)) {
      fordiag <- fordiag + 1
      fordiag[which(fordiag == k + 1)] <- 1
      diag(mtemp) <- fordiag
      res1 <- rbind(res1, mtemp)
      res2 <- c(res2, fordiag)
    }
    res2 <- cbind(rep(p, k - 1), res2)
    colnames(res2) <- NULL
  } else { #at least one group has only mgs elements
    ind.rm <- as.numeric(sapply(which(tp == mgs), function (x) { #it returns the indices of the partition corresponding to groups presenting only mgs elements. as.numeric is probably not necessary
      which(p == x)
    }))
    for (i in 1:(k - 1)) {
      fordiag <- fordiag + 1
      fordiag[which(fordiag == k + 1)] <- 1
      diag(mtemp) <- fordiag
      res1 <- rbind(res1, mtemp[-ind.rm,])
      res2 <- c(res2, fordiag[-ind.rm])
    }
    res2 <- cbind(rep(p[-ind.rm], k - 1), res2)
    colnames(res2) <- NULL
  }
  mat.neig <- list(p.list = res1, pairs = res2) #matrix of neighbouring partitions

  mat.neig.tdv <- sapply(1:nrow(mat.neig$p.list), function (x) { #like this, it is difficult to parallelize!
    pn <- mat.neig$p.list[x,]
    kc <- mat.neig$pairs[x,]
    return(tdv.neig(mt = mt, k = k, ns = ns, nr = nr, ofda = ofda, ifp =ifp, afg = afg, empty.size = empty.size, gct = gct, i.mul = i.mul, DV = DV, pn = pn, kc = kc))
  })
  return(list(tdv = tdv <- max(mat.neig.tdv), p = mat.neig$p.list[which(mat.neig.tdv == tdv)[1],]))
}

#tdv.neig: an auxiliary function for tdv calculation (of a 1-neighbour partition, pn) having as starting point a partition p (aiming at efficiency of the optimization functions)
#ofda, ifp, afg, empty.size, gct, i.mul and DV relate to partition p
#pn is the neighbouring partition
#kc gives the pair of swapping groups
tdv.neig <- function (mt, k, ns, nr, ofda, ifp, afg, empty.size, gct, i.mul, DV, pn, kc) {
  tp.n <- tabulate(pn) #size of each group (inner) in neighbour partition
  outer.size.n <- nr - tp.n #sum of the sizes of the outer groups in neighbour partition

  #updating ofda, ifp, afg, empty.size and gct (not changing their names)
  for (i in kc) { #(updates afg) no. of relevés containing the taxa, within each group (absolute frequency in each group), only for the two groups that changed a relevé!
    afg[i,] <- colSums(mt[which(pn == i),,drop = FALSE])
  }
  i.aff.tx <- afg[kc,][1,] > 0 | afg[kc,][2,] > 0 #taxa affected by the swap of groups (i.e. present in at least one of the swapping groups)
  empty.size[kc,] <- (afg == 0)[kc,] * tp.n[kc] #(updates empty.size) no. of relevés of each group, when the taxon is not present #i.aff.tx must not be used here!
  gct[i.aff.tx] <- colSums(afg > 0)[i.aff.tx] #(updates gct) no. of groups containing the taxon [e]
  i.mul[i.aff.tx] <- (gct > 1)[i.aff.tx] #(updates i.mul) indices of the taxa occurring in more than one group (they must occur in at least one)
  for (g in 1:k) { #fills matrices ofda [c/d] and ifp [a/b], only when the taxon is present in the group and only for the affected taxa!
    i.tx <- afg[g,] > 0 #indices of the taxa present in the group g
    if (sum(i.tx & i.aff.tx) > 0) { #in the case that the group g has affected taxa
      ofda[g,i.tx & !i.mul] <- 1 #ofda is 1 for the taxa occurring in one group only!
      if (sum(i.mul) > 0) { #if there are taxa occurring in more than one group
        i.tx.mul <- i.tx & i.mul #taxa of group g occurring also in other group than g
        ofda[g,i.tx.mul] <- colSums(empty.size[-g,i.tx.mul,drop = FALSE] / outer.size.n[g]) #size of outer empty groups divided by the sum of the sizes of the outer groups [c/d]
      }
      ofda[g,i.aff.tx & !i.tx] <- 0 #inserts a zero to the affected taxa that is no more present in the group!
      ifp[g,i.aff.tx] <- afg[g,i.aff.tx] / tp.n[g] #presences inside group g divided by the group g size [a/b]
    }
  }
  DV <- colSums(ifp[,i.aff.tx,drop = FALSE] * ofda[,i.aff.tx,drop = FALSE]) / gct[i.aff.tx]
  return(sum(DV)/ns)
}

#tdv.aux: an auxiliary function for direct tdv calculation (aiming at efficiency of the optimization functions)
tdv.aux <- function (mt, p, k, ns, nr) {
  tp <- tabulate(p) #size of each group (inner)
  outer.size <- nr - tp #sum of the sizes of the outer groups
  ofda <- ifp <- matrix(0, k, ns) #matrices to store [a/b], i.e. the inner frequency of presences (ifp) and [c/d], i.e. the outer frequency of differenciating absences (ofda)
  afg <- rowsum(mt, group = p) #no. of relevés containing the taxon, within each group (absolute frequency in each group)
  empty.size <- (afg == 0) * tp #no. of relevés of each group (group size), when the taxon is not present
  gct <- colSums(afg > 0) #no. of groups containing the taxon [e]
  i.mul <- gct > 1 #indices of the taxa occurring in more than one group (taxa must occur in at least one group)
  for (g in 1:k) { #fills matrices ofda [c/d] and ifp [a/b], only when the taxon is present in the group!
    i.tx <- afg[g,] > 0 #indices of the taxa present in the group g
    ofda[g, i.tx & !i.mul] <- 1 #ofda is 1 for the taxa occurring in one group only!
    if (sum(i.mul) > 0) { #if there are taxa occurring in more than one group
      i.tx.mul <- i.tx & i.mul #taxa of group g occurring also in other group than g
      ofda[g,i.tx.mul] <- colSums(empty.size[-g,i.tx.mul,drop = FALSE] / outer.size[g]) #size of outer empty groups divided by the sum of the sizes of the outer groups [c/d]
    }
    ifp[g,i.tx] <- afg[g,i.tx] / tp[g] #presences inside group g divided by the group g size [a/b]
  }
  return(sum(colSums(ifp * ofda) / gct) / ns) #for TDV1 replace gct by gct^2
}

#Two auxiliary functions to GRASP_partition_tdv

#get_DV_01: auxiliary function for GRASP efficiency
#this function uses current calculation matrix (mat_cur) to obtain, for each (usable) row (and for each group g!), the result of DiffVal by introducing a new relevé. For each row and for each group the DiffVal is presented in two columns considering that new relevé brings a 0 or a 1 to that row.
#usable.row   the lines of mat_cur that are still useful for the TDV calculation (i.e. the respective taxa is not in all groups), but the respective taxa could still be absent from all groups
#p_a_u        (present and usable) the lines of mat_cur that are not null and that are still useful for the TDV calculation (i.e. the respective taxa is not in all groups)
#ind_0        k indices to store values when adding 0 values
#ind_1        k indices to store values when adding 1 values

get_DV_01 <- function(k, mat_cur, usable.row, p_a_u, ns, ind_0, ind_1, ind_a, ind_b, ind_c, ind_d, ind_e, ind_usable) {
  res.mat.01 <- matrix(0, ns, k*2) #to store DiffVal for each group (depending if the new relevé to enter the group has a 0 or 1 in the row)
  res.mat.02 <- matrix(NA, ns, k) #to store p_a_u_new for each group
  for (g in 1:k) {

    #for the 0 column
    b_local <- mat_cur[p_a_u, ind_b, drop = FALSE]
    b_local[, g] <- b_local[, g] + 1
    ab_local <- mat_cur[p_a_u, ind_a, drop = FALSE] / b_local
    c_local <- mat_cur[p_a_u, ind_c, drop = FALSE]
    c_local[, -g] <- mapply(aux_function_c_if0, a_wg = mat_cur[p_a_u, ind_a[g]], c_og = mat_cur[p_a_u, ind_c[-g]])
    d_local <- mat_cur[p_a_u, ind_d, drop = FALSE]
    d_local[, -g] <- d_local[, -g] + 1
    cd_local <- c_local / d_local

    res.mat.01[p_a_u, ind_0[g]] <- rowSums(ab_local * cd_local) / (mat_cur[p_a_u, ind_e]) #DiffVal

    #for the 1 column
    e_parameter <- mat_cur[, ind_e] #get the e parameter (all rows)
    #auxiliary indices to update "e" and "usable"
    ind.for.e <- mat_cur[usable.row, ind_a[g]] == 0 #previously absent in group g, but will be present as is present in newcol
    ind.for.usable <- ind.for.e & (e_parameter[usable.row] == (k - 1)) #previously absent in group g, but will be present as is present in newcol, and g was the only empty group before
    #updating parameter "e" locally
    e_parameter[usable.row][ind.for.e] <- (e_parameter[usable.row][ind.for.e] + 1)
    #updating "usable" locally
    usable_col <- mat_cur[, ind_usable] #get the entire column (all rows)
    usable_col[usable.row][ind.for.usable] <- 0 #lost rows
    #updating "usable.row" and "p_a_u"
    usable.row_new <- as.logical(usable_col) #taxa not present in all groups (taxa that are present in all groups have DiffVal = 0)
    p_a_u_new <- usable.row_new

    a_local <- mat_cur[p_a_u_new,ind_a, drop = FALSE]
    a_local[, g] <- a_local[, g] + 1
    b_local <- mat_cur[p_a_u_new, ind_b, drop = FALSE]
    b_local[, g] <- b_local[,g] + 1
    ab_local <- a_local / b_local
    c_local <- mat_cur[p_a_u_new, ind_c, drop = FALSE]
    c_local[,-g] <- mapply(aux_function_c_if1, a_wg = mat_cur[p_a_u_new, ind_a[g]], c_og = mat_cur[p_a_u_new, ind_c[-g]], b_wg = mat_cur[p_a_u_new, ind_b[g]])
    d_local <- mat_cur[p_a_u_new, ind_d, drop = FALSE]
    d_local[,-g] <- d_local[,-g] + 1
    cd_local <- c_local / d_local

    res.mat.01[p_a_u_new, ind_1[g]] <- rowSums(ab_local * cd_local) / (e_parameter[p_a_u_new]) #DiffVal
    res.mat.02[,g] <- p_a_u_new
  }
  return(list(res.mat.01, res.mat.02))
}

### get_tdv_newcol: auxiliary function for GRASP efficiency (based on get_DV_01)
#this function uses current mat_cur to calculate TDV efficiently, given a newcol and DVchanges, which is the output of function get_DV_01.
#newcol           the index of the column of m.bin (the relevé) to be included in group g
#g                the group (of the partition) where the new relevé is to be included
#p_a_u            the lines of mat_cur that are not null and that are still useful for the TDV calculation (i.e. the respective taxa is not in all groups)
#ind_0            k indices to store values when adding 0 values
#ind_1            k indices to store values when adding 1 values
#ns               the number ot rows of m.bin (not used now). It would be needed to calculate TDV, but as ns is a constant it is not needed in this auxiliary function.

get_tdv_newcol <- function(m.bin, newcol, g, p_a_u, present.in.groups, currDV, DVchanges, ind_0, ind_1) {
  newrel <- as.logical(m.bin[, newcol])
  #update when 0
  currDV[!newrel & p_a_u] <- DVchanges[[1]][, ind_0[g]][!newrel & p_a_u]
  #update when 1
  p_a_u_new <- DVchanges[[2]][,g]
  p_a_u_final <- p_a_u
  p_a_u_final[newrel] <- p_a_u_new[newrel]
  currDV[newrel & p_a_u_new] <- DVchanges[[1]][, ind_1[g]][newrel & p_a_u_new]

  present.in.groups_final <- present.in.groups | newrel

  return(sum(currDV[p_a_u_final]) / sum(present.in.groups_final)) #Attention this is not divided by ns (as it should if TDV is needed). As ns is a constant it is not needed for the optimization.
}

#Auxiliary function for GRDTP efficiency
#get_TDV: this function uses current calculation matrix (mat_cur) to obtain, for each group g, the result of TDV by introducing a new relevé in to that group g.
#usable.row   the lines of mat_cur that are still useful for the TDV calculation (i.e. the respective taxa is not in all groups), but the respective taxa could still be absent from all groups
#p_a_u        (present and usable) the lines of mat_cur that are not null and that are still useful for the TDV calculation (i.e. the respective taxa is not in all groups)

get_TDV <- function(m, k, newcol, mat_cur, usable.row, p_a_u, ind_a, ind_b, ind_c, ind_d, ind_e, ind_usable) {
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

#Three auxiliary functions to decrease GRASP and GRDTP computation time
#These are functions to recalculate parameter "c" (outside group g) given a new column relevé to be included in group g
#these functions could be adapted to the use of function ifelse(); apparently the gain in time is only for smaller datasets,
#as mapply() returned faster for big datasets.
#a_wg   the "a" parameter within the group g
#newrel the new relevé to be included
#c_og   the "c" parameter, outside the group g #this vector could be longer, the other vectors are recycled by the mapply function
#b_wg   the "b" parameter, within the group g

#aux_function_c:
aux_function_c <- function(a_wg, newrel, c_og, b_wg) {
  if (a_wg == 0) {
    if (newrel == 0) {
      return (c_og + 1)
    } else {
      return (c_og - b_wg)
    }
  } else {
    return(c_og)
  }
}

#aux_function_c_if0:
aux_function_c_if0 <- function(a_wg, c_og) { #simpler function, to use when newrel is a vector of 0s
  if (a_wg == 0) {
    return (c_og + 1)
  } else {
    return(c_og)
  }
}

#aux_function_c_if1:
aux_function_c_if1 <- function(a_wg, c_og, b_wg) {  #simpler function, to use when newrel is a vector of 1s
  if (a_wg == 0) {
    return (c_og - b_wg)
  } else {
    return(c_og)
  }
}

#Auxiliary function for SimulAnne_optim_tdv
#random.neighbour.SA: selects a random neighbour (or the same partition), changing just one relevé

random.neighbour.SA <- function (p, nr, k) {
  tp <- table(p)
  k.int <- which(tp > 1) #Check behaviour when all groups have only one element...
  change.col <- sample((1:nr)[p %in% k.int], 1) #Changing just one column
  #tp[p[change.col]] <- tp[p[change.col]] - 1 #Changing just one column
  p[change.col] <- sample(1:k, length(change.col), replace=TRUE)
  return(p)
}

# Auxiliary functions to GUROBI_k2_optim_tdv
# These are two auxiliary functions to assist the preparation of all necessary objects to pass to Gurobi solver, from a binary table (e.g. a binary phytosociological table)
#table 		a binary matrix
#t			 	the size of one of the groups

#optim_tdv_gurobi_td: for the t-dependent formulation
optim_tdv_gurobi_td <- function (table, t, n, alphai) {

  m <- nrow(table) #i index

  #number of lines for each restriction
  num.restr_1 <- 1
  num.restr_2 <- sum(table)
  num.restr_3 <- num.restr_2

  #restrictions matrix
  mat <- matrix(0, num.restr_1 + num.restr_2 + num.restr_3, n + m + m)
  nlinha <- 0

  #restriction 1
  #sum x_j = t
  nlinha <- nlinha + 1
  mat[nlinha, 1:n] <- 1

  #restriction 2
  #x_j + G1_i ≤ 1 #for a_i_j = 1
  for (i in 1:m) {
    for (j in 1:n) {
      if (table[i,j] == 1) {
        nlinha <- nlinha + 1
        mat[nlinha, c(j, n+i)] <- c(1, 1)
      }
    }
  }

  #restriction 3
  #x_j - G2_i ≥ 0 #for a_i_j = 1
  for (i in 1:m) {
    for (j in 1:n) {
      if (table[i,j] == 1) {
        nlinha <- nlinha + 1
        mat[nlinha, c(j, n+m+i)] <- c(1, -1)
      }
    }
  }

  #right hand side vector
  rhs <- rep(0, num.restr_1 + num.restr_2 + num.restr_3)
  rhs[1] <- t
  rhs[1+(1:(num.restr_2))] <- 1
  rhs[(1+num.restr_2)+(1:(num.restr_3))] <- 0

  #objective function vector
  obj <- c(rep(0, n), alphai/(n-t), alphai/t)

  #directions vector/sense
  sense <- c(rep("=",1), rep("<=", num.restr_2), rep(">=", num.restr_3))

  #objective function aim
  modelsense <- "max"

  #variable types
  types <- c(rep("B", n), rep("C", m), rep("C", m))

  return(list(
    A=mat,
    obj=obj,
    modelsense=modelsense,
    rhs=rhs,
    sense=sense,
    vtype=types
  ))
}

# optim_tdv_gurobi_ti: for the t-independent formulation
optim_tdv_gurobi_ti <- function (table, n, alphai) {
  m <- nrow(table) #i index

  #rows of restrictions matrix
  n.restr1 <- 1					#sum x_j >= 1
  n.restr2 <- n.restr1			#sum x_j <= n-1
  n.restr3 <- sum(table)		#x_{j} + G1_{i} ≤ 1 #for a_i_j = 1
  n.restr4 <- n.restr3			#x_{j} - G2_{i} ≥ 0 #for a_i_j = 1
  n.restr5 <- n*m				#Y1_{i} - Z1_{ij} >= 0
  n.restr6 <- n.restr5			#x_{j} - Z1_{ij} >= 0
  #restriction 7 was abandoned
  n.restr8 <- n.restr5			#x_{j} + Y1_{i} - Z1_{ij} <= 1
  n.restr9 <- m					#alpha_{i}*G1_{i} - sum_{j}(Z1_{i,j}) = 0
  n.restr10 <- n.restr5		#Y2_{i} - Z2_{ij} >= 0
  n.restr11 <- n.restr5		#x_{j} - Z2_{ij} >= 0
  #restriction 12 was abandoned
  n.restr13 <- n.restr5		#x_{j} + Y2_{i} - Z2_{ij} <= 1
  n.restr14 <- m					#alpha_{i}*G2_{i} - n * Y2_{i} + sum_{j}(Z2_{i,j}) = 0

  total.rest.rows <- 2 * n.restr1 + 2 * n.restr3 + 6 * n.restr5 + 2 * m

  #columns of  restrictions matrix
  n.var1 <- n #x
  n.var2 <- m #G1
  n.var3 <- m #G2
  n.var4 <- m #Y1
  n.var5 <- m #Y2
  n.var6 <- n*m #Z1
  n.var7 <- n*m #Z2

  col_ini_G1 <- n.var1
  col_ini_G2 <- col_ini_G1 + n.var2
  col_ini_Y1 <- col_ini_G2 + n.var3
  col_ini_Y2 <- col_ini_Y1 + n.var4
  col_ini_Z1 <- col_ini_Y2 + n.var5
  col_ini_Z2 <- col_ini_Z1 + n.var6

  total.var.cols <- col_ini_Z2 + n.var7

  #empty restriction matrix
  mat <- matrix(0, total.rest.rows, total.var.cols)
  nlinha <- 0

  #RESTRICTION 1
  #sum x_j >= 1
  nlinha <- nlinha + 1
  mat[nlinha, 1:n.var1] <- 1

  #RESTRICTION 2
  #sum x_j <= n-1 #TO DO: experimentar com <= floor(n/2)
  nlinha <- nlinha + 1
  mat[nlinha, 1:n.var1] <- 1

  #RESTRICTION 3
  #x_{j} + G2_{i} ≤ 1 #para a_i_j = 1
  for (i in 1:m) {
    for (j in 1:n) {
      if (table[i,j] == 1) {
        nlinha <- nlinha + 1
        mat[nlinha, c(j, col_ini_G2 + i)] <- c(1, 1)
      }
    }
  }

  #RESTRICTION 4
  #x_{j} - G1_{i} ≥ 0 #para a_i_j = 1
  for (i in 1:m) {
    for (j in 1:n) {
      if (table[i,j] == 1) {
        nlinha <- nlinha + 1
        mat[nlinha, c(j, col_ini_G1 + i)] <- c(1, -1)
      }
    }
  }

  #RESTRICTION 5
  #Y1_{i} - Z1_{ij} >= 0
  for (i in 1:m) {
    for (j in 1:n) {
      nlinha <- nlinha + 1
      mat[nlinha, c(
        col_ini_Y1 + i,
        col_ini_Z1 + j + n * (i - 1)
      )] <- c(1, -1)
    }
  }

  #RESTRICTION 6
  #x_{j} - Z1_{ij} >= 0
  for (i in 1:m) {
    for (j in 1:n) {
      nlinha <- nlinha + 1
      mat[nlinha, c(
        j,
        col_ini_Z1 + j + n * (i - 1)
      )] <- c(1, -1)
    }
  }

  #RESTRICTION 8
  #x_{j} + Y1_{i} - Z1_{ij} <= 1
  for (i in 1:m) {
    for (j in 1:n) {
      nlinha <- nlinha + 1
      mat[nlinha, c(
        j,
        col_ini_Y1 + i,
        col_ini_Z1 + j + n * (i - 1)
      )] <- c(1, 1, -1)
    }
  }

  #RESTRICTION 9
  #alpha_{i}*G1_{i} - sum_{j}(Z1_{i,j}) = 0
  for (i in 1:m) {
    nlinha <- nlinha + 1
    mat[nlinha, col_ini_G1 + i] <- alphai[i]
    mat[nlinha, (col_ini_Z1 + 1 + n * (i - 1)):(col_ini_Z1 + n + n * (i - 1))] <- -1

  }

  #RESTRICTION 10
  #Y2_{i} - Z2_{ij} >= 0
  for (i in 1:m) {
    for (j in 1:n) {
      nlinha <- nlinha + 1
      mat[nlinha, c(
        col_ini_Y2 + i,
        col_ini_Z2 + j + n * (i - 1)
      )] <- c(1, -1)
    }
  }

  #RESTRICTION 11
  #x_{j} - Z2_{ij} >= 0
  #x_{j} + Z2_{ij} <= 1
  for (i in 1:m) {
    for (j in 1:n) {
      nlinha <- nlinha + 1
      mat[nlinha, c(
        j,
        col_ini_Z2 + j + n * (i - 1)
      )] <- c(1, -1)
      #)] <- c(1, 1)
    }
  }

  #RESTRICTION 13
  #x_{j} + Y2_{i} - Z2_{ij} <= 1
  for (i in 1:m) {
    for (j in 1:n) {
      nlinha <- nlinha + 1
      mat[nlinha, c(
        j,
        col_ini_Y2 + i,
        col_ini_Z2 + j + n * (i - 1)
      )] <- c(1, 1, -1)
    }
  }

  #RESTRICTION 14
  #alpha_{i}*G2_{i} - n * Y2_{i} + sum_{j}(Z2_{i,j}) = 0
  for (i in 1:m) {
    nlinha <- nlinha + 1
    mat[nlinha, c(col_ini_G2 + i, col_ini_Y2 + i)] <- c(alphai[i], -n)
    mat[nlinha, (col_ini_Z2 + 1 + n * (i - 1)):(col_ini_Z2 + n + n * (i - 1))] <- 1
  }

  #right hand side vector
  rhs <- rep(0, total.rest.rows)
  rhs[n.restr1] <- 1
  rhs[n.restr1 + n.restr2] <- n - 1 #TO DO: experimentar com floor(n/2)
  rhs[n.restr1 + n.restr2 + (1:n.restr3)] <- 1
  rhs[n.restr1 + n.restr2 + n.restr3 + (1:n.restr4)] <- 0
  rhs[n.restr1 + n.restr2 + n.restr3 + n.restr4 + (1:n.restr5)] <- 0
  rhs[n.restr1 + n.restr2 + n.restr3 + n.restr4 + n.restr5 + (1:n.restr6)] <- 0
  rhs[n.restr1 + n.restr2 + n.restr3 + n.restr4 + n.restr5 + n.restr6 + (1:n.restr8)] <- 1
  rhs[n.restr1 + n.restr2 + n.restr3 + n.restr4 + n.restr5 + n.restr6 + n.restr8 + (1:n.restr9)] <- 0
  rhs[n.restr1 + n.restr2 + n.restr3 + n.restr4 + n.restr5 + n.restr6 + n.restr8 + n.restr9 + (1:n.restr10)] <- 0
  rhs[n.restr1 + n.restr2 + n.restr3 + n.restr4 + n.restr5 + n.restr6 + n.restr8 + n.restr9 + n.restr10 + (1:n.restr11)] <- 0
  rhs[n.restr1 + n.restr2 + n.restr3 + n.restr4 + n.restr5 + n.restr6 + n.restr8 + n.restr9 + n.restr10 + n.restr11 + (1:n.restr13)] <- 1
  rhs[n.restr1 + n.restr2 + n.restr3 + n.restr4 + n.restr5 + n.restr6 + n.restr8 + n.restr9 + n.restr10 + n.restr11 + n.restr13 + (1:n.restr14)] <- 0

  #directions vector/sense
  sense <- c(
    rep(">=", n.restr1),
    rep("<=", n.restr2),
    rep("<=", n.restr3),
    rep(">=", n.restr4 + n.restr5 + n.restr6),
    rep("<=", n.restr8),
    rep("=", n.restr9),
    rep(">=", n.restr10 + n.restr11),
    rep("<=", n.restr13),
    rep("=", n.restr14)
  )

  #objective function vector
  obj <- c(rep(0, n.var1 + n.var2 + n.var3), rep(1, n.var4 + n.var5), rep(0, n.var6 + n.var7))

  #objective function aim
  modelsense <- "max"

  #variable types
  types <- c(rep("B", n.var1), rep("C", n.var2 + n.var3), rep("C", n.var4 + n.var5), rep("C", n.var6 + n.var7))

  return(list(
    A=mat,
    obj=obj,
    modelsense=modelsense,
    rhs=rhs,
    sense=sense,
    vtype=types
  ))
}
