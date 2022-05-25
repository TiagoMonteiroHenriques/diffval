#SOME AUXILIARY FUNCTIONS

#Auxiliary function to HillClimb_optim_tdv
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

#Auxiliary function to HillClimb_optim_tdv
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
    return(tdv.neig(mt = mt, tp = tp, k = k, ns = ns, nr = nr, ofda = ofda, ifp =ifp, afg = afg, empty.size = empty.size, gct = gct, i.mul = i.mul, DV = DV, pn = pn, kc = kc))
  })
  return(list(tdv = tdv <- max(mat.neig.tdv), p = mat.neig$p.list[which(mat.neig.tdv == tdv)[1],]))
}

#Auxiliary function for SimulAnne_optim_tdv
#random.neighbour.SA: selects a random neighbour (or the same partition), changing just one relevé

random.neighbour.SA <- function (p, nr, k) {
  tp <- table(p)
  k.int <- which(tp > 1) #Check behaviour when all groups have only one element...
  change.col <- sample((1:nr)[p %in% k.int], 1) #Changing just one column
  #tp[p[change.col]] <- tp[p[change.col]] - 1 #Changing just one column
  p[change.col] <- sample(1:k, length(change.col), replace=TRUE) #CF replace = TRUE (added later)
  return(p)
}


#Auxiliary functions to decrease GRASP and GRDTP computation time
#These are three auxiliary functions to recalculate parameter "c" (outside group g) given a new column relevé to be included in group g
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
