# GUROBI_k2_optim_tdv.R
#'
#' @title TotDiffVal (or TotDiffVal1) optimization using GUROBI
#'
#' @description This function finds the partition of the columns of a given matrix that maximizes the TotDiffVal (or TotDiffVal1) index.
#'
#' @param m A `matrix`, i.e. a phytosociological table of 0s (absences) and 1s (presences), where rows correspond to taxa and columns correspond to relevés.
#' @param formulation A `character`,  of integer numbers with the initial partition of the relevés (i.e. a vector with values from 1 to k, with length equal to the number of columns of m, ascribing each relevé to one of the k groups). By default, "random", generates a random initial partition.
#' @param index A `character`, selecting which index to optimize "TotDiffVal", or "TotDiffVal1". #STILL TODO
#' @param TimeLimit A `numeric` ("double") with the time limit (in seconds) to be passed as a parameter to GUROBI, Defaults to `Inf`, but see Details.
#'
#' @details Given a phytosociological table `m` (rows corresponding to taxa and columns corresponding to relevés) this function finds a 2-partition (a partition in two groups) that maximizes TotDiffVal (or TotDiffVal1) index (see http://home.isa.utl.pt/~tmh/), i.e. finds, using the GUROBI otpimizer (see \code{\link[gurobi]{gurobi}}). This partition is a global maximum of TotDiffVal (or TotDiffVal1) for any 2-partitions of the dataset.
#'
#' For medium-sized matrices the computation time might became prohibitive, thus the use of a time limit (`TimeLimit`) is very advisable.
#'
#' @return A `list` with the following components:
#'
#' \describe{
#'   \item{par}{A `vector` with the 2-partition that maximizes TotDiffVal (or TotDiffVal1). A global maximum.}
#'   \item{objval}{A `numeric` with the maximum TotDiffVal (or TotDiffVal1) found.}
#' }
#'
#' @author Jorge Orestes Cerdeira and Tiago Monteiro-Henriques. E-mail: \email{tiagomonteirohenriques@@gmail.com}.
#'
#' @export
#'
####
#### GUROBI_k2_optim_tdv
####

GUROBI_k2_optim_tdv <- function(m, formulation = c("t-independent", "t-dependent"), index = "TotDiffVal1", TimeLimit = Inf) {
  n <- ncol(m)
  ns <- nrow(m)
  alphai <- rowSums(m) #number of ones in each row
  if ("t-independent" %in% formulation) {
    params <- list(OutputFlag = 0, TimeLimit = TimeLimit)
    LISTgurobi <- optim_tdv_gurobi_ti(table=m, n=n, alphai=alphai)
    RES <- gurobi::gurobi(LISTgurobi, params)
    return(list(status = RES$status, par = RES$x[1:n] + 1, objval = RES$objval/ns))
  }
  if ("t-dependent" %in% formulation) {
    res_objval_1 <- NULL
    res_par_1 <- list()
    res_status_1 <- NULL
    params <- list(OutputFlag=0, TimeLimit = TimeLimit)
    t <- 1
    LISTgurobi <- optim_tdv_gurobi_td(table=m, t=t, n=n, alphai=alphai)
    RES <- gurobi::gurobi(LISTgurobi, params)
    res_objval_1 <- c(res_objval_1, RES$objval/ns)
    res_par_1[[t]] <- RES$x[1:n] + 1
    res_status_1 <- c(res_status_1, RES$status)
    if (floor(n/2) > 1) {
      for (t in 2:floor(n/2)) {
        LISTgurobi$rhs[1] <- t
        LISTgurobi$obj <- c(rep(0, n), alphai/(n-t), alphai/t)
        RES <- gurobi::gurobi(LISTgurobi, params)
        res_objval_1 <- c(res_objval_1, RES$objval/ns)
        res_par_1[[t]] <- RES$x[1:n]+1
        res_status_1 <- c(res_status_1, RES$status)
      }
    }
    #return(list(objval_1 = res_objval_1, par_1 = res_par_1))
    max.sol <- which.max(res_objval_1)
    return(list(status.all = res_status_1, status.max = res_status_1[[max.sol]], par = res_par_1[[max.sol]], objval = res_objval_1[[max.sol]]))
  }
  stop("In GUROBI_k2_optim_tdv, formulation must be 't-independent' or 't-dependent'.")
}

####
#### optim_tdv_gurobi_td
####

#Function to prepare all necessary objects to pass to Gurobi (t dependent) solver from a binary table (e.g. a binary phytosociological table)

#table 		a binary matrix
#t			 	the number of groups

optim_tdv_gurobi_td <- function (table, t, n, alphai) {

  m <- nrow(table) #i index
  #n <- ncol(table) #j index #this is now given as a function parameter

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

####
#### optim_tdv_gurobi_ti
####

#Function to prepare all necessary objects to pass to Gurobi (t independent) solver from a binary table (e.g. a binary phytosociological table)

#table 			a binary matrix

optim_tdv_gurobi_ti <- function (table, n, alphai) {
  m <- nrow(table) #i index
  #n <- ncol(table) #j index #this is now given as a function parameter

  #rows of restrictions matrix
  n.restr1 <- 1					#sum x_j >= 1
  n.restr2 <- n.restr1			#sum x_j <= n-1
  n.restr3 <- sum(table)		#x_{j} + G1_{i} ≤ 1 #for a_i_j = 1
  n.restr4 <- n.restr3			#x_{j} - G2_{i} ≥ 0 #for a_i_j = 1
  n.restr5 <- n*m				#Y1_{i} - Z1_{ij} >= 0
  n.restr6 <- n.restr5			#x_{j} - Z1_{ij} >= 0
  #n.restr7 <- n.restr5		#-1/(n-1)x_{j} + Z1_{ij} >= 0
  n.restr8 <- n.restr5			#x_{j} + Y1_{i} - Z1_{ij} <= 1
  n.restr9 <- m					#alpha_{i}*G1_{i} - sum_{j}(Z1_{i,j}) = 0
  n.restr10 <- n.restr5		#Y2_{i} - Z2_{ij} >= 0
  n.restr11 <- n.restr5		#x_{j} - Z2_{ij} >= 0
  #n.restr11 <- n.restr5		#x_{j} + Z2_{ij} <= 1
  #n.restr12 <- n.restr5		#-1/(n-1)x_{j} + Z2_{ij} >= 0
  n.restr13 <- n.restr5		#x_{j} + Y2_{i} - Z2_{ij} <= 1
  n.restr14 <- m					#alpha_{i}*G2_{i} - n * Y2_{i} + sum_{j}(Z2_{i,j}) = 0

  #total.rest.rows <- 2 * n.restr1 + 2 * n.restr3 + 8 * n.restr5 + 2 * m
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

  # #RESTRICTION 7
  # #-1/(n-1)x_{j} + Z1_{ij} >= 0
  # coef_r7 <- -1/(n-1)
  # for (i in 1:m) {
  # for (j in 1:n) {
  # nlinha <- nlinha + 1
  # mat[nlinha, c(
  # j,
  # col_ini_Z1 + j + n * (i - 1)
  # )] <- c(coef_r7, 1)
  # }
  # }

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

  # #RESTRICTION 12
  # #-1/(n-1)x_{j} + Z2_{ij} >= 0
  # coef_r12 <- -1/(n-1)
  # for (i in 1:m) {
  # for (j in 1:n) {
  # nlinha <- nlinha + 1
  # mat[nlinha, c(
  # j,
  # col_ini_Z2 + j + n * (i - 1)
  # )] <- c(coef_r12, 1)
  # }
  # }

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


  # #rows of restrictions matrix
  # n.restr1 <- 1			#sum x_j >= 1
  # n.restr2 <- n.restr1		#sum x_j <= n-1
  # n.restr3 <- sum(table)	#x_{j} + G1_{i} ≤ 1 #for a_i_j = 1
  # n.restr4 <- n.restr3		#x_{j} - G2_{i} ≥ 0 #for a_i_j = 1
  # n.restr5 <- n*m			#Y1_{i} - Z1_{ij} >= 0
  # n.restr6 <- n.restr5		#x_{j} - Z1_{ij} >= 0
  # n.restr7 <- n.restr5		#-1/(n-1)x_{j} + Z1_{ij} >= 0
  # n.restr8 <- n.restr5		#x_{j} + Y1_{i} - Z1_{ij} <= 1
  # n.restr9 <- m			#alpha_{i}*G1_{i} - sum_{j}(Z1_{i,j}) = 0
  # n.restr10 <- n.restr5		#Y2_{i} - Z2_{ij} >= 0
  # n.restr11 <- n.restr5		#x_{j} - Z2_{ij} >= 0
  # n.restr12 <- n.restr5		#-1/(n-1)x_{j} + Z2_{ij} >= 0
  # n.restr13 <- n.restr5		#x_{j} + Y2_{i} - Z2_{ij} <= 1
  # n.restr14 <- m			#alpha_{i}*G2_{i} - n * Y2_{i} + sum_{j}(Z2_{i,j}) = 0


  # #right hand side vector
  # rhs <- rep(0, total.rest.rows)
  # rhs[n.restr1] <- 1
  # rhs[n.restr1 + n.restr2] <- n - 1
  # rhs[n.restr1 + n.restr2 + (1:n.restr3)] <- 1
  # rhs[n.restr1 + n.restr2 + n.restr3 + (1:n.restr4)] <- 0
  # rhs[n.restr1 + n.restr2 + n.restr3 + n.restr4 + (1:n.restr5)] <- 0
  # rhs[n.restr1 + n.restr2 + n.restr3 + n.restr4 + n.restr5 + (1:n.restr6)] <- 0
  # rhs[n.restr1 + n.restr2 + n.restr3 + n.restr4 + n.restr5 + n.restr6 + (1:n.restr7)] <- 0
  # rhs[n.restr1 + n.restr2 + n.restr3 + n.restr4 + n.restr5 + n.restr6 + n.restr7 + (1:n.restr8)] <- 1
  # rhs[n.restr1 + n.restr2 + n.restr3 + n.restr4 + n.restr5 + n.restr6 + n.restr7 + n.restr8 + (1:n.restr9)] <- 0
  # rhs[n.restr1 + n.restr2 + n.restr3 + n.restr4 + n.restr5 + n.restr6 + n.restr7 + n.restr8 + n.restr9 + (1:n.restr10)] <- 0
  # rhs[n.restr1 + n.restr2 + n.restr3 + n.restr4 + n.restr5 + n.restr6 + n.restr7 + n.restr8 + n.restr9 + n.restr10 + (1:n.restr11)] <- 0
  # rhs[n.restr1 + n.restr2 + n.restr3 + n.restr4 + n.restr5 + n.restr6 + n.restr7 + n.restr8 + n.restr9 + n.restr10 + n.restr11 + (1:n.restr12)] <- 0
  # rhs[n.restr1 + n.restr2 + n.restr3 + n.restr4 + n.restr5 + n.restr6 + n.restr7 + n.restr8 + n.restr9 + n.restr10 + n.restr11 + n.restr12 + (1:n.restr13)] <- 1
  # rhs[n.restr1 + n.restr2 + n.restr3 + n.restr4 + n.restr5 + n.restr6 + n.restr7 + n.restr8 + n.restr9 + n.restr10 + n.restr11 + n.restr12 + n.restr13 + (1:n.restr14)] <- 0

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
  #rhs[n.restr1 + n.restr2 + n.restr3 + n.restr4 + n.restr5 + n.restr6 + n.restr8 + n.restr9 + n.restr10 + (1:n.restr11)] <- 1
  rhs[n.restr1 + n.restr2 + n.restr3 + n.restr4 + n.restr5 + n.restr6 + n.restr8 + n.restr9 + n.restr10 + n.restr11 + (1:n.restr13)] <- 1
  rhs[n.restr1 + n.restr2 + n.restr3 + n.restr4 + n.restr5 + n.restr6 + n.restr8 + n.restr9 + n.restr10 + n.restr11 + n.restr13 + (1:n.restr14)] <- 0



  # #directions vector/sense
  # sense <- c(
  # rep(">=", n.restr1),
  # rep("<=", n.restr2),
  # rep("<=", n.restr3),
  # rep(">=", n.restr4 + n.restr5 + n.restr6 + n.restr7),
  # rep("<=", n.restr8),
  # rep("=", n.restr9),
  # rep(">=", n.restr10 + n.restr11 + n.restr12),
  # rep("<=", n.restr13),
  # rep("=", n.restr14)
  # )

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
