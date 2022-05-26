# optim_tdv_GUROBI_k2.R
#'
#' @title Total Differential Value optimization using GUROBI
#'
#' @description Given a phytosociological matrix, this function finds the partition of its columns that maximizes the Total Differential Value (TDV).
#'
#' @param m A `matrix`, i.e. a phytosociological table of 0s (absences) and 1s (presences), where rows correspond to taxa and columns correspond to relevés.
#' @param formulation A `character` selecting wich formulation to use. Possible values are "t-dependent" (the default) or "t-independent". See Details.
#' @param TimeLimit A `numeric` ("double") with the time limit (in seconds) to be passed as a parameter to GUROBI, Defaults to 5 seconds, but see Details.
#'
#' @details Given a phytosociological table `m` (rows corresponding to taxa and columns corresponding to relevés) this function finds a 2-partition (a
#' partition in two groups) that maximizes TDV, using the GUROBI otpimizer (see \code{\link[gurobi]{gurobi}}).
#' When successful, this partition is a global maximum of TDV for any 2-partitions of the columns on `m`.
#'
#' See \code{\link{tdv}} for an explanation on the Total Differential Value of a phytosociological table.
#'
#' The function implements two different mixed-integer linear programming formulations of the problem. The formulations differ as one is independent of the size
#' of the obtained groups (t-independent), while the other formulation fixes the size of the obtained groups (t-dependent). The t-dependent formulation is
#' implemented to run GUROBI as many times as necessary to cover all possible group sizes; this approach can result in faster total computation time.
#'
#' For medium-sized matrices the computation time might become already prohibitive, thus the use of a time limit (`TimeLimit`) is advisable.
#'
#' @return For `formulation = "t-dependent"`, a `list` with the following components:
#'
#' \describe{
#'   \item{status.runs}{A `character` with GUROBI output status for all the runs.}
#'   \item{objval}{A `numeric` with the maximum TDV found by GUROBI.}
#'   \item{par}{A `vector` with the 2-partition corresponding to the the maximum TDV found by GUROBI.}
#' }
#'
#' For `formulation = "t-independent"`, a `list` with the following components:
#'
#' \describe{
#'   \item{status}{A `character` with GUROBI output status.}
#'   \item{objval}{A `numeric` with the maximum TDV found by GUROBI.}
#'   \item{par}{A `vector` with the 2-partition corresponding to the the maximum TDV found by GUROBI.}
#' }
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
#' #obtaining the 2-partition that maximizes TDV using the GUROBI solver, by
#' #mixed-integer linear programming
#' \dontrun{
#' #requires the suggested package 'gurobi'
#' optim_tdv_GUROBI_k2(taxus_bin)
#' }
#'
#' @export
#'
optim_tdv_GUROBI_k2 <- function(m, formulation = "t-dependent", TimeLimit = 5) {
  if (!requireNamespace("gurobi", quietly = TRUE)) {
    stop(
      "Package \"gurobi\" must be installed to use this function.",
      call. = FALSE
    )
  }
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
    return(list(status.runs = res_status_1, par = res_par_1[[max.sol]], objval = res_objval_1[[max.sol]]))
  }
  stop('In optim_tdv_GUROBI_k2, formulation must be "t-independent" or "t-dependent".')
}
