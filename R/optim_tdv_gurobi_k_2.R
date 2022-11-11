#' Total Differential Value optimization using Gurobi
#'
#' Given a phytosociological matrix, this function finds the partition of its
#'   columns that maximizes the Total Differential Value (TDV).
#'
#' @param m_bin A matrix. A phytosociological table of 0s (absences) and 1s
#'   (presences), where rows correspond to taxa and columns correspond to
#'   relevés.
#' @param formulation A character selecting which formulation to use. Possible
#'   values are "t-dependent" (the default) or "t-independent". See Details.
#' @param time_limit A numeric ("double") with the time limit (in seconds) to
#'   be passed as a parameter to Gurobi, Defaults to 5 seconds, but see Details.
#'
#' @details Given a phytosociological table `m_bin` (rows corresponding to taxa
#'   and columns corresponding to relevés) this function finds a 2-partition (a
#'   partition in two groups) that maximizes TDV, using the Gurobi optimizer.
#'
#'   [Gurobi](https://www.gurobi.com/) is a commercial software for which a free
#'   academic license can be obtained if you are affiliated with a recognized
#'   educational institution. Package 'prioritizr' contains a comprehensive
#'   vignette ([Gurobi Installation Guide](https://prioritizr.net/articles/gurobi_installation_guide.html)),
#'   which can guide you trough the process of obtaining a license, installing
#'   the [Gurobi optimizer](https://www.gurobi.com/products/gurobi-optimizer/),
#'   activating the license and eventually installing the R package 'gurobi'.
#'
#'   [optim_tdv_gurobi_k_2()] returns, when the optimization is successful, a
#'   2-partition which is a global maximum of TDV for any 2-partitions of the
#'   columns on `m_bin`.
#'
#'   See [tdv()] for an explanation on the Total Differential Value of a
#'   phytosociological table.
#'
#'   The function implements two different mixed-integer linear programming
#'   formulations of the problem. The formulations differ as one is independent
#'   of the size of the obtained groups (t-independent), while the other
#'   formulation fixes the size of the obtained groups (t-dependent). The
#'   t-dependent formulation is implemented to run Gurobi as many times as
#'   necessary to cover all possible group sizes; this approach can result in
#'   faster total computation time.
#'
#'   For medium-sized matrices the computation time might become already
#'   prohibitive, thus the use of a time limit (`time_limit`) is advisable.
#'
#' @return For `formulation = "t-dependent"`, a list with the following
#'   components:
#'
#'   \describe{
#'     \item{status.runs}{A character vector with Gurobi output status for all
#'     the runs.}
#'     \item{objval}{A numeric with the maximum TDV found by Gurobi.}
#'     \item{par}{A vector with the 2-partition corresponding to the the
#'     maximum TDV found by Gurobi.}
#'   }
#'
#'   For `formulation = "t-independent"`, a list with the following components:
#'
#'   \describe{
#'     \item{status}{A character with Gurobi output status.}
#'     \item{objval}{A numeric with the maximum TDV found by Gurobi.}
#'     \item{par}{A vector with the 2-partition corresponding to the the
#'     maximum TDV found by Gurobi.}
#'   }
#'
#' @author Jorge Orestes Cerdeira and Tiago Monteiro-Henriques.
#'   E-mail: \email{tmh.dev@@icloud.com}.
#'
#' @examples
#' # Getting the Taxus baccata forests data set
#' data(taxus_bin)
#'
#' # Obtaining the 2-partition that maximizes TDV using the Gurobi solver, by
#' # mixed-integer linear programming
#' \dontrun{
#' # Requires the suggested package 'gurobi'
#' optim_tdv_gurobi_k_2(taxus_bin)
#' }
#'
#' @export
optim_tdv_gurobi_k_2 <- function(m_bin,
                                 formulation = "t-dependent",
                                 time_limit = 5) {
  if (!requireNamespace("gurobi", quietly = TRUE)) {
    stop(
      "Package 'gurobi' must be installed to use this function.",
      call. = FALSE
    )
  }
  n <- ncol(m_bin)
  ns <- nrow(m_bin)
  alphai <- rowSums(m_bin) # Number of ones in each row
  if ("t-independent" %in% formulation) {
    params <- list(OutputFlag = 0, TimeLimit = time_limit)
    list_gurobi <- optim_tdv_gurobi_ti(table = m_bin, n = n, alphai = alphai)
    res_gurobi <- gurobi::gurobi(list_gurobi, params)
    return(list(
      status = res_gurobi$status,
      par    = res_gurobi$x[1:n] + 1,
      objval = res_gurobi$objval / ns
    ))
  }
  if ("t-dependent" %in% formulation) {
    res_objval_1 <- NULL
    res_par_1 <- list()
    res_status_1 <- NULL
    params <- list(OutputFlag = 0, TimeLimit = time_limit)
    t <- 1
    list_gurobi <- optim_tdv_gurobi_td(
      table = m_bin,
      t = t,
      n = n,
      alphai = alphai
    )
    res_gurobi <- gurobi::gurobi(list_gurobi, params)
    res_objval_1 <- c(res_objval_1, res_gurobi$objval / ns)
    res_par_1[[t]] <- res_gurobi$x[1:n] + 1
    res_status_1 <- c(res_status_1, res_gurobi$status)
    if (floor(n / 2) > 1) {
      for (t in 2:floor(n / 2)) {
        list_gurobi$rhs[1] <- t
        list_gurobi$obj <- c(rep(0, n), alphai / (n - t), alphai / t)
        res_gurobi <- gurobi::gurobi(list_gurobi, params)
        res_objval_1 <- c(res_objval_1, res_gurobi$objval / ns)
        res_par_1[[t]] <- res_gurobi$x[1:n] + 1
        res_status_1 <- c(res_status_1, res_gurobi$status)
      }
    }
    max_sol <- which.max(res_objval_1)
    return(list(
      status.runs = res_status_1,
      par         = res_par_1[[max_sol]],
      objval      = res_objval_1[[max_sol]]
    ))
  }
  stop('In optim_tdv_gurobi_k_2, formulation must be "t-independent" or
       "t-dependent".')
}
