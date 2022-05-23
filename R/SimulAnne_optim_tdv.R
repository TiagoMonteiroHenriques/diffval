# SimulAnne_optim_tdv.R
#'
#' @title Total Differential Value optimization using a simulated annealing (and GRASP) algorithm(s)
#'
#' @description This function searches for `k`-partitions of the columns of a given matrix (i.e. a partition of the columns in `k` groups), optimizing
#' the Total Differential Value (TDV) using a stochastic global optimization method called simulated annealing (SANN) algorithm. Optionally, a Greedy
#' Randomized Adaptive Search Procedure (GRASP) can be used to find a initial partition (seed) to be passed to the SANN algorithm.
#'
#' @param m A `matrix`, i.e. a phytosociological table of 0s (absences) and 1s (presences), where rows correspond to taxa and columns correspond to relevés.
#' @param k A `numeric` giving the number of desired groups.
#' @param p.initial A `vector` of integer numbers with the partition of the relevés (i.e. a `k`-partition, consisting in a vector with values from 1 to k, with length equal to the number of columns of m, ascribing each relevé to one of the k groups), to be used as initial partition in the simulated annealing. For a random partition use `p.initial = "random"`. This argument is ignored if `use.GRASP = TRUE`.
#' @param full.output A `logical`. If `FALSE` (the default) the best `n.sol` partitions and respective indices are returned. If `TRUE` the output will contain, additionally, data on the optimization steps, for all runs.
#' @param n.runs A `numeric` giving the number of runs. Defaults to 10.
#' @param n.sol A `numeric`, giving the number of best solutions to keep in the final output (only used if `full.output` is `FALSE`; if `full.output` is `TRUE` all runs will produce an output). Defaults to 5.
#' @param T_inic A `numeric` giving the initial temperature. Must be greater than 0 and maximum admitted value is 1. Defaults to 0.3.
#' @param T_final A `numeric` giving the initial temperature. Must be bounded between 0 and 1. Usually very low values are needed to ensure convergence. Defaults to 0.000001.
#' @param alpha A `numeric` giving the fraction of temperature drop to be used in the temperature reduction scheme (see Details). Must be bounded between 0 and 1. Defaults to 0.05.
#' @param n.iter A `numeric` giving the initial temperature. Defaults to 1000
#' @param use.GRASP A `logical`. Defaults to `TRUE`. IF `TRUE`, a GRASP is used to obtain the initial partitions for the simulated annealing. If `FALSE` the user should provide an initial partition or use or use `p.initial = "random"` for a random one.
#' @param thr A `numeric` giving a threshold value (from 0 to 1 ) with the probability used to compute the sample quantile, in order to get the best `m` columns from which to select one to be include in the GRASP solution (in each step of the procedure). Only needed if `use.GRASP` is `TRUE`.
#' @param full.output A `logical`. Defaults to `FALSE`. IF `TRUE` extra information is presented in the output. See Value.
#'
#' @details Given a phytosociological table (`m`, with rows corresponding to taxa and columns corresponding to relevés) this function searches
#' for a k-partition (`k`, defined by the user) optimizing TDV, i.e. searches, using a SANN algorithm (optionally working upon GRASP solutions), for
#' the global maximum of TDV (by rearranging the relevés into `k` groups).
#'
#' This function uses two main algorithms:
#'
#' 1) An optional GRASP, which is used to obtain initial solutions (partitions of `m`) using function \code{\link{GRASP_partition_tdv}}. Such initial
#' solutions are then submitted to the SANN algorithm.
#' 2) The (main) SANN algorithm, which is used to search for the global maximum of TDV. The initial partition for each run of SANN can be a partition
#' obtained from GRASP (if `use.GRASP = TRUE`) or, (if `use.GRASP = FALSE`), a partition given by the user (using `p.initial`) or a random partition
#' (using `p.initial = "random"`).
#'
#' The SANN algorithm decreases the temperature multiplying the current temperature by `1 - alpha` according to a predefined schedule, which is
#' automatically calculated from the given values for `T_inic`, `T_final`, `alpha` and `n.iter`.
#' Specifically, the cooling schedule is obtained calculating the number of times that the temperature has to be decreased in order to
#' approximate `T_final` starting from `T_inic`. The number of times that the temperature decreases, say `nt`, is calculated by the
#' expression `floor(n.iter/((n.iter * log(1 - alpha)) / (log((1 - alpha) * T_final / T_inic))))`.
#' Finally, these decreasing stages are scattered through the desired iterations (`n.iter`) homogeneously, by calculating the indices of the
#' interactions that will experience a decreasing temperature, using the expression `floor(n.iter / nt * (1:nt))`.
#'
#' SANN is often seen as an exploratory technique where the temperature settings are challenging and dependent on the problem. This function
#' tries to restrict temperature values taking into account that TDV is always between 0 and 1.
#' Even though, obtaining values of temperature that allow convergence can be challenging. `full.output = TRUE` allows the user to inspect the
#' behaviour of `current.tdv` and check if convergence fails.
#' Generally, convergence failure can be spotted when final SANN TDV values are similar to the initial `current.tdv`, specially when coming
#' from random partitions.
#' In such cases, as a rule of thumb, it is advisable to decrease `T_final`.
#'
#' @return If `full.output = FALSE` (the default), a `list` with the following components (the GRASP component is only returned if `use.GRASP = TRUE`):
#'
#' \describe{
#'   \item{GRASP}{A `list` with at most `n.sol` components, each one containing also a `list` with two components:
#'   \itemize{
#'   \item{'par', containing the partition obtained by GRASP;}
#'   \item{'tdv', with the value of TDV of the above partition.}
#'   }
#'   }
#'   \item{SANN}{A `list` with at most `n.sol` components, each one containing also a `list` with two components:
#'   \itemize{
#'   \item{'par', containing the partition with highest value of TDV, obtained by the SANN algorithm;}
#'   \item{'tdv', with the value of TDV of the above partition.}
#'   }
#'   }
#'}
#'
#' If `full.output = TRUE`, a `list` with the following components (the GRASP component is only returned if `use.GRASP = TRUE`):
#'
#' \describe{
#'   \item{GRASP}{A `list` with `n.runs` components, each one containing also a `list` with two components:
#'   \itemize{
#'   \item{'par', containing the partition obtained by GRASP;}
#'   \item{'tdv', with the value of TDV of the above partition.}
#'   }
#'   }
#'   \item{SANN}{A `list` with `n.runs` components, each one containing also a `list` with six components:
#'   \itemize{
#'   \item{'par', containing the partition with highest value of TDV, obtained by the SANN algorithm;}
#'   \item{'tdv', with the value of TDV of the above partition;}
#'   \item{'current.tdv', a `vector` of length `n.iter` with the current TDV of each SANN iteration;}
#'   \item{'alternative.tdv', a `vector` of length `n.iter` with the alternative TDV used in each SANN iteration;}
#'   \item{'probability', a `vector` of length `n.iter` with the probability used in each SANN iteration;}
#'   \item{'temperature', a `vector` of length `n.iter` with the temperature of each SANN iteration.}
#'   }
#'   }
#' }
#'
#' @author Jorge Orestes Cerdeira and Tiago Monteiro-Henriques. E-mail: \email{tiagomonteirohenriques@@gmail.com}.
#'
#' @export
#'
SimulAnne_optim_tdv <- function(m, k, p.initial = NULL, n.runs = 10, n.sol = 5, T_inic = 0.3, T_final = 0.000001, alpha = 0.05, n.iter = 1000, use.GRASP = TRUE, thr = 0.95, full.output = FALSE) {
  stopifnot(is.matrix(m))
  if (alpha >= 1 | alpha <= 0) {stop("Please note that 0 < alpha < 1 is mandatory.")}
  if (T_inic > 1 | T_inic <= 0) {stop("Please note that 0 < T_inic <= 1 is mandatory.")}
  if (T_final >= 1 | T_final <= 0) {stop("Please note that 0 < T_final < 1 is mandatory.")}
  if (T_final >= T_inic) {stop("Please note that T_final must be lower than T_inic.")}
  if (T_inic*(1 - alpha)^n.iter > T_final) {stop("Desired T_final is not achieved with given alpha and n.iter (alpha and/or n.iter are too small).")}
  mode(m) <- "integer"
  if (!identical(c(0L,1L), sort(unique(as.vector(m))))) {stop("Matrix m should contain only 0's and 1's")}
  if (min(rowSums(m))==0) {stop("At least one taxa is not present in any relev\u00e9")}
  if (min(colSums(m))==0) {stop("At least one relev\u00e9 contains no taxa")}
  nr <- ncol(m) # no. of relevés
  ns <- nrow(m) # no. of taxa
  if (k > nr) {stop("Given k size is too big.")}
  if (k <= 1) {stop("Given k size is too small.")}

  if (use.GRASP) {
    RES.GRASP <- list()
    RES.SANN <- list()
    for (i in 1:n.runs) {
      #GRASP initial partitions
      par.GRASP <- GRASP_partition_tdv(m, k, thr = 0.95, verify = FALSE)
      GRASP_result <- list(par = par.GRASP, tdv = tdv(m, par.GRASP, output.type = "fast"))
      #SANN step
      SANN_result <- SimulAnne_optim_tdv(m = m, k = k, p.initial = par.GRASP, n.runs = 1, n.sol = 1, T_inic = T_inic, T_final = T_final, alpha = alpha, n.iter = n.iter, use.GRASP = FALSE, full.output = full.output)$SANN[[1]]

      if (full.output) {
        RES.GRASP[[i]] <- GRASP_result
        RES.SANN[[i]] <- SANN_result
      }
      if (!full.output) {
        if (i == 1) {
          RES.GRASP[[1]] <- GRASP_result
          RES.SANN[[1]] <- SANN_result
        } else {
          already.in.bestsol <- any(sapply(RES.SANN, function (x) {
            equivalent_partition(x$par, SANN_result$par)
          }))
          if (!already.in.bestsol) {
            if (length(RES.SANN) < n.sol) {
              RES.GRASP[[length(RES.SANN) + 1]] <- GRASP_result
              RES.SANN[[length(RES.SANN) + 1]] <- SANN_result
            } else { #already n.sol components in RES.SANN
              best.sol.values <- sapply(RES.SANN, function (x) {
                x$tdv
              })
              if (SANN_result$tdv >  min(best.sol.values)) {
                worse.bestsol <- which.min(best.sol.values) #selects the first in case of ties!
                RES.GRASP[[worse.bestsol]] <- GRASP_result
                RES.SANN[[worse.bestsol]] <- SANN_result
              }
            }
          }
        }
      }
    }
    ind.order <- order(sapply(RES.SANN, function (x) {
      x$tdv
    }), decreasing = TRUE)
    return(list(GRASP = RES.GRASP[ind.order], SANN = RES.SANN[ind.order]))
  }
  #if use.GRASP == FALSE
  #simple SANN on a given initial partition

  if (is.null(p.initial)) {stop("Argument p.initial can not be NULL when use.GRASP = FALSE.")}
  if (p.initial[1] != "random") {
    if (!identical(length(p.initial), nr)) {stop("Object p.initial must be a partition of the columns of m")}
    mode(p.initial) <- "integer"
    if (!identical(sort(unique(p.initial)), 1:k)) {stop("Object p.initial is not a valid partition of the columns of m")}
    p.ini <- p.initial
  }

  #Cooling schedule (CooS)
  nt <- floor(n.iter/((n.iter * log(1 - alpha)) / (log((1 - alpha) * T_final / T_inic))))
  CooS <- floor(n.iter/nt*(1:nt))

  RES.SANN <- list()
  for (i in 1:n.runs) {
    if (p.initial[1] == "random") {
      p.ini <- sample(c(1:k, sample(k, nr-k, replace=TRUE)))
    }
    best.p <- cur.p <- p.ini
    best.tdv <- cur.tdv <- tdv.ini <- tdv(m, p.ini, output.type = "fast")
    Temp <- T_inic
    res.cur.tdv <- res.alt.tdv <- res.alt.tdv <- res.probabi <- res.tempera <- NULL

    for (iter in 1:n.iter) {
      alt.p <- random.neighbour.SA(p = cur.p, nr = nr, k = k)
      alt.tdv <- tdv(m, alt.p, output.type = "fast")
      res.alt.tdv <- c(res.alt.tdv, alt.tdv)
      res.cur.tdv <- c(res.cur.tdv, cur.tdv)
      if (alt.tdv >= cur.tdv) {
        cur.p <- alt.p
        cur.tdv <- alt.tdv
        if (cur.tdv > best.tdv) {
          best.p <- cur.p
          best.tdv <- cur.tdv
        }
        probabi <- 1
      } else {
        probabi <- exp((alt.tdv-cur.tdv)/Temp)
        if (stats::runif(1) < probabi) {
          cur.p <- alt.p
          cur.tdv <- alt.tdv
        }
      }
      if (iter %in% CooS) {Temp <- (1 - alpha) * Temp} #Cooling
      res.probabi <- c(res.probabi, probabi)
      res.tempera <- c(res.tempera, Temp)
    }
    if (full.output) {
      RES.SANN[[i]] <- list(par = best.p, tdv = best.tdv, current.tdv = res.cur.tdv, alternative.tdv = res.alt.tdv, probability = res.probabi, temperature = res.tempera)
    }
    if (!full.output) {
      if (i == 1) {
        RES.SANN[[1]] <- list(par = best.p, tdv = best.tdv)
      } else {
        already.in.bestsol <- any(sapply(RES.SANN, function (x) {equivalent_partition(x$par, best.p)}))
        if (!already.in.bestsol) {
          if (length(RES.SANN) < n.sol) {
            RES.SANN[[length(RES.SANN)+1]] <- list(par = best.p, tdv = best.tdv)
          } else { #already n.sol components in RES.SANN
            best.sol.values <- sapply(RES.SANN, function (x) {x$tdv})
            if (best.tdv >  min(best.sol.values)) {
              worse.bestsol <- which.min(best.sol.values) #selects the first in case of ties!
              RES.SANN[[worse.bestsol]] <- list(par = best.p, tdv = best.tdv)
            }
          }
        }
      }
    }
  }
  ind.order <- order(sapply(RES.SANN, function (x) {
    x$tdv
  }), decreasing = TRUE)
  return(list(SANN = RES.SANN[ind.order]))
}
