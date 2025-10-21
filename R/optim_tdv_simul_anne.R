#' Total Differential Value optimization using a Simulated Annealing (and GRASP)
#'   algorithm(s)
#'
#' This function searches for `k`-partitions of the columns of a given matrix
#'   (i.e., partitions of the columns into `k` groups), optimizing the Total
#'   Differential Value (TDV) using a stochastic global optimization method
#'   known as the Simulated Annealing (SANN) algorithm. Optionally, a Greedy
#'   Randomized Adaptive Search Procedure (GRASP) can be used to find an initial
#'   partition (seed) to be passed to the SANN algorithm.
#'
#' @param m_bin A matrix. A phytosociological table of 0s (absences) and
#'   1s (presences), where rows correspond to taxa and columns correspond to
#'   relevés.
#' @param k A numeric giving the number of desired groups.
#' @param p_initial A vector of integer numbers with the partition of the
#'   relevés (i.e., a `k`-partition, consisting in a vector with values from 1
#'   to `k`, with length equal to the number of columns of `m_bin`, ascribing
#'   each relevé to one of the `k` groups), to be used as initial partition in
#'   the Simulated Annealing. For a random partition use `p_initial = "random"`.
#'   This argument is ignored if `use_grasp = TRUE`.
#' @param full_output A logical. If `FALSE` (the default) the best `n_sol`
#'   partitions and respective indices are returned. If `TRUE` the output will
#'   contain, additionally, data on the optimization steps, for all runs.
#' @param n_runs A numeric giving the number of runs. Defaults to 10.
#' @param n_sol A numeric giving the number of best solutions to keep in the
#'   final output (only used if `full_output` is `FALSE`; if `full_output` is
#'   `TRUE` all runs will produce an output). Defaults to 1.
#' @param t_inic A numeric giving the initial temperature. Must be greater
#'   than 0 and the maximum admitted value is 1. Defaults to 0.3.
#' @param t_final A numeric giving the final temperature. Must be bounded
#'   between 0 and 1. Usually very low values are needed to ensure convergence.
#'   Defaults to 0.000001.
#' @param alpha A numeric giving the fraction of temperature drop to be used
#'   in the temperature reduction scheme (see Details). Must be bounded between
#'   0 and 1. Defaults to 0.05.
#' @param n_iter A numeric giving the number of iterations. Defaults to 1000.
#' @param use_grasp A logical. Defaults to `TRUE`. IF `TRUE`, a GRASP is used
#'   to obtain the initial partitions for the Simulated Annealing. If `FALSE`
#'   the user should provide an initial partition or use or use
#'   `p_initial = "random"` for a random one.
#' @param thr A numeric giving a threshold value (from 0 to 1 ) with the
#'   probability used to compute the sample quantile, in order to get the best
#'   `m_bin` columns from which to select one to be include in the GRASP
#'   solution (in each step of the procedure). Only needed if `use_grasp` is
#'   `TRUE`.
#' @param full_output A logical. Defaults to `FALSE`. If `TRUE` extra
#'   information is presented in the output. See Value.
#'
#' @details Given a phytosociological table (`m_bin`, with rows corresponding to
#'   taxa and columns corresponding to relevés) this function searches for a
#'   `k`-partition (`k`, defined by the user) optimizing the TDV, i.e.,
#'   searches, using a SANN algorithm (optionally working upon GRASP solutions),
#'   for a global maximum of TDV (by rearranging the relevés into `k` groups).
#'
#'   In the terminology of cluster analysis, taxa correspond to features,
#'   variables, or attributes, while relevés correspond to objects or samples.
#'
#'   This function uses two main algorithms:
#'
#'   1) An optional GRASP, which is used to obtain initial solutions
#'   (partitions of `m_bin`) using function [partition_tdv_grasp()].
#'   Such initial solutions are then submitted to the SANN algorithm.
#'   2) The (main) SANN algorithm, which is used to search for a global
#'   maximum of TDV. The initial partition for each run of SANN can be a
#'   partition obtained from GRASP (if `use_grasp = TRUE`) or, (if
#'   `use_grasp = FALSE`), a partition given by the user (using `p_initial`) or
#'   a random partition (using `p_initial = "random"`).
#'
#'   The SANN algorithm decreases the temperature multiplying the current
#'   temperature by `1 - alpha` according to a predefined schedule, which is
#'   automatically calculated from the given values for `t_inic`, `t_final`,
#'   `alpha` and `n_iter`.
#'   Specifically, the cooling schedule is obtained calculating the number of
#'   times that the temperature has to be decreased in order to approximate
#'   `t_final` starting from `t_inic`. The number of times that the temperature
#'   decreases, say `nt`, is calculated by the expression:
#'
#'   `floor(log(t_final / t_inic) / log(1 - alpha))`.
#'
#'   Finally, these decreasing stages are scattered through the desired
#'   iterations (`n_iter`) homogeneously, by calculating the indices of the
#'   iterations that will experience a decrease in temperature using
#'   `floor(n_iter / nt * (1:nt))`.
#'
#'   SANN is often seen as an exploratory technique where the temperature
#'   settings are challenging and dependent on the problem. This function tries
#'   to restrict temperature values taking into account that TDV is always
#'   between 0 and 1. Even though, obtaining values of temperature that allow
#'   convergence can be challenging. `full_output = TRUE` allows the user to
#'   inspect the behaviour of `current.tdv` and check if convergence fails.
#'   Generally, convergence failure can be spotted when final SANN TDV values
#'   are similar to the initial `current.tdv`, specially when coming from random
#'   partitions. In such cases, as a rule of thumb, it is advisable to decrease
#'   `t_final` by a factor of 10.
#'
#' @return If `full_output = FALSE` (the default), a list with the following
#'   components (the GRASP component is only returned if `use_grasp = TRUE`):
#'
#'   \describe{
#'     \item{GRASP}{A list with at most `n_sol` components, each one
#'     containing also a list with two components:
#'     \describe{
#'       \item{par}{A vector with the partition of highest TDV obtained by
#'       GRASP;}
#'       \item{tdv}{A numeric with the TDV of `par`.}
#'       }
#'       }
#'     \item{SANN}{A list with at most `n_sol` components, each one containing
#'       also a list with two components:
#'     \describe{
#'       \item{par}{A vector with the partition of highest TDV obtained by the
#'       (GRASP +) SANN algorithm(s);}
#'       \item{tdv}{A numeric with the TDV of `par`.}
#'     }
#'     }
#'   }
#'
#'   If `full_output = TRUE`, a list with the following components (the GRASP
#'   component is only returned if `use_grasp = TRUE`):
#'
#'   \describe{
#'     \item{GRASP}{A list with `n_runs` components, each one containing also a
#'     list with two components:
#'     \describe{
#'       \item{par}{A vector with the partition of highest TDV obtained by
#'       GRASP.}
#'       \item{tdv}{A numeric with the TDV of `par`.}
#'     }
#'     }
#'     \item{SANN}{A list with `n_runs` components, each one containing also a
#'     list with six components:
#'     \describe{
#'       \item{current.tdv}{A vector of length `n_iter` with the current TDV of
#'       each SANN iteration.}
#'       \item{alternative.tdv}{A vector of length `n_iter` with the alternative
#'       TDV used in each SANN iteration.}
#'       \item{probability}{A vector of length `n_iter` with the probability
#'       used in each SANN iteration.}
#'       \item{temperature}{A vector of length `n_iter` with the temperature of
#'       each SANN iteration.}
#'       \item{par}{A vector with the partition of highest TDV obtained by the
#'       (GRASP +) SANN algorithm(s).}
#'       \item{tdv}{A numeric with the TDV of `par`.}
#'     }
#'     }
#'   }
#'
#' @author Jorge Orestes Cerdeira and Tiago Monteiro-Henriques.
#'   E-mail: \email{tmh.dev@@icloud.com}.
#'
#' @examples
#' # Getting the Taxus baccata forests data set
#' data(taxus_bin)
#'
#' # Removing taxa occurring in only one relevé in order to
#' # reproduce the example in the original article of the data set
#' taxus_bin_wmt <- taxus_bin[rowSums(taxus_bin) > 1, ]
#'
#' # Obtaining a partition that maximizes TDV using the Simulated Annealing
#' # algorithm
#' result <- optim_tdv_simul_anne(
#'   m_bin = taxus_bin_wmt,
#'   k = 3,
#'   p_initial = "random",
#'   n_runs = 5,
#'   n_sol = 5,
#'   use_grasp = FALSE,
#'   full_output = TRUE
#' )
#'
#' # Inspect the result
#' # The TDV of each run
#' sapply(result[["SANN"]], function(x) x$tdv)
#' # The best partition that was found (i.e., with highest TDV)
#' result[["SANN"]][[1]]$par
#'
#' # A TDV of 0.1958471 indicates you are probably reproducing the three
#' # groups (Estrela, Gerês and Galicia) from the original article. A solution
#' # with TDV = 0.2005789 might also occur, but note that one group has only two
#' # elements. For now, a minimum group size is not implemented in function
#' # optim_tdv_simul_anne() as it is in the function optim_tdv_hill_climb().
#'
#' # Inspect how the optimization progressed (should increase towards the right)
#' plot(
#'   result[["SANN"]][[1]]$current.tdv,
#'   type = "l",
#'   xlab = "Iteration number",
#'   ylab = "TDV of the currently accepted solution"
#' )
#' for (run in 2:length(result[["SANN"]])) {
#'   lines(result[["SANN"]][[run]]$current.tdv)
#' }
#'
#' # Plot the sorted (or tabulated) phytosociological table, using the best
#' # partition that was found
#' tabul <- tabulation(
#'   m_bin = taxus_bin_wmt,
#'   p = result[["SANN"]][[1]]$par,
#'   taxa_names = rownames(taxus_bin_wmt),
#'   plot_im = "normal"
#' )
#'
#' @export
optim_tdv_simul_anne <- function(m_bin,
                                 k,
                                 p_initial = NULL,
                                 n_runs = 10,
                                 n_sol = 1,
                                 t_inic = 0.3,
                                 t_final = 0.000001,
                                 alpha = 0.05,
                                 n_iter = 1000,
                                 use_grasp = TRUE,
                                 thr = 0.95,
                                 full_output = FALSE) {
  stopifnot(is.matrix(m_bin))
  if (alpha >= 1 || alpha <= 0) {
    stop("Please note that 0 < `alpha` < 1 is mandatory.")
  }
  if (t_inic > 1 || t_inic <= 0) {
    stop("Please note that 0 < `t_inic` <= 1 is mandatory.")
  }
  if (t_final >= 1 || t_final <= 0) {
    stop("Please note that 0 < `t_final` < 1 is mandatory.")
  }
  if (t_final >= t_inic) {
    stop("Please note that `t_final` must be lower than `t_inic`.")
  }
  if (t_inic * (1 - alpha)^n_iter > t_final) {
    stop("Desired `t_final` is not reached with given `alpha` and `n_iter`
    (i.e.,`alpha` and/or `n_iter` are too small).")
  }
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
  nr <- ncol(m_bin) # no. of relevés
  if (k > nr) {
    stop("Given `k` size is too big.")
  }
  if (k <= 1) {
    stop("Given `k` size is too small.")
  }

  if (use_grasp) {
    res_grasp <- list()
    res_sann <- list()
    for (i in 1:n_runs) {
      # GRASP initial partitions
      par_grasp <- partition_tdv_grasp(m_bin, k, thr = 0.95, verify = FALSE)
      grasp_result <- list(
        par = par_grasp,
        tdv = tdv(m_bin, par_grasp, output_type = "fast")
      )
      # SANN step
      sann_result <- optim_tdv_simul_anne(
        m_bin = m_bin,
        k = k,
        p_initial = par_grasp,
        n_runs = 1,
        n_sol = 1,
        t_inic = t_inic,
        t_final = t_final,
        alpha = alpha,
        n_iter = n_iter,
        use_grasp = FALSE,
        full_output = full_output
      )$SANN[[1]]

      if (full_output) {
        res_grasp[[i]] <- grasp_result
        res_sann[[i]] <- sann_result
      }
      if (!full_output) {
        if (i == 1) {
          res_grasp[[1]] <- grasp_result
          res_sann[[1]] <- sann_result
        } else {
          already_in_bestsol <- any(sapply(res_sann, function(x) {
            identical_partition(x$par, sann_result$par)
          }))
          if (!already_in_bestsol) {
            if (length(res_sann) < n_sol) {
              res_grasp[[length(res_sann) + 1]] <- grasp_result
              res_sann[[length(res_sann) + 1]] <- sann_result
            } else { # Already n_sol components in res_sann
              bestsol_values <- sapply(res_sann, function(x) {
                x$tdv
              })
              if (sann_result$tdv > min(bestsol_values)) {
                worse_bestsol <- which.min(bestsol_values) # Selects the first
                # in case of ties!
                res_grasp[[worse_bestsol]] <- grasp_result
                res_sann[[worse_bestsol]] <- sann_result
              }
            }
          }
        }
      }
    }
    ind_order <- order(sapply(res_sann, function(x) {
      x$tdv
    }), decreasing = TRUE)
    return(list(GRASP = res_grasp[ind_order], SANN = res_sann[ind_order]))
  }
  # If use_grasp == FALSE
  # Simple SANN on a given initial partition

  if (is.null(p_initial)) {
    stop("Argument p_initial can not be NULL when use_grasp = FALSE.")
  }
  if (p_initial[1] != "random") {
    if (!identical(length(p_initial), nr)) {
      stop("Object `p_initial` must be a partition of the columns of `m_bin`")
    }
    mode(p_initial) <- "integer"
    if (!identical(sort(unique(p_initial)), 1:k)) {
      stop(
        "Object `p_initial` is not a valid partition of the columns of `m_bin`"
      )
    }
    p_ini <- p_initial
  }

  # Cooling schedule (cool_sched)
  # nt <- floor(n_iter / ((n_iter * log(1 - alpha)) / (log((1 - alpha) *
  #   t_final / t_inic))))
  # nt <- floor(log((1 - alpha) * t_final / t_inic) / log(1 - alpha))
  nt <- floor(log(t_final / t_inic) / log(1 - alpha))
  cool_sched <- floor(n_iter / nt * (1:nt))

  res_sann <- list()
  for (i in 1:n_runs) {
    if (p_initial[1] == "random") {
      p_ini <- sample(c(1:k, sample(k, nr - k, replace = TRUE)))
    }
    best_p <- cur_p <- p_ini
    best_tdv <- cur_tdv <- tdv(m_bin, p_ini, output_type = "fast")
    tempera <- t_inic
    res_cur_tdv <- res_alt_tdv <- res_alt_tdv <- NULL
    res_probabi <- res_tempera <- NULL

    for (iter in 1:n_iter) {
      alt_p <- random_neighbour_sa(p = cur_p, nr = nr, k = k)
      alt_tdv <- tdv(m_bin, alt_p, output_type = "fast")
      res_alt_tdv <- c(res_alt_tdv, alt_tdv)
      res_cur_tdv <- c(res_cur_tdv, cur_tdv)
      if (alt_tdv >= cur_tdv) {
        cur_p <- alt_p
        cur_tdv <- alt_tdv
        if (cur_tdv > best_tdv) {
          best_p <- cur_p
          best_tdv <- cur_tdv
        }
        probabi <- 1
      } else {
        probabi <- exp((alt_tdv - cur_tdv) / tempera)
        if (stats::runif(1) < probabi) {
          cur_p <- alt_p
          cur_tdv <- alt_tdv
        }
      }
      if (iter %in% cool_sched) {
        tempera <- (1 - alpha) * tempera # Cooling
      }
      res_probabi <- c(res_probabi, probabi)
      res_tempera <- c(res_tempera, tempera)
    }
    if (full_output) {
      res_sann[[i]] <- list(
        current.tdv     = res_cur_tdv,
        alternative.tdv = res_alt_tdv,
        probability     = res_probabi,
        temperature     = res_tempera,
        par             = best_p,
        tdv             = best_tdv
      )
    }
    if (!full_output) {
      if (i == 1) {
        res_sann[[1]] <- list(par = best_p, tdv = best_tdv)
      } else {
        already_in_bestsol <- any(sapply(res_sann, function(x) {
          identical_partition(x$par, best_p)
        }))
        if (!already_in_bestsol) {
          if (length(res_sann) < n_sol) {
            res_sann[[length(res_sann) + 1]] <- list(
              par = best_p,
              tdv = best_tdv
            )
          } else { # Already n_sol components in res_sann
            bestsol_values <- sapply(res_sann, function(x) {
              x$tdv
            })
            if (best_tdv > min(bestsol_values)) {
              worse_bestsol <- which.min(bestsol_values) # Selects the first in
              # case of ties!
              res_sann[[worse_bestsol]] <- list(par = best_p, tdv = best_tdv)
            }
          }
        }
      }
    }
  }
  ind_order <- order(sapply(res_sann, function(x) {
    x$tdv
  }), decreasing = TRUE)
  list(SANN = res_sann[ind_order])
}
