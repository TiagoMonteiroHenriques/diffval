#' Total Differential Value optimization using Hill-climbing algorithms
#'
#' This function searches for partitions of the columns of a given matrix,
#'   optimizing the Total Differential Value (TDV).
#'
#' @param m_bin A matrix. A phytosociological table of 0s (absences) and 1s
#'   (presences), where rows correspond to taxa and columns correspond to
#'   relevés.
#' @param k A numeric giving the number of desired groups.
#' @param p_initial A vector or a character. A vector of integer numbers
#'   with the initial partition of the relevés (i.e., a vector with values from
#'   1 to `k`, with length equal to the number of columns of `m_bin`, ascribing
#'   each relevé to one of the `k` groups). By default, `p_initial = "random"`,
#'   generates a random initial partition.
#' @param n_runs A numeric giving the number of runs to perform.
#' @param n_sol A numeric giving the number of best solutions to keep
#'   in the final output. Defaults to 1.
#' @param maxit A numeric giving the number of iterations of the Hill-climbing
#'   optimization.
#' @param min_g_size A numeric. The minimum number of relevés that a group can
#'   contain (must be 1 or higher).
#' @param stoch_first A logical. `FALSE` (the default), performs only
#'   Hill-climbing on the 1-neighbours; `TRUE` first, performs a Stochastic
#'   Hill-climbing on n-neighbours (n is defined by the parameter
#'   `stoch_neigh_size`), and only after runs the Hill-climbing search on the
#'   1-neighbours; see description above.
#' @param stoch_neigh_size A numeric giving the size (n) of the
#'   n-neighbours for the Stochastic Hill-climbing; only used if
#'   `stoch_first = TRUE`. Defaults to 1.
#' @param stoch_maxit A numeric giving the number of iterations of the
#'   Stochastic Hill-climbing optimization; only used if `stoch_first = TRUE`.
#'   Defaults to 100.
#' @param full_output A logical. If `FALSE` (the default) the best `n_sol`
#'   partitions and respective indices are returned. If `TRUE` (only available
#'   for `n_sol = 1`) the output will also contain information on the
#'   optimization steps (see below).
#' @param verbose A logical. If `FALSE` nothing is printed during the runs.
#'   If `TRUE`, after each run, the run number is printed as well as and
#'   indication if the found partition is a 1-neighbour local maximum.
#'
#' @details Given a phytosociological table (`m_bin`, rows corresponding to
#'   taxa and columns corresponding to relevés) this function searches for
#'   a `k`-partition (`k` defined by the user) optimizing TDV, i.e., searches,
#'   using a Hill-climbing algorithm, for patterns of differential taxa by
#'   rearranging the relevés into `k` groups.
#'
#'   Optimization can start from a random partition (`p_ini = "random"`), or
#'   from a given partition (`p_ini`, defined by the user or produced by any
#'   clustering method, or even a manual classification of the relevés).
#'
#'   Each iteration searches for a TDV improvement screening all 1-neighbours,
#'   until the given number of maximum iterations (`maxit`) is reached. A
#'   1-neighbour of a given partition is another partition obtained by changing
#'   1 relevé (of the original partition) to a different group. A n-neighbour
#'   is obtained, equivalently, ascribing n relevés to different groups.
#'
#'   Optionally, a faster search (Stochastic Hill-climbing) can be performed in
#'   a first step (`stoch_first = TRUE`), consisting on searching for TDV
#'   improvements, by randomly selecting, in each iteration, one n-neighbour (n
#'   defined by the user in the parameter `stoch_neigh_size`), accepting that
#'   n-neighbour partition as a better solution if it improves TDV. This is
#'   repeated until a given number of maximum iterations (`stoch_maxit`) is
#'   reached. Stochastic Hill-climbing might be helpful for big tables (where
#'   the screening of all 1-neighbours might be too time consuming).
#'
#'   Several runs of this function (i.e., multiple starts) should be
#'   tried out, as several local maxima are usually present and the
#'   Hill-climbing algorithm converges easily to local maxima.
#'
#'   Trimming your table by a 'constancy' range or using the result of other
#'   cluster methodologies as input, might help finding interesting partitions.
#'   Specially after trimming the table by a 'constancy' range, getting a random
#'   initial partition with TDV greater than zero might be unlikely; on such
#'   cases using a initial partition from [partition_tdv_grasp()] or
#'   [partition_tdv_grdtp()] (or even the result of other clustering
#'   strategies) as an input partition might be useful.
#'
#' @return If `full_output = FALSE`, a list with (at most) `n_sol` best
#'   solutions (equivalent solutions are removed). Each best solution is also
#'   a list with the following components:
#'
#'   \describe{
#'     \item{local_maximum}{A logical indicating if `par` is a 1-neighbour
#'     local maximum.}
#'     \item{par}{A vector with the partition of highest TDV obtained by the
#'     Hill-climbing algorithm(s).}
#'     \item{tdv}{A numeric with the TDV of `par`.}
#'   }
#'
#'   If `full_output = TRUE`, a list with just one component (one run only),
#'   containing also a list with the following components:
#'
#'   \describe{
#'     \item{res.stoch}{A matrix with the iteration number (of the Stochastic
#'     Hill-climbing phase), the maximum TDV found until that iteration, and the
#'     TDV of the randomly selected n-neighbour in that iteration.}
#'     \item{par.stoch}{A vector with the best partition found in the Stochastic
#'     Hill-climbing phase.}
#'     \item{tdv.stoch}{A numeric showing the maximum TDV found in the
#'     Stochastic Hill-climbing phase (if selected).}
#'     \item{res}{A matrix with the iteration number (of the Hill-climbing), the
#'     maximum TDV found until that iteration, and the highest TDV among all
#'     1-neighbours.}
#'     \item{local_maximum}{A logical indicating if `par` is a 1-neighbour local
#'     maximum.}
#'     \item{par}{A vector with the partition of highest TDV obtained by the
#'     Hill-climbing algorithm(s).}
#'     \item{tdv}{A numeric with the TDV of `par`.}
#'   }
#'
#' @author Tiago Monteiro-Henriques. E-mail: \email{tmh.dev@@icloud.com}.
#'
#' @examples
#' # Getting the Taxus baccata forests data set
#' data(taxus_bin)
#'
#' # Removing taxa occurring in only one relevé in order to
#' # reproduce the example in the original article of the data set
#' taxus_bin_wmt <- taxus_bin[rowSums(taxus_bin) > 1, ]
#'
#' # Obtaining a partition that maximizes TDV using the Stochastic Hill-climbing
#' # and the Hill-climbing algorithms
#'
#' result <- optim_tdv_hill_climb(
#'   m_bin = taxus_bin_wmt,
#'   k = 3,
#'   n_runs = 7,
#'   n_sol = 2,
#'   min_g_size = 3,
#'   stoch_first = TRUE,
#'   stoch_maxit = 500,
#'   verbose = TRUE
#' )
#'
#' # Inspect the result. The highest TDV found in the runs.
#' result[[1]]$tdv
#' # If result[[1]]$tdv is 0.1958471 you are probably reproducing the three
#' # groups (Estrela, Gerês and Galicia) from the original article. If not
#' # try again the optim_tdv_hill_climb function (maybe increasing n_runs).
#'
#' # Plot the sorted (or tabulated) phytosociological table
#' tabul1 <- tabulation(
#'   m_bin = taxus_bin_wmt,
#'   p = result[[1]]$par,
#'   taxa_names = rownames(taxus_bin_wmt),
#'   plot_im = "normal"
#' )
#'
#' # Plot the sorted (or tabulated) phytosociological table, also including
#' # taxa occurring just once in the matrix
#' tabul2 <- tabulation(
#'   m_bin = taxus_bin,
#'   p = result[[1]]$par,
#'   taxa_names = rownames(taxus_bin),
#'   plot_im = "normal"
#' )
#'
#' @export
optim_tdv_hill_climb <- function(m_bin,
                                 k,
                                 p_initial = "random",
                                 n_runs = 1,
                                 n_sol = 1,
                                 maxit = 10,
                                 min_g_size = 1,
                                 stoch_first = FALSE,
                                 stoch_neigh_size = 1,
                                 stoch_maxit = 100,
                                 full_output = FALSE,
                                 verbose = FALSE) {
  stopifnot("`m_bin` must be a matrix." = is.matrix(m_bin))
  mode(m_bin) <- "integer"
  if (!identical(c(0L, 1L), sort(unique(as.vector(m_bin))))) {
    stop("Matrix `m_bin` must contain only 0's and 1's.")
  }
  if (min(rowSums(m_bin)) == 0) {
    stop("At least one taxa is not present in any relev\u00e9.")
  }
  if (min(colSums(m_bin)) == 0) {
    stop("At least one relev\u00e9 contains no taxa.")
  }
  mgs <- as.integer(min_g_size)
  if (mgs < 1) {
    stop("Object `min_g_size` must be greater than or equal to 1.")
  }
  mt <- t(m_bin)
  nr <- nrow(mt) # No. of relevés
  ns <- ncol(mt) # No. of taxa
  if (nr <= mgs * k) {
    stop(paste0(
      "Random partition cannot guarantee at least ",
      mgs,
      " relev\u00e9s per group!"
    ))
  }

  if (p_initial[1] != "random") {
    p_ini <- p_initial
    if (!identical(as.integer(sort(unique(p_ini))), 1:k)) { # maybe when
      # p_initial is given k could be ignored
      stop("Object `p` is not a valid partition of the columns of m.")
    }
    if (!identical(length(p_ini), nr)) {
      stop(
        "Object `p_ini` must be a partition of the columns of matrix `m_bin`."
      )
    }
    tp <- tabulate(p_ini) # size of each group (inner)
    if (min(tp) < mgs) {
      stop(paste0(
        "At least one group of the provided partition has less than ",
        mgs,
        " elements."
      ))
    }
    if (max(tp) == mgs) {
      stop(paste0(
        "At least one group of the provided partition has to have more than ",
        mgs,
        " elements."
      ))
    }
  }
  if (n_sol > n_runs) {
    stop("The number of runs (`n_runs`) should not be lower than the desired
         number of best solutions (`n_sol`).")
  }
  if ((n_sol != 1 || n_runs != 1) && full_output == TRUE) {
    stop("The option `full_output = TRUE` is only available for `n_runs == 1`
         and `n_sol = 1`.")
  }

  # Special case of n_sol == 1 & n_runs == 1 & full_output == TRUE
  if (!full_output) {
    res_list <- list()
  }
  for (n.run in 1:n_runs) {
    if (p_initial[1] == "random") {
      p_ini <- sample(c(rep(1:k, mgs), sample(k, nr - mgs * k, replace = TRUE)))
      tp <- tabulate(p_ini) # Size of each group (inner)
    } else {
      p_ini <- p_initial # Needed only for n.run > 1... (can be improved,
      # additionally for p_ini not random, should not repeat
      # the first calculation!)
      tp <- tabulate(p_ini)
    }

    # First calculation of TDV (i.e., current value, curr_val, for p_ini),
    # keeping the intermediate steps of the tdv calculation
    outer_size <- nr - tp # Sum of the sizes of the outer groups
    ofda <- ifp <- matrix(0, k, ns) # Matrices to store [a/b], i.e., the inner
    # frequency of presences (ifp) and [c/d], i.e., the outer frequency of
    # differenciating absences (ofda)
    afg <- rowsum(mt, group = p_ini) # No. of relevés containing the taxon,
    # within each group (absolute frequency in each group)
    empty_size <- (afg == 0) * tp # No. of relevés of each group (group size),
    # when the taxon is not present
    gct <- colSums(afg > 0) # No. of groups containing the taxon [e]
    i_mul <- gct > 1 # Indices of the taxa occurring in more than one group
    # (taxa must occur in at least one group)
    for (g in 1:k) { # Fills matrices ofda [c/d] and ifp [a/b], only when the
      # taxon is present in the group!
      i_tx <- afg[g, ] > 0 # Indices of the taxa present in the group g
      ofda[g, i_tx & !i_mul] <- 1 # ofda is 1 for the taxa occurring in one
      # group only!
      if (sum(i_mul) > 0) { # If there are taxa occurring in more than one group
        i_tx_mul <- i_tx & i_mul # Taxa of group g occurring also in other
        # group than g
        ofda[g, i_tx_mul] <- colSums(
          empty_size[-g, i_tx_mul, drop = FALSE] / outer_size[g]
        ) # Size of outer empty groups divided by the sum of the sizes of the
        # outer groups [c/d]
      }
      ifp[g, i_tx] <- afg[g, i_tx] / tp[g] # Presences inside group g divided by
      # the group g size [a/b]
    }
    dv <- colSums(ifp * ofda) / gct
    curr_val <- sum(dv) / ns

    res_stoch <- NULL
    par_stoch <- NULL
    curr_val_stoch <- NULL

    if (stoch_first == TRUE) { # Stochastic Hill-climbing (stoch_first = TRUE)
      p_curr <- p_ini
      if (full_output == TRUE) {
        res_stoch <- matrix(0, stoch_maxit, 3)
      }
      for (iter in 1:stoch_maxit) {
        p_neig <- random_neighbour_hc(
          p = p_curr,
          k = k,
          mgs = mgs,
          stoch_neigh_size = stoch_neigh_size
        )
        neig_val <- tdv_aux(mt = mt, p = p_neig, k = k, ns = ns, nr = nr)
        if (neig_val >= curr_val) {
          p_curr <- p_neig
          curr_val <- neig_val
        }
        if (full_output == TRUE) {
          res_stoch[iter, ] <- c(iter, curr_val, neig_val)
        }
      }
      p_ini <- p_curr
      if (full_output == TRUE) {
        par_stoch <- p_curr
        curr_val_stoch <- curr_val
      }

      # Replacing the intermediate steps of the tdv calculation, for
      # p_ini/p_curr coming from the Stochastic Hill-climbing (it is
      # calculating again, actually, but it is preferable this way)
      tp <- tabulate(p_ini) # Size of each group (inner)
      outer_size <- nr - tp # Sum of the sizes of the outer groups
      ofda <- ifp <- matrix(0, k, ns) # Matrices to store [a/b], i.e., the inner
      # frequency of presences (ifp) and [c/d], i.e., the outer frequency of
      # differenciating absences (ofda)
      afg <- rowsum(mt, group = p_ini) # No. of relevés containing the taxon,
      # within each group (absolute frequency in each group)
      empty_size <- (afg == 0) * tp # No. of relevés of each group (group size),
      # when the taxon is not present
      gct <- colSums(afg > 0) # No. of groups containing the taxon [e]
      i_mul <- gct > 1 # Indices of the taxa occurring in more than one
      # group (taxa must occur in at least one group)
      for (g in 1:k) { # Fills matrices ofda [c/d] and ifp [a/b],
        # only when the taxon is present in the group!
        i_tx <- afg[g, ] > 0 # Indices of the taxa present in the group g
        ofda[g, i_tx & !i_mul] <- 1 # ofda is 1 for the taxa occurring in one
        # group only!
        if (sum(i_mul) > 0) { # If there are taxa occurring in more than
          # one group
          i_tx_mul <- i_tx & i_mul # Taxa of group g occurring also in other
          # group than g
          ofda[g, i_tx_mul] <- colSums(
            empty_size[-g, i_tx_mul, drop = FALSE] / outer_size[g]
          ) # Size of outer empty groups divided by the sum of the sizes of the
          # outer groups [c/d]
        }
        ifp[g, i_tx] <- afg[g, i_tx] / tp[g] # Presences inside group g divided
        # by the group g size [a/b]
      }
      dv <- colSums(ifp * ofda) / gct
      # curr_val was already calculated in the Stochastic Hill-climbing
    }

    # Hill-climbing
    if (maxit == 0) {
      if (full_output == TRUE) {
        res <- NULL
      }
      p_curr <- p_ini
      loc_max <- FALSE
    } else {
      loc_max <- FALSE
      p_curr <- p_ini

      # curr_val is already calculated (during or before the Stochastic
      # Hill-climbing)
      if (full_output == TRUE) {
        res <- matrix(0, maxit, 3)
      }
      for (iter in 1:maxit) {
        temp <- max_tdv_neighbour(
          mt = mt,
          p = p_curr,
          k = k,
          ns = ns,
          nr = nr,
          mgs = mgs,
          ofda = ofda,
          ifp = ifp,
          afg = afg,
          empty_size = empty_size,
          gct = gct,
          i_mul = i_mul,
          dv = dv
        )
        p_neig <- temp$p
        neig_val <- temp$tdv
        if (full_output == TRUE) {
          res[iter, ] <- c(iter, curr_val, neig_val)
        }
        if (neig_val > curr_val) {
          p_curr <- p_neig
          curr_val <- neig_val
          ofda <- temp$ofda
          ifp <- temp$ifp
          afg <- temp$afg
          empty_size <- temp$empty_size
          gct <- temp$gct
          i_mul <- temp$i_mul
          dv <- temp$dv
        } else {
          loc_max <- TRUE
          if (full_output == TRUE) {
            res[iter, 3] <- curr_val
          }
          break
        }
      }
      if (verbose) {
        cat("Run number:", n.run, "Confirmed local maximum:", loc_max, "\n")
      }
    }
    if (full_output == TRUE) {
      if (!is.null(res)) {
        rows_not_zero <- apply(res, 1, function(x) {
          any(as.logical(x))
        })
        res <- res[rows_not_zero, , drop = FALSE]
      }
      res_list <- list()
      res_list[[1]] <- list(
        res.stoch     = res_stoch,
        par.stoch     = par_stoch,
        tdv.stoch     = curr_val_stoch,
        res           = res,
        local_maximum = loc_max,
        par           = p_curr,
        tdv           = curr_val
      )
      return(res_list)
    } else { # I.e., full_output is FALSE keeping only (at most) n_sol best
      # solutions (deleting repeated ones)
      if (n.run == 1) {
        res_list[[1]] <- list(
          local_maximum = loc_max,
          par           = p_curr,
          tdv           = curr_val
        )
      } else {
        already_in_bestsol <- any(sapply(res_list, function(x) {
          identical_partition(x$par, p_curr)
        }))
        if (!already_in_bestsol) {
          if (length(res_list) < n_sol) {
            res_list[[length(res_list) + 1]] <- list(
              local_maximum = loc_max,
              par           = p_curr,
              tdv           = curr_val
            )
          } else { # Already n_sol components in res_list
            bestsol_values <- sapply(res_list, function(x) {
              x$tdv
            })
            if (curr_val > min(bestsol_values)) {
              worse_bestsol <- which(
                bestsol_values == min(bestsol_values)
              )[1] # [1] selects one, in case of ties!
              res_list[[worse_bestsol]] <- list(
                local_maximum = loc_max,
                par           = p_curr,
                tdv           = curr_val
              )
            }
          }
        }
      }
    }
  }
  ind_order <- order(sapply(res_list, function(x) {
    x$tdv
  }), decreasing = TRUE)
  res_list[ind_order]
}
