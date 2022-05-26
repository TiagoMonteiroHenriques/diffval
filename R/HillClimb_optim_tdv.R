# HillClimb_optim_tdv.R
#'
#' @title Total Differential Value optimization using Hill-climbing algorithms
#'
#' @description This function searches for partitions of the columns of a given matrix, optimizing the Total Differential Value (TDV).
#'
#' @param m.bin A `matrix`, i.e. a phytosociological table of 0s (absences) and 1s (presences), where rows correspond to taxa and columns correspond to relevés.
#' @param p.initial A `vector` or a `character`. A `vector` of integer numbers with the initial partition of the relevés (i.e. a vector with values from 1 to k, with length equal to the number of columns of `m.bin`, ascribing each relevé to one of the k groups). By default, `p.initial = "random"`, generates a random initial partition.
#' @param k A `numeric`, giving the number of desired groups.
#' @param n.runs A `numeric`, giving the number of runs to perform.
#' @param n.sol A `numeric`, giving the number of best solutions to keep in the final output. Defaults to 1.
#' @param maxit A `numeric` giving the number of iterations of the Hill-climbing optimization.
#' @param min.g.size A `numeric` The minimum number of relevés that a group can contain (must be 1 or higher).
#' @param stoch.first A `logical`. `FALSE` (the default), performs only Hill-climbing on the 1-neighbours; `TRUE` first, performs a Stochastic Hill-climbing on `n`-neighbours (`n` is definded by the parameter `stoch.neigh.size`), and only after runs the Hill-climbing search on the 1-neighbours; see description above.
#' @param stoch.neigh.size	A `numeric`, giving the size (`n`) of the `n`-neighbour for the Stochastic Hill-climbing; only used if `stoch.first` = `TRUE`. Defaults to 1.
#' @param stoch.maxit A `numeric`, giving the number of iterations of the Hill-climbing optimization; only used if `stoch.first` = `TRUE`. Defaults to 100.
#' @param full.output A `logical`. If `FALSE` (the default) the best `n.sol` partitions and respective indices are returned. If `TRUE` (only available for `n.sols = 1`) the output will also contain information on the optimization steps (see below).
#' @param verbose A `logical`. If `FALSE` nothing is printed during the runs. If `TRUE`, after each run, the run number is printed as well as and indication if the found partition is a 1-neighbour local maximum.
#'
#' @details Given a phytosociological table (`m.bin`, rows corresponding to taxa and columns corresponding to relevés) this function searches for
#' a k-partition (`k` defined by the user) optimizing TDV, i.e. searches, using a Hill-climbing algorithm, for patterns of differential taxa by
#' rearranging the relevés into k groups.
#'
#' Optimization can start from a random partition (`p.ini = "random"`), or from a given partition (`p.ini`, defined by the user or produced by any
#' clustering method, or even a manual classification).
#'
#' Each iteration searches for a TDV improvement screening all 1-neighbours, until the given number of maximum iterations (`maxit`) is reached. A
#' 1-neighbour of a given partition is another partition obtained by changing 1 relevé (of the original partition) to a different group. A n-neighbour
#' is obtained, equivalently, ascribing n relevés to different groups.
#'
#' Optionally, a faster search (Stochastic Hill-climbing) can be performed in a first step (`stoch.first` = `TRUE`), consisting on searching for
#' TDV improvements, by randomly selecting n-neighbours (n defined by the user with the parameter `stoch.neigh.size`), accepting that neighour partition
#' as a better solution if it improves TDV. This is repeated until a given number of maximum iterations (`stoch.maxit`) is reached. Stochastic Hill-climbing
#' might be helpful for big tables (where the simple screening of all 1-neighbours might be too time consuming).
#'
#' Several runs of HillClimb_optim_tdv (i.e. multiple starts) should be tried out, as several local maxima are usually present and the Hill-climbing
#' algorithm converges easily to local maxima.
#'
#' Trimming your table by a 'constancy' range or using the result of other cluster methodologies as input, might help finding interesting partitions.
#' Specially after trimming the table by a 'constancy' range, getting a random initial partition with TDV greater than zero might be very unlikely; on
#' such cases using a initial partition from \code{\link{GRASP_partition_tdv}} or \code{\link{GRDTP_partition_tdv}} (or even the result of other clustering
#' strategies) as input might be useful.
#'
#' @return If `full.output = FALSE`, a `list` with (at most) `n.sol` best solutions (equivalent solutions are removed). Each best solution is also
#' a `list` with the following components:
#'
#' \describe{
#'   \item{local_maximum}{A `logical` indicating if `par` is a 1-neighbour local maximum.}
#'   \item{par}{A `vector` with the partition of highest TDV obtained by the Hill-climbing algorithm(s).}
#'   \item{tdv}{A `numeric` with the TDV of `par`.}
#' }
#'
#' If `full.output = TRUE`, a `list` with just one component (one run only), containing also a list with the following components:
#'
#' \describe{
#'   \item{res.stoch}{A `matrix` with the iteration number (of the Stochastic Hill-climbing phase), the maximum TDV found until that iteration, and the higher TDV among all 1-neighbours.}
#'   \item{par.stoch}{A `vector` with the best partition found in the Stochastic Hill-climbing phase.}
#'   \item{tdv.stoch}{A `numeric` showing the maximum TDV found in the Stochastic Hill-climbing phase (if selected).}
#'   \item{res}{A `matrix` with the iteration number (of the Hill-climbing), the maximum TDV found until that iteration, and the higher TDV among all 1-neighbours.}
#'   \item{local_maximum}{A `logical` indicating if `par` is a 1-neighbour local maximum.}
#'   \item{par}{A `vector` with the partition of highest TDV obtained by the Hill-climbing algorithm(s).}
#'   \item{tdv}{A `numeric` with the TDV of `par`.}
#' }
#'
#' @author Tiago Monteiro-Henriques. E-mail: \email{tmh.dev@@icloud.com}.
#'
#' @examples
#'
#' #getting the Taxus baccata forests data set
#' data(taxus_bin)
#'
#' #removing taxa occurring in only one relevé in order to
#' #try reproducing the example in the original article of the data set
#' taxus_bin_wmt <- taxus_bin[rowSums(taxus_bin) > 1,]
#'
#' #obtaining a partition that maximizes TDV using the Stochastic Hill-climbing
#' #and the Hill-climbing algorithms
#'
#' result <- HillClimb_optim_tdv(taxus_bin_wmt, k = 3, n.runs = 7, n.sol = 2,
#'   min.g.size = 3, stoch.first = TRUE, stoch.maxit = 500, verbose = TRUE)
#'
#' #inspect the result
#' #the highest TDV found in the runs
#' result[[1]]$tdv
#' #if result[[1]]$tdv is 0.1958471 you are probably reproducing the three
#' #groups (Estrela, Gerês and Galicia) from the original article; if not
#' #try again the HillClimb_optim_tdv function (maybe increasing n.runs)
#'
#' #plot the sorted (or tabulated) phytosociological table
#' tabul1 <- tabulation(taxus_bin_wmt, result[[1]]$par, rownames(taxus_bin_wmt), "normal")
#'
#' #plot the sorted (or tabulated) phytosociological table, also including
#' #taxa with occurring just once in the matrix
#' tabul2 <- tabulation(taxus_bin, result[[1]]$par, rownames(taxus_bin), "normal")
#'
#' @export
#'
HillClimb_optim_tdv <- function(m.bin, p.initial = "random", k, n.runs = 1, n.sol = 1, maxit = 10, min.g.size = 1, stoch.first = FALSE, stoch.neigh.size = 1, stoch.maxit = 100, full.output = FALSE, verbose = FALSE) {
  stopifnot(is.matrix(m.bin))
  mode(m.bin) <- "integer"
  if (!identical(c(0L,1L), sort(unique(as.vector(m.bin))))) {stop("Matrix m.bin should contain only 0's and 1's.")}
  if (min(rowSums(m.bin)) == 0) {stop("At least one taxa is not present in any relev\u00e9.")}
  if (min(colSums(m.bin)) == 0) {stop("At least one relev\u00e9 contains no taxa.")}
  mgs <- as.integer(min.g.size)
  if (mgs < 1) {stop("Object min.g.size must be greater than or equal to 1.")}
  mt <- t(m.bin)
  nr <- nrow(mt) # no. of relevés
  ns <- ncol(mt) #no. of taxa
  if (nr <= mgs * k) {stop(paste0("Random partition cannot guarantee at least ", mgs, " relev\u00e9s per group!"))}

  if (p.initial[1] != "random") {
    p.ini <- p.initial
    if (!identical(as.integer(sort(unique(p.ini))), 1:k)) { #maybe when p.ini is given k could be ignored
      stop("Object p is not a valid partition of the columns of m")
    }
    if (!identical(length(p.ini), nr)) {
      stop("Object p.ini must be a partition of the columns of matrix m.bin.")
    }
    tp <- tabulate(p.ini) #size of each group (inner)
    if (min(tp) < mgs) {
      stop(paste0("At least one group of the provided partition has less than ", mgs, " elements"))
    }
    if (max(tp) == mgs) {
      stop(paste0("At least one group of the provided partition has to have more than ", mgs, " elements"))
    }
  }
  if (n.sol > n.runs) {
    stop("The number of runs ('n.runs') should not be lower than the desired number of best solutions ('n.sol').")
  }
  if ((n.sol != 1 | n.runs != 1) & full.output == TRUE) {
    stop("The option 'full.output = TRUE' is only available for 'n.runs == 1' and 'n.sol = 1'.")
  }

  #SPECIAL CASE OF n.sol == 1 & n.runs == 1 & full.output == TRUE
  if (!full.output) {
    res.list <- list()
  }
  for (n.run in 1:n.runs) {
    if (p.initial[1] == "random") {
      p.ini <- sample(c(rep(1:k, mgs), sample(k, nr - mgs * k, replace = TRUE)))
      tp <- tabulate(p.ini) #size of each group (inner)
    }  else {
      p.ini <- p.initial #needed only for n.run > 1... (can be improved, aditionally for p.ini not random should not repeat the "first calculation!")
      tp <- tabulate(p.ini)
    }

    #first calculation of TDV (i.e. current value, curr.val, for p.ini), keeping the intermediate steps of the tdv calculation
    outer.size <- nr - tp #sum of the sizes of the outer groups
    ofda <- ifp <- matrix(0, k, ns) #matrices to store [a/b], i.e. the inner frequency of presences (ifp) and [c/d], i.e. the outer frequency of differenciating absences (ofda)
    afg <- rowsum(mt, group = p.ini) #no. of relevés containing the taxon, within each group (absolute frequency in each group)
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
    DV <- colSums(ifp * ofda) / gct
    curr.val <- sum(DV) / ns
    #curr.val <- tdv.aux(mt = mt, p = p.ini, tp = tp, k = k, ns = ns, nr = nr)

    res.stoch <- NULL
    par.stoch <- NULL
    curr.val.stoch <- NULL

    if (stoch.first == TRUE) { #STOCHASTIC HILL CLIMBING (stoch.first = TRUE)
      p.curr <- p.ini
      if (full.output == TRUE) {
        res.stoch <- matrix(0, stoch.maxit, 3)
      }
      for (iter in 1:stoch.maxit) {
        p.neig <- random.neighbour(p = p.curr, k = k, mgs = mgs, stoch.neigh.size = stoch.neigh.size)
        neig.val <- tdv.aux(mt = mt, p = p.neig, k = k, ns = ns, nr = nr)
        if (neig.val >= curr.val) {
          p.curr <- p.neig
          curr.val <- neig.val
        }
        if (full.output == TRUE) {
          res.stoch[iter,] <- c(iter, curr.val, neig.val)
        }
      }
      p.ini <- p.curr
      if (full.output == TRUE) {
        par.stoch <- p.curr
        curr.val.stoch <- curr.val
      }

      #replacing the intermediate steps of the tdv calculation, for p.ini/p.curr coming from Stochastic Hill-climbing (it is calculating again, actually, but it is preferable this way)
      tp <- tabulate(p.ini) #size of each group (inner)
      outer.size <- nr - tp #sum of the sizes of the outer groups
      ofda <- ifp <- matrix(0, k, ns) #matrices to store [a/b], i.e. the inner frequency of presences (ifp) and [c/d], i.e. the outer frequency of differenciating absences (ofda)
      afg <- rowsum(mt, group = p.ini) #no. of relevés containing the taxon, within each group (absolute frequency in each group)
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
      DV <- colSums(ifp * ofda) / gct
      #curr.val was already calculated in the Stochastic Hill-climbing
    }

    #HILL CLIMBING
    if (maxit == 0) {
      if (full.output == TRUE) {
        res <- NULL
      }
      p.curr <- p.ini
      loc_max <- FALSE
    } else {
      loc_max <- FALSE
      p.curr <- p.ini
      #curr.val is already calculated (during or before the Stochastic Hill-climbing)
      if (full.output == TRUE) {
        res <- matrix(0, maxit, 3)
      }
      for (iter in 1:maxit) {
        temp <- max.tdv.neighbour(mt = mt, p = p.curr, k = k, ns = ns, nr = nr, mgs = mgs, ofda = ofda, ifp = ifp, afg = afg, empty.size = empty.size, gct = gct, i.mul = i.mul, DV = DV)
        p.neig <- temp$p
        neig.val <- temp$tdv
        if (full.output == TRUE) {
          res[iter,] <- c(iter, curr.val, neig.val)
        }
        if (neig.val > curr.val) {
          p.curr <- p.neig
          curr.val <- neig.val
        } else {
          loc_max <- TRUE
          if (full.output == TRUE) {
            res[iter, 3] <- curr.val
          }
          break
        }
      }
      if (verbose) {
        cat("Run number:", n.run,"Confirmed local maximum:", loc_max, "\n")
      }
    }
    if (full.output == TRUE) {
      if (!is.null(res)) {
        rows.not.zero <- apply(res, 1, function (x) {any(as.logical(x))})
        res <- res[rows.not.zero,,drop = FALSE]
      }
      res.list <- list()
      res.list[[1]] <- list(res.stoch = res.stoch, par.stoch = par.stoch, tdv.stoch = curr.val.stoch, res = res, local_maximum = loc_max, par = p.curr, tdv = curr.val)
      return(res.list)
    } else { # full.output == FALSE
      #keeping only (at most) n.sol best solutions (deleting repeated ones)
      if (n.run == 1) {
        res.list[[1]] <- list(local_maximum = loc_max, par = p.curr, tdv = curr.val)
      } else {
        already.in.bestsol <- any(sapply(res.list, function (x) {
          equivalent_partition(x$par, p.curr)
        }))
        if (!already.in.bestsol) {
          if (length(res.list) < n.sol) {
            res.list[[length(res.list) + 1]] <- list(local_maximum = loc_max, par = p.curr, tdv = curr.val)
          } else { #already n.sol components in res.list
            best.sol.values <- sapply(res.list, function (x) {x$tdv})
            if (curr.val >  min(best.sol.values)) {
              worse.bestsol <- which(best.sol.values == min(best.sol.values))[1] #selects one, in case of ties!
              res.list[[worse.bestsol]] <- list(local_maximum = loc_max, par = p.curr, tdv = curr.val)
            }
          }
        }
      }
    }
  }
  ind.order <- order(sapply(res.list, function (x) {x$tdv}), decreasing = TRUE)
  return(res.list[ind.order])
}
