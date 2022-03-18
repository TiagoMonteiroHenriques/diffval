# SimulAnne_optim_tdv.R
#'
#' @title TotDiffVal (or TotDiffVal1) optimization using the simulated annealing (and GRASP) algorithm(s)
#'
#' @description This function searches for partitions of the columns of a given matrix, optimizing the TotDiffVal (or TotDiffVal1) index using simulated annealing (and GRASP) algorithm(s).
#'
#' @param m A `matrix`, i.e. a phytosociological table of 0s (absences) and 1s (presences), where rows correspond to taxa and columns correspond to relevés.
#' @param k A `numeric` giving the number of desired groups.
#' @param index A `character` selecting which index is calculated: "TotDiffVal" or "TotDiffVal1".
#' @param p.initial A `vector` of integer numbers with the partition of the relevés (i.e. a k-partition, consisting in a vector with values from 1 to k, with length equal to the number of columns of m, ascribing each relevé to one of the k groups), to be used as initial partition in the simulated annealing, when `use.GRASP = FALSE`.
#' @param full.output A `logical`. If `FALSE` (the default) the best `n.sol` partitions and respective indices are returned. If `TRUE` (only available for `n.sols = 1`) the output will also contain information on the optimization steps (see below).
#' @param n.runs A `numeric` giving the number of runs (it will define the number of seeds). Defaults to 10.
#' @param n.sol A `numeric`, giving the number of best solutions to keep in the final output. Defaults to 5.
#' @param T_inic A `numeric` giving the initial temperature. Defaults to 1.
#' @param T_final A `numeric` giving the initial temperature. Defaults to 0.001.
#' @param alpha A `numeric` giving the fraction of temperature drop to be used in the temperature reduction scheme (see Details). Defaults to 0.05.
#' @param n.iter A `numeric` giving the initial temperature. Defaults to 1000
#' @param use.GRASP A `logical`. Defaults to `TRUE`. IF `TRUE`, a GRASP procedure is performed on the seeds and then passed to a new simulated annealing phase. If `FALSE` the previously obtained seeds are passed to subsequent phase. If `TRUE` only that first phase of simulated annealing is performed returning a collection of possible seeds for further analysis.
#' @param full.output A `logical`. Defaults to `FALSE`. IF `TRUE` extra information is presented in the output. See Value..
#' @param consider.zeros A `logical. TEMPORARY ARGUMENT TO DO SOME TESTING. ELIMINATE AFTER.
#' @param improve.seeds A `logical. TEMPORARY ARGUMENT TO DO SOME TESTING. ELIMINATE AFTER.
#'
#' @details This function uses two main algorithms:
#' 1) A greedy randomized adaptive search procedure (GRASP), which is used to obtain first solutions (partitions of `m`) that can then be submitted to the next algorithm.
#' 2) A stochastic global optimization method, simulated annealing (SANN), which is used to search for the global maximum of TotDiffVal.
#'
#' Given a phytosociological table (`m`, with rows corresponding to taxa and columns corresponding to relevés) this function searches for a k-partition (`k`, defined by the user) optimizing TotDiffVal (or TotDiffVal1) index (see http://home.isa.utl.pt/~tmh/), i.e. searches, using a SANN algorithm (possible working on GRASP solutions), for the global maximum of the selected index (by rearranging the relevés into k groups).
#'
#' The default procedure (`use.GRASP = TRUE`, which ignores `p.ini`) performs the following steps:
#'
#' 1) A first step where: i) an initial set of 'group seeds' is obtained randomly (i.e. `n.runs` sets of `k` columns of `m`); ii) these 'group seeds' are submitted to a first SANN optimization, obtaining sets of 'improved group seeds' (i.e., `k` columns of `m` with high value of TotDiffVal (or TotDiffVal1)); iii) The 'improved group seeds' are submitted to the GRASP to obtain partitions of the columns of `m` matrix.
#' 2) In a second step, the partitions obtained by GRASP are used as initial partitions in a SANN optimization.
#'
#' The temperature reduction scheme (CooS) used in the SANN procedure is...
#'
#' Argument `alpha`must be bounded between 0 and 1.
#'
#' @return If `full.output = FALSE` (the default), a `list` with the following component:
#'
#' \describe{
#'   \item{SANN}{A `list` with `n.runs` components, each one containing also a `list` with two components:
#'   1) 'par', containing the partition found in that run, corresponding to the highest value of TotDiffVal (or TotDiffVal1);
#'   2) 'tdv', with the value of TotDiffVal (or TotDiffVal1) of the referred partition.}
#'}
#'
#' If `full.output = TRUE`, a `list` with the following components (the first two components are only given if `use.GRASP = TRUE`):
#'
#' \describe{
#'   \item{improved.seeds}{A `list` with `n.runs` components, each one containing also a `list` with two components: 'seed', containing the `k` columns of `m` used as improved seed in that run, and 'partial.tdv', which corresponds to the TotDiffVal (or TotDiffVal1) of only those `k` columns.}
#'   \item{GRASP}{A `list` with `n.runs` components, each one containing also a `list` with two components: 'par', containing the partition obtained by GRASP, and 'tdv', which corresponds to the TotDiffVal (or TotDiffVal1) of that partition.}
#'   \item{SANN}{A `list` with `n.runs` components, each one containing also a `list` with six components:
#'   1) 'par', containing the partition found in that run, corresponding to the highest value of TotDiffVal (or TotDiffVal1);
#'   2) 'tdv', with the value of TotDiffVal (or TotDiffVal1) of the referred partition;
#'   3) 'cur.tdv', a `vector` of length `n.iter` with the current TotDiffVal (or TotDiffVal1) in each SANN iteration;
#'   4) 'alt.tdv', a `vector` of length `n.iter` with the alternative TotDiffVal (or TotDiffVal1) in each SANN iteration;
#'   5) 'probability', a `vector` of length `n.iter` with the probability TotDiffVal (or TotDiffVal1) in each SANN iteration;
#'   6) 'temperature', a `vector` of length `n.iter` with the temperature TotDiffVal (or TotDiffVal1) in each SANN iteration.
#'   }
#' }
#'
#' @author Jorge Orestes Cerdeira and Tiago Monteiro-Henriques. E-mail: \email{tiagomonteirohenriques@@gmail.com}.
#'
#' @export
#'
SimulAnne_optim_tdv <- function(m, k, index = "TotDiffVal", p.initial = NULL, n.runs = 10, n.sol = 5, T_inic = 1, T_final = 0.001, alpha = 0.05, n.iter = 1000, use.GRASP = TRUE, full.output = FALSE, consider.zeros = TRUE, improve.seeds = TRUE) {
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

  if (use.GRASP) { #First step
    mt <- as.matrix(t(m))

    #SANN for improved seeds + GRASP for initial partitions
    RES.seeds <- list()
    RES.GRASP <- list()
    RES.SANN <- list()
    for (i in 1:n.runs) {

      if (improve.seeds) {#SANN improved seeds
        RES <- SANN_seeds_index(m = m, k = k, nr = nr, ns = ns, T_inic = T_inic, T_final = T_final, alpha = alpha, n.iter = n.iter)
        RES.seeds[[i]] <- list(seed = RES$par, partial.tdv = RES$partial.tdv)
      } else {
        #simple random seeds
        seed.random <- sample(1:nr, k)
        RES <- list()
        RES$partial.tdv <- tdv.constructive(m[, seed.random], 1:k, k, index, consider.zeros)
        RES$par <- seed.random
        RES.seeds[[i]] <- list(seed = RES$par, partial.tdv = RES$partial.tdv)
      }

      #GRASP initial partitions
      thr <- 0.95
      par.GRASP <- rep(0, nr)
      seed <- RES.seeds[[i]]$seed
      par.GRASP[seed] <- 1:k

      #preparing mat_cur (the matrix to assist DV and TDV calculation)
      mat_cur <- matrix(0, ns, 6 * k + 2)
      ind_a <- 1 + 0:(k-1) * 6
      ind_b <- ind_a + 1
      ind_c <- ind_b + 1
      ind_d <- ind_c + 1
      ind_ab <- ind_d + 1
      ind_cd <- ind_ab + 1
      ind_e <- k * 6 + 1
      ind_usable <- ind_e + 1

      presences.in.groups <- rowSums(m[,seed])
      absences.in.groups <- k - presences.in.groups

      #for the special case of seed (i.e. only one relevé in each of the groups):
      mat_cur[,ind_a] <- m[,seed]
      mat_cur[,ind_b] <- 1
      mat_cur[,ind_c] <- absences.in.groups - !m[,seed] #CF. THIS
      mat_cur[,ind_d] <- k-1
      mat_cur[,ind_ab] <- mat_cur[,ind_a] #at this stage mat_cur[,ind_b] is always 1
      mat_cur[,ind_cd] <- mat_cur[,ind_c]/mat_cur[,ind_d]
      mat_cur[,ind_e] <- presences.in.groups
      mat_cur[,ind_usable] <- as.numeric(mat_cur[,ind_e] != k)

      present.in.groups <- presences.in.groups != 0 #present in at least one group (i.e. not absent)
      usable.row <- as.logical(mat_cur[,ind_usable]) #taxa not present in all groups (as taxa that is in all groups have DiffVal = 0)
      p_a_u <- present.in.groups & usable.row #(present and usable) present taxa in at least in one group but not in all of them

      ind_0 <- (1:k-1)*2+1 #k indices to store values when adding 0 values
      ind_1 <- (1:k-1)*2+2 #k indices to store values when adding 1 values

      #while (0 %in% par.GRASP & any(as.logical(mat_cur[,ind_usable]))) { #do while there is an empty column in the partition and at least one line is "usable" (CF!)
      while (0 %in% par.GRASP) { #do while there is an empty column in the partition (what happen if all lines are not usable? EXPERIMENTAR COM MATRIZ SÓ DE 1s?)
        rel.ind <- which(par.GRASP==0)
        res.mat <- matrix(0, length(rel.ind), k)
        currDV <- rep(0, ns)
        if (all(!usable.row)) {warning("Attention, there are no usable taxa for TDV calculation!")}
        currDV[p_a_u] <- rowSums(mat_cur[p_a_u, ind_ab, drop = FALSE] * mat_cur[p_a_u, ind_cd, drop = FALSE] / mat_cur[p_a_u, ind_e]^2)
        DVchanges <- get_DV_01(m = m, k = k, mat_cur = mat_cur, usable.row = usable.row, p_a_u = p_a_u, ns = ns, ind_0 = ind_0, ind_1 = ind_1, ind_a =ind_a, ind_b = ind_b, ind_c = ind_c, ind_d = ind_d, ind_e = ind_e, ind_usable = ind_usable)
        for (rel in 1:length(rel.ind)) {
          for (g in 1:k) {
            res.mat[rel, g] <- get_tdv_newcol(m = m, mat_cur = mat_cur, newcol = rel.ind[rel], g = g, p_a_u = p_a_u, present.in.groups = present.in.groups, currDV = currDV, DVchanges = DVchanges, ind_0 = ind_0, ind_1 = ind_1, ns = ns)
          }
        }
        #res.mat <- round(res.mat, 10)
        res.mat[res.mat < stats::quantile(res.mat, thr)] <- 0
        res.mat.cum <- res.mat
        res.mat.cum[] <- cumsum(as.vector(res.mat))
        #set.seed(1)
        prob <- stats::runif(1, 0, res.mat.cum[length(rel.ind), k]) #CHECK IF '/ns' could be removed previously
        ind.aux <- which(res.mat.cum >= prob, arr.ind=TRUE)[1,]
        par.GRASP[rel.ind[ind.aux[1]]] <- ind.aux[2]

        #update mat_cur
        newcol <- rel.ind[ind.aux[1]]
        g <- ind.aux[2]
        #auxiliary indices to update "e" and "usable"
        ind.for.e <- mat_cur[usable.row, ind_a[g]] == 0 & m[usable.row, newcol] == 1 #previously absent in group g, but will be present as is present in newcol
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
        mat_cur[usable.row, ind_c[-g]] <- mapply(aux_function_c, a_wg = mat_cur[usable.row, ind_a[g]], newrel = m[usable.row, newcol], c_og = mat_cur[usable.row, ind_c[-g]], b_wg = mat_cur[usable.row, ind_b[g]])
        #updating parameter "d"
        mat_cur[usable.row, ind_d[-g]] <- mat_cur[usable.row, ind_d[-g]] + 1
        #updating "c/d"
        mat_cur[p_a_u, ind_cd[-g]] <- mat_cur[p_a_u, ind_c[-g]] / mat_cur[p_a_u,ind_d[-g]]
        #within group g
        #updating parameter "a"
        mat_cur[p_a_u, ind_a[g]] <- mat_cur[p_a_u, ind_a[g]] + m[p_a_u, newcol]
        #updating parameter "b"
        mat_cur[usable.row, ind_b[g]] <- mat_cur[usable.row, ind_b[g]] + 1
        #updating "a/b"
        mat_cur[p_a_u, ind_ab[g]] <- mat_cur[p_a_u,ind_a[g]] / mat_cur[p_a_u,ind_b[g]]
      }
      #print(sum(p_a_u))
      #print(sum(mat_cur[p_a_u, ind_ab] * mat_cur[p_a_u, ind_cd] / mat_cur[p_a_u, ind_e]^2)/ns)

      #OLD SLOWER APPROACH TO DELETE
      # while (0 %in% par.GRASP) {
      #   rel.ind <- which(par.GRASP==0)
      #   res.mat <- matrix(0, length(rel.ind), k)
      #   for (rel in 1:length(rel.ind)) {
      #     for (g in 1:k) {
      #       p.temp <- par.GRASP
      #       p.temp[rel.ind[rel]] <- g
      #       res.mat[rel, g] <- tdv.constructive(m[, p.temp>0], p.temp[p.temp>0], k, index, consider.zeros)
      #     }
      #   }
      #   res.mat[res.mat < stats::quantile(res.mat, thr)] <- 0
      #   res.mat.cum <- res.mat
      #   res.mat.cum[] <- cumsum(as.vector(res.mat))
      #   prob <- stats::runif(1, 0, res.mat.cum[length(rel.ind), k])
      #   ind.aux <- which(res.mat.cum > prob, arr.ind=TRUE)[1,] #CHEK IF THIS SHOULD BE >= INSTEAD
      #   par.GRASP[rel.ind[ind.aux[1]]] <- ind.aux[2]
      # }

      RES.GRASP[[i]] <- list(par=par.GRASP, tdv=tdv(m, par.GRASP, index, output.type = "fast"))

      #Final SANN step
      RES.SANN[[i]] <- SimulAnne_optim_tdv(m = m, k = k, index = index, p.initial = RES.GRASP[[i]]$par, n.runs = 1, T_inic = T_inic, T_final = T_final, alpha = alpha, n.iter = n.iter, use.GRASP = FALSE, full.output = full.output)$SANN[[1]]
    }
    if (full.output) {
      return(list(improved.seeds = RES.seeds, GRASP = RES.GRASP, SANN = RES.SANN))
    }
    if (!full.output) {
      return(list(SANN = RES.SANN))
    }
  } #if use.GRASP == FALSE

  #simple SANN on an initial partition
  RES.SANN <- list()
  if (is.null(p.initial)) {stop("Argument p.initial can not be NULL when use.GRASP = FALSE.")}

  if (p.initial[1] != "random") {
    if (!identical(length(p.initial), nr)) {stop("Object p.initial must be a partition of the columns of m")}
    mode(p.initial) <- "integer"
    if (!identical(sort(unique(p.initial)), 1:k)) {stop("Object p.initial is not a valid partition of the columns of m")}
    p.ini <- p.initial
  }

  #Cooling schedule (CooS)
  #l <-  floor(log(T_final)/log(T_inic * (1-alpha)))
  #CooS <- floor(n.iter/l*(1:l))

  #T_final = (1 - alpha)^(n.iter/nt - 1) * T_inic
  #nt =  (n.iter * log(1 - alpha)) / (log((1 - alpha) * T_final / T_inic))
  nt <-  (n.iter * log(1 - alpha)) / (log((1 - alpha) * T_final / T_inic))
  l <- floor(n.iter/nt)
  CooS <- floor(n.iter/l*(1:l))

  for (i in 1:n.runs) {
    if (p.initial[1] == "random") {
      p.ini <- sample(c(1:k, sample(k, nr-k, replace=TRUE)))
    }
    best.p <- cur.p <- p.ini
    best.tdv <- cur.tdv <- tdv.ini <- tdv(m, p.ini, index, output.type = "fast")
    Temp <- T_inic
    res.cur.tdv <- NULL
    res.alt.tdv <- NULL
    res.probabi <- NULL
    res.tempera <- NULL

    for (iter in 1:n.iter) {
      alt.p <- random.neighbour(p = cur.p, nr = nr, k = k)
      alt.tdv <- tdv(m, alt.p, index, output.type = "fast")
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
        #alt.tdv <- alt.tdv * 100 #TEST
        #cur.tdv <- cur.tdv * 100 #TEST

        exp((3-3.5)/0.1)

        probabi <- exp((alt.tdv-cur.tdv)/Temp)
        if (stats::runif(1) < probabi) {
          cur.p <- alt.p
          cur.tdv <- alt.tdv
        }
        #alt.tdv <- alt.tdv / 100 #TEST
        #cur.tdv <- cur.tdv / 100 #TEST
      }
      if (iter %in% CooS) {Temp <- (1-alpha)*Temp}
      res.probabi <- c(res.probabi, probabi)
      res.tempera <- c(res.tempera, Temp)
    }
    if (full.output) {
      RES.SANN[[i]] <- list(par = best.p, tdv = best.tdv, cur.tdv = res.cur.tdv, alt.tdv = res.alt.tdv, probability = res.probabi, temperature = res.tempera)
    }
    if (!full.output) {
      #RES.SANN[[i]] <- list(par = best.p, tdv = best.tdv)
      if (i == 1) {RES.SANN[[1]] <- list(par = best.p, tdv = best.tdv)} else {
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
  ind.order <- order(sapply(RES.SANN, function (x) {x$tdv}))
  return(list(SANN = RES.SANN))
}

###
### random.neighbour
###

#Auxiliary function for simulated annealing
#CONFIRMAR SE FUNÇÃO FALHA COM min.g.size = 1!!!!

random.neighbour <- function (p, nr, k) {
  neigh.size <- 1 #eliminar (ou implementar poder ser variável)
  min.g.size <- 1 #se não falhar, eliminar
  tp <- table(p)
  col.troca <- NULL
  for (i in 1:neigh.size) { #eliminar (ou implementar poder ser variável)
    k.int <- which(tp > min.g.size)
    col.troca <- c(col.troca, sample((1:nr)[p %in% k.int], 1))
    tp[p[col.troca]] <- tp[p[col.troca]] - 1 #eliminar (ou implementar poder ser variável)
  }
  p[col.troca] <- sample(1:k, length(col.troca), replace=TRUE) #CF ESTE REPLACE QUE PUS MAIS TARDE!
  return(p)
}

####
#### SANN_seeds_index
####

#Special function to improve random seeds using SANN

SANN_seeds_index <- function(m, k, nr, ns, T_inic, T_final, alpha, n.iter) {
  p.ini <- sample(1:nr, k)
  best.p <- cur.p <- p.ini
  best.tdv <- cur.tdv <- tdv.ini <- objective_f_to_tdv(objective_f(m[,p.ini], k=k), k=k)/ns

  Temp <- T_inic
  res.cur.tdv <- NULL
  res.alt.tdv <- NULL
  res.probabi <- NULL
  res.tempera <- NULL

  #Temperature reduction scheme (CooS)
  #l <-  floor(log(T_final)/log(T_inic * (1-alpha)))
  #CooS <- floor(n.iter/l*(1:l))

  #T_final = (1 - alpha)^(n.iter/nt - 1) * T_inic
  #nt =  (n.iter * log(1 - alpha)) / (log((1 - alpha) * T_final / T_inic))
  nt <-  (n.iter * log(1 - alpha)) / (log((1 - alpha) * T_final / T_inic))
  l <- floor(n.iter/nt)
  CooS <- floor(n.iter/l*(1:l))

  for (iter in 1:n.iter) {
    alt.p <- random.neighbour_seeds_index(p = cur.p, nr = nr)
    alt.tdv <- objective_f_to_tdv(objective_f(m[,alt.p], k=k), k=k)/ns
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
    if (iter %in% CooS) {Temp <- (1-alpha)*Temp}

    res.probabi <- c(res.probabi, probabi)
    res.tempera <- c(res.tempera, Temp)
  }
  return(list(par = best.p, partial.tdv = best.tdv, res.cur.partial.tdv = res.cur.tdv, res.alt.partial.tdv = res.alt.tdv, res.probabi = res.probabi, res.tempera = res.tempera))
}

####
#### random.neighbour_seeds_index
####

#Auxiliary function
random.neighbour_seeds_index <- function (p, nr) {
  col_out <- sample(p, 1)
  col_in <- sample((1:nr)[-p], 1)
  p.new <- p[-which(p %in% col_out)]
  p.new <- c(p.new, col_in)
  return(p.new)
}

######
###### objective_f
######

#simplified form of tdv, for k relevés as if each one is in a different group

#kcol <- DATAMATRIX[,c(7,2,6)]

#objective_f <- function(kcol, k) { #in this function k is not used (if is the final function delete k and k=k in SANN_seeds_index)
#	p_i <- rowSums(kcol)
#	p_i <- p_i[p_i!=0]
#	goal <- sum(1/p_i)
#	return(goal)
#}

objective_f <- function(kcol, k) {
  p_i <- rowSums(kcol)
  p_i <- p_i[p_i!=0]
  goal <- sum(k/p_i-1)
  return(goal)
}

######
###### objective_f_to_tdv
######

#CONFIRMAR O CÁLCULO DO tdv #FUNDIR COM A DE CIMA!

objective_f_to_tdv <- function(goal, k) {
  tdv <- (1/(k-1))*goal
  return(tdv)
}


######
###### tdv.constructive
######

#function to calculate tdv taking into account lines of zeros (not expected in the tdv() function), which is needed for constructive optimization procedures such as GRASP

tdv.constructive <- function (m, p, k, index = "TotDiffVal", consider.zeros = TRUE) {
  m.clean <- maniphyt::clean_empty_rows(m)
  mt <- t(m.clean)
  ns <- ncol(mt) #no. of taxa
  tp <- tabulate(p)
  tpb <- length(p)-tp
  ad <- af <- matrix(0L, k, ns)
  afg <- rowsum(mt, group=p) #no. of relevés containing the taxon, within each group (absolute frequency in each group)
  t1 <- (afg == 0)*tp #no. of relevés of each group, when the taxon is not present
  t2 <- colSums(afg>0) #no. of groups containing the taxon [e]
  t3 <- t2 > 1 #indices of the taxa occurring in more than one group (they must occur in at least one)
  for (i in 1:k) { #fills matrices ad [c/d] and af [a/b], only when the taxon is present in the group!
    spe.i <- afg[i,] > 0 #indices of the taxa present in the group (in any of its columns)
    ad[i, spe.i & !t3] <- 1 #ad is 1 for the taxa occurring in one group only!
    if (sum(t3) > 0) { #if taxa occurring in more than one group do exist
      t4 <- spe.i & t3 #taxa of group k occurring also in other group ≠k
      ad[i, t4] <- colSums(t1[-i, t4, drop = FALSE]/tpb[i])
    }
    af[i, spe.i] <- afg[i, spe.i]/tp[i]
  }

  if (consider.zeros) { #if consider.zeros = TRUE, ns is changed by this line, then index calculation considers the possible lines of zeros (yet it implicitly assumes that TDV or TDV1 is zero for that line, yet it is in practice NaN)
    ns <- nrow(m)
  }  #if consider.zeros = FALSE ns is unchanged it doesn't considers the possible lines of zeros

  if (index == "TotDiffVal") {
      return(sum(colSums(af*ad)/t2)/ns)
  }
  if (index == "TotDiffVal1") {
    return(sum(colSums(af*ad)/t2^2)/ns)
  }
}

###
### Auxiliary functions to decrease GRASP computation time
###

### aux_function_c; aux_function_c_if0; aux_function_c_if1

#auxiliary function to recalculate parameter "c" (outside group g) given a new column relevé to be included in group g
#these functions could be adapted to the use of function ifelse(); apparently the gain in time is only for smaller datasets, as mapply() returned faster for very big datasets.

#a_wg   the "a" parameter within the group g
#newrel the new relevé to be included
#c_og   the "c" parameter, outside the group g #this vector could be longer, the other vectors are recycled by the mapply function
#b_wg   the "b" parameter, within the group g

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

aux_function_c_if0 <- function(a_wg, c_og) { #simpler function, to use when newrel is a vector of 0s
  if (a_wg == 0) {
    return (c_og + 1)
  } else {
    return(c_og)
  }
}

aux_function_c_if1 <- function(a_wg, c_og, b_wg) {  #simpler function, to use when newrel is a vector of 1s
  if (a_wg == 0) {
    return (c_og - b_wg)
  } else {
    return(c_og)
  }
}

### get_DV_01

#this function uses current calculation matrix (mat_cur) to obtain, for each (usable) row (and for each group g!), the result of DiffVal by introducing a new relevé. For each row and for each group the DiffVal is presented in two columns considering that new relevé brings a 0 or a 1 to that row.

#usable.row      the lines of mat_cur that are still useful for the TDV calculation (i.e. the respective taxa is not in all groups), but the respective taxa could still be absent from all groups
#p_a_u            (present and usable) the lines of mat_cur that are not null and that are still useful for the TDV calculation (i.e. the respective taxa is not in all groups)
#ns

#ind_0 <- (1:k-1)*k+1 #k indices to store values when adding 0 values
#ind_1 <- (1:k-1)*k+2 #k indices to store values when adding 1 values

get_DV_01 <- function(m, k, mat_cur, usable.row, p_a_u, ns, ind_0, ind_1, ind_a, ind_b, ind_c, ind_d, ind_e, ind_usable) {
  res.mat.01 <- matrix(0, ns, k*2) #to keep DiffVal for each group (depending if the new relevé to enter the group has a 0 or 1 in the row)
  res.mat.02 <- matrix(NA, ns, k) #to keep p_a_u_new for each group
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
    #present.in.groups_new <- rep(TRUE, ns) #all will have at least presence in group g (i.e. all taxa not absent from all groups)
    p_a_u_new <- usable.row_new #as present.in.groups_new is all TRUE it is the same as usable.row_new (present and usable) taxa present in at least one group but not present in all of them

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

#auxiliary function for GRASP efficiency (based on get_DV_01)

#this function uses current mat_cur to calculate TDV efficiently (given a newcol)

#newcol           the index of the column of m (the relevé) to be included in group g
#g                the group (of the partition) where the new relevé is to be included
#usable.row      the lines of mat_cur that are still useful for the TDV calculation (i.e. the respective taxa is not in all groups), but the respective taxa could still be absent from all groups
#p_a_u            the lines of mat_cur that are not null and that are still useful for the TDV calculation (i.e. the respective taxa is not in all groups)
#ns               the number ot rows of m (only needed to calculate TDV "with zeros", i.e. taking into account the empty lines)

#ind_0 <- (1:k-1)*k+1 #k indices to store values when adding 0 values
#ind_1 <- (1:k-1)*k+2 #k indices to store values when adding 1 values

#currDV <- rep(0, ns)
#currDV[p_a_u] <- sum(mat_cur[p_a_u, ind_ab] * mat_cur[p_a_u, ind_cd] / mat_cur[p_a_u, ind_e])
#DVchanges <- get_DV_01(m, mat_cur, usable.row, p_a_u, ns, ind_0, ind_1)

get_tdv_newcol <- function(m, mat_cur, newcol, g, p_a_u, present.in.groups, currDV, DVchanges, ind_0, ind_1, ns) {
  newrel <- as.logical(m[, newcol])
  #update when 0
  currDV[!newrel & p_a_u] <- DVchanges[[1]][, ind_0[g]][!newrel & p_a_u]
  #update when 1
  p_a_u_new <- DVchanges[[2]][,g]
  p_a_u_final <- p_a_u
  p_a_u_final[newrel] <- p_a_u_new[newrel]
  currDV[newrel & p_a_u_new] <- DVchanges[[1]][, ind_1[g]][newrel & p_a_u_new]

  present.in.groups_final <- present.in.groups | newrel

  return(sum(currDV[p_a_u_final]) / sum(present.in.groups_final)) #Attention this is not divided by ns (as it should if TotDiffVal is needed)
}






######
###### New approach to the constructive phase of GRASP, based on mat_cur (a matrix to store and improve the calculation of necessary tdvs)
######


#
# ns <- nrow(m)
# nr <- ncol(m)
# k <- 3
# seeds <- 1:k
# seed <- 1:k
#
# mat_cur <- matrix(0, ns, 6 * k + 2)
# ind_a <- 1 + 0:(k-1) * 6
# ind_b <- ind_a + 1
# ind_c <- ind_b + 1
# ind_d <- ind_c + 1
# ind_ab <- ind_d + 1
# ind_cd <- ind_ab + 1
# ind_e <- k * 6 + 1
# ind_usable <- ind_e + 1
#
# presences.in.groups <- rowSums(m[,seeds])
# absences.in.groups <- k - presences.in.groups
#
# #for the special case of the seed (i.e. only one relevé in each of the groups):
# mat_cur[,ind_a] <- m[,seeds]
# mat_cur[,ind_b] <- 1
# mat_cur[,ind_c] <- absences.in.groups - !m[,seeds] #CF. THIS
# mat_cur[,ind_d] <- k-1
# mat_cur[,ind_ab] <- mat_cur[,ind_a] #at this stage mat_cur[,ind_b] is always 1
# mat_cur[,ind_cd] <- mat_cur[,ind_c]/mat_cur[,ind_d]
# mat_cur[,ind_e] <- presences.in.groups
# mat_cur[,ind_usable] <- as.numeric(mat_cur[,ind_e] != k)
#
# present.in.groups <- presences.in.groups != 0 #present in at least one group (i.e. not absent)
# usable.row <- as.logical(mat_cur[,ind_usable]) #taxa not present in all groups (as taxa that is in all groups have DiffVal = 0)
# p_a_u <- present.in.groups & usable.row #(present and usable) present taxa in at least in one group but not in all of them
#
# thr <- 0.95
# par.GRASP <- rep(0, nr)
# #seed <- RES.seeds[[i]]$seed
# par.GRASP[seed] <- 1:k
#
#
# new.f2 <- function(m, k, par.GRASP, mat_cur, usable.row, p_a_u, present.in.groups, thr, ns) {
#   #while (0 %in% par.GRASP & any(as.logical(mat_cur[,ind_usable]))) { #do while there is an empty column in the partition and at least one line is "usable" (CF!)
#   ind_0 <- (1:k-1)*2+1 #k indices to store values when adding 0 values
#   ind_1 <- (1:k-1)*2+2 #k indices to store values when adding 1 values
#
#   while (0 %in% par.GRASP) { #do while there is an empty column in the partition (what happen if all lines are not usable? EXPERIMENTAR COM MATRIZ SÓ DE 1s?)
#     rel.ind <- which(par.GRASP==0)
#     res.mat <- matrix(0, length(rel.ind), k)
#     currDV <- rep(0, ns)
#     if (all(!usable.row)) {warning("Attention, there are no usable taxa for TDV calculation!")}
#     currDV[p_a_u] <- rowSums(mat_cur[p_a_u, ind_ab, drop = FALSE] * mat_cur[p_a_u, ind_cd, drop = FALSE] / mat_cur[p_a_u, ind_e]^2)
#     DVchanges <- get_DV_01(m, k, mat_cur, usable.row, p_a_u, ns, ind_0, ind_1)
#     for (rel in 1:length(rel.ind)) {
#       for (g in 1:k) {
#         res.mat[rel, g] <- get_tdv_newcol(m = m, mat_cur = mat_cur, newcol = rel.ind[rel], g = g, p_a_u = p_a_u, present.in.groups = present.in.groups, currDV = currDV, DVchanges = DVchanges, ind_0 = ind_0, ind_1 = ind_1, ns = ns)
#       }
#     }
#     #res.mat <- round(res.mat, 10)
#     res.mat[res.mat < stats::quantile(res.mat, thr)] <- 0
#     res.mat.cum <- res.mat
#     res.mat.cum[] <- cumsum(as.vector(res.mat))
#     #set.seed(1)
#     prob <- stats::runif(1, 0, res.mat.cum[length(rel.ind), k])
#     ind.aux <- which(res.mat.cum > prob, arr.ind=TRUE)[1,] #CHEK IF THIS SHOULD BE >= INSTEAD
#     par.GRASP[rel.ind[ind.aux[1]]] <- ind.aux[2]
#
#     #update mat_cur
#     newcol <- rel.ind[ind.aux[1]]
#     g <- ind.aux[2]
#     #auxiliary indices to update "e" and "usable"
#     ind.for.e <- mat_cur[usable.row, ind_a[g]] == 0 & m[usable.row, newcol] == 1 #previously absent in group g, but will be present as is present in newcol
#     ind.for.usable <- ind.for.e & mat_cur[usable.row, ind_e] == k - 1 #previously absent in group g, but will be present as is present in newcol, and g was the only empty group before
#     #updating parameter "e"
#     mat_cur[usable.row, ind_e][ind.for.e] <- mat_cur[usable.row, ind_e][ind.for.e] + 1
#     #updating "usable"
#     mat_cur[usable.row, ind_usable][ind.for.usable] <- 0
#     #updating "usable.row" and "p_a_u"
#     usable.row <- as.logical(mat_cur[,ind_usable]) #taxa not present in all groups (taxa that are present in all groups have DiffVal = 0)
#     present.in.groups <- mat_cur[, ind_e] != 0 #taxa not absent from all groups
#     p_a_u <- present.in.groups & usable.row #taxa present in at least one group but not present in all of them
#     #outside group g
#     #updating parameter "c"
#     mat_cur[usable.row, ind_c[-g]] <- mapply(aux_function_c, a_wg = mat_cur[usable.row, ind_a[g]], newrel = m[usable.row, newcol], c_og = mat_cur[usable.row, ind_c[-g]], b_wg = mat_cur[usable.row, ind_b[g]])
#     #updating parameter "d"
#     mat_cur[usable.row, ind_d[-g]] <- mat_cur[usable.row, ind_d[-g]] + 1
#     #updating "c/d"
#     mat_cur[p_a_u, ind_cd[-g]] <- mat_cur[p_a_u, ind_c[-g]] / mat_cur[p_a_u,ind_d[-g]]
#     #within group g
#     #updating parameter "a"
#     mat_cur[p_a_u, ind_a[g]] <- mat_cur[p_a_u, ind_a[g]] + m[p_a_u, newcol]
#     #updating parameter "b"
#     mat_cur[usable.row, ind_b[g]] <- mat_cur[usable.row, ind_b[g]] + 1
#     #updating "a/b"
#     mat_cur[p_a_u, ind_ab[g]] <- mat_cur[p_a_u,ind_a[g]] / mat_cur[p_a_u,ind_b[g]]
#   }
#   print(sum(p_a_u))
#   print(sum(mat_cur[p_a_u, ind_ab] * mat_cur[p_a_u, ind_cd] / mat_cur[p_a_u, ind_e]^2)/ns)
#   return(list(par.GRASP, mat_cur))
# }


















#
#
#
# #COM CICLO WHILE
#
# old.f <- function(m, k, par.GRASP, res.mat, thr) {
#   while (0 %in% par.GRASP) { #do while there is an empty column in the partition
#     rel.ind <- which(par.GRASP==0)
#     res.mat <- matrix(0, length(rel.ind), k)
#     for (rel in 1:length(rel.ind)) {
#       for (g in 1:k) {
#         p.temp <- par.GRASP
#         p.temp[rel.ind[rel]] <- g
#         res.mat[rel, g] <- tdv.constructive(m[, p.temp>0], p.temp[p.temp>0], k, index = "TotDiffVal", consider.zeros = FALSE)
#       }
#     }
#     res.mat <- round(res.mat, 10)
#     res.mat[res.mat < stats::quantile(res.mat, thr)] <- 0
#     res.mat.cum <- res.mat
#     res.mat.cum[] <- cumsum(as.vector(res.mat))
#     prob <- stats::runif(1, 0, res.mat.cum[length(rel.ind), k])
#     ind.aux <- which(res.mat.cum > prob, arr.ind=TRUE)[1,] #CHEK IF THIS SHOULD BE >= INSTEAD
#     par.GRASP[rel.ind[ind.aux[1]]] <- ind.aux[2]
#   }
#   return(par.GRASP)
# }
#
# new.f1 <- function(m, k, par.GRASP, mat_cur, usable.row, p_a_u, thr, ns) {
#   #while (0 %in% par.GRASP & any(as.logical(mat_cur[,ind_usable]))) { #do while there is an empty column in the partition and at least one line is "usable" (CF!)
#   while (0 %in% par.GRASP) { #do while there is an empty column in the partition (what happen if all lines are not usable? EXPERIMENTAR COM MATRIZ SÓ DE 1s?)
#     rel.ind <- which(par.GRASP==0)
#     res.mat <- matrix(0, length(rel.ind), k)
#     for (rel in 1:length(rel.ind)) {
#       for (g in 1:k) {
#         res.mat[rel, g] <- get_tdv_newcol1(m = m, mat_cur = mat_cur, newcol = rel.ind[rel], g = g, usable.row = usable.row, p_a_u = p_a_u, ns = ns)
#       }
#     }
#     #res.mat <- round(res.mat, 10)
#     res.mat[res.mat < stats::quantile(res.mat, thr)] <- 0
#     res.mat.cum <- res.mat
#     res.mat.cum[] <- cumsum(as.vector(res.mat))
#     #set.seed(1)
#     prob <- stats::runif(1, 0, res.mat.cum[length(rel.ind), k])
#     ind.aux <- which(res.mat.cum > prob, arr.ind=TRUE)[1,] #CHEK IF THIS SHOULD BE >= INSTEAD
#     par.GRASP[rel.ind[ind.aux[1]]] <- ind.aux[2]
#
#     #update mat_cur
#     newcol <- rel.ind[ind.aux[1]]
#     g <- ind.aux[2]
#     #auxiliary indices to update "e" and "usable"
#     ind.for.e <- mat_cur[usable.row, ind_a[g]] == 0 & m[usable.row, newcol] == 1 #previously absent in group g, but will be present as is present in newcol
#     ind.for.usable <- ind.for.e & mat_cur[usable.row, ind_e] == k - 1 #previously absent in group g, but will be present as is present in newcol, and g was the only empty group before
#     #updating parameter "e"
#     mat_cur[usable.row, ind_e][ind.for.e] <- (mat_cur[usable.row, ind_e][ind.for.e] + 1)
#     #updating "usable"
#     mat_cur[usable.row, ind_usable][ind.for.usable] <- 0
#     #updating "usable.row" and "p_a_u"
#     usable.row <- as.logical(mat_cur[,ind_usable]) #taxa not present in all groups (taxa that are present in all groups have DiffVal = 0)
#     present.in.groups <- mat_cur[, ind_e] != 0 #taxa not absent from all groups
#     p_a_u <- present.in.groups & usable.row #taxa present in at least one group but not present in all of them
#     #outside group g
#     #updating parameter "c"
#     mat_cur[p_a_u, ind_c[-g]] <- mapply(aux_function_c, a_wg = mat_cur[p_a_u, ind_a[g]], newrel = m[p_a_u, newcol], c_og = mat_cur[p_a_u, ind_c[-g]], b_wg = mat_cur[p_a_u, ind_b[g]])
#     #updating parameter "d"
#     mat_cur[p_a_u, ind_d[-g]] <- (mat_cur[p_a_u, ind_d[-g]] + 1)
#     #updating "c/d"
#     mat_cur[p_a_u, ind_cd[-g]] <- (mat_cur[p_a_u,ind_c[-g]] / mat_cur[p_a_u,ind_d[-g]])
#     #within group g
#     #updating parameter "a"
#     mat_cur[p_a_u, ind_a[g]] <- (mat_cur[p_a_u, ind_a[g]] + m[p_a_u, newcol])
#     #updating parameter "b"
#     mat_cur[p_a_u, ind_b[g]] <- (mat_cur[p_a_u, ind_b[g]] + 1)
#     #updating "a/b"
#     mat_cur[p_a_u, ind_ab[g]] <- (mat_cur[p_a_u,ind_a[g]] / mat_cur[p_a_u,ind_b[g]])
#   }
#   print(sum(mat_cur[p_a_u, ind_ab] * mat_cur[p_a_u, ind_cd] / mat_cur[p_a_u, ind_e]^2)/ns)
#   return(par.GRASP)
# }


#this function uses current mat_cur to calculate TDV efficiently (given a newcol)

#newcol           the index of the column of m (the relevé) to be included in group g
#g                the group (of the partition) where the new relevé is to be included
#usable.row      the lines of mat_cur that are still useful for the TDV calculation (i.e. the respective taxa is not in all groups), but the respective taxa could still be absent from all groups
#p_a_u            the lines of mat_cur that are not null and that are still useful for the TDV calculation (i.e. the respective taxa is not in all groups)
#ns               the number ot rows of m (only needed to calculate TDV "with zeros", i.e. taking into account the empty lines)

# get_tdv_newcol1 <- function(m, mat_cur, newcol, g, usable.row, p_a_u, ns) {
#   #auxiliary indices to update "e" and "usable"
#   ind.for.e <- mat_cur[usable.row, ind_a[g]] == 0 & m[usable.row, newcol] == 1 #previously absent in group g, but will be present as is present in newcol
#   ind.for.usable <- ind.for.e & mat_cur[usable.row, ind_e] == k - 1 #previously absent in group g, but will be present as is present in newcol, and g was the only empty group before
#
#   #updating parameter "e"
#   mat_cur[usable.row, ind_e][ind.for.e] <- (mat_cur[usable.row, ind_e][ind.for.e] + 1)
#
#   #updating "usable"
#   mat_cur[usable.row, ind_usable][ind.for.usable] <- 0
#
#   #updating "usable.row" and "p_a_u"
#   usable.row <- as.logical(mat_cur[,ind_usable]) #taxa not present in all groups (taxa that are present in all groups have DiffVal = 0)
#   present.in.groups <- mat_cur[, ind_e] != 0 #taxa not absent from all groups
#   p_a_u <- present.in.groups & usable.row #taxa present in at least one group but not present in all of them
#
#   #outside group g
#   #updating parameter "c"
#   mat_cur[p_a_u, ind_c[-g]] <- mapply(aux_function_c, a_wg = mat_cur[p_a_u, ind_a[g]], newrel = m[p_a_u, newcol], c_og = mat_cur[p_a_u, ind_c[-g]], b_wg = mat_cur[p_a_u, ind_b[g]])
#   #updating parameter "d"
#   mat_cur[p_a_u, ind_d[-g]] <- (mat_cur[p_a_u, ind_d[-g]] + 1)
#   #updating "c/d"
#   mat_cur[p_a_u, ind_cd[-g]] <- (mat_cur[p_a_u,ind_c[-g]] / mat_cur[p_a_u,ind_d[-g]])
#
#   #within group g
#   #updating parameter "a"
#   mat_cur[p_a_u, ind_a[g]] <- (mat_cur[p_a_u, ind_a[g]] + m[p_a_u, newcol])
#   #updating parameter "b"
#   mat_cur[p_a_u, ind_b[g]] <- (mat_cur[p_a_u, ind_b[g]] + 1)
#   #updating "a/b"
#   mat_cur[p_a_u, ind_ab[g]] <- (mat_cur[p_a_u,ind_a[g]] / mat_cur[p_a_u,ind_b[g]])
#
#   #tdv
#
#   return(sum(mat_cur[p_a_u, ind_ab] * mat_cur[p_a_u, ind_cd] / mat_cur[p_a_u, ind_e]) / sum(present.in.groups)) #without zeros
#   #return(sum(mat_cur[p_a_u, ind_ab] * mat_cur[p_a_u, ind_cd] / mat_cur[p_a_u, ind_e]) / ns) #with zeros
# }
