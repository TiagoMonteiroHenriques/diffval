# HillClimb_optim_tdv.R
#'
#' @title TotDiffVal (or TotDiffVal1) optimization using hill climbing algorithms
#'
#' @description This function searches for partitions of the columns of a given matrix, optimizing the TotDiffVal (or TotDiffVal1) index.
#'
#' @param m A `matrix`, i.e. a phytosociological table of 0s (absences) and 1s (presences), where rows correspond to taxa and columns correspond to relevés.
#' @param p.initial A `vector` of integer numbers with the initial partition of the relevés (i.e. a vector with values from 1 to k, with length equal to the number of columns of m, ascribing each relevé to one of the k groups). By default, "random", generates a random initial partition.
#' @param k A `numeric`, giving the number of desired groups.
#' @param n.starts A `numeric`, giving the number of starts to perform.
#' @param n.sol A `numeric`, giving the number of best solutions to keep in the final output.
#' @param index A `character`, selecting wich index to optimize "TotDiffVal", or "TotDiffVal1". #STILL TODO
#' @param maxit A `numeric` giving the number of iterations of the hill-climbing optimization.
#' @param min.g.size A `numeric` The minimum number of relevés that a group can contain (must be 2 or higher).
#' @param random.first A `logical`. `FALSE` (the default), performs only hill-climbing on the 1-neighbours; `TRUE` first, performs a stochastic hill-climbing on random `n`-neighbours (`n` is definded by the parameter `rf.neigh.size`), and only after runs the hill-climbing search on the 1-neighbours; see description above.
#' @param rf.neigh.size	A `numeric`, giving the size (`n`) of the `n`-neighbour for the stochastic hill-climbing; only used if `random.first` = `TRUE`.
#' @param rf.maxit A `numeric`, giving the number of iterations of the hill-climbing optimization; only used if `random.first` = `TRUE`.
#' @param full.output A `logical`. If `FALSE` (the default) the best `n.sol` partitions and respective indices are returned. If `TRUE` (only available for `n.sols = 1`) the output will also contain information on the optimization steps (see below).
#'
#' @details Given a phytosociological table (`m`, rows corresponding to taxa and columns corresponding to relevés) this function searches for a k-partition (k, defined by the user) optimizing TotDiffVal (or TotDiffVal1) index (see http://home.isa.utl.pt/~tmh/), i.e. searches, using a hill-climbing algorithm, for patterns of differential taxa by rearranging the relevés into k groups.
#'
#' Optimization can start from a random partition (`p.ini` = "random"), or from a given partition (`p.ini`, defined by the user or produced by any clustering method, or even a manual classification).
#'
#' Each iteration searches for a TotDiffVal1 improvement screening all 1-neighbours, until the given number of maximum iterations (`maxit`) is reached. A 1-neighbour of a given partition is another partition obtained by ascribing 1 relevé (of the original partition) to a different group.  A n-neighbour is obtained, equivalently, ascribing n relevés to different groups.
#'
#' Optionally, a faster search (stochastic hill-climbing) can be performed in a first moment (`random.first` = `TRUE`), consisting on searching for TotDiffVal1 improvements, by randomly selecting n-neighbours (n defined by the user with the parameter `rf.neigh.size`), until a given number of maximum iterations (`rf.maxit`) is reached. Stochastic hill-climbing might be helpful for big tables (where the simple screening of all 1-neighbours might be too time consuming).
#'
#' Several runs of HillClimb_optim_tdv (multi-starts) should be tried out, as several local maxima are usually present and the hill-climbing algorithm converges easily to local maxima. Sometimes, converging to a known high-valued partition is very unlikely, being dependent on data structure and on the initial partition.
#'
#' Trimming your table by a 'constancy' range (see \code{\link{select_taxa}} function) or using the result of other cluster methodologies as input, might help finding interesting partitions. Specially after trimming the table by a 'constancy' range, getting a random initial partition with TotDiffVal1 greater than zero might be very unlike; on such cases using the result of other clustering strategies as input is useful.
#'
#' @return A `list` with the following components:
#'
#' \describe{
#'   \item{res.rf}{A `matrix` with the iteration number (of the stochastic hill-climbing phase), the maximum TotDiffVal1 found until that iteration, and the higher TotDiffVal1 among all 1-neighbours; a first line of zeros is being added at the beginning of the matrix, which should be removed in future.}
#'   \item{par.rf}{A `vector` with the best partition found in the stochastic hill-climbing phase.}
#'   \item{max.TotDiffVal1.rf}{A `numeric` showing the maximum TotDifVal1 found in the stochastic hill-climbing phase (if selected).}
#'   \item{res}{A `matrix` with the iteration number (of the hill-climbing), the maximum TotDiffVal1 found until that iteration, and the higher TotDiffVal1 among all 1-neighbours; a first line of zeros is being added at the beginning of the matrix, which should be removed in future.}
#'   \item{par}{A `vector` with the best partition found in the hill-climbing phase.}
#'   \item{local_maximum}{A `logical` indicating if `par` is a 1-neighbour local maximum.}
#'   \item{time}{The total amount of time of the run (from function \code{\link[base]{Sys.time}}).}
#'   \item{max.TotDiffVal1}{A `numeric` with the maximum TotDifVal1 found in the run.}
#' }
#'
#' @author Tiago Monteiro-Henriques. E-mail: \email{tiagomonteirohenriques@@gmail.com}.
#'
#' @export
#'
HillClimb_optim_tdv <- function(m, p.initial="random", k, n.starts = 1, n.sol = 1, index = "TotDiffVal1", maxit = 10, min.g.size = 2, random.first = FALSE, rf.neigh.size = 1, rf.maxit = 500, full.output = FALSE) {
  stopifnot(is.matrix(m))
  mode(m) <- "integer"
  if (!identical(c(0L,1L), sort(unique(as.vector(m))))) {stop("Matrix m should contain only 0's and 1's.")}
  if (min(rowSums(m))==0) {stop("At least one taxa is not present in any relev\u00e9.")}
  if (min(colSums(m))==0) {stop("At least one relev\u00e9 contains no taxa.")}
  if ((mgs <- as.integer(min.g.size))<2) {stop("Object min.g.size must be greater than 2.")} #TODO: check if this can be eliminated
  mt <- t(m)
  nr <- nrow(mt) # no. of relevés
  ns <- ncol(mt) #no. of taxa
  if (nr <= mgs*k) {stop(paste0("Random partition cannot guarantee at least ", mgs, " relev\u00e9s per group!"))}

  if (p.initial[1] != "random") {
    p.ini <- p.initial
    if (!identical(as.integer(sort(unique(p.ini))), 1:k)) {stop("Your partition doesn't have k groups!")} #maybe when p.ini is given k could be ignored
    if (!identical(length(p.ini), nr)) {stop("Object p.ini must be a partition of the columns of matrix m.")}
    if (min(tp)<mgs) {stop(paste0("At least one group of the provided partition has less than ", mgs, " elements"))}
    if (max(tp)==mgs) {stop(paste0("At least one group of the provided partition has to have more than ", mgs, " elements"))}
    tp <- table(p.ini)
  }

  if (n.sol > n.starts) {stop("The number of starts ('n.starts') should not be lower than the desired number of best solutions ('n.sol').")}

  if ((n.sol > 1 | n.starts > 1) & full.output == TRUE) {stop("The option 'full.output = TRUE' is only available for 'n.sol == 1' and 'n.starts = 1'.")}

  if (n.sol == 1 & n.starts == 1 & full.output == TRUE) {
    t <- Sys.time()

    if (p.initial[1] == "random") {
      p.ini <- sample(c(rep(1:k, mgs), sample(k, nr-mgs*k, replace=TRUE)))
      tp <- table(p.ini)
    }

    arf <- matrix(0, k, ns); colnames(arf) <- colnames(mt); rownames(arf) <- 1:k; adp <- arf #arf is the adjusted relative frequency in each group and adp is the adjusted differential proportion
    afg <- rowsum(mt, group=p.ini) #no. of relevés containing the taxon, within each group (absolute frequency in each group)
    t1 <- (afg==0)*as.vector(tp) #no. of relevés of each group, when the taxon is not present
    t2 <- colSums(afg>0) #no. of groups containing the taxon
    t3 <- t2 > 1 #indices of the taxa occurring in more than one group (they must occur in at least one)
    for (i in 1:k) { #fills matrices adp and arf, only when the taxon is present in the group!
      rel.i <- which(p.ini==i) #indices of the relevés of the group
      spe.i <- afg[i,]>0 #indices of the taxa present in the group
      adp[i,spe.i & !t3] <- 1 #adp is 1 for the taxa occurring in one group only!
      if (sum(t3) > 0) { #for taxa occurring in more than one group
        t4 <- spe.i & t3
        adp[i, t4] <- colSums(t1[-i, t4, drop=FALSE]/sum(tp[-i]))/(t2[t4])
      }
      arf[i, spe.i] <- afg[i, spe.i]/as.vector(tp[i])/(t2[spe.i])
    }
    cuscor <- sum(colSums(arf*adp))/ns

    res.rf = NULL
    par.rf = NULL
    cuscor.rf=NULL
    if (random.first==TRUE) {

      #random.neighbour: this function returns randomly one of the (k-1)*nr 1-neighbour partitions, assuring that the minimum group size (mgs) is respected
      random.neighbour <- function (p) {
        tp <- table(p)
        k.int <- which(tp>mgs) #k of interest to sample
        k.sam <- tp[k.int]-mgs #k samplable
        troca <- sample(rep(k.int,k.sam), min(rf.neigh.size, sum(k.sam))) #neighbourhood.size cannot be greater than sum(k.sam)
        niter <- table(troca)
        gn.tot <- NULL
        in.tot <- NULL
        for (gv in sort(unique(troca))) {#it assumes that table function also sorts!
          if(length(c(1:k)[-gv])!=1) {gn.tot <- c(gn.tot,sample(c(1:k)[-gv],niter[paste(gv)],replace=TRUE))} else {gn.tot <- c(gn.tot,rep(c(1:k)[-gv],niter[paste(gv)]))}
          in.tot <- c(in.tot,sample(which(p==gv),niter[paste(gv)]))
        }
        p[in.tot] <- gn.tot
        return(p)
      }

      #STOCHASTIC HILL CLIMBING (random.first=TRUE)
      parcor <- p.ini
      res.rf <- c(0,0,0)
      for (iter in 1:rf.maxit) {
        parviz <- random.neighbour(parcor)
        tp <- table(parviz)
        arf <- matrix(0, k, ns); colnames(arf) <- colnames(mt); rownames(arf) <- 1:k; adp <- arf #arf is the adjusted relative frequency in each group and adp is the adjusted differential proportion
        afg <- rowsum(mt, group=parviz) #no. of relevés containing the taxon, within each group (absolute frequency in each group)
        t1 <- (afg==0)*as.vector(tp) #no. of relevés of each group, when the taxon is not present
        t2 <- colSums(afg>0) #no. of groups containing the taxon
        t3 <- t2>1 #indices of the taxa occurring in more than one group (they must occur in at least one)
        for (i in 1:k) { #fills matrices adp and arf, only when the taxon is present in the group!
          rel.i <- which(parviz==i) #indices of the relevés of the group
          spe.i <- afg[i,]>0 #indices of the taxa present in the group
          adp[i,spe.i & !t3] <- 1 #adp is 1 for the taxa occurring in one group only!
          if (sum(t3)>0) { #for taxa occurring in more than one group
            t4 <- spe.i & t3
            adp[i, t4] <- colSums(t1[-i, t4, drop=FALSE]/sum(tp[-i]))/(t2[t4])
          }
          arf[i, spe.i] <- afg[i,spe.i]/as.vector(tp[i])/(t2[spe.i])
        }
        cusviz <- sum(colSums(arf*adp))/ns
        if (cusviz >= cuscor) {parcor <- parviz; cuscor <- cusviz}
        res.rf <- rbind(res.rf,c(iter,cuscor,cusviz))
      }
      p.ini <- par.rf <- parcor
      cuscor.rf <- cuscor
      tp <- table(p.ini)
      #afg, t1, t2, t3, arf and adp must be recalculated for parcor (as a worse neighbour might have been used to calculate them before)!
      arf <- matrix(0, k, ns); colnames(arf) <- colnames(mt); rownames(arf) <- 1:k; adp <- arf #arf is the adjusted relative frequency in each group and adp is the adjusted differential proportion
      afg <- rowsum(mt, group=p.ini) #no. of relevés containing the taxon, within each group (absolute frequency in each group)
      t1 <- (afg==0)*as.vector(tp) #no. of relevés of each group, when the taxon is not present
      t2 <- colSums(afg>0) #no. of groups containing the taxon
      t3 <- t2>1 #indices of the taxa occurring in more than one group (they must occur in at least one)
      for (i in 1:k) { #fills matrices adp and arf, only when the taxon is present in the group!
        rel.i <- which(p.ini==i) #indices of the relevés of the group
        spe.i <- afg[i,]>0 #indices of the taxa present in the group
        adp[i,spe.i & !t3] <- 1 #adp is 1 for the taxa occurring in one group only!
        if (sum(t3)>0) { #for taxa occurring in more than one group
          t4 <- spe.i & t3
          adp[i, t4] <- colSums(t1[-i, t4, drop=FALSE]/sum(tp[-i]))/(t2[t4])
        }
        arf[i, spe.i] <- afg[i,spe.i]/as.vector(tp[i])/(t2[spe.i])
      }
    }

    #max.tdv.neighbour: this function returns the(one of the) 1-neighbouring partition(s) presenting greater TotDiffVal1 (tdv)
    max.tdv.neighbour <- function (p) {
      tp <- table(p); mtemp <- matrix(p, nr, nr, byrow=TRUE); fordiag <- p; res1 <- NULL; res2 <- NULL
      if (min(tp)>mgs) { #all groups have more than mgs elements
        for (i in 1:(k-1)) {
          fordiag <- fordiag+1
          fordiag[which(fordiag==k+1)] <- 1
          diag(mtemp) <- fordiag
          res1 <- rbind(res1,mtemp)
          res2 <- c(res2,fordiag)
        }
        res2 <- cbind(rep(p, k-1),res2)
        colnames(res2) <- NULL
      } else { #at least one group has only mgs elements
        ind.rm <- as.numeric(sapply(which(tp==mgs), function (x) {which(p==x)})) #it returns the indices of the partition corresponding to groups presenting only mgs elements #as.numeric is probably not necessary, it is possibly done by default in []
        for (i in 1:(k-1)) {
          fordiag <- fordiag+1
          fordiag[which(fordiag==k+1)] <- 1
          diag(mtemp) <- fordiag
          res1 <- rbind(res1,mtemp[-ind.rm,])
          res2 <- c(res2,fordiag[-ind.rm])
        }
        res2 <- cbind(rep(p[-ind.rm], k-1),res2)
        colnames(res2) <- NULL
      }
      mat.neig <- list(p.list=res1,pairs=res2)

      mat.neig.tdv <- sapply(1:nrow(mat.neig$p.list), function (x) { #like this, it is difficult to parallelize!
        pn <- mat.neig$p.list[x,]
        kc <- mat.neig$pairs[x,]
        tpn <- table(pn) #neighbour partition
        for (i in kc) {afg[i,] <- colSums(mt[which(pn==i),])} #(updates afg) no. of relevés containing the taxa, within each group, only for the two groups that changed a relevé!
        t5 <- afg[kc,][1,]>0 | afg[kc,][2,]>0 #taxa affected by the change of relevés
        t1[kc,] <- (afg==0)[kc,]*as.vector(tpn)[kc] #(updates t1) no. of relevés of each group, when the taxon is not present #t5 must not be used here!
        t2[t5] <- colSums(afg>0)[t5] #(new t2) no. of groups containing the taxon
        t3[t5] <- (t2>1)[t5] #indices of the taxa occurring in more than one group (they must occur in at least one)
        for (i in 1:k) { #changes matrices adp and arf, only when the taxon is present in the group!
          rel.i <- which(pn==i) #indices of the relevés of the group
          spe.i <- afg[i,]>0 #indices of the taxa present in the group
          if (sum(spe.i & t5)>0) { #caso o grupo em causa contenha espécies afetadas
            adp[i, spe.i & !t3] <- 1 #adp is 1 for the taxa occurring in one group only!
            if (sum(t3)>0) { #for taxa occurring in more than one group
              t4 <- spe.i & t3
              adp[i, t4] <- colSums(t1[-i, t4, drop=FALSE]/sum(tpn[-i]))/(t2[t4])
            }
            adp[i, t5 & !spe.i] <- 0 #inserts a zero to the affected taxa that is no more present in the group!
            arf[i, t5] <- afg[i, t5]/as.vector(tpn[i])/(t2[t5])
          }
        }
        return(sum(colSums(arf*adp))/ns)
      })
      return(list(tdv=tdv <- max(mat.neig.tdv), p=mat.neig$p.list[which(mat.neig.tdv==tdv)[1],]))
    }

    #HILL CLIMBING
    if (maxit==0) {res <- NULL; parcor=p.ini; loc_max = NA} else {
      loc_max <- NA
      parcor <- p.ini
      #cuscor is already calculated
      res <- c(0,0,0)
      for (iter in 1:maxit) {
        temp <- max.tdv.neighbour(parcor)
        parviz <- temp$p
        cusviz <- temp$tdv
        if (cusviz > cuscor) {parcor <- parviz; cuscor <- cusviz} else {loc_max <- TRUE; print("Local maximum reached"); break}
        res <- rbind(res,c(iter,cuscor,cusviz))
      }
    }
    res.t <- Sys.time() - t
    return(list(res.rf=res.rf, par.rf=par.rf, max.TotDiffVal1.rf=cuscor.rf,res=res, par=parcor, local_maximum=loc_max, time=res.t, max.TotDiffVal1=cuscor))
  }

  ###### THE COMMON CASE
  if (full.output == FALSE) {
    res.list <- list()
    for (n.run in 1:n.starts) {

      if (p.initial[1] == "random") {
        p.ini <- sample(c(rep(1:k, mgs), sample(k, nr-mgs*k, replace=TRUE)))
        tp <- table(p.ini)
      } else {
        p.ini <- p.initial
        tp <- table(p.ini)
      }

      arf <- matrix(0, k, ns); colnames(arf) <- colnames(mt); rownames(arf) <- 1:k; adp <- arf #arf is the adjusted relative frequency in each group and adp is the adjusted differential proportion
      afg <- rowsum(mt, group=p.ini) #no. of relevés containing the taxon, within each group (absolute frequency in each group)
      t1 <- (afg==0)*as.vector(tp) #no. of relevés of each group, when the taxon is not present
      t2 <- colSums(afg>0) #no. of groups containing the taxon
      t3 <- t2 > 1 #indices of the taxa occurring in more than one group (they must occur in at least one)
      for (i in 1:k) { #fills matrices adp and arf, only when the taxon is present in the group!
        rel.i <- which(p.ini==i) #indices of the relevés of the group
        spe.i <- afg[i,]>0 #indices of the taxa present in the group
        adp[i,spe.i & !t3] <- 1 #adp is 1 for the taxa occurring in one group only!
        if (sum(t3) > 0) { #for taxa occurring in more than one group
          t4 <- spe.i & t3
          adp[i, t4] <- colSums(t1[-i, t4, drop=FALSE]/sum(tp[-i]))/(t2[t4])
        }
        arf[i, spe.i] <- afg[i, spe.i]/as.vector(tp[i])/(t2[spe.i])
      }
      cuscor <- sum(colSums(arf*adp))/ns

      if (random.first==TRUE) {

        #random.neighbour: this function returns randomly one of the (k-1)*nr 1-neighbour partitions, assuring that the minimum group size (mgs) is respected
        random.neighbour <- function (p) {
          tp <- table(p)
          k.int <- which(tp>mgs) #k of interest to sample
          k.sam <- tp[k.int]-mgs #k samplable
          troca <- sample(rep(k.int,k.sam), min(rf.neigh.size, sum(k.sam))) #neighbourhood.size cannot be greater than sum(k.sam)
          niter <- table(troca)
          gn.tot <- NULL
          in.tot <- NULL
          for (gv in sort(unique(troca))) {#it assumes that table function also sorts!
            if(length(c(1:k)[-gv])!=1) {gn.tot <- c(gn.tot,sample(c(1:k)[-gv],niter[paste(gv)],replace=TRUE))} else {gn.tot <- c(gn.tot,rep(c(1:k)[-gv],niter[paste(gv)]))}
            in.tot <- c(in.tot,sample(which(p==gv),niter[paste(gv)]))
          }
          p[in.tot] <- gn.tot
          return(p)
        }

        #STOCHASTIC HILL CLIMBING (random.first=TRUE)
        parcor <- p.ini
        for (iter in 1:rf.maxit) {
          parviz <- random.neighbour(parcor)
          tp <- table(parviz)
          arf <- matrix(0, k, ns); colnames(arf) <- colnames(mt); rownames(arf) <- 1:k; adp <- arf #arf is the adjusted relative frequency in each group and adp is the adjusted differential proportion
          afg <- rowsum(mt, group=parviz) #no. of relevés containing the taxon, within each group (absolute frequency in each group)
          t1 <- (afg==0)*as.vector(tp) #no. of relevés of each group, when the taxon is not present
          t2 <- colSums(afg>0) #no. of groups containing the taxon
          t3 <- t2>1 #indices of the taxa occurring in more than one group (they must occur in at least one)
          for (i in 1:k) { #fills matrices adp and arf, only when the taxon is present in the group!
            rel.i <- which(parviz==i) #indices of the relevés of the group
            spe.i <- afg[i,]>0 #indices of the taxa present in the group
            adp[i,spe.i & !t3] <- 1 #adp is 1 for the taxa occurring in one group only!
            if (sum(t3)>0) { #for taxa occurring in more than one group
              t4 <- spe.i & t3
              adp[i, t4] <- colSums(t1[-i, t4, drop=FALSE]/sum(tp[-i]))/(t2[t4])
            }
            arf[i, spe.i] <- afg[i,spe.i]/as.vector(tp[i])/(t2[spe.i])
          }
          cusviz <- sum(colSums(arf*adp))/ns
          if (cusviz >= cuscor) {parcor <- parviz; cuscor <- cusviz}
        }
        p.ini <- parcor
        tp <- table(p.ini)
        #afg, t1, t2, t3, arf and adp must be recalculated for parcor (as a worse neighbour might have been used to calculate them before)!
        arf <- matrix(0, k, ns); colnames(arf) <- colnames(mt); rownames(arf) <- 1:k; adp <- arf #arf is the adjusted relative frequency in each group and adp is the adjusted differential proportion
        afg <- rowsum(mt, group=p.ini) #no. of relevés containing the taxon, within each group (absolute frequency in each group)
        t1 <- (afg==0)*as.vector(tp) #no. of relevés of each group, when the taxon is not present
        t2 <- colSums(afg>0) #no. of groups containing the taxon
        t3 <- t2>1 #indices of the taxa occurring in more than one group (they must occur in at least one)
        for (i in 1:k) { #fills matrices adp and arf, only when the taxon is present in the group!
          rel.i <- which(p.ini==i) #indices of the relevés of the group
          spe.i <- afg[i,]>0 #indices of the taxa present in the group
          adp[i,spe.i & !t3] <- 1 #adp is 1 for the taxa occurring in one group only!
          if (sum(t3)>0) { #for taxa occurring in more than one group
            t4 <- spe.i & t3
            adp[i, t4] <- colSums(t1[-i, t4, drop=FALSE]/sum(tp[-i]))/(t2[t4])
          }
          arf[i, spe.i] <- afg[i,spe.i]/as.vector(tp[i])/(t2[spe.i])
        }
      }

      #max.tdv.neighbour: this function returns the(one of the) 1-neighbouring partition(s) presenting greater TotDiffVal1 (tdv)
      max.tdv.neighbour <- function (p) {
        tp <- table(p); mtemp <- matrix(p, nr, nr, byrow=TRUE); fordiag <- p; res1 <- NULL; res2 <- NULL
        if (min(tp)>mgs) { #all groups have more than mgs elements
          for (i in 1:(k-1)) {
            fordiag <- fordiag+1
            fordiag[which(fordiag==k+1)] <- 1
            diag(mtemp) <- fordiag
            res1 <- rbind(res1, mtemp)
            res2 <- c(res2,fordiag)
          }
          res2 <- cbind(rep(p, k-1),res2)
          colnames(res2) <- NULL
        } else { #at least one group has only mgs elements
          ind.rm <- as.numeric(sapply(which(tp==mgs), function (x) {which(p==x)})) #it returns the indices of the partition corresponding to groups presenting only mgs elements #as.numeric is probably not necessary, it is possibly done by default in []
          for (i in 1:(k-1)) {
            fordiag <- fordiag+1
            fordiag[which(fordiag==k+1)] <- 1
            diag(mtemp) <- fordiag
            res1 <- rbind(res1, mtemp[-ind.rm,])
            res2 <- c(res2, fordiag[-ind.rm])
          }
          res2 <- cbind(rep(p[-ind.rm], k-1),res2)
          colnames(res2) <- NULL
        }
        mat.neig <- list(p.list=res1, pairs=res2)

        mat.neig.tdv <- sapply(1:nrow(mat.neig$p.list), function (x) { #like this, it is difficult to parallelize!
          pn <- mat.neig$p.list[x,]
          kc <- mat.neig$pairs[x,]
          tpn <- table(pn) #neighbour partition
          for (i in kc) {afg[i,] <- colSums(mt[which(pn==i),])} #(updates afg) no. of relevés containing the taxa, within each group, only for the two groups that changed a relevé!
          t5 <- afg[kc,][1,]>0 | afg[kc,][2,]>0 #taxa affected by the change of relevés
          t1[kc,] <- (afg==0)[kc,]*as.vector(tpn)[kc] #(updates t1) no. of relevés of each group, when the taxon is not present #t5 must not be used here!
          t2[t5] <- colSums(afg>0)[t5] #(new t2) no. of groups containing the taxon
          t3[t5] <- (t2>1)[t5] #indices of the taxa occurring in more than one group (they must occur in at least one)
          for (i in 1:k) { #changes matrices adp and arf, only when the taxon is present in the group!
            rel.i <- which(pn==i) #indices of the relevés of the group
            spe.i <- afg[i,]>0 #indices of the taxa present in the group
            if (sum(spe.i & t5)>0) { #caso o grupo em causa contenha espécies afetadas
              adp[i, spe.i & !t3] <- 1 #adp is 1 for the taxa occurring in one group only!
              if (sum(t3)>0) { #for taxa occurring in more than one group
                t4 <- spe.i & t3
                adp[i, t4] <- colSums(t1[-i, t4, drop=FALSE]/sum(tpn[-i]))/(t2[t4])
              }
              adp[i, t5 & !spe.i] <- 0 #inserts a zero to the affected taxa that is no more present in the group!
              arf[i, t5] <- afg[i, t5]/as.vector(tpn[i])/(t2[t5])
            }
          }
          return(sum(colSums(arf*adp))/ns)
        })
        return(list(tdv=tdv <- max(mat.neig.tdv), p=mat.neig$p.list[which(mat.neig.tdv==tdv)[1],]))
      }

      #HILL CLIMBING
      if (maxit==0) {parcor=p.ini; loc_max = NA} else {
        loc_max <- NA
        parcor <- p.ini
        #cuscor is already calculated
        for (iter in 1:maxit) {
          temp <- max.tdv.neighbour(parcor)
          parviz <- temp$p
          cusviz <- temp$tdv
          if (cusviz > cuscor) {parcor <- parviz; cuscor <- cusviz} else {loc_max <- TRUE; print("Local maximum reached"); break}
        }
      }

      if (n.run == 1) {res.list[[1]] <- list(local_maximum = loc_max, par = parcor, max.TotDiffVal1 = cuscor)} else {
        already.in.bestsol <- any(sapply(res.list, function (x) {equivalent_partition(x$par, parcor)}))
        if (!already.in.bestsol) {
          if (length(res.list) < n.sol) {
            res.list[[length(res.list)+1]] <- list(local_maximum=loc_max, par=parcor, max.TotDiffVal1=cuscor)
          } else { #already n.sol components in res.list
            best.sol.values <- sapply(res.list, function (x) {x$max.TotDiffVal1})
            if (cuscor >  min(best.sol.values)) {
              worse.bestsol <- which(best.sol.values == min(best.sol.values))[1] #selects one, in case of ties!
              res.list[[worse.bestsol]] <- list(local_maximum=loc_max, par=parcor, max.TotDiffVal1=cuscor)
            }
          }
        }
      }
    }
    ind.order <- order(sapply(res.list, function (x) {x$max.TotDiffVal1}))
    return(res.list[ind.order])
  }
}
