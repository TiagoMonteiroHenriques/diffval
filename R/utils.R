# SOME AUXILIARY FUNCTIONS

# Auxiliary function to HillClimb_optim_tdv
#random.neighbour: this function returns randomly one of the rf.neigh.size-neighbour partitions, assuring that the minimum group size (mgs) is respected
random.neighbour <- function (p, k, mgs, rf.neigh.size) {
  tp <- tabulate(p)
  k.int <- which(tp > mgs) #k of interest to sample
  k.sam <- tp[k.int] - mgs #k samplable
  swap <- sample(rep(k.int, k.sam), min(rf.neigh.size, sum(k.sam))) #neighbourhood.size cannot be greater than sum(k.sam)
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

# Auxiliary function to HillClimb_optim_tdv
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
#random.neighbour.SA selects a random neighbour (or the same partition), changing just one relevé

random.neighbour.SA <- function (p, nr, k) {
  tp <- table(p)
  k.int <- which(tp > 1) #Check behaviour when all groups have only one element...
  change.col <- sample((1:nr)[p %in% k.int], 1) #Changing just one column
  #tp[p[change.col]] <- tp[p[change.col]] - 1 #Changing just one column
  p[change.col] <- sample(1:k, length(change.col), replace=TRUE) #CF replace = TRUE (added later)
  return(p)
}


# Auxiliary functions to decrease GRASP and GRDTP computation time
# These are three auxiliary functions to recalculate parameter "c" (outside group g) given a new column relevé to be included in group g
#these functions could be adapted to the use of function ifelse(); apparently the gain in time is only for smaller datasets,
# as mapply() returned faster for big datasets.

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
