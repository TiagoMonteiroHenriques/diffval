# tdv.R
#'
#' @title The total differential value of a phytosociological table
#'
#' @description Calculates TotDiffVal index (or TotDiffVal1) given a phytosociological table (see http://home.isa.utl.pt/~tmh/).
#'
#' @param m A `matrix`, i.e. a phytosociological table of 0s (absences) and 1s (presences), where rows correspond to taxa and columns correspond to relevés.
#' @param p A `vector` of integer numbers with the partition of the relevés (i.e. a k-partition, consisting in a vector with values from 1 to k, with length equal to the number of columns of m, ascribing each relevé to one of the k groups).
#' @param index A `character` selecting which index is calculated: "TotDiffVal", "TotDiffVal1" or "both". If `output.type = "fast"`, "both" can not be selected.
#' @param output.type A `character`, determining the amount of information returned by the function and also pre-validations (defaults to "normal").
#'
#' @details The function accepts a phytosociological table (`m`) and a k-partition of its columns (`p`), returning the respective TotDiffVal (or TotDiffVal1) index.
#' TotDiffVal was proposed by Monteiro-Henriques and Bellu (2014). Monteiro-Henriques (2016) proposed TotDiffVal1, modifying the previous index slightly with the objective of ensuring the index to be bounded between 0 and 1. Yet, TotDiffVal is already bounded between 0 and 1.
#' In practice, both indices are bounded between 0 and 1, but TotDiffVal1 reduces further the contribution of differential taxa present in more than one group.
#'
#' @return If `output.type` = "normal" (the default) pre-validations are done and a `list` is returned, with the following components:
#'
#' \describe{
#'   \item{arf}{A `matrix` with the adjusted relative frequency (a/b/e) of each taxon in each group (see http://home.isa.utl.pt/~tmh/).}
#'   \item{adp}{A `matrix` with the adjusted differential proportion (c/d/e) of each taxon in each group (see http://home.isa.utl.pt/~tmh/).}
#'   \item{DiffVal1}{A `matrix` with the DiffVal1 of each taxon (see http://home.isa.utl.pt/~tmh/).}
#'   \item{TotDiffVal1}{A `numeric` with the TotDiffVal1 of matrix `m` given the partition `p` (see http://home.isa.utl.pt/~tmh/).}
#' }
#'
#' If `output.type` = "full", some extra components are added to the output: afg, t1, t2, t3. These are intermediate matrices used in the computation of TotDiffVal (or TotDiffVal1).
#'
#' If `output.type` = "fast", only the value of the selected index is returned and no pre-validations are done.
#'
#' @author Tiago Monteiro-Henriques. E-mail: \email{tiagomonteirohenriques@@gmail.com}.
#'
#' @export
#'
tdv <- function (m, p, index = "TotDiffVal", output.type = "normal") {
  if (output.type == "fast") {
    k <- max(p)
  } else {
    stopifnot(is.matrix(m))
    if (!identical(length(p), ncol(m))) {stop("Object p must be a partition of the columns of m")}
    k <- max(p)
    mode(p) <- "integer"
    if (!identical(sort(unique(p)), 1:k)) {stop("Object p is not a valid partition of the columns of m")}
    mode(m) <- "integer"
    if (!identical(c(0L, 1L), sort(unique(as.vector(m))))) {stop("Matrix m should contain only 0's and 1's")}
    if (min(rowSums(m))==0) {stop("At least one taxa is not present in any relev\u00e9")}
    if (min(colSums(m)) == 0) {stop("At least one relev\u00e9 contains no taxa")}
  }
  mt <- t(m)
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

  if (output.type == "fast") {
    if (index == "TotDiffVal") {
      return(sum(colSums(af*ad)/t2)/ns)
    }
    if (index == "TotDiffVal1") {
      return(sum(colSums(af*ad)/t2^2)/ns)
    }
    if (index == "both") {
      stop('Arguments combination output.type = "fast" and index = "both" is not implemented.')
    }
  }

  if (output.type == "normal") {
    colnames(ad) <- colnames(af) <- colnames(mt)
    rownames(ad) <- rownames(af) <- 1:k
    arf <- af/t2
    adp <- ad/t2
    if (index == "TotDiffVal") {
      DV <- colSums(af*ad)/t2
      return(list(arf = t(arf), adp = t(adp), DiffVal = matrix(DV, ns, 1, dimnames = list(colnames(af), c("diff.val"))), TotDiffVal = sum(DV)/ns))
    }
    if (index == "TotDiffVal1") {
      DV1 <- colSums(arf*adp)
      return(list(arf = t(arf), adp = t(adp), DiffVal1=matrix(DV1, ns, 1, dimnames = list(colnames(arf), c("diff.val1"))), TotDiffVal1 = sum(DV1)/ns))
    }
    if (index == "both") {
      DV <- colSums(af*ad)/t2
      DV1 <- colSums(arf*adp)
      return(list(arf = t(arf), adp = t(adp), DiffVal1=matrix(DV1, ns, 1, dimnames = list(colnames(arf), c("diff.val1"))), TotDiffVal1 = sum(DV1)/ns))
    }
    return(list(arf=t(arf),adp=t(adp), DiffVal = matrix(DV, ns, 1, dimnames = list(colnames(af), c("diff.val"))), DiffVal1=matrix(DV1, ns, 1, dimnames=list(colnames(arf),c("diff.val1"))), TotDiffVal = sum(DV)/ns, TotDiffVal1 = sum(DV1)/ns))
  }

  if (output.type == "full") {
    colnames(ad) <- colnames(af) <- colnames(mt)
    rownames(ad) <- rownames(af) <- 1:k
    arf <- af/t2
    adp <- ad/t2
    if (index == "TotDiffVal") {
      DV <- colSums(af*ad)/t2
      return(list(afg = afg, t1 = t1, t2 = t2, t3 = t3, arf = arf, adp = adp, DiffVal = matrix(DV, ns, 1, dimnames = list(colnames(af), c("diff.val"))), TotDiffVal = sum(DV)/ns))
    }
    if (index == "TotDiffVal1") {
      DV1 <- colSums(arf*adp)
      return(list(afg = afg, t1 = t1, t2 = t2, t3 = t3, arf = arf, adp = adp, DiffVal1 = matrix(DV1, ns, 1, dimnames = list(colnames(arf), c("diff.val1"))), TotDiffVal1=sum(DV1)/ns))
    }
    if (index == "both") {
      DV <- colSums(af*ad)/t2
      DV1 <- colSums(arf*adp)
      return(list(afg = afg, t1 = t1, t2 = t2, t3 = t3, arf = arf, adp = adp, DiffVal = matrix(DV, ns, 1, dimnames = list(colnames(af), c("diff.val"))), DiffVal1 = matrix(DV1, ns, 1, dimnames = list(colnames(arf), c("diff.val1"))), TotDiffVal = sum(DV)/ns, TotDiffVal1=sum(DV1)/ns))
    }
  }
  stop("Maybe arguments index and/or output.type were not well defined?")
}
