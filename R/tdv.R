# tdv.R
#'
#' @title The total differential value (TotDiffVal1) of a phytosociologic table
#'
#' @description Calculates TotDiffVal1 index (see http://home.isa.utl.pt/~tmh/ for the calculation formula of TotDiffVal1).
#'
#' @param m A `matrix`, i.e. a phytosociological table of 0s (absences) and 1s (presences), where rows correspond to taxa and columns correspond to relevés.
#' @param p A `vector` of integer numbers with the partition of the relevés (i.e. a k-partition, consisting in a vector with values from 1 to k, with length equal to the number of columns of m, ascribing each relevé to one of the k groups).
#' @param full.output A `logical`, determining the amount of information returned by the function (defaults to `FALSE`).
#'
#' @details The function accepts a phytosociological table (`m`) and a k-partition of its columns (`p`), returning the respective TotDiffVal1 index.
#'
#' @return If `full.output` = `FALSE`, a `list` with the following components:
#'
#' \describe{
#'   \item{arf}{A `matrix` with the adjusted relative frequency (a/b/e) of each taxon in each group (see http://home.isa.utl.pt/~tmh/).}
#'   \item{adp}{A `matrix` with the adjusted differential proportion (c/d/e) of each taxon in each group (see http://home.isa.utl.pt/~tmh/).}
#'   \item{DiffVal1}{A `matrix` with the DiffVal1 of each taxon (see http://home.isa.utl.pt/~tmh/).}
#'   \item{TotDiffVal1}{A `numeric` with the TotDiffVal1 of matrix `m` given the partition `p` (see http://home.isa.utl.pt/~tmh/).}
#' }
#'
#' If `full.output` = `TRUE`, these extra components are added: afg, t1, t2, t3. These are intermediate matrices used in the computation of TotDiffVal1.
#'
#' @author Tiago Monteiro-Henriques. E-mail: \email{tiagomonteirohenriques@@gmail.com}.
#'
#' @export
#'
tdv <- function (m, p, full.output=FALSE) {
  if (!identical(length(p), ncol(m))) {stop("Object p must be a partition of the columns of m")}
  k <- max(p)
  mode(m) <- "integer"
  if (!identical(c(0L,1L), sort(unique(as.vector(m))))) {stop("Matrix m should contain only 0's and 1's")}
  if (min(rowSums(m))==0) {stop("At least one taxa is not present in any relev\u00e9")}
  if (min(colSums(m))==0) {stop("At least one relev\u00e9 contains no taxa")}
  mt <- t(m)
  ns <- ncol(mt) #no. of taxa
  tp <- table(p)
  arf <- matrix(0, k, ns); colnames(arf) <- colnames(mt); rownames(arf) <- 1:k; adp <- arf
  afg <- rowsum(mt, group=p) #no. of relevés containing the taxon, within each group (absolute frequency in each group)
  t1 <- (afg==0)*as.vector(tp) #no. of relevés of each group, when the taxon is not present
  t2 <- colSums(afg>0) #no. of groups containing the taxon
  t3 <- t2>1 #indices of the taxa occurring in more than one group (they must occur in at least one)

  for (i in 1:k) { #fills matrices adp and arf, only when the taxon is present in the group!
    rel.i <- which(p==i) #indices of the relevés of the group
    spe.i <- afg[i,]>0 #indices of the taxa present in the group
    adp[i,spe.i & !t3] <- 1 #adp is 1 for the taxa occurring in one group only!
    if (sum(t3)>0) { #for taxa occurring in more than one group
      t4 <- spe.i & t3
      adp[i, t4] <- colSums(t1[-i, t4, drop=FALSE]/sum(tp[-i]))/(t2[t4])
    }
    arf[i, spe.i] <- afg[i,spe.i]/as.vector(tp[i])/(t2[spe.i])
  }
  if (full.output==FALSE) {return(list(arf=t(arf),adp=t(adp),DiffVal1=matrix(colSums(arf*adp), ns, 1, dimnames=list(colnames(arf),c("diff.val1"))),TotDiffVal1=sum(colSums(arf*adp))/ns))} else {
    return(list(afg=afg,t1=t1,t2=t2,t3=t3,arf=arf,adp=adp,DiffVal1=matrix(colSums(arf*adp),ns,1),TotDiffVal1=sum(colSums(arf*adp))/ns))}
}
