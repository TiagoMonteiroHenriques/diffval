# tdv.R
#'
#' @title The Total Differential Value of a phytosociological table
#'
#' @description Given a phytosociological table and a partition of its columns, this function calculates the respective Total Differential Value (TDV).
#'
#' @param m.bin A `matrix`, i.e. a phytosociological table of 0s (absences) and 1s (presences), where rows correspond to taxa and columns correspond to relevés.
#' @param p A `vector` of integer numbers with the partition of the relevés (i.e. a k-partition, consisting in a vector with values from 1 to k, with length equal to the number of columns of `m.bin`, ascribing each relevé to one of the k groups).
#' @param output.type A `character`, determining the amount of information returned by the function and also the amount of pre-validations. Possible values are "normal" (the default), "fast" and "full".
#'
#' @details The function accepts a phytosociological table (`m.bin`) and a k-partition of its columns (`p`), returning the corresponding TDV.
#' TDV was proposed by Monteiro-Henriques and Bellu (2014). Monteiro-Henriques (2016) proposed TDV1, modifying TDV slightly with the objective of ensuring a value from 0 to 1. Yet, TDV is always within that range.
#' In practice, both TDV and TDV1 have 0 as possible minimum value and 1 as possible maximum value, but TDV1 reduces further the contribution of differential taxa present in more than one group. TDV is then implemented here, for parsimony.
#'
#' TDV is calculated using the \eqn{DiffVal} index for each (and all) of the taxa present in a tabulated phytosociological table (also called sorted table).
#' \eqn{DiffVal} index aims at characterizing how well a taxon works as a differential taxon in a tabulated phytosociological table (for more information on differential taxa
#' see Mueller-Dombois & Ellenberg, 1974).
#' An archetypal differential taxon of a certain group \eqn{g} of the partition \eqn{p} is the one present in all relevés of group \eqn{g}, and absent from all the other
#' groups of that partition. Therefore, \eqn{DiffVal} has two components, an inner one (\eqn{\frac{a}{b}}{a/b}), which measures the presence of the taxon inside each of the groups,
#' and an outer one (\eqn{\frac{c}{d}}{c/d}), which measures the relevant absences of the taxon outside of each of the groups. Specifically, given a partition \eqn{p} with \eqn{k} groups, \eqn{DiffVal} is
#' calculated for each taxon \eqn{s} as:
#'
#' \deqn{DiffVal_{s,p} = \frac{1}{e}\sum_{g=1}^k{\frac{a}{b}\frac{c}{d}}}{DiffVal of taxon s, given the partition p = the summation of [(a/b)*(c/d)]/e, obtained for each group g, from g = 1 to k}
#' where,
#'
#' \itemize{
#' \item{\eqn{a}, is the total number of presences of taxon \eqn{s} within group \eqn{g}}
#' \item{\eqn{b}, is the total number of relevés of group \eqn{g}}
#' \item{\eqn{c}, is the total number of differentiating absences of taxon \eqn{s}, i.e. absences coming from the groups other than \eqn{g} from which the taxon \eqn{s} is completely absent}
#' \item{\eqn{d}, is the total number of relevés of all groups but \eqn{g} (i.e., the total number of relevés in the table - \eqn{b})}
#' \item{\eqn{e}, is the total number of groups in which the taxon \eqn{s} occurs at least once}
#' }
#'
#' Therefore, for each taxon \eqn{s} and for each group \eqn{g}, the \eqn{DiffVal} index evaluates:
#' \itemize{
#' \item{\eqn{\frac{a}{b}}{a/b}, i.e., the frequency of the presences of taxon \eqn{s}, relative to the size of group \eqn{g}; commonly called 'relative frequency.' \eqn{\frac{a}{b}}{a/b} is only 1 if and only if taxon \eqn{s} occurs in all the relevés of group \eqn{g}.}
#' \item{\eqn{\frac{c}{d}}{c/d}, i.e., the frequency of the differentiating absences of taxon \eqn{s} outside group \eqn{g}, relative to the sum of sizes of all groups but \eqn{g}. _Nota bene_: absences in \eqn{c} are counted outside the group \eqn{g} but only in the groups from which taxon \eqn{s} is completely absent (these are the relevant absences, which produce differentiation among groups); in practice \eqn{c} corresponds to the sum of the sizes of all groups other than \eqn{g} that are empty. \eqn{\frac{c}{d}}{c/d} is 1 if and only if the taxon \eqn{s} is absent from all groups but \eqn{g}.}
#' }
#'
#' Finally, \eqn{\frac{1}{e}}{1/e} ensures that \eqn{DiffVal} is a value from 0 to 1.
#'
#' The Total Differential Value (TDV or \eqn{TotDiffVal}) of a phytosociological table \eqn{m.bin} tabulated/sorted by the partition \eqn{p} is:
#'
#'  \deqn{TDV_{m.bin,p} = \frac{1}{n}\sum_{i=1}^n{Diffval_{i,p}}}{TDV of table m.bin, given the partition p = the summation of the DiffVal of all taxa in the table, divided by n}
#' where,
#'
#' \itemize{
#' \item{\eqn{n}, is the number of taxa in table \eqn{m.bin}}
#' }
#'
#' The division by the number of taxa present in \eqn{m.bin} ensures that TDV remains in the \[0,1\] interval (as \eqn{DiffVal} is also in the same interval).
#'
#' @return If `output.type` = "normal" (the default) pre-validations are done and a `list` is returned, with the following components:
#'
#' \describe{
#'   \item{ifp}{A `matrix` with the \eqn{\frac{a}{b}}{a/b} values for each taxon in each group, for short called the 'inner frequency of presences'.}
#'   \item{ofda}{A `matrix` with the \eqn{\frac{c}{d}}{c/d} values for each taxon in each group, for short called the 'outer frequency of differentiating absences'.}
#'   \item{e}{A `vector` with the \eqn{e} values for each taxon, i.e. the number of groups containing that taxon.}
#'   \item{diffval}{A `matrix` with the \eqn{DiffVal} for each taxon.}
#'   \item{tdv}{A `numeric` with the TDV of matrix `m.bin,` given the partition `p`.}
#' }
#'
#' If `output.type` = "full", some extra components are added to the output: afg, empty.size, gct (= \eqn{e}), i.mul. These are intermediate matrices used in the computation of TDV.
#'
#' If `output.type` = "fast", only TDV is returned and no pre-validations are done.
#'
#' @references
#' Monteiro-Henriques T. & Bellu A. 2014. *An optimization approach to the production of differentiated tables based on new differentiability measures*. 23rd EVS European Vegetation Survey. Presented orally. Ljubljana, Slovenia.
#'
#' Monteiro-Henriques T. 2016. *A bunch of R functions to assist phytosociological tabulation*. 25th Meeting of European Vegetation Survey. Presented in poster. Rome. Italy.
#'
#' Mueller-Dombois D. & Ellenberg H. 1974. *Aims and Methods of Vegetation Ecology*. New York: John Wiley & Sons.
#'
#' @author Tiago Monteiro-Henriques. E-mail: \email{tmh.dev@@icloud.com}.
#'
#' @examples
#'
#' #getting the Taxus baccata forests data set
#' data(taxus_bin)
#' #creating a group partition, as presented in original article of the data set
#' groups <- rep(c(1,2,3), c(3,11,19))
#'
#' #removing taxa occurring in only one relevé in order to
#' #reproduce exactly the example in the original article of the data set
#' taxus_bin_wmt <- taxus_bin[rowSums(taxus_bin) > 1,]
#'
#' #calculating TDV
#' result <- tdv(taxus_bin_wmt, groups)
#'
#' #this is the TDV
#' result$tdv
#' #this is TDV1, reproducing exactly the value from the original article
#' sum(result$diffval/result$e)/nrow(taxus_bin_wmt)
#'
#' @export
#'
tdv <- function (m.bin, p, output.type = "normal") {
  if (output.type == "fast") {
    k <- max(p)
  } else {
    stopifnot(is.matrix(m.bin))
    if (!identical(length(p), ncol(m.bin))) {stop("Object p must be a partition of the columns of m.bin")}
    k <- max(p)
    mode(p) <- "integer"
    if (!identical(sort(unique(p)), 1:k)) {stop("Object p is not a valid partition of the columns of m.bin")}
    mode(m.bin) <- "integer"
    if (!identical(c(0L, 1L), sort(unique(as.vector(m.bin))))) {stop("Matrix m.bin should contain only 0's and 1's")}
    if (min(rowSums(m.bin))==0) {stop("At least one taxa is not present in any relev\u00e9")}
    if (min(colSums(m.bin)) == 0) {stop("At least one relev\u00e9 contains no taxa")}
  }
  mt <- t(m.bin)
  ns <- ncol(mt) #no. of taxa
  tp <- tabulate(p) #size of each group (inner)
  outer.size <- length(p) - tp #sum of the sizes of the outer groups
  ofda <- ifp <- matrix(0, k, ns) #matrices to store [a/b], i.e. the inner frequency of presences (ifp) and [c/d], i.e. the outer frequency of differenciating absences (ofda)
  afg <- rowsum(mt, group = p) #no. of relevés containing the taxon, within each group (absolute frequency in each group)
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
  if (output.type == "fast") {
    return(tdv = sum(colSums(ifp * ofda) / gct) / ns) #for TDV1 replace gct by gct^2
  }
  colnames(ofda) <- colnames(ifp) <- colnames(mt)
  rownames(ofda) <- rownames(ifp) <- 1:k
  DV <- colSums(ifp * ofda) / gct
  if (output.type == "normal") {
    return(list(ifp = t(ifp), ofda = t(ofda), e = gct, diffval = matrix(DV, ns, 1, dimnames = list(colnames(ifp), c("DiffVal"))), tdv = sum(DV) / ns))
  }
  if (output.type == "full") {
    return(list(ifp = ifp, ofda = ofda, e = gct, afg = afg, empty.size = empty.size, gct = gct, i.mul = i.mul, diffval = matrix(DV, ns, 1, dimnames = list(colnames(ifp), c("DiffVal"))), tdv = sum(DV) / ns))
  }
  stop('Argument output.type must be "fast", "normal" or "full".')
}
