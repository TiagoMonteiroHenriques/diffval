# equivalent_partition.R
#'
#' @title Are the two k-partitions equivalent?
#'
#' @description Checks if two k-partitions are equivalent.
#'
#' @param p1 A `vector` of integers representing a k-partition (taking values from 1 to k), of the same length of `p2`.
#' @param p2 A `vector` of integers representing a k-partition (taking values from 1 to k), of the same length of `p1`.
#'
#' @details The function checks if the two given k-partitions are equivalent, i.e. if their groupings are the same (useful
#' in the case that the partitions have a divergent group numbering).
#'
#' @return `TRUE` if the k-partitions are equivalent; `FALSE` otherwise.
#'
#' @author Tiago Monteiro-Henriques. E-mail: \email{tmh.dev@@icloud.com}.
#'
#' @examples
#'
#' #creating three 2-partitions
#' par1 <- c(1,1,2,2,2)
#' par2 <- c(2,2,1,1,1)
#' par3 <- c(1,1,1,2,2)
#'
#' #Are they equivalent?
#' equivalent_partition(par1, par2) #TRUE
#' equivalent_partition(par1, par3) #FALSE
#' equivalent_partition(par2, par3) #FALSE
#'
#' @export
#'
equivalent_partition <- function(p1,p2) {
    if (!identical(length(p1), length(p2))) {stop("Partitions to compare must have the same length")}
    if (!identical(as.numeric(sort(unique(p1))), as.numeric(sort(unique(p2))))) {stop("Partitions to compare must use the same group names (or numbers) and must have the same total number of groups")}
    if (nrow(unique(cbind(p1,p2))) == length(unique(p1))) {
      res <- TRUE
    } else {
      res <- FALSE
    }
    return(res)
  }
