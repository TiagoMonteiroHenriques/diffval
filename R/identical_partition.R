#' Do the vectors represent the same k-partition?
#'
#' Checks if two vectors represent the same k-partition.
#'
#' @param p1 A vector of integers representing a k-partition (taking values
#'   from 1 to k), of the same length of `p2`.
#' @param p2 A vector of integers representing a k-partition (taking values
#'   from 1 to k), of the same length of `p1`.
#'
#' @details Parameters `p1`and `p2`are vectors indicating group membership.
#'   In this package context, these vectors have as many elements as the columns
#'   of a phytosociological table, indicating the group membership of each
#'   relevé to one of k groups (i.e., a k-partition).
#'   This function checks if the two given vectors `p1`and `p2` correspond, in
#'   practice, to the same k-partition, i.e., if the relevé groups are actually
#'   the same, but the group numbers are somehow swapped.
#'
#' @return `TRUE` if `p1`and `p2` represent the same k-partitions; `FALSE`
#'   otherwise.
#'
#' @author  Tiago Monteiro-Henriques and Jorge Orestes Cerdeira.
#'   E-mail: \email{tmh.dev@@icloud.com}.
#'
#' @examples
#' # Creating three 2-partitions
#' par1 <- c(1, 1, 2, 2, 2)
#' par2 <- c(2, 2, 1, 1, 1)
#' par3 <- c(1, 1, 1, 2, 2)
#'
#' # Is it the same partition?
#' identical_partition(par1, par2) # TRUE
#' identical_partition(par1, par3) # FALSE
#' identical_partition(par2, par3) # FALSE
#'
#' @export
identical_partition <- function(p1, p2) {
  if (!identical(length(p1), length(p2))) {
    stop("Partitions to compare must have the same length.")
  }
  if (!identical(
    as.numeric(sort(unique(p1))),
    as.numeric(sort(unique(p2)))
  )
  ) {
    stop(
      "Partitions to compare must use the same group names (or numbers) and must
      have the same total number of groups."
    )
  }
  if (nrow(unique(cbind(p1, p2))) == length(unique(p1))) {
    res <- TRUE
  } else {
    res <- FALSE
  }
  res
}
