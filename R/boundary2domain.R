#' Convert boundary to domain
#'
#' @description From data with each boundaries stored by line (chromosome, start and end), return GRanges with domains.
#'
#' @details Start or end position of domains are the middle of a boundary.
#' @param boundaries data frame (chromosome, start and end) or GRanges with the boundaries
#'
#' @return GRanges object with boundaries
#' @import GenomicRanges
#' @export
#'
boundary2domain <- function(boundaries) {

  if (class(boundaries) == "GRnages") {
    boundaries = as.data.frame(boundaries)[,1:3]}

  boundaries$mean = apply(boundaries[,2:3], 1, mean)
  list = split(boundaries, boundaries[,1])

  df.lst = lapply(list, function(x) {
    data.frame(
      seqnames = x[1:(nrow(x) - 1), 1],
      start = x$mean[1:(nrow(x) - 1)],
      end = x$mean[2:nrow(x)])
    })
  df = do.call("rbind", df.lst)

  return(GenomicRanges::makeGRangesFromDataFrame(df))
}
