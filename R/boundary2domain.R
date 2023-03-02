#' Convert boundary to domain
#'
#' @description From datas with each boundaries stored by line (chromosome, start and end), `boundary2domain()` return `GRanges` with domains.
#'
#' @details Start and end of domains are the middle of boundaries.
#'
#' @param boundaries `dataframe` (chromosome, start and end) or `GRanges` with the boundaries.
#'
#' @return `GRanges` object with domains
#'
#' @import GenomicRanges
#'
#' @export
#'
boundary2domain <- function(boundaries) {

  if (class(boundaries) == "GRanges") {
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
  rownames(df) = paste0(df$seqnames, "_", df$start)
  return(GenomicRanges::makeGRangesFromDataFrame(df))
}
