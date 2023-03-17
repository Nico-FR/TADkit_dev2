#' Convert bigwig file to `GRanges`
#'
#' @description From the full path of a bigwig files, return a `GRanges` object. The bin width of the `GRanges` can be specify with `bin.width` parameter.
#'
#'
#' @param bigwig.path path of bigwig file.
#' @param bin.width Bin width. Default is NULL to return the bin size as it is in the bigwig file.
#'
#' @return `GRanges`
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import rtracklayer
#'
#' @export
#'
bw2grange <- function(bigwig.path, bin.width = NULL) {

  gr = rtracklayer::import.bw(bigwig.path, as = "GRanges")

  if (is.null(bin.width)) {return(gr)
  } else {

    bin.gr <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(gr),
                                      tilewidth = bin.width, cut.last.tile.in.chrom = TRUE)

    return(rtracklayer::import.bw(bigwig.path, which = bin.gr, as = "GRanges"))
  }
}
