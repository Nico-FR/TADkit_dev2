#' Convert bigwig file to `GRanges`
#'
#' @description From the full path of a bigwig files, return a `GRanges` object. The bin width of the `GRanges` can be specify with `bin.width` parameter. These function does not work on Windows.
#'
#'
#' @param bigwig.path path of bigwig file.
#' @param bin.width Bin width. Default is NULL to return the bin size as it is in the bigwig file.
#'
#' @return `GRanges`
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import rtracklayer
#' @examples
#' library("rtracklayer")
#' if (.Platform$OS.type != "windows") { #do no not work on window
#' #translate bedgraph to bigwig file
#' export.bw(rna_seq_chr25_13to16mb.bedgraph, "./rna_seq_chr25_13to16mb.bw")
#' }
#' bw2grange("./rna_seq_chr25_13to16mb.bw")
#' file.remove("./rna_seq_chr25_13to16mb.bw")
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
