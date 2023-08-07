#' Convert bigwig file to `GRanges`
#'
#' @description From the full path of a bigwig files, return a `GRanges` object. The bin width of the `GRanges` can be specify with `bin.width` parameter. These function does not work on Windows.
#'
#'
#' @param bigwig.path path of bigwig file.
#' @param bin.width bin width. Default is NULL to return the bin size as it is in the bigwig file.
#' @param transformation.method method used to aggregate the coverage within bin.width. Three methods are available: "sum" (default), "mean" and "median".
#'
#' @return `GRanges`
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import rtracklayer
#' @examples
#' library("rtracklayer")
#' if (.Platform$OS.type != "windows") { #does no not work on window
#' #translate bedgraph to bigwig file
#' export.bw(rna_seq_chr25_13to16mb.bedgraph, "./rna_seq_chr25_13to16mb.bw")
#' print(bw2grange("./rna_seq_chr25_13to16mb.bw"))
#' file.remove("./rna_seq_chr25_13to16mb.bw")
#' }
#'
#' @export
#'
bw2grange <- function(bigwig.path, bin.width = NULL, transformation.method = "sum") {

  if (.Platform$OS.type == "windows") {stop("Bigwig file cannot be read on Windows systems")}

  gr = rtracklayer::import.bw(bigwig.path, as = "GRanges")

  if (is.null(bin.width)) {return(gr)} else {

    #bins genome
    bin.gr <- GenomicRanges::tileGenome(
      GenomeInfoDb::seqlengths(gr)[GenomeInfoDb::seqlengths(gr) > bin.width],
      tilewidth = bin.width, cut.last.tile.in.chrom = TRUE)

    #add missing gaps to input
    gaps.gr = GenomicRanges::gaps(gr)[BiocGenerics::strand(GenomicRanges::gaps(gr))=="*"]
    if (length(gaps.gr) > 0) {gaps.gr$score = 0}
    gr = c(gr, gaps.gr)

    # coverage for each bins
    hits <- GenomicRanges::findOverlaps(bin.gr, gr)

    if (transformation.method == "median") {
      agg <- S4Vectors::aggregate(gr, hits, score = median(score))
    }
    if (transformation.method == "mean") {
      agg <- S4Vectors::aggregate(gr, hits, score = mean(score))
    }
    if (transformation.method == "sum") {
      agg <- S4Vectors::aggregate(gr, hits, score = sum(score))
    }

    bin.gr$score = agg$score #bin.gr$score = agg$score
    names(S4Vectors::mcols(bin.gr)) <- paste0(transformation.method, "_score")
    return(bin.gr)
  }


}
