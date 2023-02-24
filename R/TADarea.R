#' @title Annotations analysis around TAD boundaries
#'
#' @description For each TAD boundary (start or end), `TADarea` return all annotations that are in the desired window of the boundary (i.e boundary +/- `window.size`).
#' The main use of this function is to be combined with `areaHist()` or `areaCov()` functions.
#' It is possible to use the `"start"`, `"end"` or `"center"` of the TAD to analyze the distribution of annotations.
#'
#' @inheritParams TADplot
#' @param annot.gr `GRanges` object with genomic annotations.
#' @param window.size Window in bp surrounding each side of the TAD boundaries.
#' @param tad.boundary TAD border to be analyzed; `"start"`, `"end"` or `"center"` denoting what to use as an anchor. Default is the `"start"` of each TAD.
#'
#' @return `GRanges` object with all annotations surrounding all TAD boundaries. TAD boundaries are at the window.size position.
#'
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import methods
#' @export
#'
#' @examples
#' # tad.gr with 2 TADs
#' tad.gr <- dataframes2grange(
#'   data.frame(chr = rep(1, 2), start = c(2e6, 4e6), end = c(3e6, 5e6)),
#'   data.frame(chr = "1", size = 6e6)
#' )
#'
#'
#' # annot.gr with 3 genes
#' annot.gr <- dataframes2grange(
#'   data.frame(
#'     chr = rep(1, 3),
#'     start = c(1.5e6, 3.5e6, 5.5e6),
#'     end = c(1.6e6, 3.6e6, 5.6e6),
#'     type = c("a", "b", "c")
#'   ),
#'   data.frame(chr = "1", size = 6e6),
#'   metadata.mcols = 4
#' )
#'
#' TADarea(tad.gr = tad.gr, annot.gr = annot.gr, window.size = 2e6)
TADarea <- function(tad.gr, annot.gr, window.size = 50e3, tad.boundary = "start") {
  GenomeInfoDb::seqlengths(annot.gr) <- NA # remove seqlengths to prevent out-of-bound ranges
  ##############################
  # function to use: shift.R
  ##############################
  shift <- function(x, tad.gr, annot.gr, window.size) {
    tad_1 <- tad.gr[x] # for each TAD boundary

    # Annot.gr filtering within TAD boundary +/- window.size
    data.gr <- annot.gr[
      GenomeInfoDb::seqnames(annot.gr) == as.character(GenomeInfoDb::seqnames(tad_1)) # filter chr
      & GenomicRanges::end(annot.gr) >= GenomicRanges::start(tad_1) - window.size # upstream limit
      & GenomicRanges::start(annot.gr) <= GenomicRanges::start(tad_1) + window.size # downstream limit
    ]

    # shift annot.gr (so TADboundary is at window.size distance)
    GenomicRanges::shift(data.gr, -GenomicRanges::start(tad_1) + window.size)
  }

  ##############################
  # sapply function for each TAD border
  ##############################
  # resize the TADs according to start or stop
  tad.gr <- GenomicRanges::resize(tad.gr, 1, fix = tad.boundary)
  annot.gr$window.size <- window.size # add the window.size used
  annot.gr$nb_boundary <- length(tad.gr) # add nb domains
  annot.gr$boundary_type <- tad.boundary

  # sapply function for each TAD border
  output <- unlist(methods::as(
    sapply(1:length(tad.gr), shift,
      tad.gr = tad.gr,
      annot.gr = annot.gr,
      window.size = window.size
    ),
    "GRangesList"
  ))
  return(output)
}
