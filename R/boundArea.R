#' @title Annotations analysis around boundaries
#'
#' @description For each domain boundary (start or end), `boundArea()` return all annotations that are in the desired window of the boundary (i.e boundary +/- `window.size`).
#' The main use of this function is to be combined with `areaHist()` or `areaCov()` functions.
#' It is possible to use the `"start"`, `"end"` or `"center"` of the domains to analyze the distribution of annotations.
#'
#' @param domain.gr `GRanges` with domains (TADs, compartments...).
#' @param annot.gr `GRanges` with genomic annotations (genes, repeat elements...).
#' @param window.size Window in base pair surrounding each side of the TAD boundaries.
#' @param domain.boundary Domain boundary to be analyzed: `"start"` or `"end"`. Default is the `"start"` of each domain to analyzed annotations around it. Note that it is also possible to take the center of the domains with `"center"`.
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
#' TADarea(domain.gr = tad.gr, annot.gr = annot.gr, window.size = 2e6)
boundArea <- function(domain.gr, annot.gr, window.size = 50e3, domain.boundary = "start") {
  GenomeInfoDb::seqlengths(annot.gr) <- NA # remove seqlengths to prevent out-of-bound ranges
  ##############################
  # function to use: shift.R
  ##############################
  shift <- function(x, domain.gr, annot.gr, window.size) {
    tad_1 <- domain.gr[x] # for each TAD boundary

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
  domain.gr <- GenomicRanges::resize(domain.gr, 1, fix = domain.boundary)
  annot.gr$window.size <- window.size # add the window.size used
  annot.gr$nb_boundary <- length(domain.gr) # add nb domains
  annot.gr$boundary_type <- domain.boundary

  # sapply function for each TAD border
  output <- unlist(methods::as(
    sapply(1:length(domain.gr), shift,
      domain.gr = domain.gr,
      annot.gr = annot.gr,
      window.size = window.size
    ),
    "GRangesList"
  ))
  return(output)
}
