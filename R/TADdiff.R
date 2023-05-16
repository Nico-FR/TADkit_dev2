#' @title TAD boundaries differences between 2 or more samples.
#'
#' @description `TADdiff` check if for each TAD boundary from one sample there is an other boundary in the vicinity on all other samples provided.
#' It also return the score for each boundary and the delta between each sample.
#'
#' @details `TADdiff()` takes the score used to find TAD boundaries (e.g insulation score) and the TAD boundaries as a `GRanges` object (see `dataframes2grange()` function).
#' For each TAD boundary :
#'     -checks if there is an other boundary in the other sample +/- `window.size` (e.g +/- 20kb if `bin.size = 10e3` and `window.size = NULL`),
#'     -exports the score at the boundary and the delta score at he boundary compared to the other sample. Those delta scores can be used to rank and find the boundary with the highest differences.
#' Differences analysis (new boundary and delta score at boundary) are performed between each samples (2 samples or more).
#' Because the scores are often not accurate at the extremities of the chromosomes (beginning and the end of each chromosome), these regions can be ignored using `restrict` parameters.
#'
#' @param boundaries.lst List of TAD boundaries as `GRanges`. Each `GRanges` in the list must have names and seqlengths metadata. The width of each boundary must correspond to the resolution of the matrix used and must be the same between samples.
#' @param score.lst List of insulation score as a `GRanges`. Each `GRanges` in the list must have names and seqlengths metadata.
#' @param bin.width Default is `NULL` to estimate the bin.width from `boundaries.lst` files.
#' @param window.size Defaults is `NULL` to take a value 2 times higher than `bin.width`, otherwise set the `window.size` in base pairs.
#' @param restrict Default is `2e6` (i.e 2Mb) to remove boundary at the extremities of chromosomes.
#'
#' @return list of GRanges
#' @importFrom methods as
#' @importFrom IRanges %outside%
#' @import GenomicRanges
#' @import GenomeInfoDb
#'
#' @export
#'
#'
TADdiff <- function(boundaries.lst, score.lst, bin.width = NULL, window.size = NULL, restrict = 2.e6) {

  #local variables:
  i <- j <- x <- e <- s <- chr <- e2 <- s2 <-start1 <- end1 <- start2 <- end2 <- chr1 <- chr2 <- NULL

  #parameters
  if (is.null(bin.width)) {
    bin.width <- round(IRanges::median(BiocGenerics::width(unlist(methods::as(boundaries.lst, "GRangesList")))) /2)  * 2
  }

  if (is.null(window.size)) {
    window.size <- round(bin.width/2) * 4
  }

  message("######################################################")
  message(paste0("Bin width used is: ", bin.width, "bp."))
  message(paste0("New TAD boundary is FALSE if there is another TAD boundary in the vicinity +/-", window.size, "bp."))
  message(paste0("TAD boundaries for the first and last ", restrict, "bp of each chromosomes are removed from the analysis."))
  message("######################################################")

  output = NULL
  output.names = NULL
  for (ind1 in names(boundaries.lst)) {

    message("######################################################")
    message(paste0("TADdiff: ", ind1, " processing..."))
    start_time <- Sys.time()

    # insulation score from itself
    ind1_ins.gr = score.lst[[ind1]]
    ind1_boundaries.gr = GenomicRanges::restrict(boundaries.lst[[ind1]],
                                  start = GenomeInfoDb::seqlengths(boundaries.lst[[ind1]]) - GenomeInfoDb::seqlengths(boundaries.lst[[ind1]]) + restrict,
                                  end = GenomeInfoDb::seqlengths(boundaries.lst[[ind1]]) - restrict)

    #correct width(boundaries) < bin.width (i.e cut by the restric function)
    width(ind1_boundaries.gr) = ifelse(width(ind1_boundaries.gr) < bin.width + 1, bin.width + 1, width(ind1_boundaries.gr))

    #findoverlap between boundary and bedgraph
    bin_hit = GenomicRanges::findOverlaps(subject = ind1_ins.gr, query = ind1_boundaries.gr, ignore.strand=TRUE, select="arbitrary", minoverlap = bin.width/2)
    GenomicRanges::mcols(ind1_boundaries.gr)[[paste0("score_",ind1)]] = GenomicRanges::mcols(ind1_ins.gr[bin_hit])[,1]


    for (ind2 in names(boundaries.lst)) {
      if (ind1 == ind2) {next} #skeep analysis between same sample


      ind2_ins.gr = score.lst[[ind2]]
      ind2_boundaries.gr = GenomicRanges::restrict(boundaries.lst[[ind2]],
                                                   start = GenomeInfoDb::seqlengths(boundaries.lst[[ind2]]) - GenomeInfoDb::seqlengths(boundaries.lst[[ind2]]) + restrict,
                                                   end = GenomeInfoDb::seqlengths(boundaries.lst[[ind2]]) - restrict)

      #findoverlap between boundary and bedgraph
      bin_hit = GenomicRanges::findOverlaps(subject = ind2_ins.gr, query = ind1_boundaries.gr, ignore.strand=TRUE, select="arbitrary", minoverlap = bin.width/2)
      GenomicRanges::mcols(ind1_boundaries.gr)[[paste0("score_",ind2)]] = GenomicRanges::mcols(ind2_ins.gr[bin_hit])[,1]

      # #add column to specify if there is a TAD in the other sample (+/- window.size). "newTAD" means no TAD boundaries in the vicinity i.e. +/- window.size_value.
      GenomicRanges::mcols(ind1_boundaries.gr)[[paste0("newTAD_",ind2)]] = ind1_boundaries.gr %outside% (ind2_boundaries.gr + window.size - 1) # "minus 1" to count TAD with start(indiv1) = end(indiv2)+window.size or end(indiv1) = start(indiv2)-window.size as new TAD

      # Delta score
      GenomicRanges::mcols(ind1_boundaries.gr)[[paste0("delta_score_",ind2)]] =  GenomicRanges::mcols(ind1_boundaries.gr)[[paste0("score_",ind1)]] - GenomicRanges::mcols(ind1_boundaries.gr)[[paste0("score_",ind2)]]

      nb_newTAD = length(ind1_boundaries.gr[ GenomicRanges::mcols(ind1_boundaries.gr)[[paste0("newTAD_",ind2)]]==TRUE])
      nb_TAD = length(ind1_boundaries.gr)
      message(paste0(nb_newTAD , "/", nb_TAD," (", round(nb_newTAD / nb_TAD * 100), "%)"," new boundaries compared to ", ind2))
    }

    output = append(output,list(ind1_boundaries.gr))
    output.names = append(output.names, paste0("TADdiff_",ind1))
    message(paste0("TADdiff: ", ind1, " ended succefully (duration: ",  round(difftime(Sys.time(), start_time, units = 'secs')), "sec)."))
    message("######################################################")
  }

  names(output) = output.names

  return(output)
}
