#' Create GRange with sequence lengths
#'
#' This function create a GRange object from the first dataframe and add the size of the chromosomes from the second dataframe.
#' A name of each line is added to the GRange object.
#'
#' @param annotation.table Data frame with genomic annotations (at least 3 columns: chromosome, start and stop).
#' @param chromsize Data frame with chromosome names and lengths. Chromosomes are also filtered according to this list.
#' @param chr.col Column number for chromosome names, default=1.
#' @param start.col Column number for start, default=2.
#' @param end.col Column number for end, default=3.
#' @param strand.col Column number for strands, default=NULL to set strands as "*" (instead of "+" or "-").
#' @param name.col Column number for annotation names, default=NULL to create names with chr.names and start.
#' @param metadata.mcols Column number(s) for metadata, default=NULL (ie no metadata).
#'
#'
#' @return GRange file that contains chromosome size.
#' @import GenomicRanges
#' @import S4Vectors Rle
#' @import IRanges
#' @import GenomeInfoDb
#' @importFrom BiocGenerics strand end start width
#' @export


dataframes2grange <- function(annotation.table, chromsize, chr.col = 1, start.col = 2, end.col = 3, strand.col = NULL, name.col = NULL, metadata.mcols = NULL) {
  data = annotation.table[annotation.table[, chr.col] %in% chromsize[,1] , ]
  temp <- GenomicRanges::GRanges(
    seqnames <- S4Vectors::Rle(data[, chr.col]),
    ranges <- IRanges::IRanges(
      start = data[, start.col],
      end = data[, end.col],
      names =
        if (is.null(name.col)) {
          paste0(data[, chr.col], "_", data[, start.col])
        } else {
          data[, name.col]
        }
    ),
    strand =
      if (is.null(strand.col)) {
        S4Vectors::Rle(BiocGenerics::strand("*"))
      } else {
        data[, strand.col]
      },
    mcols =
      if (is.null(metadata.mcols)) {
        NULL
      } else {
        data[, c(metadata.mcols)]
      }
  )

  names(mcols(temp)) <- names(data)[c(metadata.mcols)]

  for (i in 1:length(GenomeInfoDb::seqlengths(temp))) {
    GenomeInfoDb::seqlengths(temp)[i] <- as.numeric(chromsize[, 2][chromsize[, 1] == names(GenomeInfoDb::seqlengths(temp))[i]])
  }
  return(sort(temp))
}
