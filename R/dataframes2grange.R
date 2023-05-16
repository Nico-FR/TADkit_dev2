#' Create GRanges with metadatas
#'
#' This function create a GRange object from the first `dataframe` (annotation.table) and add the size of the chromosomes on the second `dataframe`.
#' The chromosomes are filtered according to the list of chromosomes contained in the second `dataframe` (chromsize).
#'
#' @param annotation.table `dataframe` with genomic annotations (at least 3 columns: chromosome, start and stop).
#' @param chromsize `dataframe` with chromosome names and sizes. Chromosomes are also filtered according to this list.
#' @param chr.col Column number for chromosome names, default is 1.
#' @param start.col Column number for start, default is 2.
#' @param end.col Column number for end, default is 3.
#' @param strand.col Column number for strands, default is `NULL` to set strands as "*" (i.e both strands or unknown stand) otherwise the column must contain "+" or "-" characters.
#' @param name.col Column number for annotation names, default is `NULL` to create names with chr.names and start.
#' @param metadata.mcols Column number(s) for metadata, default is `NULL` (i.e no metadata).
#'
#'
#' @return GRange file that contains chromosome size.
#' @import GenomicRanges
#' @importFrom S4Vectors Rle
#' @importFrom IRanges IRanges
#' @import GenomeInfoDb
#'
#' @export
#'
#' @examples
#' bed.df = data.frame(chr = 1, start = 198e3, end = 290e3, strand = "+", names = "geneID")
#' chromsize.df = data.frame(chr = "1", size = 400e3)
#' dataframes2grange(bed.df, chromsize.df, strand.col = 4,name.col = 5)
#'
#' dataframes2grange(TADkit::tad_1_10kb.bed, TADkit::chromsize)


dataframes2grange <- function(annotation.table, chromsize, chr.col = 1, start.col = 2, end.col = 3, strand.col = NULL, name.col = NULL, metadata.mcols = NULL) {
  data = annotation.table[annotation.table[, chr.col] %in% chromsize[,1] , ]
  temp <- GenomicRanges::GRanges(
    seqnames <- S4Vectors::Rle(as.factor(data[, chr.col])) %>% droplevels(),
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
