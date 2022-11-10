#' @title Compartments calling (A or B)
#'
#' @description From a score for each bin (i.e principal component 1), the function output a grange with compartment A for score >= 0 and compartment B for score < 0
#'
#' @details From bedgraph as input (GRange object or path) output a GRange with the compartments.
#' NA score are considered as part of a compartment A if they are between two compartments A (the other way around for B compartments).
#' Otherwise (if NA are between compartment A and B) they are not called, thus leaving a gap between compartments.
#'
#'
#' @param bedgraph GRange file or bedgraph path for PC1 values.
#'
#' @return grange
#' @import GenomeInfoDb
#' @importFrom BiocGenerics sort
#' @importFrom dplyr mutate case_when
#' @importFrom S4Vectors mcols split
#' @import GenomicRanges
#' @export
#' @examples
#' # output <- geneTADtopo(tad.gr, gene.gr)
#' # plot(output)
#'

PC1calling <- function(bedgraph) {
  if (is.character(bedgraph)) {
    bedgraph = read.table(bedgraph, h=F, sep="\t") %>% dplyr::mutate(comp = dplyr::case_when(V4 < 0 ~ 'B', V4 >= 0 ~ 'A', is.na(V4) ~ "AB"))

    grange = GenomicRanges::makeGRangesFromDataFrame(bedgraph, start.field="V2", end.field = "V3", seqnames.field = "V1", keep.extra.columns = T)
  }

  if (class(bedgraph)=="GRanges") {
    grange = bedgraph %>% as.data.frame %>% dplyr::mutate(comp = dplyr::case_when(mcols(bedgraph)[,1] < 0 ~ 'B', S4Vectors::mcols(bedgraph)[,1] >= 0 ~ 'A', is.na(S4Vectors::mcols(bedgraph)[,1]) ~ "AB")) %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  }

  data1 = S4Vectors::split(grange, grange$comp)

  if (length(data1$AB) > 0){
    B = GenomicRanges::reduce(c(data1$B, GenomicRanges::reduce(data1$AB)[countOverlaps(GenomicRanges::reduce(data1$AB), data1$B) == 2]))

    A = GenomicRanges::reduce(c(data1$A, GenomicRanges::reduce(data1$AB)[countOverlaps(GenomicRanges::reduce(data1$AB), data1$A) == 2]))
  }

  if (length(data1$AB) == 0){
    B = GenomicRanges::reduce(data1$B)

    A = GenomicRanges::reduce(data1$A)
  }


  A$comp = "A"
  B$comp = "B"

  return(BiocGenerics::sort(c(A,B)))
}
