#' @title Distribution of annotation density within domains
#'
#' @description Graph of the distribution of density of annotations (like genes) within domains (like TADs, compartments...)
#' calculates the density of annotations per bin and and then measure the relative position of each bin among domain.

#'
#'
#' @details As an exemple, this function take all TAD domains and count the relative position of all annotation features (start, end or even the center).
#' For example:
#' 1- calculates the density of genes per bin of 10kb,
#' 2- calculate the relative positions of each bin within is corresponding TAD,
#' 3- plot the bin density with relative positions.
#'
#'
#' @param domain.gr GRange file with domains (TADs, compartments...).
#' @param annot.gr GRange file with genomic annotations (genes...).
#' @param annot.type Default is NULL. Otherwise specify a column number used to separate annotations. Use "strand" to separate annotation acording to their strands.
#' @param bin.width Size of the bin in bp. This should match the size of the bins used to determine the domains.
#' @param output Default is "plot" to return a ggplot. Use "data" to return the datas used to produce the plot.
#'
#'
#' @return Return a ggplot graph
#' @import GenomeInfoDb
#' @importFrom BiocGenerics strand end start width
#' @importFrom plyr ddply
#' @import ggplot2
#' @import scales
#' @import GenomicRanges
#' @export
#'
#' @examples
#'
domain_density <- function(domain.gr, annot.gr, annot.type = NULL, bin.width = 10e3, output = "plot") {

}

