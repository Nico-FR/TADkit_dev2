#' @title Plot TADs for one individual
#'
#' @description TADplot use Gviz package to plot at least 2 tracks. The first one includes all the TADs from the specified chromosome. And the second one is a zoom in the defined area.
#' Other tracks can be added like gene annotation or read depth data.
#'
#' @details This create a plot with at least 2 tracks from TADs annotation (tad.gr).
#'
#' On the first track the alternating TADs are represented in grey and black while the interTAD regions are represented in white.
#' The second track is a zoom in the area of interest. This area is extend to the border of the TADs on both sides.
#'
#' Another track from a bigwig file can be added like read depth sequencing (RNAseq...) as an histogram.
#' The bin size of the histogram is 1Kb by default, the read depth is smoothed using the median value of each bin.
#' If the chromosome name is different (between the bigwig and the chr parameter), it can be fix using the bigwig.chr parameter ("1" versus "chr1").
#'
#' Another track with any annotations can be added.
#' This track can be group using factors in a specified column of the annot.gr files (metadata) otherwise the names of each annotation is used.
#'
#' @param matrix matrix object (data frame or matrix) or matrix path at matrix format (dense) for only one chromosome. The path can be gzip (ending by .gz)
#' @param start start in bp of the area of interest.
#' @param stop end in bp of the area of interest.
#' @param bin.width bin width of the matrix in bp.
#' @param matrix.colname logical. Does your matrix file (ie path) have column names (ie header)? Default = TRUE
#' @param matrix.rowname logical. Does your matrix file (ie path) have row names? Default = TRUE
#' @param matrix.sep the field separator character. Values on each line of the matrix file are separated by this character. Default = "\t" (ie tabulation).
#' @param matrix.diag logical. Weather or not to plot diagonal of the matrix. Default = TRUE
#' @param log2 logical. Weather or not to plot the log2 of the matrix values.
#' @param scale.colors A character string indicating the color map option to use. Eight options are available (see viridis package:
#' "magma" (or "A")
#' "inferno" (or "B")
#' "plasma" (or "C")
#' "viridis" (or "D")
#' "cividis" (or "E")
#' "rocket" (or "F")
#' "mako" (or "G")
#' "turbo" (or "H")
#'
#' @return ggplot of the matrix.
#' @import reshape2
#' @import viridis
#' @import scales
#' @import ggplot2
#' @export
#'
#' @examples


MATplot <- function(matrix, start, stop, bin.width, matrix.colname = T, matrix.rowname = T, matrix.sep = "\t", matrix.diag = T, log2 = T, scale.colors = "H") {

  ##############################"
  matrix="/home/nmary/Downloads/Bovin-669.ARS-UCD1.2.mapq_10.10000.chr1.matrix.gz"
  bin.width = 10e3
  start = 1e6 ; stop = 4e6
  matrix.colname = T ; matrix.rowname = T ; matrix.sep = "\t"
  log2 = T ; scale.colors = "H"; matrix.diag = T
  ##############################


  #sanity check

  #read matrix
  from = start %/% bin.width + 1 ; to = stop %/% bin.width #nb bin
  matrix.row.skip <- ifelse(matrix.colname == T,  1, 0)
  matrix.col.skip <- ifelse(matrix.rowname == T,  1, 0)

  if (is.character(matrix) & rev(strsplit(matrix, split = "\\.")[[1]])[1] == "gz")  {
    df = read.table(gzfile(matrix), sep = matrix.sep, h = F, row.names = matrix.col.skip, skip = matrix.row.skip + from - 1)[1:(to - from + 1), from:to]
  }

  if (is.character(matrix) & !rev(strsplit(matrix, split = "\\.")[[1]])[1] == "gz")  {
    df = read.table(matrix, sep = matrix.sep, h = F, row.names = matrix.col.skip, skip = matrix.row.skip + from - 1)[1:(to - from + 1), from:to]
  }

  if (is.data.frame(matrix))  {df = matrix[from:to, from:to]}

  if(!is.matrix(matrix)) {mat = matrix(as.matrix(df), nrow = length(df))} else {mat = matrix[from:to, from:to]}
  if(isFALSE(matrix.diag)) {diag(mat) <- NA}

  #plot
  if (log2 == T) {mat = log2(mat)}
  melted_mat <- reshape2::melt(mat)
  melted_mat$Var2 = (melted_mat$Var2 + from - 1) * -bin.width + bin.width / 2
  melted_mat$Var1 = (melted_mat$Var1 + from - 1) * bin.width - bin.width / 2
  ggplot2::ggplot(data = melted_mat, ggplot2::aes(x = Var1, y = Var2, fill = value))+
    ggplot2::geom_tile()+viridis::scale_fill_viridis(na.value = "black", option = scale.colors)+
    ggplot2::scale_x_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6))+
    ggplot2::scale_y_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6))+
    ggplot2::coord_fixed()+ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(), legend.title = ggplot2::element_blank())

}


















