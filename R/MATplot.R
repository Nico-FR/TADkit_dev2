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
#' @param matrix.sep the field separator character. Values on each line of the matrix file are separated by this character. Default = '\\t' (ie tabulation).
#' @param matrix.diag logical. Weather or not to plot diagonal of the matrix. Default = TRUE
#' @param log2 logical. Weather or not to plot the log2 of the matrix values.
#' @param tad.upper.tri bed files path, data frame or grange object with the TAD to plot as triangle in the upper part of the matrix. Default is NULL
#' @param tad.lower.tri bed files path, data frame or grange object with the TAD to plot as triangle in the upper part of the matrix. Default is NULL
#' @param tad.upper.line bed files path, data frame or grange object with the TAD to plot as line in the upper parts of the matrix. Default is NULL
#' @param tad.lower.line bed files path, data frame or grange object with the TAD to plot as line in the lower parts of the matrix. Default is NULL
#' @param tad.line.col col number use to color tad.upper.line and tad.lower.line. Default is NULL
#' @param loop.upper.bedpe bedpe files path or data frame plot on both parts of the matrix. Six columns table (chr1, start1, end1, chr2, start2, end2) that gives loops between 2 areas. Default is NULL
#' @param tad.chr chromosome name to filter in TAD files. Default is NULL
#' @param annotations.color color for loop and tri annotations. Default is "red".
#' @param line.color colors for
#' @param scale.colors A character string indicating the color map option to use. Eight options are available (see viridis package), default is "H":
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


MATplot <- function(matrix, start, stop, bin.width, matrix.colname = T, matrix.rowname = T, matrix.sep = "\t", matrix.diag = T, log2 = T,
                    tad.upper.tri = NULL, tad.lower.tri = NULL, loop.bedpe = NULL,
                    tad.chr = NULL, annotations.color = "red", line.colors = c("red", "bleu"), scale.colors = "H") {

  #bin to read
  from = start %/% bin.width + 1 ; to = stop %/% bin.width #nb bin

  #bin to read for matrix.path
  matrix.row.skip <- ifelse(matrix.colname == T,  from, from - 1)
  if(isTRUE(matrix.rowname)) matrix.col.skip <- 1 else matrix.col.skip <- NULL

  #read matrix
  if (is.character(matrix) & rev(strsplit(matrix, split = "\\.")[[1]])[1] == "gz")  {
    df = read.table(gzfile(matrix), sep = matrix.sep, h = F, row.names = matrix.col.skip, skip = matrix.row.skip)[1:(to - from + 1), from:to]
  }

  if (is.character(matrix) & !rev(strsplit(matrix, split = "\\.")[[1]])[1] == "gz")  {
    df = read.table(matrix, sep = matrix.sep, h = F, row.names = matrix.col.skip, skip = matrix.row.skip)[1:(to - from + 1), from:to]
  }

  if (is.data.frame(matrix))  {df = matrix[from:to, from:to]}

  if(!is.matrix(matrix)) {mat = matrix(as.matrix(df), nrow = length(df))} else {mat = matrix[from:to, from:to]}
  if(isFALSE(matrix.diag)) {diag(mat) <- NA}

  # matrix plot
  if (log2 == T) {mat = log2(mat)}
  melted_mat <- reshape2::melt(mat)
  melted_mat$Var2 = (melted_mat$Var2 + from - 1) * - bin.width + bin.width / 2
  melted_mat$Var1 = (melted_mat$Var1 + from - 1) * bin.width - bin.width / 2
  p <- ggplot2::ggplot()+ggplot2::geom_tile(data = melted_mat, ggplot2::aes(x = Var1, y = Var2, fill = value))+
    viridis::scale_fill_viridis(na.value = "black", option = scale.colors)+
    ggplot2::scale_x_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6), limits = c(start, stop))+
    ggplot2::scale_y_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6), limits = c(-stop, -start))+
    ggplot2::coord_fixed()+ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(), legend.title = ggplot2::element_blank())


  #upper tri
  if (!is.null(tad.upper.tri)) {

    if (is.character(tad.upper.tri)) {
      tad = read.table(tad.upper.tri, h = F, sep = "\t")[,1:3]
    }

    if (class(tad.upper.tri) == "GRanges") {
      tad = as.data.frame(tad.upper.tri)[,1:3]
    }

    if (is.data.frame(tad.upper.tri)) {
      tad = tad.upper.tri[,1:3]
    }

    names(tad) = c("chr", "s", "e")
    if (is.null(tad.chr)) {
      tad <- dplyr::filter(tad, e > start, s < stop)} else {
        tad <- dplyr::filter(tad, chr == tad.chr, e > start, s < stop)}

    tad$e2 <- ifelse(tad$e >= stop, stop, tad$e)
    tad$s2 <- ifelse(tad$s <= start, start, tad$s)

    p = p + ggplot2::geom_segment(data = tad, ggplot2::aes(x = s, y = -s, xend = e2, yend = -s2), color = annotations.color, size = 0.3)+
      ggplot2::geom_segment(data = tad, ggplot2::aes(x = e2, y = -s2, xend = e, yend = -e), color = annotations.color, size = 0.3)
  }

  #lower tri
  if (!is.null(tad.lower.tri)) {

    if (is.character(tad.lower.tri)) {
      tad = read.table(tad.lower.tri, h = F, sep = "\t")[,1:3]
    }

    if (class(tad.lower.tri) == "GRanges") {
      tad = as.data.frame(tad.lower.tri)[,1:3]
    }

    if (is.data.frame(tad.lower.tri)) {
      tad = tad.lower.tri[,1:3]
    }

    names(tad) = c("chr", "s", "e")
    if (is.null(tad.chr)) {
      tad <- dplyr::filter(tad, e > start, s < stop)} else {
        tad <- dplyr::filter(tad, chr == tad.chr, e > start, s < stop)}

    tad$e2 <- ifelse(tad$e >= stop, stop, tad$e)
    tad$s2 <- ifelse(tad$s <= start, start, tad$s)

    p = p + ggplot2::geom_segment(data = tad, ggplot2::aes(x = s2, y = -e2, xend = e, yend = -e), color = annotations.color, size = 0.3)+
      ggplot2::geom_segment(data = tad, ggplot2::aes(x = s, y = -s, xend = s2, yend = -e2), color = annotations.color, size = 0.3)
  }


  #upper line
  if (!is.null(tad.upper.line)) {

    if (is.character(tad.upper.line)) {
      tad = read.table(tad.upper.line, h = F, sep = "\t")[, c(1:3, tad.line.col)]
    }

    if (class(tad.upper.line) == "GRanges") {
      tad = as.data.frame(tad.upper.line)[, c(1:3, tad.line.col + 5)]
    }

    if (is.data.frame(tad.upper.line)) {
      tad = tad.upper.line[, c(1:3, tad.line.col)]
    }

    if (is.null(tad.line.col)) {
      tad$col = rep(c("1", "2"), length.out=nrow(tad))
    }

    names(tad) = c("chr", "s", "e", "col")
    if (is.null(tad.chr)) {
      tad <- dplyr::filter(tad, e > start, s < stop)} else {
      tad <- dplyr::filter(tad, chr == tad.chr, e > start, s < stop)
      }

    tad$e2 <- ifelse(tad$e >= stop, stop, tad$e)
    tad$s2 <- ifelse(tad$s <= start, start, tad$s)

    p = p + ggplot2::geom_segment(data = tad, ggplot2::aes(x = s2, y = -start, xend = e2, yend = -start, col = col), size = 1.5)+
      ggplot2::geom_segment(data = tad, ggplot2::aes(x = stop, y = -s2, xend = stop, yend = -e2, col = col), size = 1.5)+
      ggplot2::scale_colour_manual(values = line.colors)
  }

  #lower line
  if (!is.null(tad.lower.line)) {

    if (is.character(tad.lower.line)) {
      tad = read.table(tad.lower.line, h = F, sep = "\t")[, c(1:3, tad.line.col)]
    }

    if (class(tad.lower.line) == "GRanges") {
      tad = as.data.frame(tad.lower.line)[, c(1:3, tad.line.col + 5)]
    }

    if (is.data.frame(tad.lower.line)) {
      tad = tad.lower.line[, c(1:3, tad.line.col)]
    }

    if (is.null(tad.line.col)) {
      tad$col = rep(c("1", "2"), length.out=nrow(tad))
    }

    names(tad) = c("chr", "s", "e", "col")
    if (is.null(tad.chr)) {
      tad <- dplyr::filter(tad, e > start, s < stop)} else {
        tad <- dplyr::filter(tad, chr == tad.chr, e > start, s < stop)
      }

    tad$e2 <- ifelse(tad$e >= stop, stop, tad$e)
    tad$s2 <- ifelse(tad$s <= start, start, tad$s)

    if (is.null(tad.chr)) {
      tad <- dplyr::filter(tad, e > start, s < stop)} else {
        tad <- dplyr::filter(tad, chr == tad.chr, e > start, s < stop)}

    p = p + ggplot2::geom_segment(data = tad, ggplot2::aes(x = start, y = -s2, xend = start, yend = -e2, col = col), size = 1)+
      ggplot2::geom_segment(data = tad, ggplot2::aes(x = s2, y = -stop, xend = e2, yend = -stop, col = col), size = 1)+
      ggplot2::scale_colour_manual(values = line.colors)
  }

  #loop
  if (!is.null(loop.bedpe)) {

    if (is.character(loop.bedpe)) {
      loop = read.table(loop.bedpe, h = F, sep = "\t")[,1:6]
    }

    if (is.data.frame(loop.bedpe)) {
      loop = loop.bedpe[,1:6]
    }

    names(loop) = c("chr1", "start1", "end1", "chr2", "start2", "end2")
    if (is.null(tad.chr)) {
      loop <- dplyr::filter(loop, start1 >= start, end1 < stop, start2 > start, end2 < stop)} else {
        loop <- dplyr::filter(loop, chr1 == tad.chr, chr2 == tad.chr, start1 >= start, end1 < stop, start2 > start, end2 < stop)}


    p = p + ggplot2::geom_rect(data = loop, ggplot2::aes(xmin = start2, xmax = end2, ymin = -start1, ymax = -end1), fill = NA, color = annotations.color, size = 0.3)+
      ggplot2::geom_rect(data = loop, ggplot2::aes(xmin = start1, xmax = end1, ymin = -start2, ymax = -end2), fill = NA, color = annotations.color, size = 0.3)
  }

  return(p)
}





















