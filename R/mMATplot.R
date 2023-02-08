#' @title Plot 2 matrix with annotations
#'
#' @description mMATplot allow to plot 2 matrix files on the upper or lower part of the matrix.
#' It is also possibles to add 3 types of annotations:
#' -domains (like TADs or compartments) plot as triangles or lines on the upper or/and lower part of the matrix.
#' -interaction between bins (loop) plot as a square on the upper and lower part of the matrix.
#'
#' @details The matrix datas can be R object (dataframe, matrix or array) or the path of the files.
#' In that case you can specify is the files as column and row names usely used to write the bin numbers or coordonates.
#' All bed files can be R object (dataframe or GRange) or the path of the files.
#' All bedpe files can be a dataframe or the path of the files.
#' Chromosome of bedpe and bed files can be filter using tad.chr parameter.
#' Upper and lower line can be used to represent compartments A and B and colored accordingly using the tad.line.col parameter.
#'
#' @param matrix.upper matrix object (data frame or matrix) or matrix path at matrix format (dense) for only one chromosome. The path can be gzip (ending by .gz)
#' @param matrix.lower matrix object (data frame or matrix) or matrix path at matrix format (dense) for only one chromosome. The path can be gzip (ending by .gz)
#' @param start start in bp of the area of interest.
#' @param stop end in bp of the area of interest.
#' @param bin.width bin width of the matrix in bp.
#' @param matrix.colname logical. Does your matrix file (ie path) have column names (ie header)? Default = TRUE
#' @param matrix.rowname logical. Does your matrix file (ie path) have row names? Default = FALSE
#' @param matrix.upper.txt text to write on the upper part of the plot
#' @param matrix.lower.txt text to write on the lower part of the plot
#' @param matrix.sep the field separator character. Values on each line of the matrix file are separated by this character. Default = '\\t' (ie tabulation).
#' @param log2 logical. Use the log2 of the matrix values. Default is TRUE
#' @param scale.colors A character string indicating the color map option to use. Eight colors palettes are available from viridis package. Another palette "OE" is made for data centered on 0 (ie log2(observed/expected) matrix). Default is "H":
#' "ObsExp" (or "OE")
#' "magma" (or "A")
#' "inferno" (or "B")
#' "plasma" (or "C")
#' "viridis" (or "D")
#' "cividis" (or "E")
#' "rocket" (or "F")
#' "mako" (or "G")
#' "turbo" (or "H")
#' @param tad.upper.tri bed files path, data frame or grange object with the TAD to plot as triangle in the upper part of the matrix. Default is NULL
#' @param tad.lower.tri bed files path, data frame or grange object with the TAD to plot as triangle in the upper part of the matrix. Default is NULL
#' @param tad.upper.line bed files path, data frame or grange object with the TAD to plot as line in the upper parts of the matrix. Default is NULL
#' @param tad.lower.line bed files path, data frame or grange object with the TAD to plot as line in the lower parts of the matrix. Default is NULL
#' @param tad.line.col col number of the tad.line files that contain factors used to color tad.upper.line and tad.lower.line. Default is NULL
#' @param loop.upper.bedpe bedpe files path or data frame to plot on both parts of the matrix. Six columns table (chr1, start1, end1, chr2, start2, end2) that gives loops between 2 areas. Default is NULL
#' @param tad.chr chromosome name to filter bed and bedpe files. Default is NULL
#' @param annotations.color color for loop and tri annotations. Default is "red".
#' @param line.color colors for upper and lower lines.

#'
#' @return ggplot of the matrix.
#' @import reshape2
#' @import viridis
#' @import scales
#' @import ggplot2
#' @export
#'
#' @examples


mMATplot <- function(matrix.upper, matrix.lower, start, stop, bin.width, log2 = T, scale.colors = "H",
                    matrix.colname = T, matrix.rowname = F, matrix.sep = "\t",
                    matrix.upper.txt = NULL, matrix.lower.txt = NULL,
                    tad.upper.tri = NULL, tad.lower.tri = NULL, loop.bedpe = NULL, tad.chr = NULL, annotations.color = "red",
                    tad.upper.line = NULL, tad.lower.line = NULL, tad.line.col = NULL, line.colors = c("red", "blue")) {

  #bin to read
  from = start %/% bin.width + 1 ; to = stop %/% bin.width #nb bin

  ####################
  #mat1
  #bin to read for matrix.path
  matrix.row.skip <- ifelse(matrix.colname == T,  from, from - 1)
  if(isTRUE(matrix.rowname)) matrix.col.skip <- 1 else matrix.col.skip <- NULL

  #read matrix path in data.frame
  if (is.character(matrix.upper))  {
    if (substr(matrix.upper, nchar(matrix.upper) - 2, nchar(matrix.upper)) == ".gz") {
      df = read.table(gzfile(matrix.upper), sep = matrix.sep, h = F, row.names = matrix.col.skip, skip = matrix.row.skip)[1:(to - from + 1), from:to]
    } else {
      df = read.table(matrix.upper, sep = matrix.sep, h = F, row.names = matrix.col.skip, skip = matrix.row.skip)[1:(to - from + 1), from:to]
    }
  }

  #read data frame
  if (is.data.frame(matrix.upper))  {df = matrix[from:to, from:to]}

  #read matrix or create one
  if(!is.matrix(matrix.upper)) {mat.upper = matrix(as.matrix(df), nrow = length(df))} else {mat.upper = matrix[from:to, from:to]}
  diag(mat.upper) <- NA
  ####################

  #create empty matrix
  mat = matrix(ncol=ncol(mat.upper),nrow=ncol(mat.upper))

  #write upper matrix
  mat[upper.tri(mat)] = mat.upper[upper.tri(mat.upper)]

  ####################
  #mat2
  #read matrix path in data.frame
  if (is.character(matrix.lower))  {
    if (substr(matrix.lower, nchar(matrix.lower) - 2, nchar(matrix.lower)) == ".gz") {
      df = read.table(gzfile(matrix.lower), sep = matrix.sep, h = F, row.names = matrix.col.skip, skip = matrix.row.skip)[1:(to - from + 1), from:to]
    } else {
      df = read.table(matrix.lower, sep = matrix.sep, h = F, row.names = matrix.col.skip, skip = matrix.row.skip)[1:(to - from + 1), from:to]
    }
  }

  #read data frame
  if (is.data.frame(matrix.lower))  {df = matrix[from:to, from:to]}

  #read matrix or create one
  if(!is.matrix(matrix.lower)) {mat.lower = matrix(as.matrix(df), nrow = length(df))} else {mat.lower = matrix[from:to, from:to]}
  diag(mat.lower) <- NA
  ####################

  #write lower matrix
  mat[lower.tri(mat)] = mat.lower[lower.tri(mat.lower)]

  # matrix plot
  if (log2 == T) {mat = log2(mat)}
  melted_mat <- reshape2::melt(t(mat))
  melted_mat$Var2 = (melted_mat$Var2 + from - 1) * - bin.width + bin.width / 2
  melted_mat$Var1 = (melted_mat$Var1 + from - 1) * bin.width - bin.width / 2

  if (scale.colors == "OE" | scale.colors == "ObsExp") {
    p <- ggplot2::ggplot()+ggplot2::geom_tile(data = melted_mat, ggplot2::aes(x = Var1, y = Var2, fill = value))+
      ggplot2::scale_fill_gradient2(low = "blue", high = "red",midpoint = 0, mid="white", na.value = "white")+
      ggplot2::scale_x_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6), limits = c(start, stop))+
      ggplot2::scale_y_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6), limits = c(-stop, -start))+
      ggplot2::coord_fixed()+ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(), legend.title = ggplot2::element_blank())

  } else {
    p <- ggplot2::ggplot()+ggplot2::geom_tile(data = melted_mat, ggplot2::aes(x = Var1, y = Var2, fill = value))+
      viridis::scale_fill_viridis(na.value = "black", option = scale.colors)+
      ggplot2::scale_x_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6), limits = c(start, stop))+
      ggplot2::scale_y_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6), limits = c(-stop, -start))+
      ggplot2::coord_fixed()+ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(), legend.title = ggplot2::element_blank())
  }


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
    tad$col = as.factor(tad$col)
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
    tad$col = as.factor(tad$col)
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


  if (!is.null(matrix.upper.txt)) {
    p = p + geom_text(x = to, y = -from, label = matrix.upper.txt, col="red", size = 10, hjust = "right", vjust = "top")
  }

  if (!is.null(matrix.lower.txt)) {
    p = p + geom_text(x = from, y = -to, label = matrix.lower.txt, col="red", size = 10, hjust="right", vjust="top")
  }

  p
}





















