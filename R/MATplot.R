#' @title Plot interaction matrix
#'
#' @description `MATplot` allow to plot matrix with 2 types of annotations:
#'     * domains (e.g. TADs or compartments) are plot as triangles and/or lines on the upper or/and lower part of the matrix.
#'     * interactions between domains/bins (loop) are plot as squares on the upper and lower part of the matrix.
#'
#' @details The matrix input must be a `dgCMatrix` or a `matrix` object for only one chromosome (see `coolFetch` function to read cool files).
#' All domains (TADs or compartments) are bed files (3 columns: chr, start and end) and can be R object (`dataframe` or `GRanges`) or the path of the files.
#' For `tad.lines`, another column can be used to specify different classes of domains (e.g compartment A or B). To use those domain classes, specify the column number (of the `tad.upper.line` and `tad.lower.line` inputs) with `tad.line.col` parameter and a custom set of colors with `line.colors` parameter.
#' Loop are stored in bedpe files (6 columns: chr1, start1, end, chr2, start2 and end2) and can be a `dataframe` object or the path of the file.
#' Chromosome domains and loops can be filter using `tad.chr` parameter.
#'
#' @inheritParams matObsExp
#' @inheritParams TADplot
#' @param matrix `dgCMatrix` or `matrix` object for only one chromosome.
#' @param bin.width Bin width of the matrix in base pair.
#' @param matrix.diag logical. Weather or not to plot diagonal values of the matrix. Default = `TRUE`.
#' @param log2 logical. Use the log2 of the matrix values. Default is `TRUE`.
#' @param scale.colors A character string indicating the color map option to use. Eight colors palettes are available from `viridis` package. Another palette `"OE"` is made for data centered on 0 (e.g log2 of observed/expected matrix). Default is `"H"`:
#' * `"ObsExp"` (or `"OE"`),
#' * `"magma"` (or `"A"`),
#' * `"inferno"` (or `"B"`),
#' * `"plasma"` (or `"C"`),
#' * `"viridis"` (or `"D"`),
#' * `"cividis"` (or `"E"`),
#' * `"rocket"` (or `"F"`),
#' * `"mako"` (or `"G"`),
#' * `"turbo"` (or `"H"`).
#' @param tad.upper.tri,tad.lower.tri `data.frame`, `GRanges` or the bed files path with the TAD to plot as triangle in the upper or lower part of the matrix. Default is `NULL`.
#' @param tad.upper.line,tad.lower.line `data.frame`, `GRanges` or the bed files path with the TAD to plot as line in the upper or lower parts of the matrix. Default is `NULL`.
#' @param tad.line.col Column number of `tad.upper.line` and `tad.lower.line` files that contain factors used to color upper and lower lines. Default is `NULL`.
#' @param loop.bedpe `data.frame` or bedpe files path to plot on both parts of the matrix. Six columns table (chr1, start1, end1, chr2, start2, end2) that gives loops between 2 areas. Default is `NULL`
#' @param tad.chr Chromosome name to filter annotations (domains and loop). Default is `NULL`
#' @param annotations.color Color for loop and triangular annotations. Default is `"red"`.
#' @param line.colors Colors for `tad.upper.line` and `tad.lower.line`. Default is `c("red", "blue")`.

#'
#' @return `ggplot`
#'
#' @importFrom Matrix triu summary diag tril
#' @importFrom viridis scale_fill_viridis
#' @importFrom scales unit_format
#' @importFrom dplyr filter
#' @importFrom BiocGenerics t
#' @import ggplot2
#' @export
#'
#' @examples
#' MATplot(matrix_1_chr25_50kb,
#'     start = 10e6, stop = 30e6,
#'     bin.width = 50e3, log2 = TRUE,
#'     scale.colors = "H", #color of matrix, try "D" or "H"
#'     annotations.color = "red")

MATplot <- function(matrix, start, stop, bin.width, log2 = T, scale.colors = "H", matrix.diag = T,
                    tad.upper.tri = NULL, tad.lower.tri = NULL, loop.bedpe = NULL, tad.chr = NULL, annotations.color = "red",
                    tad.upper.line = NULL, tad.lower.line = NULL, tad.line.col = NULL, line.colors = c("red", "blue")) {

  #sanity check
  if(!inherits(matrix, c("dgCMatrix", "matrix"))) {
    stop("input matrix is not a matrix or dgCMatrix object")}

  #bin to read
  from = start %/% bin.width + 1 ; to = stop %/% bin.width #nb bin

  #filter matrix area
  mat = Matrix::triu(matrix[from:to, from:to])
  mat[Matrix::triu(mat == 0)] <- NA

  #get log2
  if (log2 == T) {mat@x = log2(mat@x)}

  #melt matrix
  upper_mat = Matrix::summary(Matrix::triu(mat, 1))
  diag_mat = data.frame(i = 1:nrow(mat), j = 1:nrow(mat), x = Matrix::diag(mat))
  lower_mat = Matrix::summary(Matrix::tril(BiocGenerics::t(mat), -1))

  #diag
  if(matrix.diag) {melted_mat = rbind(upper_mat, lower_mat, diag_mat)} else {melted_mat = rbind(upper_mat, lower_mat)}

  #add genomic coordinates
  melted_mat$j = (melted_mat$j + from - 1) * bin.width - bin.width / 2
  melted_mat$i = (melted_mat$i + from - 1) * - bin.width + bin.width / 2

  if (scale.colors == "OE" | scale.colors == "ObsExp") {
    p <- ggplot2::ggplot()+ggplot2::geom_tile(data = melted_mat, ggplot2::aes(y = i, x = j, fill = x))+
      ggplot2::scale_fill_gradient2(low = "blue", high ="red",midpoint = 0, mid="white", na.value = "white")+
      ggplot2::scale_x_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6), limits = c(start, stop))+
      ggplot2::scale_y_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6), limits = c(-stop, -start))+
      ggplot2::coord_fixed()+ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(), legend.title = ggplot2::element_blank())

  } else {
    p <- ggplot2::ggplot()+ggplot2::geom_tile(data = melted_mat, ggplot2::aes(y = i, x = j, fill = x))+
      viridis::scale_fill_viridis(na.value = "black", option = scale.colors)+
      ggplot2::scale_x_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6), limits = c(start, stop))+
      ggplot2::scale_y_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6), limits = c(-stop, -start))+
      ggplot2::coord_fixed()+ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(), legend.title = ggplot2::element_blank())
  }


  #upper tri
  if (!is.null(tad.upper.tri)) {

    if (is.character(tad.upper.tri)) {
      tad = read.table(tad.upper.tri, header = F, sep = "\t")[,1:3]
    }

    if (inherits(tad.upper.tri, "GRanges")) {
      tad = as.data.frame(tad.upper.tri, row.names = NULL)[,1:3]
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
      tad = read.table(tad.lower.tri, header = F, sep = "\t")[,1:3]
    }

    if (inherits(tad.lower.tri, "GRanges")) {
      tad = as.data.frame(tad.lower.tri, row.names = NULL)[,1:3]
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
      tad = read.table(tad.upper.line, header = F, sep = "\t")[, c(1:3, tad.line.col)]
    }

    if (inherits(tad.upper.line, "GRanges")) {

      if (is.null(tad.line.col)) {
        tad = as.data.frame(tad.upper.line, row.names = NULL)[, c(1:3)]
      } else {
        temp = as.data.frame(tad.upper.line, row.names = NULL)

        if (length(temp) < (tad.line.col + 5)) {
          stop(paste0("There is only ", length(temp) - 5, " column(s) with metadata in tad.upper.line input but you select column ", tad.line.col))
        }

        tad = temp[, c(1:3, tad.line.col + 5)]
      }
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
      tad = read.table(tad.lower.line, header = F, sep = "\t")[, c(1:3, tad.line.col)]
    }

    if (inherits(tad.lower.line, "GRanges")) {

      if (is.null(tad.line.col)) {
        tad = as.data.frame(tad.lower.line, row.names = NULL)[, c(1:3)]
      } else {
        temp = as.data.frame(tad.lower.line, row.names = NULL)

        if (length(temp) < (tad.line.col + 5)) {
          stop(paste0("There is only ", length(temp) - 5, " column(s) with metadata in tad.lower.line input but you select column ", tad.line.col))
        }

        tad = temp[, c(1:3, tad.line.col + 5)]
      }
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
      loop = read.table(loop.bedpe, header = F, sep = "\t")[,1:6]
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

  p
}





















