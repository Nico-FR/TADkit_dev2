#' @title Import Orca matrix into full chromosome matrix
#'
#' @description This function create an empty matrix of the full chromosome which is then filled with Orca matrix (i.e. 250x250 bins).
#'
#' @details From the matrices produced by orca_predict.genomepredict function.
#'
#' @param df_prediction.path The full path of the Orca matrix. Data frame without header which contains orca predictions, i.e. log(Obseved/Expected).
#' @param scale Size of the scale (i.e. window) of the predictions in base pair (i.e. 32Mb, 16Mb, 8Mb, 4Mb, 2Mb and 1Mb).
#' @param mpos The coordinate to zoom into for multiscale prediction. By default mpos is the center of the 32Mb sequence provided (i.e. start + 16Mb).
#' @param chromsize Chromosome size in base pair of the related chromosome.
#' @param output Character: Default is "OE" to return observed / expected counts. "Obs" to return observed counts.
#' @param df_normmats.path Optional. Background distance-based expected balanced contact matrices. Only needed to return observed counts.
#'
#' @return A `Matrix` class object: upper triangular and sparse Matrix
#'
#' @importFrom Matrix triu sparseMatrix
#' @importFrom methods as
#' @importFrom magrittr %>%
#'
#' @export
#'
#'
orca2matrix <- function(df_prediction.path, sep = "\t", mpos, scale, chromsize, output = "OE", df_normmats.path = NULL) {

  #matrix specifications
  bin.width = scale / 250
  nbins = ifelse(chromsize %/% bin.width == chromsize / bin.width, chromsize %/% bin.width, chromsize %/% bin.width + 1)  #nb bins of the final matrix
  mpos2 = ifelse(scale == 32e6, mpos, mpos - bin.width)
  bin_start = ifelse(mpos2 %/% bin.width <= 124,
                     1, mpos2 %/% bin.width - 124) #bin number of the first bin of orca matrix
  bin_end = bin_start + 249 #position of the last bin
  #orca_start_bp = (mpos - 126 * bin.width)
  #orca_start_bin = (mpos - 126 * bin.width) / bin.width + 1
  #orca_end_bp = mpos + 124 * bin.width
  message("Creating matrix of ", nbins, "x", nbins, " bins with ", output, " counts of bins ", bin_start, " to ", bin_end,
          " at ", bin.width / 1e3, "kb resolution (i.e. from ", (bin_start - 1) * bin.width / 1e6, "Mb to ", (bin_end) * bin.width / 1e6, "Mb).")

  #create empty matrix
  mat = Matrix::Matrix(0, nrow = nbins, ncol = nbins, sparse = TRUE)

  #read orca output
  if (output == "OE") {
    orca = read.table(df_prediction.path, h=F, sep=sep) %>% as.matrix() %>% exp
  }
  if (output == "Obs") {
    predictions = read.table(df_prediction.path, h=F, sep=sep) %>% as.matrix() %>% exp
    normats = read.table(df_normmats.path, h=F, sep=sep) %>% as.matrix()
    orca = predictions * normats
  }

  #add upper input in the full upper matrix
  mat[bin_start:bin_end, bin_start:bin_end][upper.tri(mat[bin_start:bin_end, bin_start:bin_end], diag = TRUE)] <- orca[upper.tri(orca, diag = TRUE)]

  return(methods::as(mat, "CsparseMatrix"))
  }



