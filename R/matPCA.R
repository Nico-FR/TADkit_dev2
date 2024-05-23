#' @title perform principal component analysis f
#'
#' @description
#'
#' @details
#'
#' @param matrix `Matrix` or `matrix` object.
#' @param intput input matrix: "Obs" for observed count or "OE" (default) for Observed / Expected counts.
#' @param bin.width Bin width of the matrix in base pair.
#' @param seqname chromosome names as character, default = "1".
#' @param nb_PC number of principal component to compute, default = 1.
#'
#' @return `GRanges` bedgraph with PC score.
#'
#' @importFrom stats prcomp
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#'
#' @export
#'
#' @examples
#' matPCA(matrix = mat_HCT116_chr19_50kb,
#'     bin.width = 50e3,
#'     input = "Obs",
#'     seqname = "chr1")
#'

matPCA <- function(matrix, bin.width, input = "OE", seqname = "1", nb_PC = 1) {

  if(!inherits(matrix, c("Matrix", "matrix"))) {
    stop("input matrix is not a matrix or dgCMatrix object")}

  if (input == "Obs") {
    message("computing oberved / expected matrix...")
    matrix = matObsExp(matrix, output = "OE")
  }

  matrix = as.matrix(matrix)

  if (!isSymmetric(matrix)) {
    matrix[lower.tri(matrix)] <- t(matrix)[lower.tri(matrix)]
  }

  gap = colMeans(matrix, na.rm = TRUE) %in% c(0, "NaN", "NA") #bin without data

  matrix[matrix == 0] <- 1
  matrix[is.na(matrix)] <- 1

  matCorr = suppressWarnings(cor(matrix, method = "pearson"))
  diag(matCorr) <- NA
  matCorr[is.na(matCorr)] <- mean(matCorr, na.rm = T)

  eigen_result <- stats::prcomp(matCorr, rank. = nb_PC)
  df = data.frame(
    seqnames = seqname,
    start = (1:length(eigen_result$rotation[,1]) - 1) * bin.width,
    stop = 1:length(eigen_result$rotation[,1]) * bin.width)

  eigen_result$rotation[gap] <- NA #filter bin without data
  df = cbind(df, eigen_result$rotation)

  gr = GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)

  return(gr)
}












