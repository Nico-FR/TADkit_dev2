#' @title observed / expected matrix
#'
#' @description For each bin of the matrix (interaction count observed) `matObsExp()` return the ratio: observed / expected.
#' The input matrix must be a `dgCMatrix` or `matrix` object (i.e matrix with as many rows and columns than the number of bin).
#' The output can be plot with `MATplot(output, log2 = TRUE, scale.colors = "ObsExp")`.
#'
#' @details The expected number of interaction corresponds to the average of the interaction counts according to the bin distances.
#' Only the upper part of the matrix is used.
#'
#' @param matrix `dgCMatrix` or `matrix` object.
#'
#' @return `dgCMatrix` object: upper triangular and sparse Matrix
#'
#' @importFrom stats toeplitz
#' @importFrom Matrix triu
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom methods as
#'
#' @export
#'
#' @examples
#' mat_obsexp = matObsExp(matrix_1_chr25_50kb)
#' MATplot(matrix = mat_obsexp,
#'     start = 10e6, stop = 30e6,
#'     bin.width = 50e3,
#'     log2 = TRUE, scale.colors = "OE")
#'

matObsExp <- function(matrix) {

  if(!inherits(matrix, c("dgCMatrix", "matrix"))) {
    stop("input matrix is not a matrix or dgCMatrix object")}

  #using dgcmatrix object
  #mat = Matrix::triu(matrix)
  #mat[Matrix::triu(mat == 0)] <- NA

  #using matrix (i.e faster)
  mat =  matrix(matrix, nrow = nrow(matrix), ncol = ncol(matrix))
  #mat[mat == 0] <- NA

  #mean diag
  mean_diag = sapply(1:(ncol(mat) - 1),
                     function (x){mean(diag(mat[,x:ncol(mat)]), na.rm = T)})

  mat_expected = stats::toeplitz(c(mean_diag, mat[1,ncol(mat)])) #create expected matrix

  output = matrix / mat_expected

  return(
    if(inherits(output, "dgCMatrix")) {output} else {methods::as(output, "dgCMatrix")}
    )
  }












