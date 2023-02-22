#' @title observed / expected matrix
#'
#' @description For each bin of the matrix (interaction count observed) this function return a ratio: observed / expected.
#' The input matrix must be a dgCMatrix or matrix object (i.e array with as many rows and columns than the number of bin). Only the upper part of the matrix is used.
#' The output can be plot with MATplot function with log2=T and scale.colors = "ObsExp" parameters.
#'
#' @details The expected number of interaction corresponds to the average of the interaction counts according to the bin distances (NB: bins with zero counts are not considered).
#'
#' @param matrix a dgCMatrix or matrix object
#' @return a dgCMatrix object: upper triangular and sparse Matrix
#' @importFrom stats toeplitz
#' @import Matrix
#' @importFrom methods as
#' @export
#'

matObsExp <- function(matrix) {

  if(isFALSE(class(matrix)[1] %in% c("dgCMatrix", "matrix"))) {
    stop("input matrix is not a matrix or dgCMatrix object")}

  #using dgcmatrix object
  #mat = Matrix::triu(matrix)
  #mat[Matrix::triu(mat == 0)] <- NA

  #using matrix
  mat =  matrix(matrix, nrow = nrow(matrix), ncol = ncol(matrix))
  mat[mat == 0] <- NA

  #mean diag
  mean_diag = sapply(1:(ncol(mat) - 1),
                     function (x){mean(diag(mat[,x:ncol(mat)]), na.rm = T)})

  mat_expected = stats::toeplitz(c(mean_diag, mat[1,ncol(mat)])) #create expected matrix

  output = matrix / mat_expected

  return(
    if(class(output) == "matrix") {as(output, "dgCMatrix")} else {output}
    )
  }












