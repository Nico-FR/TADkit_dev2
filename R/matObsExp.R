#' @title observed / expected matrix
#'
#' @description For each bin of the matrix (interaction coutn observed) this function return a ratio: observed / expected.
#' The input file path must be in a matrix format (ie. array with as many rows and columns than the number of bin). Only the upper part of the matrix is used.
#' The output can be plot with MATplot function with log2=T
#'
#' @details The expected interaction corresponds to the average of the interactions according to the bin distances.
#'
#' @param matrix R object (data frame or matrix) or the matrix file path (in matrix format) for only one chromosome. The path can be gzip file (ending by .gz)
#' @param matrix.colname logical. Does your matrix file (ie path) have column names (ie header)? Default = TRUE
#' @param matrix.rowname logical. Does your matrix file (ie path) have row names? Default = TRUE
#' @param matrix.sep the field separator character. Values on each line of the matrix file are separated by this character. Default = '\\t' (ie tabulation).
#' @return matrix array
#' @export
#'

matObsExp <- function(matrix, matrix.colname = T, matrix.rowname = T, matrix.sep = "\t") {

  #read matrix
  #read matrix.path
  if(isTRUE(matrix.rowname)) matrix.col.skip <- 1 else matrix.col.skip <- NULL

  if (is.character(matrix) & rev(strsplit(matrix, split = "\\.")[[1]])[1] == "gz")  {
    df = read.table(gzfile(matrix), sep = matrix.sep, header = matrix.colname, row.names = matrix.col.skip)
  }

  if (is.character(matrix) & !rev(strsplit(matrix, split = "\\.")[[1]])[1] == "gz")  {
    df = read.table(matrix, sep = matrix.sep, header = matrix.colname, row.names = matrix.col.skip)
  }

  #read data frame matrix
  if (is.data.frame(matrix))  {df = matrix}

  #create matrix
  if(!is.matrix(matrix)) {mat = matrix(as.matrix(df), nrow = length(df))} else {mat = matrix}

  #mean diag
  mean_diag = sapply(1:(ncol(mat) - 1),
                     function (x){mean(diag(mat[,x:ncol(mat)]), na.rm = T)})

  mat_obsexp = matrix(ncol = ncol(mat), nrow = ncol(mat)) #create new matrix

  #obs/exp for each diag
  for (x in (1:(ncol(mat)-1))) {
    diag(mat_obsexp[,x:ncol(mat_obsexp)]) <- diag(mat[,x:ncol(mat)])/mean_diag[x]}

  #add lower tri
  mat_obsexp[lower.tri(mat_obsexp)] = t(mat_obsexp)[lower.tri(mat_obsexp)]

  return(mat_obsexp)
  }












