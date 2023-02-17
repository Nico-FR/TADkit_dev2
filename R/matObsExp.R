#' @title observed / expected matrix
#'
#' @description For each bin of the matrix (interaction count observed) this function return a ratio: observed / expected.
#' The input file path must be in a matrix format (i.e array with as many rows and columns than the number of bin). Only the upper part of the matrix is used.
#' The output can be plot with MATplot function with log2=T and scale.colors = "ObsExp" parameters.
#'
#' @details The expected number of interaction corresponds to the average of the interaction counts according to the bin distances.
#'
#' @param matrix R object (data frame or matrix) or the matrix file path (in full matrix format) for only one chromosome. The path can be gzip file (ending by .gz)
#' @param matrix.colname logical. Does your matrix file path have column names (i.e header)? Default = TRUE
#' @param matrix.rowname logical. Does your matrix file path have row names? Default = TRUE
#' @param matrix.sep the field separator character of the matrix file path. Default = '\\t' (i.e tabulation).
#' @return matrix array
#' @export
#'

matObsExp <- function(matrix, matrix.colname = T, matrix.rowname = T, matrix.sep = "\t") {

  if(isTRUE(matrix.rowname)) matrix.col.skip <- 1 else matrix.col.skip <- NULL

  #read matrix path in data.frame
  if (is.character(matrix))  {
    if (substr(matrix, nchar(matrix) - 2, nchar(matrix)) == ".gz") {
      df = read.table(gzfile(matrix), sep = matrix.sep, header = matrix.colname, row.names = matrix.col.skip)
    } else {
      df = read.table(matrix, sep = matrix.sep, header = matrix.colname, row.names = matrix.col.skip)
    }
  }

  #read data frame
  if (is.data.frame(matrix))  {df = matrix}

  #read matrix or create one
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












