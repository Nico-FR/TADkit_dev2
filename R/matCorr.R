#' @title Correlation between matrices
#'
#' @description Computes the correlations between bins (e.g. between counts) of 2 or more matrices.
#'
#'
#' @param matrice.lst list of `dgCMatrix` or `matrix` object for only one chromosome.
#' @param log2 logical. Use the log2 of the matrix values. Default is `TRUE`.
#' @param method a character string indicating which correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#' @param max.distance Maximum distance between bins. Default is NULL for no filtering. For example, if distance = 1e6 , bins with distances greater than 1Mb are not taken into account.
#' @param bin.width Default is NULL but mandatory if max.distance is set.
#' @param self_interaction logical. Whether or not to add interactions within bins (e.g. interactions within bin 1).
#' @param output default is "corr" to output correlation values. "data" to return values (e.g. count between each bin) and "plot" to return correlation graph of the first two matrices (only 2e5 points are randomly selected to produce plot).
#'
#' @return matrix with correlation values
#'
#' @importFrom Matrix triu summary diag tril
#' @importFrom dplyr full_join filter
#' @import ggplot2
#' @examples
#' matrice.lst = list(ind1 = mat_HCT116_chr19_50kb, ind2 = mat_HCT116_chr19_50kb)
#' matCorr(matrice.lst)
#'
#' @export
#'
matCorr <- function(matrice.lst, log2 = TRUE, output = "corr", self_interaction = FALSE, max.distance = NULL, bin.width = NULL, method = "pearson") {


  . <- NULL
  ########################################"
  #Sanity check
  ########################################"
  if (!is.list(matrice.lst)) {
    stop("matrice.lst must be a list")
  }

  #stop if not a matrix
  if(!inherits(matrice.lst[[1]], c("Matrix", "matrix"))) {
    stop("input matrix is not a matrix or Matrix object")}

  #warning if matrices with diff size
  if(!lapply(matrice.lst, nrow) %>% unname %>% do.call(identical, .)) {
    warning("matrices do not have the same size!")
  }

  #set to CsparseMatrix
  matrice.lst = base::lapply(matrice.lst, function(mat) {
    if (!inherits(mat, "CsparseMatrix")) {
      as(Matrix::triu(mat), "CsparseMatrix")} else {mat}
  })

  #melt matrices and filter diagonal
  melted.lst = base::lapply(matrice.lst, function(mat){
    if(isFALSE(self_interaction)) {
      Matrix::summary(Matrix::triu(mat, 1))
    } else {Matrix::summary(mat)}
  })

  #filter max.distance
  if(!is.null(max.distance)) {
    if (is.null(bin.width)) {
      stop("bin.width value is needed!")
    }
    melted.lst = base::lapply(melted.lst, function(melted){
      melted %>% dplyr::filter(abs(i-j) <= max.distance %/% bin.width)
    })
  }

  #merge melted.lst
  merge = melted.lst %>% Reduce(function(...) dplyr::full_join(..., by = c("i", "j")), .) %>%
    `colnames<-`(c("bin1", "bin2", names(matrice.lst)))

  #log2
  if(isTRUE(log2)) {
    merge[,3:ncol(merge)] <- log2(merge[,3:ncol(merge)])
  }

  ########################################
  #output correlation matrix
  ########################################
  corr = stats::cor(merge %>% select(-bin1, -bin2), use = "complete.obs", method = method)

  if (output == "data") {
    return(merge)
  }

  if (output == "corr") {
    return(corr)
  }

  if (output == "plot") {

    for (i in 3:ncol(merge)) {
      j = i + 1
      if (j > ncol(merge)) {next}
      df = merge[,c(i,j)]  %>% tidyr::drop_na()
      if(nrow(df) > 2e5) {df = dplyr::sample_n(df, 2e5)}
      names(df) = c("a", "b")
        p =  ggplot2::ggplot(data = df, aes(x=a, y=b))+
          ggplot2::geom_point(alpha=0.1)+
          ggplot2::ylab(names(matrice.lst)[j - 2])+ggplot2::xlab(names(matrice.lst)[i - 2])+
          ggplot2::annotate("text", x = min(df[,1], na.rm = TRUE), y =  max(df[,2], na.rm = TRUE),
                   label = corr[names(matrice.lst)[i - 2], names(matrice.lst)[j - 2]],
                   col = "red", vjust = "inward", hjust = "inward")+ggplot2::theme_minimal()
        return(p)
    }
  }
}
