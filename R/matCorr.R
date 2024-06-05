#' @title Correlation between matrices
#'
#' @description Computes the correlations between bins (e.g. between counts) of 2 or more matrices.
#'
#' @param matrix.lst list of `dgCMatrix` or `matrix` object for only one chromosome.
#' @param log2 logical. Use the log2 of the matrix values. Default is `TRUE`.
#' @param method a character string indicating which correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#' @param max.distance Maximum distance between bins. Default is NULL for no filtering. For example, if distance = 1e6 , bins with distances greater than 1Mb are not taken into account.
#' @param bin.width Default is NULL but mandatory if max.distance is set.
#' @param self_interaction logical. Whether or not to add interactions within bins (e.g. interactions within bin 1). If FALSE, interactions along the diagonal will be removed from the matrices before computation.
#' @param output default is "corr" to output correlation values. "data" to return values (e.g. count between each bin) and "plot" to return correlation graph of the first two matrices (only 2e5 points are randomly selected to produce plot).
#' @param start,stop Area in bp of the chromosome to be analyzed. Default is NULL to use the entire chromosome (i.e. entire matrix).
#' @return matrix with correlation values
#'
#' @importFrom Matrix triu summary diag tril
#' @importFrom dplyr full_join filter
#' @import ggplot2
#' @examples
#' matrix.lst = list(ind1 = mat_HCT116_chr19_50kb, ind2 = mat_HCT116_chr19_50kb)
#' matCorr(matrix.lst)
#'
#' @export
#'
matCorr <- function(matrix.lst, log2 = TRUE, output = "corr", self_interaction = FALSE, max.distance = NULL, bin.width = NULL, start = NULL, stop = NULL, method = "pearson") {

  . <- i <- j <- bin1 <- bin2 <- a <- b <- NULL
  ########################################"
  #Sanity check
  ########################################"
  if (!is.list(matrix.lst)) {
    stop("matrix.lst must be a list")
  }

  #stop if not a matrix
  if(!inherits(matrix.lst[[1]], c("Matrix", "matrix"))) {
    stop("input matrix is not a matrix or Matrix object")}

  #warning if matrices with diff size
  if(!lapply(matrix.lst, nrow) %>% unname %>% unlist() %>% unique %>% length() == 1) {
    warning("matrices do not have the same size!")
  }
  if(!lapply(matrix.lst, ncol) %>% unname %>% unlist() %>% unique %>% length() == 1) {
    warning("matrices do not have the same size!")
  }

  if (is.null(names(matrix.lst))) {names(matrix.lst) = paste0("mat", 1:length(matrix.lst))}

  #set to CsparseMatrix
  matrix.lst = base::lapply(matrix.lst, function(mat) {
    if (!inherits(mat, "CsparseMatrix")) {
      as(Matrix::triu(mat), "CsparseMatrix")} else {mat}
  })

  #melt matrices and filter diagonal
  melted.lst = base::lapply(matrix.lst, function(mat){
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

  #filter area (start and stop)
  if(!(is.null(start) | is.null(stop))) { #if start & stop != NULL

    if (is.null(bin.width)) {
      stop("bin.width value is needed!")
    }

    #bin indexes
    if ((start %/% bin.width == start / bin.width) & (stop %/% bin.width == stop / bin.width)) {
      start = start / bin.width
      stop = stop / bin.width} else {
        start = round(start / bin.width) + 1 # +1 because matrix coordinates start to 1. So first nucleotide (i.e. vp.start = 1bp) is bin nb 1.
        stop = ifelse((round(stop / bin.width) + 1) <= start, start, round(stop / bin.width))
        message(paste0("start/stop are not multiples of bin.width, round to ", (start - 1) * bin.width,
                       " and ", stop * bin.width, " (i.e. ", stop + 1 - start, " bins)."))
      }

    #filter area
    melted.lst = base::lapply(melted.lst, function(melted){
      melted %>% dplyr::filter(i <= stop, j >= start)
    })
  }

  #merge melted.lst
  merge = melted.lst %>% Reduce(function(...) dplyr::full_join(..., by = c("i", "j")), .) %>%
    `colnames<-`(c("bin1", "bin2", names(matrix.lst)))

  #log2
  if(isTRUE(log2)) {
    merge[,3:ncol(merge)] <- log2(merge[,3:ncol(merge)])
  }

  ########################################
  #output correlation matrix
  ########################################
  corr = stats::cor(merge %>% select(-bin1, -bin2), use = "pairwise.complete.obs", method = method)

  if (output == "data") {
    return(merge)
  }

  if (output == "corr") {
    return(corr)
  }

  if (output == "plot") {
    df = merge[,3:4]  %>% tidyr::drop_na()
    if(nrow(df) > 2e5) {df = dplyr::sample_n(df, 2e5)}
    names(df) = c("a", "b")
    p =  ggplot2::ggplot(data = df, aes(x=a, y=b))+
      ggplot2::geom_point(alpha=0.1)+
      ggplot2::ylab(names(matrix.lst)[2])+ggplot2::xlab(names(matrix.lst)[1])+
      ggplot2::annotate("text", x = min(df[,1], na.rm = TRUE), y =  max(df[,2], na.rm = TRUE),
                        label = corr[names(matrix.lst)[1], names(matrix.lst)[2]],
                        col = "red", vjust = "inward", hjust = "inward")+ggplot2::theme_minimal()
    return(p)
  }
}

