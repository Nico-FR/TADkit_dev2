#' @title MD plot between matrices
#'
#' @description Plot of log2(mat1 / mat2) between interactions values of 2 matrices according to bin distances. Note that only 1e6 interaction pairs are randomly selected to produce the plot.
#'
#' @param matrix.lst list of `dgCMatrix` or `matrix` object for only one chromosome.
#' @param bin.width Default is NULL, but if bin size is set, real distances between interactions is used on the x axis.
#' @param mean_centering Logical, default is FALSE. Whether or to subtracting the mean.
#' @return ggplot
#'
#' @importFrom Matrix triu summary
#' @importFrom dplyr filter mutate sample_n
#' @importFrom viridis scale_fill_viridis
#' @importFrom tidyr drop_na
#' @import ggplot2
#' @examples
#' #create matrices with rnadom number of interactions:
#' nrow = nrow(mat_HCT116_chr19_50kb)
#' mat2 = matrix(abs(rnorm(nrow*nrow)), nrow, nrow)
#'
#' #create list of matrices
#' matrix.lst = list(ind1 = mat_HCT116_chr19_50kb, ind2 = mat2)
#'
#' MDplot(matrix.lst, bin.width = 50e3, mean_centering = TRUE)
#'
#' @export
#'
MDplot <- function(matrix.lst, bin.width = NULL, mean_centering = FALSE) {

  i <- j <- x <- y <- NULL
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

  #melt matrices
  mat = matrix.lst[[1]] / matrix.lst[[2]]
  mat = if (!inherits(mat, "CsparseMatrix")) {
    as(Matrix::triu(mat), "CsparseMatrix")} else {
      Matrix::triu(mat)}

  #melted
  df = Matrix::summary(mat) %>%
    dplyr::mutate(distance = if(is.null(bin.width)) {abs(i - j)} else {abs(i - j) * bin.width}) %>%
    dplyr::mutate(y = log2(x)) %>%
    tidyr::drop_na(y) %>%
    dplyr::filter(!y == Inf)

  ########################################
  #plot
  ########################################
  ind1 = names(matrix.lst)[1]
  ind2 = names(matrix.lst)[2]

  if (isTRUE(mean_centering)) {
    m = mean(df$y)
    df$y = df$y - m
  }

  if(nrow(df) > 1e6) {df = dplyr::sample_n(df, 1e6)}

  p <- ggplot2::ggplot(df, aes(y = y, x = distance))+stat_bin_2d(bins = 50)+
    ggplot2::geom_hline(yintercept = mean(df$y), col = "red")+
    ggplot2::stat_smooth(se=FALSE, col="white", span = 0.5)+
    viridis::scale_fill_viridis(option = "D")+
    ggplot2::ylab(paste0("log2(", ind1," / ", ind2, ")"))

  if(!is.null(bin.width)) {
    p <- p+ggplot2::scale_x_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6))
  }

  return(p)

}
