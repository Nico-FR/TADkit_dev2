#' expected interactions versus genetic distances
#'
#' @description This function calculates the average number of interactions (i.e. expected number of interaction) as a function of the distance between bins.
#'
#' @param matrix.lst List of `dgCMatrix` or `matrix` object for only one chromosome.
#' @param bin.width Bin width of the matrices in base pair.
#' @param colors.lst Set of 8 colors used for plot.
#' @return ggplot
#'
#' @importFrom scales unit_format
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' matExpDist(matrix.lst = list(HCT116 = mat_HCT116_chr19_50kb), bin.width = 50e3)
#'
#'
matExpDist <- function(matrix.lst, bin.width,
                       colors.lst = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")) {

  #sanity check
  if (!is.list(matrix.lst)) {
    stop("matrix.lst must be a list")
  }

  message(paste0("Analysis of ", length(matrix.lst), " matrix/matrices..."))

  output = lapply(matrix.lst, function(m) {
    sapply(1:(ncol(m) - 1), function (x){
      mean(diag(matrix(m, nrow = nrow(m), ncol = ncol(m))[,x:ncol(m)]), na.rm = T)}
      )
    })

  df = do.call(data.frame, output)
  df$distance = 1:(ncol(matrix.lst[[1]]) - 1) * bin.width / 2

  if (is.null(names(matrix.lst))) {names(df) = c(1:length(matrix.lst), "distance")}

  p = ggplot2::ggplot(data = tidyr::gather(df, "matrix", "expected", 1:length(matrix.lst)),
                      aes(y = expected, x = distance, color = matrix))+
    ggplot2::geom_line()+scale_color_manual(values = colors.lst)+
    scale_x_continuous(labels = scales::unit_format(unit = "Mb", scale = 1e-6), trans = "log10")+
    scale_y_continuous(trans = "log10")

  return(p)
}






































