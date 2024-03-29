% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matCorr.R
\name{matCorr}
\alias{matCorr}
\title{Correlation between matrices}
\usage{
matCorr(
  matrice.lst,
  log2 = TRUE,
  output = "corr",
  self_interaction = FALSE,
  max.distance = NULL,
  bin.width = NULL,
  method = "pearson"
)
}
\arguments{
\item{matrice.lst}{list of \code{dgCMatrix} or \code{matrix} object for only one chromosome.}

\item{log2}{logical. Use the log2 of the matrix values. Default is \code{TRUE}.}

\item{output}{default is "corr" to output correlation values. "data" to return values (e.g. count between each bin) and "plot" to return correlation graph of the first two matrices (only 2e5 points are randomly selected to produce plot).}

\item{self_interaction}{logical. Whether or not to add interactions within bins (e.g. interactions within bin 1).}

\item{max.distance}{Maximum distance between bins. Default is NULL for no filtering. For example, if distance = 1e6 , bins with distances greater than 1Mb are not taken into account.}

\item{bin.width}{Default is NULL but mandatory if max.distance is set.}

\item{method}{a character string indicating which correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.}
}
\value{
matrix with correlation values
}
\description{
Computes the correlations between bins (e.g. between counts) of 2 or more matrices.
}
\examples{
matrice.lst = list(ind1 = mat_HCT116_chr19_50kb, ind2 = mat_HCT116_chr19_50kb)
matCorr(matrice.lst)

}
