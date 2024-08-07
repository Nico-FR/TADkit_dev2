% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/viewPointInteract.R
\name{viewPointInteract}
\alias{viewPointInteract}
\title{plot view point interaction}
\usage{
viewPointInteract(
  matrix.lst,
  bin.width,
  vp.start,
  vp.stop,
  output = "GRanges",
  start = NULL,
  stop = NULL,
  self_interaction = FALSE,
  Qnorm = FALSE,
  log2 = FALSE,
  seqname = "1",
  colors.lst = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
    "#E5C494", "#B3B3B3")
)
}
\arguments{
\item{matrix.lst}{List of \code{dgCMatrix} or \code{matrix} object for only one chromosome.}

\item{bin.width}{Bin width of the matrix in base pair.}

\item{vp.start}{Start of the view point in base pair.}

\item{vp.stop}{Stop/end of the view point in base pair.}

\item{output}{"GRanges" (default) to return a bedgraph with the average interaction of the view point with the area analyzed and "plot" to return a graph.}

\item{start, stop}{Area in bp of the chromosome to be analyzed. Default is NULL to use the entire chromosome (i.e. entire matrix).}

\item{self_interaction}{logical. Whether or not to add interactions within the view-point.}

\item{Qnorm}{Quantile normalization of the average interaction. Default is FALSE.
This parameter can be useful for the graph, as it allows a better comparison of distributions between matrices. Note that normalization is best performed on the entire chromosome than anly on an area.}

\item{log2}{logical. Use the log2 of the matrix values. Default is \code{FALSE}. Note that if TRUE, interaction with 0 count are removed from the analysis.}

\item{seqname}{chromosome names as character, default = "1".}

\item{colors.lst}{Set of 8 colors used for plot.}
}
\value{
\code{GRanges} bedgraph with average interaction.
}
\description{
This function (also known as virtual4C) calculates the average number of interactions of a region of interest (e.g. a TAD) with respect to the rest of the matrix.
}
\examples{
# get domains from boundaries:
boundaries.gr = dataframes2grange(tad_HCT116_5kb.bed, human_chromsize)
domains.gr = boundary2domain(boundaries.gr)

MATplot(mat_HCT116_chr19_50kb,
    start = 7e6, stop = 15e6,
    bin.width = 50e3, log2 = TRUE,
    tad.upper.tri = domains.gr["chr19_11597500",])

#plot interactions of the domain above, between 7Mb and 15Mb
viewPointInteract(matrix.lst = list(HCT116 = mat_HCT116_chr19_50kb), bin.width = 50e3,
  vp.start = 11597500, vp.stop = 12672500,
  start = 7e6, stop = 15e6,
  self_interaction = FALSE, #FALSE to not plot intra-TAD interactions
  output = "plot")


}
