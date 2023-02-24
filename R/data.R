#' Matrix of individual 1 and chromosome 25.
#' @format ## `matrix_1_chr25_50kb`
#' A dgCMatrix object with 847 rows and columns containing the raw counts (number of interactions).
#' For individual 1, chromosome 25 with a bin width of 50e3 bases pairs.
#' Only the upper matrix and non zero count are stored.
"matrix_1_chr25_50kb"
#'
#' Matrix of individual 1 and chromosome 26.
#' @format ## `matrix_1_chr26_50kb`
#' A dgCMatrix object with 1040 rows and columns containing the raw counts (number of interactions).
#' For individual 1, chromosome 26 with a bin width of 50e3 bases pairs.
#' Only the upper matrix and non zero count are stored.
"matrix_1_chr26_50kb"
#'
#' Matrix of individual 2 and chromosome 25.
#' @format ## `matrix_2_chr25_50kb`
#' A dgCMatrix object with 847 rows and columns containing the raw counts (number of interactions).
#' For individual 2, chromosome 25 with a bin width of 50e3 bases pairs.
#' Only the upper matrix and non zero count are stored.
"matrix_2_chr25_50kb"
#'
#' Matrix of individual 2 and chromosome 26.
#' @format ## `matrix_2_chr26_50kb`
#' A dgCMatrix object with 1040 rows and columns containing the raw counts (number of interactions).
#' For individual 2, chromosome 26 with a bin width of 50e3 bases pairs.
#' Only the upper matrix and non zero count are stored.
"matrix_2_chr26_50kb"
#'
#' Topological Associated Domain of individual 1 and chromosome 25/26.
#' @format ## `tad_1_10kb.bed`
#' A data frame containing TADs (domains) predicted by HiCExplorer on 10kb matrix.
#' For individual 1, chromosome 25 and 26.
#' \describe{
#'   \item{chr}{chromosome name}
#'   \item{start}{TAD start}
#'   \item{end}{TAD end}
#' }
"tad_1_10kb.bed"
#'
#' Topological Associated Domain of individual 2 and chromosome 25/26.
#' @format ## `tad_2_10kb.bed`
#' A data frame containing TADs (domains) predicted by HiCExplorer on 10kb matrix.
#' For individual 2, chromosome 25 and 26.
#' \describe{
#'   \item{chr}{chromosome name}
#'   \item{start}{TAD start}
#'   \item{end}{TAD end}
#' }
"tad_2_10kb.bed"
#'
#' Insulation score of individual 1 and chromosome 25/26.
#' @format ## `IS_1_10kb.bedgraph`
#' A data frame containing insulation score calculated by HiCExplorer on 10kb matrix.
#' For individual 1, chromosome 25 and 26.
#' \describe{
#'   \item{chr}{chromosome name}
#'   \item{start}{bin start}
#'   \item{end}{bin end}
#'   \item{IS}{insulation score}
#' }
"IS_1_10kb.bedgraph"
#'
#' Insulation score of individual 2 and chromosome 25/26.
#' @format ## `IS_2_10kb.bedgraph`
#' A data frame containing insulation score calculated by HiCExplorer on 10kb matrix.
#' For individual 2, chromosome 25 and 26.
#' \describe{
#'   \item{chr}{chromosome name}
#'   \item{start}{bin start}
#'   \item{end}{bin end}
#'   \item{IS}{insulation score}
#' }
"IS_2_10kb.bedgraph"
#'
#' First principal component score of individual 1 and chromosome 25/26.
#' @format ## `PC1_1_50kb.bedgraph`
#' A data frame containing PC1 score calculated by dcHiC on 50kb matrix.
#' For individual 1, chromosome 25 and 26.
#' \describe{
#'   \item{chr}{chromosome name}
#'   \item{start}{bin start}
#'   \item{end}{bin end}
#'   \item{PC1}{principal component score}
#' }
"PC1_1_50kb.bedgraph"
#'
#' First principal component score of individual 2 and chromosome 25/26.
#' @format ## `PC1_2_50kb.bedgraph`
#' A data frame containing PC1 score calculated by dcHiC on 50kb matrix.
#' For individual 2, chromosome 25 and 26.
#' \describe{
#'   \item{chr}{chromosome name}
#'   \item{start}{bin start}
#'   \item{end}{bin end}
#'   \item{PC1}{principal component score}
#' }
"PC1_2_50kb.bedgraph"
