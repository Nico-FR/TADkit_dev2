#' Matrix of bovine 1 and chromosome 25.
#' @format ## `matrix_1_chr25_50kb`
#' A dgCMatrix object with 847 rows and columns containing the raw counts (number of interactions).
#' For bovine 1, chromosome 25 with a bin width of 50kb.
#' Only the upper matrix and non zero count are stored.
"matrix_1_chr25_50kb"
#'
#' Matrix of bovine 2 and chromosome 25.
#' @format ## `matrix_2_chr25_50kb`
#' A dgCMatrix object with 847 rows and columns containing the raw counts (number of interactions).
#' For bovine 2, chromosome 25 with a bin width of 50kb.
#' Only the upper matrix and non zero count are stored.
"matrix_2_chr25_50kb"
#'
#'TAD of bovine 1.
#' @format ## `tad_1_10kb.bed`
#' A data frame containing TADs (domains) predicted by HiCExplorer on 10kb matrix.
#' \describe{
#'   \item{chr}{chromosome name}
#'   \item{start}{TAD start}
#'   \item{end}{TAD end}
#' }
"tad_1_10kb.bed"
#'
#'TAD of bovine 2.
#' @format ## `tad_2_10kb.bed`
#' A data frame containing TADs (domains) predicted by HiCExplorer on 10kb matrix.
#' \describe{
#'   \item{chr}{chromosome name}
#'   \item{start}{TAD start}
#'   \item{end}{TAD end}
#' }
"tad_2_10kb.bed"
#'
#' Insulation score of bovine 1 and chromosome 25.
#' @format ## `IS_1_10kb.bedgraph`
#' A data frame containing insulation score calculated by HiCExplorer on 10kb matrix.
#' \describe{
#'   \item{chr}{chromosome name}
#'   \item{start}{bin start}
#'   \item{end}{bin end}
#'   \item{IS}{insulation score}
#' }
"IS_1_10kb.bedgraph"
#'
#' Insulation score of bovine 2 and chromosome 25.
#' @format ## `IS_2_10kb.bedgraph`
#' A data frame containing insulation score calculated by HiCExplorer on 10kb matrix.
#' \describe{
#'   \item{chr}{chromosome name}
#'   \item{start}{bin start}
#'   \item{end}{bin end}
#'   \item{IS}{insulation score}
#' }
"IS_2_10kb.bedgraph"
#'
#' First principal component score of bovine 1.
#' @format ## `PC1_1_50kb.bedgraph`
#' A data frame containing PC1 score calculated by dcHiC on 50kb matrix.
#' \describe{
#'   \item{chr}{chromosome name}
#'   \item{start}{bin start}
#'   \item{end}{bin end}
#'   \item{PC1}{principal component score}
#' }
"PC1_1_50kb.bedgraph"
#'
#' First principal component score of bovine 2.
#' @format ## `PC1_2_50kb.bedgraph`
#' A data frame containing PC1 score calculated by dcHiC on 50kb matrix.
#' \describe{
#'   \item{chr}{chromosome name}
#'   \item{start}{bin start}
#'   \item{end}{bin end}
#'   \item{PC1}{principal component score}
#' }
"PC1_2_50kb.bedgraph"
#'
#' Size of cow autosomes.
#' @format ## `chromsize`
#' A data frame containing the size in base pair of bovine chromosomes
#' \describe{
#'   \item{chr}{chromosome name}
#'   \item{size}{sine in base pair}
#' }
"chromsize"
#'
#' RNA seq datas for chr 25 download on Ensembl
#' "http://ftp.ensembl.org/pub/release-104/bamcov/bos_taurus/genebuild/ARS-UCD1.2.ENA.heart.1.bam.bw).
#' @format ## `rna_seq_chr25_13to16mb.bedgraph`
#' Grange with mean coverage.
"rna_seq_chr25_13to16mb.bedgraph"
