#' HiC matrix of Human HCT116 cell.
#' @format ## `mat_HCT116_chr19_50kb`
#' A dgCMatrix object with 1173 rows and columns containing the raw counts (number of interactions).
#' For chromosome 19 from Human HCT116 cells with a bin width of 50kb.
#' Downloaded from the 4DN portal with accession number 4DNESNSTBMBY.
"mat_HCT116_chr19_50kb"
#'
#'TADs of Human HCT116 cell.
#' @format ## `tad_HCT116_5kb.bed`
#' A data frame containing TAD bounaries.
#' For autosomes from Human HCT116 cells with a bin width of 50kb.
#' Downloaded from the 4DN portal with accession number 4DNESNSTBMBY.
#' \describe{
#'   \item{chr}{chromosome name}
#'   \item{start}{boundary start}
#'   \item{end}{boundary end}
#' }
"tad_HCT116_5kb.bed"
#'
#' Size of Human autosomes.
#' @format ## `human_chromsize`
#' A data frame containing the size in base pair of human autosomes
#' \describe{
#'   \item{chr}{chromosome name}
#'   \item{size}{chromosome size in base pair}
#' }
"human_chromsize"
#'
#' Insulation score of Human chromosome 19.
#' @format ## `IS_HCT116_chr19_5kb.bedgraph`
#' A data frame containing insulation score.
#' For chromosome 19 from Human HCT116 cells with a bin width of 5kb.
#' Downloaded from the 4DN portal with accession number 4DNESNSTBMBY.
#' \describe{
#'   \item{chr}{chromosome name}
#'   \item{start}{bin start}
#'   \item{end}{bin end}
#'   \item{IS}{insulation score}
#' }
"IS_HCT116_chr19_5kb.bedgraph"
#'
#' RNA seq coverage for Human chr 19 download on Ensembl
#' "http://ftp.ensembl.org/pub/release-104/bamcov/homo_sapiens/genebuild/GRCh38.illumina.colon.1.bam.bw.
#' @format ## `rna_seq_chr19_10.1to10.6mb.gr`
#' Grange with RNA coverage.
"rna_seq_chr19_10.1to10.6mb.gr"
#'
#' First principal of Human chromosomes.
#' @format ## `PC1_1_50kb.bedgraph`
#' A GRanges containing PC1 score from Human HCT116 cells with a bin width of 250kb.
#' Downloaded from the 4DN portal with accession number 4DNESNSTBMBY.
"PC1_250kb.gr"
