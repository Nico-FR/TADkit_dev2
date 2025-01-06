#' @title randomize domains order
#'
#' @description From an ordered set of TADs, this function return un random set of TADs.
#'
#' @details `TADshuffling` shuffle the TADs order for each chromosomes (i.e keep the TADs sizes).
#' The TADs input should not have gaps (region between two TADs) excepted for the extremities of chromosomes.
#'
#' @inheritParams TADplot
#'
#' @return `GRanges`
#'
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#'
#' @export
#'
#' @examples
#' #get domains from boundaries:
#' boundaries.gr = dataframes2grange(tad_HCT116_5kb.bed, human_chromsize)
#' domains.gr = boundary2domain(boundaries.gr)
#'
#' TADshuffling(domains.gr)
#'
TADshuffling <- function(tad.gr) {

  output = NULL

  for (c in as.character(levels(GenomeInfoDb::seqnames(tad.gr)))) {
    tad <- as.data.frame(tad.gr) %>% dplyr::filter(seqnames == c)
    random_rows <- sample(nrow(tad))
    random_tad <- tad[random_rows,]

    l <- length(random_tad[,1]) #number of TAD on the chr
    if (l == 0) {next} #if no TAD skeep this chr

    start.tmp = c(min(random_tad$start), cumsum(random_tad$width-1)[1:length(random_tad$width)-1]+min(random_tad$start))

    data <- data.frame(seqnames = rep(c, l),
                         start = start.tmp,
                         end = start.tmp + round(random_tad$width,-1),
                         width = round(random_tad$width,-1))

    if (length(random_tad)>5) {data <- cbind(data, random_tad[,6:length(random_tad)])}

    output = rbind(output, data)
  }

output.gr = TADkitdev2::dataframes2grange(output, data.frame(chr = as.character(levels(seqnames(tad.gr))),
                                                 size = as.numeric(seqlengths(tad.gr))),
                              metadata.mcols = if (length(output)>4) {5:length(output)} else {NULL})

  return(output.gr)
}



