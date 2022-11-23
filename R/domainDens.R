#' @title Distribution of annotation density within domains
#'
#' @description Distribution of annotations density (e.g. genes) within domains (e.g. TADs, compartments...).
#' This function calculates the density per bin of annotations and measures the relative position of each bin within each domain
#'
#' @details As an exemple, this function take all TAD domains and calculate the density of annotations .
#' For example:
#' 1- calculates the density of genes per bin of 10kb,
#' 2- calculate the relative positions of each bin within his corresponding TAD,
#' 3- plot the smoothed bin density (and smoothed zscore bin density) according to relative positions.
#'
#'
#' @param domain.gr GRange file with domains (TADs, compartments...).
#' @param annot.gr GRange file with genomic annotations (genes, repeat elements...).
#' @param annot.col Column number (metadata columns) used to separate annotations. By default, annot.col = NULL. Use "strand" in order to separate annotations according to their strands.
#' @param domain.col Column number (metadata columns) used to separate domains (e.g. column to differentiate between compartments A and B). By default, domain.col = NULL.
#' @param bin.width Size of the bin in bp. This should match the size of the bins used to determine the domains.
#' @param output Default is "plot" to return a ggplot. Use "data" to return the datas used to produce the plot.
#'
#'
#' @return Return a ggplot graph
#' @import GenomeInfoDb
#' @importFrom BiocGenerics strand end start width
#' @importFrom plyr ddply
#' @import ggplot2
#' @import scales
#' @import GenomicRanges
#' @export
#'
#' @examples
#'
domainDens <- function(domain.gr, annot.gr, domain.col = NULL, annot.col = NULL, bin.width = 50e3, output = "plot") {

  domain.gr = domain.gr; annot.gr = gene.gr; domain.col = 1; annot.col = "strand"; bin.width = 50e3

  #annot.col parameter
  nb_metadatacolumns = length(GenomicRanges::mcols(annot.gr))
  if (is.null(annot.col)) {
    GenomicRanges::mcols(annot.gr)[,nb_metadatacolumns + 1] = "all"
    names(GenomicRanges::mcols(annot.gr))[nb_metadatacolumns + 1] = "type"
    annot.col = nb_metadatacolumns + 1
  }

  if (annot.col == "strand") {
    GenomicRanges::mcols(annot.gr)[,nb_metadatacolumns + 1] = GenomicRanges::strand(annot.gr)
    names(GenomicRanges::mcols(annot.gr))[nb_metadatacolumns + 1] = "strand"
    annot.col = nb_metadatacolumns + 1
  }

  if (length(GenomicRanges::mcols(annot.gr)) < annot.col) {
    stop(paste0("Wrong annot.col number. There is only ", length(GenomicRanges::mcols(annot.gr)), " column(s) with metadata."))
  }

  #genome binning
  if (is.na(mean(GenomeInfoDb::seqlengths(annot.gr)))) {
    stop("There is no seqlenths in domain.gr, see dataframes2grange function")
  }
  bin.gr = GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(domain.gr), tilewidth = bin.width, cut.last.tile.in.chrom = TRUE) #50k bin

  # add density for each annot.type (annot.col)
  #GenomicRanges::mcols(annot.gr)[, annot.col] <- gsub("-", "M", fixed=T, GenomicRanges::mcols(annot.gr)[, annot.col]) #replace dashes
  #GenomicRanges::mcols(annot.gr)[, annot.col] <- gsub("+", "P", fixed=T, GenomicRanges::mcols(annot.gr)[, annot.col]) #replace dashes
  annot.type = GenomicRanges::mcols(annot.gr)[, annot.col] %>% unique
  for (c in annot.type) {
    bin.gr = GenomicRanges::binnedAverage(bins = bin.gr,
                                          numvar = GenomicRanges::coverage(annot.gr[GenomicRanges::mcols(annot.gr)[, annot.col] == c]),
                                          na.rm = FALSE, varname = c)
  }

  bin.gr <- GenomicRanges::resize(bin.gr, 1, fix = "center") #keep middle of each been

  #add domain names (ie row number) overlapped for each bin
  bin.gr$domainHit <- GenomicRanges::findOverlaps(bin.gr, domain.gr, select = "arbitrary") # arbitrary to randomly select the domain when the center of a bin have the same position than a boundary (this should not happen)

  bin.gr = bin.gr[!is.na(bin.gr$domainHit)] #remove bin that do not overlap a domain

  # add domains datas needed (start, end, domain.col)
  bin.gr$domain_start = GenomicRanges::start(domain.gr[bin.gr$domainHit])
  bin.gr$domain_end = GenomicRanges::end(domain.gr[bin.gr$domainHit])
  if (is.numeric(domain.col)) {
    bin.gr$domain = GenomicRanges::mcols(domain.gr)[bin.gr$domainHit, domain.col]
  }


  #add relative position
  bin.gr$relative_position <- (GenomicRanges::start(bin.gr) - bin.gr$domain_start) / (bin.gr$domain_end - bin.gr$domain_start)

  #add zscore
  for (c in annot.type) {
    GenomicRanges::mcols(bin.gr)[, paste0("zscore_", c)] = scale(GenomicRanges::mcols(bin.gr)[, c])
  }


  if (output == "data") {
    return(bin.gr)
  }

  if (output == "plot") {

    if (is.null(domain.col)) {
      data = (bin.gr %>% as.data.frame())[, c(paste0("dens_", annot.type), "relative_position")] %>% tidyr::gather(annot.col, "density", 1:length(annot.type))
      #names(data) = c("relative_position", names(GenomicRanges::mcols(annot.gr[,annot.col])), "density")
      p1 = ggplot2::ggplot(data, ggplot2::aes(y = density, x = relative_position, color = annot.col))+
        ggplot2::geom_smooth()+
        ggplot2::labs(color = names(GenomicRanges::mcols(annot.gr[,annot.col])))
      print(p1)
    }

    if (is.numeric(domain.col)) {
      #data <- (as.data.frame(bin.gr))[,c(paste0("dens_", annot.type), "relative_position","domain")] %>% tidyr::gather(annot.col, "density", 1:length(annot.type))
      data <- as.data.frame(bin.gr)
      names(data) = c("seqnames", "start", "end", "width", "strand", names(mcols(bin.gr)))
      data2 = data %>% select(as.character(annot.type), "relative_position","domain") %>% tidyr::gather(annot.col, "density", 1:length(annot.type))
      p1 = ggplot2::ggplot(data2, ggplot2::aes(y = density, x = relative_position, color = domain))+
        ggplot2::geom_smooth()+ggplot2::facet_wrap(.~annot.col, scales = "free_y")+
        ggplot2::labs(color = names(GenomicRanges::mcols(domain.gr[, domain.col])))
      print(p1)

      for (c in annot.type) {
        data3 <-  data2 %>% filter(annot.col == c)
        p2 = ggplot2::ggplot(data3, ggplot2::aes(y = density, x = relative_position))+
          ggplot2::geom_smooth()+ggplot2::facet_wrap(.~domain)+ggplot2::ggtitle(c)
        print(p2)
      }

      data4 <- data %>% select(paste0("zscore_", annot.type), "relative_position","domain")
      names(data4) <- c(as.character(annot.type), names(data4)[(length(annot.type)+1):length(data4)])
      data5 <- data4 %>% tidyr::gather(annot.col, "zscore", 1:length(annot.type))
      p3 = ggplot2::ggplot(data5, ggplot2::aes(y = zscore, x = relative_position, color = annot.col))+
        ggplot2::geom_smooth()+ggplot2::facet_wrap(.~domain)+
      ggplot2::labs(color = names(GenomicRanges::mcols(annot.gr[,annot.col])))
      print(p3)
    }
  }
}


























