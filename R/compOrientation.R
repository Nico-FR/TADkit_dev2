#' @title Orientation of PC1 score
#'
#' @description Considering that the A compartments have a stronger expression than the B compartments.
#' This function allows the score (i.e PC1 score) to be oriented (positif or negatif).
#' Then the compartments A (actif) and B (inactif) can be called accurately.
#'
#' @details This function is divided into 3 steps:
#' * Compartments calling from the score (see PC1calling function).
#' * Measures median gene expression in compartments A and B per chromosome.
#' * Reverse the score if median expression of A is smaller than B compartments.
#'
#' Return a list with 3 `dataframe`.
#' * bedgraph as `GRanges` object with the score oriented.
#' * `dataframe` with expression of each gene (before orientation).
#' * `dataframe` with median expression of each compartment (before orientation).
#'
#' @param bedgraph `GRanges` file with the score to be used for compartment calling (i.e PC1).
#' @param gene.gr `GRanges` file with gene annotations.
#' @param expression.data.frame `dataframe` with 3 columns:
#' * 1: gene IDs, must be identical IDs than `names(gene.gr)`,
#' * 2: common name for genes, could be the gene ID,
#' * 3: expression of each gene, such as raw count or `log(raw_count + 1)`.
#'
#' @return S3 class object with 3 `dataframes`.
#'
#' @import GenomeInfoDb
#' @import dplyr
#' @importFrom tidyr spread
#' @import GenomicRanges
#'
#' @export
#'
#' @examples
#' # output <- compOrientation(bedgraph, gene.gr, expression.data.frame)
#'

compOrientation <- function(bedgraph, gene.gr, expression.data.frame) {

  #sanity check
  if (is.null(names(gene.gr))) {
    stop("gene.gr object must have names (i.e gene IDs)")
  }
  if (length(gene.gr[duplicated(names(gene.gr))]) > 0) {
    stop("gene.gr object must have unique names (i.e gene IDs)")
  }
  if (length(expression.data.frame[duplicated(expression.data.frame[,1]),1]) > 0) {
    stop("expression.data.frame[,1] must have unique IDs (i.e gene IDs)")
  }
  if (length(expression.data.frame) < 3) {
    stop("expression.data.frame must have 3 columns")
  }
  if (length(expression.data.frame) > 3) {
    warning("expression.data.frame have more than 3 columns, the third column is used for expression values")
  }
  if (class(bedgraph)!="GRanges") {
    stop("bedgraph must be a GRange object")
  }
  if (is.na(mean(seqlengths(bedgraph), na.rm=T))) {
    stop("bedgraph must have seqlengths datas (see dataframes2grange function)")
  }
  if (length(S4Vectors::mcols(bedgraph)) != 1) {
    warning("bedgraph object should have 1 metadata column, additional columns are ignored!")
  }

  #compartment calling
  undirected_comp.gr = TADkit::PC1calling(bedgraph)

  # add chromosomes lengths
  seqlengths(undirected_comp.gr)[sort(names(seqlengths(undirected_comp.gr)))] = seqlengths(bedgraph)[sort(names(seqlengths(bedgraph)))]

  # add unique names to each domain
  names(undirected_comp.gr) = 1:length(undirected_comp.gr)

  # position of the genes in relation to each compartment
  geneTopo <- TADkit::geneTADtopo(tad.gr = undirected_comp.gr, gene.gr = gene.gr, expression.data.frame = expression.data.frame, ifoverlap = "remove")

  # merge gene, domain and expression datas
  undirected_comp.df = undirected_comp.gr %>% as.data.frame() %>% dplyr::select(seqnames, comp)
  undirected_comp.df$TAD_id = rownames(undirected_comp.df)
  rownames(undirected_comp.df) = NULL

  data = merge(geneTopo$geneTADtopo, undirected_comp.df)

  # add comp and chr to expression.df
  for (i in 1:length(data[,1])) {
    data$gene_exp.lst[[i]]$comp = data$comp[i]
    data$gene_exp.lst[[i]]$chr = data$seqnames[i]
    }

  # merge expression.df list
  data = do.call(rbind, data$gene_exp.lst)

  names(data) = c(names(data)[1:2], "exp", names(data)[4:5]) # add "comp" to the column with expression values

  # stats
  medianExp = data %>% dplyr::group_by(comp, chr) %>% dplyr::summarise(med_exp = median(exp)) %>% dplyr::arrange(chr) %>% tidyr::spread(., comp, med_exp) %>% as.data.frame()
  medianExp$to_be_inverted = ifelse(medianExp$A < medianExp$B, -1, 1)

  message(paste0(length(medianExp[medianExp$to_be_inverted == -1,1]), "/", length(medianExp[,1]), " chromosomes have been reversed"))

  # if exp A == exp B ==> warning
  if (length((medianExp %>% filter(A == B))[,1])) {
    message(paste0("median expression of chromosome(s) ", ((medianExp %>% filter(A == B))[,1]), " are identical between compartments A and B"))
  }

  #reverse score
  bedgraph_oriented = NULL
  for (c in unique(medianExp$chr)) {
    tmp = as.data.frame(bedgraph[GenomeInfoDb::seqnames(bedgraph) == c])
    tmp[,6] = tmp[,6] * medianExp$to_be_inverted[medianExp$chr == c]
    bedgraph_oriented = rbind(bedgraph_oriented, tmp)
  }
  bedgraph_oriented.gr = GenomicRanges::makeGRangesFromDataFrame(bedgraph_oriented, keep.extra.columns = T)

  # add chromosomes lengths...
  seqlengths(bedgraph_oriented.gr)[sort(names(seqlengths(bedgraph_oriented.gr)))] = seqlengths(bedgraph)[sort(names(seqlengths(bedgraph)))]

  output = list()
  output$bedgraph_oriented = bedgraph_oriented.gr
  output$expression = data
  output$medianExp = medianExp[,1:3]
  return(output)
}








