#' @title Genes topology according to TAD
#'
#' @description For each domain (i.e TAD) the function return a data frame with the numbers, names and strands of genes. It also return some statistics which can be plotted using plot(output).
#'
#' @details Return a S3 object with 3 dataframes. The second and the third dataframe can be plotted using the plot function
#'
#'
#' @param tad.gr GRange file with TADs.
#' @param gene.gr GRange file with genomic annotations.
#' @param ifoverlap In case of annotation overlap a TAD boundary, few options are available:
#'    (1)"remove" to remove all annotations that overlaps a TAD boundary,
#'    (2)"real" to select the TAD in which the start of annotations are located,
#'    (3)"best" to select the TAD in which annotations have the best overlay.
#' @param expression.data.frame Data frame with expression data (raw counts...). The first two columns should be used to identify the genes. The 1st one must use the same gene id than names(gene.gr), ie unique ID. The 2nd can be used with usual gene names. Others columns give expression count, one column per sample/experiment.
#'
#' @return Return a S3 class object with 3 data frames.
#' @import GenomeInfoDb
#' @importFrom BiocGenerics strand end start width
#' @import GenomicRanges
#' @import dplyr
#' @export
#' @examples
#' # output <- geneTADtopo(tad.gr, gene.gr)
#' # plot(output)
#'
geneTADtopo <- function(tad.gr, gene.gr, ifoverlap = "best", expression.data.frame = NULL) {

  if (length(tad.gr) != length(unique(names(tad.gr)))) {
    stop("domain names must be unique!")
  }

  data1 <- TADhist(tad.gr, gene.gr, output = "data", annot.border = "start", ifoverlap = ifoverlap)

  # annot.gr ID [ie names(annot.gr)] must be unique
  if (length(gene.gr[duplicated(names(gene.gr))]) > 0) {
    stop("gene.gr must have unique names!")
  }

  # Count number of gene of each TAD
  nbgene_TAD = data1 %>% as.data.frame() %>% dplyr::count(nameHit, name="nb_genes")
  names(nbgene_TAD)=c("TAD_id", "nb_genes")

  # for each TAD get genes strands
  nbgene_TAD$gene_strands = sapply(nbgene_TAD[, 1],
                                   function(x){
                                     BiocGenerics::strand(data1[data1$nameHit == x]) %>%
                                       as.character() %>% paste(collapse = ' ')
                                   })

  #for each TAD get genes IDs
  nbgene_TAD$gene_id.lst = sapply(nbgene_TAD[, 1],
                                  function(x){
                                    names(
                                      data1[data1$nameHit == x]) %>%
                                      as.character()
                                  })

  #for each TAD get gene expressions
  if (is.data.frame(expression.data.frame)) {
    message(paste0(
      expression.data.frame[expression.data.frame[,1] %in% (nbgene_TAD$gene_id.lst %>% unlist), 1] %>%
        length(),"/", nbgene_TAD$gene_id.lst %>% unlist %>% length(),
      " genes have data expression"))
    nbgene_TAD$gene_exp.lst = lapply(nbgene_TAD$gene_id.lst,
                                     function(x){expression.data.frame %>%
                                         filter(expression.data.frame[,1] %in% x)})
  }

  # numbers of TAD according to number of genes
  nbTAD_nbgene = nbgene_TAD %>% dplyr::count(nb_genes, name="nb_TADs")
  nbTAD_nbgene = rbind(c(0, length(tad.gr) - length(nbgene_TAD$TAD_id)), nbTAD_nbgene) # add number of TAD with 0 genes

  # numbers of TAD according to number of genes strand order
  nbTAD_nbgenestrand = nbgene_TAD %>% dplyr::count(gene_strands, nb_genes, name="nb_TADs") %>% arrange(desc(nb_TADs))

  output <- list( geneTADtopo = nbgene_TAD,
                  nbTAD_nbgene = nbTAD_nbgene,
                  nbTAD_nbgenestrand = nbTAD_nbgenestrand)

  class(output) <- "geneTopology"
  return(output)
}



