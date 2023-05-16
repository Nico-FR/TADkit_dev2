#' @title Genes topology according to TAD
#'
#' @description From the annotation of TADs and genes, this function return 3 data frames:
#' 1- for each TAD: numbers of genes, their strands and IDs (and the expression data if it is provided).
#' It also return some statistics:
#' 2- the number of TAD according to the number of genes,
#' 3- the number of TAD according to the number of genes and strand.
#'
#'
#' @details Return a list of 3 dataframes.
#'
#' @inheritParams domainHist
#' @param annot.gr `GRanges` with genomic annotations.
#' @param expression.data.frame `dataframe` with expression data (raw counts...).
#' The first two columns should be used to identify the genes. The 1st one must use the same gene id than `names(annot.gr)` (i.e unique ID). The 2nd can be used with usual gene names.
#' Others columns give expression count, one column per sample/experiment.
#'
#' @return Return a S3 class object with 3 data frames.
#'
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @importFrom dplyr count desc
#'
#' @export
#'
#' @examples
#' # not run
#'
#' # get gene annotation
#' # txdb <- makeTxDbFromBiomart(biomart = "ensembl", dataset = "btaurus_gene_ensembl")
#' # genomic.gr = genes(txdb)
#' #
#' # create expression data set with random raw count
#' # expression.data.frame = data.frame(ID = names(genes.gr), Name = names(genes.gr),
#' #   raw_count = sample(1:100, length(genes.gr), replace = TRUE))
#'
#' # output <- geneTADtopo(dataframes2grange(tad_1_10kb.bed, chromsize),
#' #   genomic.gr, expression.data.frame = expression.data.frame)
#'
geneTADtopo <- function(domain.gr, annot.gr, ifoverlap = "best", expression.data.frame = NULL) {

  #local variables:
  nameHit <- nb_genes <- gene_strands <- nb_TADs <- NULL

  if (length(domain.gr) != length(unique(names(domain.gr)))) {
    stop("domain names must be unique!")
  }

  data1 <- domainHist(domain.gr, annot.gr, output = "data", annot.boundary = "start", ifoverlap = ifoverlap)

  # annot.gr ID [ie names(annot.gr)] must be unique
  if (length(annot.gr[duplicated(names(annot.gr))]) > 0) {
    stop("annot.gr must have unique names!")
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
  nbTAD_nbgene = rbind(c(0, length(domain.gr) - length(nbgene_TAD$TAD_id)), nbTAD_nbgene) # add number of TAD with 0 genes

  # numbers of TAD according to number of genes strand order
  nbTAD_nbgenestrand = nbgene_TAD %>% dplyr::count(gene_strands, nb_genes, name="nb_TADs") %>% arrange(desc(nb_TADs))

  output <- list(geneTADtopo = nbgene_TAD,
                  nbTAD_nbgene = nbTAD_nbgene,
                  nbTAD_nbgenestrand = nbTAD_nbgenestrand)

  return(output)
}



