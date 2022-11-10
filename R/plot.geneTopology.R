#' @title plot method for "geneTopology" class
#' @param geneTopology_obj A geneTopology object (ie output of genTADtopo function).
#' @param TAD_id TAD id as character
#' @return Return two ggplot graph
#' @import ggplot2
#' @import dplyr
#' @importFrom tidyr gather
#' @export plot.geneTopology
#' @export
plot.geneTopology <- function(geneTopology_obj, TAD_id = NULL){
  if (is.null(TAD_id)) {

  p1 <- ggplot(geneTopology_obj$nbTAD_nbgene %>% dplyr::filter(nb_TADs >= stats::median(geneTopology_obj$nbTAD_nbgene$nb_TADs)), aes(y = nb_TADs, x = nb_genes))+
    geom_line()+geom_point()+
    scale_x_continuous(minor_breaks = seq(0, max(geneTopology_obj$nbTAD_nbgene$nb_genes),1))+
    ggtitle(paste0("Frequency of the number of genes per TAD"))

  data = geneTopology_obj$nbTAD_nbgenestrand %>% dplyr::filter( nb_genes <= 4) %>% dplyr::arrange(desc(nb_TADs))
  p2 <- ggplot(data, aes(x = gene_strands, y = nb_TADs, fill = as.factor(nb_genes)))+
    geom_bar(stat = "identity")+scale_x_discrete(limits = c(data$gene_strands))+
    theme(axis.text.x = element_text(face = "bold", color = "black", size = 10, angle = -60, vjust = 1, hjust = 0))+
    ggtitle(paste0("Frequency of gene strand order per TAD"))+
    labs(fill = "nb gene\nper TAD")

  print(p2)
  print(p1)
  }

  if (is.character(TAD_id) == TRUE & is.null(geneTopology_obj$geneTADtopo$gene_exp.lst)) {
    stop("no expression data in the input file")}

  if (is.character(TAD_id) == TRUE & !is.null(geneTopology_obj$geneTADtopo$gene_exp.lst)) {
    if (is.null(geneTopology_obj$geneTADtopo$gene_exp.lst[[TAD_id]])) {
      stop("no expression data for ", TAD_id)
    } else {
      data = geneTopology_obj$geneTADtopo$gene_exp.lst[[TAD_id]]
      p1 <- ggplot( tidyr::gather(data, "ech", "exp", -gene_id, -gene_name))+
      geom_bar(aes(y=exp,x=gene_name,fill=ech), stat = "identity", position = position_dodge(), color="black")+
      theme(axis.text.x = element_text(color = "black", size = 10, angle = -60, vjust = 1, hjust = 0))
      print(p1)
      }
  }}


