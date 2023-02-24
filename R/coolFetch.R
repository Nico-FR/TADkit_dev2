#' @title Fetch a chromosome of a cool or mcool file
#'
#' @description From cool or mcool files, `coolFetch` import the interaction matrix for one chromosome.
#' If balance = `TRUE` returns the normalized counts
#'
#' @details Using the Indexes provided by the cooler data model :
#' https://cooler.readthedocs.io/en/latest/schema.html,
#' extract the intra-chromosomal counts for a given chromosome
#'
#' @param cool.path The full path of the cool or mcool file.
#' @param bin.width Bin width (i.e resolution) of the mcool matrix in bp. It must not be used for cool file.
#' @param chr The selected chromosome.
#' @param balance Logical. Weather or not to use balanced counts instead of raw counts. Default = `FALSE`.
#'
#' @return A `dgCMatrix` object: upper triangular and sparse Matrix
#'
#' @import Matrix
#' @importFrom rhdf5 h5read
#' @importFrom methods as
#'
#' @export
#'
#'
coolFetch <- function(cool.path, chr, bin.width = NA, balance = FALSE) {
    message("\nParsing '", cool.path, "'.")

  #mcool path
    uri <- function(cool.path) {
        if (!is.numeric(bin.width)) return(cool.path)
        return(
            paste(
                "resolutions",
                format(bin.width, scientific = FALSE),
                cool.path,
                sep = "/"
            )
        )
    }
    #TODO: add message for unavailable resolution and available resolution
    #command to get available resolutions (it must be possible to make it simpler):
    #(rhdf5::h5ls(cool.path) %>% filter(group == "/resolutions"))$name %>% as.numeric %>% sort %>% format(scientific = FALSE) %>% paste0(collapse = ", ")

    # The list of avalaible chromosomes
    chromosomes = rhdf5::h5read(file = cool.path, name = uri("chroms/name"))

    if (!(chr %in% chromosomes)) {
       stop("\n '", chr, "' is not a valid chromosome.", "\nChromosomes available are: ", paste0(chromosomes, collapse = ", "), ".")
   }

    # Index of chrom in the chromosome list
    cid = match(chr, chromosomes)

    # chrom_offset: which row in the bin table each chromosome first appears (0-based)
    chrom_offset = rhdf5::h5read(file = cool.path, name = uri("indexes/chrom_offset"))
    # chrom_lo: first occurence of chrom in the bin table 0-based
    chrom_lo = chrom_offset[cid] + 1   # 1-based
    # chrom_hi : last occurence of next chrom in the bin table 0-based)
    # chrom_offset[cid + 1] first occurence of next chrom (0-based)
    # in 1-based notation => last occurence of chrom in the bin table
    chrom_hi = chrom_offset[cid + 1]

    # bin1_offset:  which row in the pixel table each bin1 ID first appears.
    bin1_offset = rhdf5::h5read(file = cool.path, name = uri("indexes/bin1_offset"))
    # lo first occurence of bin number chrom_lo in the pixel table (0-based)
    lo = bin1_offset[chrom_lo] + 1  # 1-based
    # hi first occurence of bin number chrom_hi + 1 in the pixels table (0-based)
    hi = bin1_offset[chrom_hi + 1]  # 1-based => last line with chrom_hi

    # chrom intra-chromosomal interactions
    id1 = rhdf5::h5read(file = cool.path, name = uri("pixels/bin1_id"), index = list(lo:hi))
    id2 = rhdf5::h5read(file = cool.path, name = uri("pixels/bin2_id"), index = list(lo:hi))
    interactions = h5read(file = cool.path, name = uri("pixels/count"), index = list(lo:hi))

    # Now we limit the interactions to intrachromosome interactions
    # hence removing bins starting from chrom_hi
    # chrom_hi corresponds to next chrom in 0-based
    # Let the indexes start at 0 (hence remove chrom_lo)
    i = id1[which(id2 < chrom_hi)] - chrom_lo + 1
    j = id2[which(id2 < chrom_hi)] - chrom_lo + 1
    x = interactions[which(id2 < chrom_hi)]

    m = Matrix::sparseMatrix(i = i + 1, j = j + 1, x  =  as.numeric(x))

    if (balance) {
      message("\nBalancing")
      # Fetch the weights corresponding to the chromosome
      bins <- data.frame(
         chromosome = rhdf5::h5read(file = cool.path, name = uri("bins/chrom")),
         start = rhdf5::h5read(file = cool.path, name = uri("bins/start")),
         end = rhdf5::h5read(file = cool.path, name = uri("bins/end")),
         weight = rhdf5::h5read(file = cool.path, name = uri("bins/weight"))
      )
      # all indexes, id1, id2 are 0-based, hence we set index as 0-based
      bins$index = 0:(nrow(bins) - 1)
      # restricting weights to the actual bins under consideration
      min_id1 = min(id1[which(id2 < chrom_hi)])
      max_id2 = max(id2[which(id2 < chrom_hi)])
      w = bins$weight[bins$index >= min_id1 & bins$index <= max_id2] #use instead??: w = (bins %>% filter(chromosome == chr))$weight
      #upper matrix weight for balancing
      mat_weight = Matrix::triu((w %*% t(w)))
      mat_weight[is.na(mat_weight)] <- 0 #remove NaN
      #cell by cell multiplication by the matrix weight
      return(m * mat_weight)
    } else {
      return(m)
    }
}
