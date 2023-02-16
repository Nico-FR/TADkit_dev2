#' @title Fetch a chromosome of a mcool file
#'
#' @description Using the Indexes provided byt the cooler data model :
#'       https://cooler.readthedocs.io/en/latest/schema.html,
#'     extract the intra-chromosomal counts for a given chromosome
#'
#' @details

#' @param path  the full path of the mcool file
#' @param resolution the desired resolution
#' @param chrom the selected chromosome
#'
#' @return a sparse matrix
#' @importFrom methods as
#' @import Matrix
#' @import rhdf5
#' @export
#'
#'
coolFetch <- function(path, resolution, chrom) {
    message("\nParsing '", path, "'.")

    uri <- function(path) {
        return(
            paste(
                "resolutions",
                format(resolution, scientific = FALSE),
                path,
                sep = "/"
            )
        )
    }

    # The list of avalaible chromosomes
    chromosomes = h5read(file = path, name = uri("chroms/name"))

    if (!(chrom %in% chromosomes)) {
       stop("\n '", chrom, "' is not a valid chromosome")
   }

    # Index of chrom in the chromosome list
    cid = match(chrom, chromosomes)

    # chrom_offset: which row in the bin table each chromosome first appears (0-based)
    chrom_offset = h5read(file = path, name = uri("indexes/chrom_offset"))
    # chrom_lo: first occurence of chrom in the bin table 0-based
    chrom_lo = chrom_offset[cid] + 1   # 1-based
    # chrom_hi : last occurence of next chrom in the bin table 0-based)
    # chrom_offset[cid + 1] first occurence of next chrom (0-based)
    # in 1-based notation => last occurence of chrom in the bin table
    chrom_hi = chrom_offset[cid + 1]

    # bin1_offset:  which row in the pixel table each bin1 ID first appears.
    bin1_offset = h5read(file = path, name = uri("indexes/bin1_offset"))
    # lo first occurence of bin number chrom_lo in the pixel table (0-based)
    lo = bin1_offset[chrom_lo] + 1  # 1-based
    # hi first occurence of bin number chrom_hi + 1 in the pixels table (0-based)
    hi = bin1_offset[chrom_hi + 1]  # 1-based => last line with chrom_hi

    # chrom intra-chromosomal interactions
    id1 = h5read(file = path, name = uri("pixels/bin1_id"), index=list(lo:hi))
    id2 = h5read(file = path, name = uri("pixels/bin2_id"), index=list(lo:hi))
    interactions = h5read(file = path, name = uri("pixels/count"), index=list(lo:hi))

    # Now we limit the interactions to intrachromosome interactions
    # hence removing bins starting from chrom_hi
    # chrom_hi corresponds to next chrom in 0-based
    # Let the indexes start at 0 (hence remove chrom_lo)
    i = id1[which(id2 < chrom_hi)] - chrom_lo + 1
    j = id2[which(id2 < chrom_hi)] - chrom_lo + 1
    x = interactions[which(id2 < chrom_hi)]

    # sparseMatrix requires 1-based indices
    return(sparseMatrix(i=i+1, j=j+1, x = x))
}
