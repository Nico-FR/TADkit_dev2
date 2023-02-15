#' @title Fetch a chromosome of a mcool file
#'
#' @description Using the Indexes provided byt the cooler data model : https://cooler.readthedocs.io/en/latest/schema.html, extract the intra-chromosome counts for a given chromosome
#'
#' @details

#' @param path  the full path of the mcool file
#' @param resolution the desired resolution
#' @param chrom
#' @param extend Defaults is NULL to take a value 2 times higher than bin.width.
#' @param restrict Default is 2e6 base pair to remove boundary at the extremities of chromosomes.
#'
#' @return a sparse matrix
#' @importFrom methods as
#' @import
#' @import
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

    # Index of chrom in the chromosome list
    cid = match(chrom, chromosomes)

    # chrom_offset which row in the bin table each chromosome first appears (0-based)
    chrom_offset = h5read(file = path, name = uri("indexes/chrom_offset"))
    # chrom_lo: first occurence of chrom in the bin table
    # chrom_hi : last occurence of next chrom in the bin table
    chrom_lo = chrom_offset[cid] + 1   # 1-based
    chrom_hi = chrom_offset[cid + 1]   # 1-based, hence last occurence of chrom in the bin table

    # bin1_offset  which row in the pixel table each bin1 ID first appears.
    bin1_offset = h5read(file = path, name = uri("indexes/bin1_offset"))
    # lo first occurence of bin number chrom_lo in the pixel table (0-based)
    # hi first occurence of bin number chrom_hi + 1 in the pixels table (0-based)
    lo = bin1_offset[chrom_lo] + 1  # 1-based
    hi = bin1_offset[chrom_hi + 1]  # 1-based

    # the chromosome intra-chromosomal interactions
    id1 = h5read(file = path, name = uri("pixels/bin1_id"), index=list(lo:hi))
    id2 = h5read(file = path, name = uri("pixels/bin2_id"), index=list(lo:hi))
    interactions = h5read(file = path, name = uri("pixels/count"), index=list(lo:hi))

    # Now we limit the interactions to intrachromosome interactions
    # hence removing bins starting from chrom_hi (corresponding to next chrom)
    i = id1[which(id2 < chrom_hi)] - chrom_lo + 2
    j = id2[which(id2 < chrom_hi)] - chrom_lo + 2
    x = interactions[which(id2 < chrom_hi)]

    return(sparseMatrix(i=i, j=j, x = interactions))
}
