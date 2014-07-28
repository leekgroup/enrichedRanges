#' Create genome tiles for different bin sizes
#'
#' Given a set of chromosomes, tile the genome for a set of bin sizes.
#'
#' @param seqlengths See \link[GenomicRanges]{tileGenome}.
#' @param gaps A GRanges object with the set of gaps to remove from 
#' consideration.
#' @param binSize A vector with the length of the bins to use.
#'
#' @import GenomicRanges
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @examples
#' ## Create some input ranges
#' library('GenomicRanges')
#'
#' ## Get chr lengths
#' data(hg19Ideogram, package = 'biovizBase', envir = environment())
#' seqlengths <- seqlengths(hg19Ideogram)[c("chr1", "chrX")]
#' 
#' ## Define some gaps to remove from tiling
#' gaps <- GRanges(c('chr1', 'chrX'), IRanges(c(2e8, 1e8), width = 1e6))
#'
#' ## Create the tiles
#' tiles <- createTiles(seqlengths, gaps, binSize = c(1e4, 1e5))
#' tiles

createTiles <- function(seqlengths, gaps, binSize = c(1e3, 1e4, 1e5)) {
    bins <- lapply(binSize, function(L) tileGenome(seqlengths, tilewidth = L, cut.last.tile.in.chrom=TRUE))
    
    # drop gap regions from bins
    bins <- lapply(bins, function(x) {
    	o  <- findOverlaps(x, gaps)
    	x[-queryHits(o)]
    })
    names(bins) <- binSize
    
    return(bins)
}