#' Region enrichment
#'
#' Determine the relative enrichment of a set of ranges versus random intervals #' across a set of tiles.
#'
#' @param query A GRanges object with the regions of interest to check for 
#' enrichment.
#' @param space The result from using \link{generateSpace}. You can use the 
#' query ranges width as shown in the example.
#' @param tiles The result from using \link{createTiles}.
#'
#' @return A list with two components. The first one is a vector with the number
#' of overlaps for each of the query ranges versus the set space ranges for the 
#' number of permutations used to create the space.
#'
#' The second component is a list with a vector for each of the tiles
#' \code{binSize}. It has the enrichment of the query ranges versus the space
#' ranges for each of the permutations.
#'
#' @import GenomicRanges
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @examples
#' ## Create some input ranges
#' library('GenomicRanges')
#' I <- GRanges(c('chr1', 'chrX'), IRanges(c(1e6, 3e6), width = 1e5), strand = 
#'     c('+', '-'))
#' 
#' ## Intervals of interest to resample from
#' gr <- GRanges(c('chr1', 'chrX'), IRanges(c(1e6, 3e6), width = 1e3), strand = 
#'     c('+', '-'))
#'
#' ## Set the seed if you want to reproduce the random results
#' set.seed(20140728)
#' 
#' ## Obtain 4 random intervals
#' space <- generateSpace(I, gr, nPermute = 1000)
#'
#' ## Get chr lengths
#' data(hg19Ideogram, package = 'biovizBase', envir = environment())
#' seqlengths <- seqlengths(hg19Ideogram)[c("chr1", "chrX")]
#' 
#' ## Define some gaps to remove from tiling
#' gaps <- GRanges(c('chr1', 'chrX'), IRanges(c(2e8, 1e8), width = 1e6))
#'
#' ## Create the tiles
#' tiles <- createTiles(seqlengths, gaps, binSize = c(1e3, 1e4))
#' 
#' ## Build query
#' query <- GRanges(c('chr1', 'chrX'), IRanges(c(1e6, 3e6), width = c(1234, 
#'     4321)), strand = c('+', '-'))
#'
#' ## Get the results
#' enrich <- enrichRanges(query, space, tiles)
#'
#' ## Number of overlaps between query and space ranges for all permutations
#' ## per query range
#' head(enrich["Query.vs.Space"])
#'
#' ## Enrichment result for first tiles binSize comparing query overlaps with
#' ## the tiles and space overlaps with the tiles for each permutation
#' head(enrich[["Query.and.Space.vs.Tiles"]][[1]])
#' 
#' ## Should be equal to the number of permutations used
#' length(enrich[["Query.and.Space.vs.Tiles"]][[1]])
 
enrichRanges <- function(query, space, tiles) {
    
    stopifnot(is(query, "GRanges"))
    stopifnot(is(space, "GRanges"))
    stopifnot('permutation' %in% colnames(values(space)))
    
    ov.query <- findOverlaps(query, space)
    querySpace <- table(queryHits(ov.query))
    
    
    tilesRes <- lapply(tiles, .TileStep, query = query, space = space)
    
    return(list("Query.vs.Space" = querySpace, "Query.and.Space.vs.Tiles" = tilesRes))
}


.TileStep <- function(tile, query, space) {
    ## Overlaps vs tile
    o.tile.query <- findOverlaps(query, tile)
    o.tile.space <- findOverlaps(space, tile)
    
    ## Find unique tile hits
    query.hits <- unique(subjectHits(o.tile.query))
    space.hits <- lapply(split(subjectHits(o.tile.space), space$permutation[queryHits(o.tile.space)]), unique)
    
    ## Calculate pieces
    topleft <- sapply(space.hits, function(x) {
        sum(query.hits %in% x)
    })
    topright <- length(query.hits) - topleft
    bottomleft <- sapply(space.hits, length) - topleft
    bottomright <- length(tile) - topleft - topright - bottomleft
    
    ## Finish
    res <- topleft/topright/bottomleft * bottomright
    return(res)
}