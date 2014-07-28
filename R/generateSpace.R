#' Obtain the null distribution of a set of intervals
#'
#' Given a set of regions of the genome, get random intervals for a number of 
#' permutations which will be used to determine the null distribution.
#'
#' @param I A GRanges object giving intervals from which to sample from.
#' @param nPermute Number of permutations to run.
#' @param n Number of intervals to generate.
#' @param ms Length of intervals to generate.
#' @param strand.specific Whether the random intervals should have be strand
#' specific or not.
#' @param randomize.strand When \code{strand.specific=TRUE}, whether the strand
#' should be randomized or depend on the strand from the regions in \code{I}.
#' @param seed Seed to be used for reproducing the results.
#'
#' @import GenomicRanges
#'
#' @author Leonardo Collado-Torres
#' @export
#'
#' @examples
#' ## Create some input ranges
#' library('GenomicRanges')
#' gr <- GRanges(c('chr1', 'chrX'), IRanges(c(1e6, 3e6), width = 1e5), strand = 
#'     c('+', '-'))
#' gr
#' 
#' ## Set the seed if you want to reproduce the random results
#' set.seed(20140728)
#' 
#' ## Obtain space of random intervals
#' space <- generateSpace(gr, nPermute = 1000, n = 2, ms = c(1234, 4321))
#' space

generateSpace <- function(I, nPermute = 1000, n = 1, ms = 1e4, 
    strand.specific = FALSE, randomize.strand = FALSE, seed = NULL) {
    stopifnot(is(I, "GRanges"))
    
    if(!is.null(seed)) {
        set.seed(seed)
    }
    space <- randomInterval(I, n = n * nPermute, ms = ms, strand.specific = strand.specific, randomize.strand = randomize.strand)
    values(space)$permutation <- rep(seq_len(nPermute), n)
    return(space)
}