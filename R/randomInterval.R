#' Get random intervals based on a set of regions
#'
#' Given a set of regions, generate genomic intervals uniformly at 
#' random from those regions. These can be strand specific or not.
#'
#' @param I A GRanges object giving intervals from which to sample from.
#' @param n Number of intervals to generate.
#' @param ms Length of intervals to generate. It can be a vector but 
#' recycling rules will be used if \code{ms} is of a length smaller than 
#' \code{n}.
#' @param strand.specific Whether the random intervals should have be strand
#' specific or not.
#' @param randomize.strand When \code{strand.specific=TRUE}, whether the strand
#' should be randomized or depend on the strand from the regions in \code{I}.
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
#' ## Obtain 4 randomg intervals
#' randomInterval(gr, 4)
#'
#' ## Strand specific without randomizing the strand
#' randomInterval(gr, 4, strand.specific = TRUE)
#'
#' ## This function is a modified version of seqbias::random.intervals
#' ## Check the original code below:
#' library('seqbias')
#' random.intervals
#' ## Some key differences are that randomInterval works with repeated
#' ## seqnames and starts at the beginning of the input ranges, not at base 1.

randomInterval <- function(I, n = 1, ms = 10000,
    strand.specific = FALSE, randomize.strand = FALSE) {
    stopifnot(is(I, "GRanges"))
    
    
    seqs <- width(I)
    sample_sequence <- function(m, seqs = seqs) {
        ps <- pmax(0, as.numeric(seqs - m + 1))
        psum <- sum(ps)
        if (psum == 0) {
            stop(paste("no sequence is long enough to sample an interval of length ", 
                m, sep = ""))
        }
        prob <- ps / psum
        sample(seq_len(length(seqs)), size = 1, prob = prob)
    }
    ## vector ms
    xs <- sapply(cbind(1:n, ms)[, 2], sample_sequence, seqs = seqs)
    
    starts <- mapply(function(min, max) {
        as.integer(round(runif(n = 1, min = min, max = max)))
        }, start(I)[xs], end(I)[xs])
    
    if(strand.specific) {
        if(randomize.strand) {
            strand <- sample(c("+", "-"), size = length(starts), replace = TRUE)
        } else {
            strand <- strand(I)[xs]
        }
        
    } else {
        strand <- "*"
    }
    
    seqnames <- seqnames(I)[xs]
    gr <- GRanges(seqnames = seqnames, ranges = IRanges(starts, starts + 
        ms - 1), strand = strand, seqlengths = seqlengths(I)[unique(seqnames)] )
    
    return(gr)
}