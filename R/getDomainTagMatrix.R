#' This script is adapted from the ChIPseeker Bioconudctor package. It mimics the
#' getTagMatrix() function but with better resolution of the domain coverage:
#' The read counts are generated from the bam files, not bed files with score.
#' 
.getDomianTagMatrix <- function(domains, coverage_regions, ga_granges, cores=2L) {
    #' must keep width of all coverage regions must be the same
    tag_vector <- mclapply(domains, function(x) {
        ov_gr <- subsetByOverlaps(ga_granges, x)
        cov <- coverage(ov_gr, weight=1L)
        ov <- subsetByOverlaps(coverage_regions, ov_gr, ignore.strand=TRUE)
        if (identical(0, length(ov))) {
            cov_cnt <- rep(0, width(coverage_regions[1]))
        } 
        if (length(ov) > 0) {
            cov_list <- lapply(ov, function(gr) {
                cov_chrom <- cov[[as.character(seqnames(gr))]] 
                region_cov <- decode(cov_chrom[start(gr):end(gr)])
                if (runValue(strand(gr)) == "-")
                    region_cov <- rev(region_cov)
                return(region_cov)
            })
            cov_cnt <- colSums(do.call(rbind, cov_list))
        }
        return(cov_cnt)
    }, mc.cores=cores, mc.preschedule=FALSE)
    
    tag_matrix <- do.call(rbind, tag_vector)
}

getDomainTagMatrixPerSample <- function(file_bam, domains,
                                        coverage_regions,
                                        ga_granges, cores=2L) {
    #' get the GAlignmentPairs, only works for PE
    require(GenomicAlignments)
    require(parallel)
    flag <- scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE)
    what <- c("qname", "flag", "qwidth", "isize")
    param <- ScanBamParam(flag=flag, what=what, tag="XS",
                      which=coverage_regions)
    ga <- readGAlignmentPairs(file=file_bam,
                          param=param)
    ga_granges <- granges(ga)
    tag_matrix <- .getDomainTagMatrix(domains, coverage_regions=keep_TSS,
                               ga_granges, cores=cores)
}
