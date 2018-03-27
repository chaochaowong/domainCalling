removeDuplicates <- function(ga) {
    #' This is not recommended for CAR-seq
    #' input: GAlignemtPairs, output: GRanges
    gr <- granges(ga)
    f <- first(ga)
    s <- second(ga)
    strand_location <- paste0(as.character(granges(f)), "::",
                              as.character(granges(s)))
    dup <- duplicated(strand_location)
    message("Total ", length(ga), " paired reads. Remove ", sum(dup), " duplicated reads.")
    gr[!dup]
}


makeFragCoverageTrack <- function(sample_info, 
                                  norm.factors=rep(1, nrow(sample_info)),
                                  format=c("bedGraph", "BigWig", "both"),
                                  paired.end=TRUE,
                                  destination=".",
                                  filter.by.fragment=FALSE,
                                  max.fragment=700, min.fragment=0, cores=1L) {
    require(parallel)
    
    format <- match.arg(format)
    #' sanity check: sample_info (file.exists?) norm.factor, valid length?
    if (!file.exists(destination)) stop(destination, " does not exists.")
    if (length(norm.factors) != nrow(sample_info))
        stop("The length of norm.factor must equal to the number of sample.")
    if (any(!is.numeric(norm.factors) | norm.factors < 0))
        stop("norm.factors muct be positive numeric")
    if (!all(file.exists(sample_info$file_bam)))
        stop("Not all bam files exist.")
   
    mcmapply(makeFragCovPerSample, sample_info$file_bam,
             norm.factors,
             MoreArgs=list(format, paired.end,
                           destination, filter.by.fragment,
                           min.fragment, max.fragment),
             mc.cores=cores, mc.preschedule=FALSE)

    return(invisible())
}

.filterByFragmentLength <- function(ga, min.fragment=0, max.fragment=700) {
    gr <- granges(ga)
    idx <- width(gr) <= max.fragment 
    if (min.fragment > 0)
        idx <- idx & width(gr) >=min.fragment

    ga <- ga[idx]
}

makeFragCovPerSample <- function(bam_file, 
                                 norm.factor=1,
                                 format=c("bedGraph", "BigWig", "both"),
                                 paired.end=TRUE,
                                 destination=".",
                                 filter.by.fragment=FALSE,
                                 min.fragment=0, max.fragment=800) {
    #' this funciton makes spike-normalized pileup bedGraph (bg) and bigWig (bw)
    #' files
    require(GenomicAlignments)
    require(rtracklayer)

    format <- match.arg(format)
    message(basename(bam_file))

    #' checking inputs
    if (!file.exists(bam_file)) stop(bam_file, " does not exists")
    if (!file.exists(destination)) stop(destination, " does not exists.")
    if (!is.numeric(norm.factor) | norm.factor < 0) 
        stop("norm.factor must be positive numeric")

    flag <- scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE)
    what <- c("qname", "flag", "qwidth", "isize")
    param <- ScanBamParam(flag = flag, what = what, tag = "XS")

    #' if pair-ended, use readGAlignmetnPairs. if single-ended, use readGalignmetnPair
    if (paired.end) { #' for pair-ended
        flag <- scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE,
                            isPaired=TRUE, isProperPair=TRUE)
        bamFlag(param) <- flag
        ga <- readGAlignmentPairs(bam_file, param=param)
    }
    else {#' for single-ended
        ga <-  readGAlignmentPairs(bam_file, param=param)
    }
    

    #' filter by fragment length
    if (filter.by.fragment)
        ga <- .filterByFragmentLength(ga, min.fragment=min.fragment,
                                      max.fragment=max.fragment)


    #' remove duplicated fraqment
    #if (remove.duplicate) 
    #    gr <- removeDuplicates(ga)
    #else 
    #    gr <- granges(ga)

    weight <- 1/norm.factor
    prefix <- sub(".bam", "", basename(bam_file))
    prefix <- paste0(prefix, "_scalledBy_", format(weight, digit=3))
    message("Getting coverage ...")
    cov <- GenomicAlignments::coverage(ga, weight=weight)
    
    #' output bedGraph files
    if (format %in% c("bedGraph", "both")) {
        output <- file.path(destination, paste0(prefix, ".bedGraph"))
         message("Exporting ", output)
        rtracklayer::export(cov , con=output, format="bedGraph")
    }
    #' output bigWig files
    if (format %in% c("BigWig", "both")) {
        output <- file.path(destination, paste0(prefix, ".bw"))
        message("Exporting ", output) 
        rtracklayer::export(cov, con=output, format="BigWig")
    } 
  
    return(invisible()) 
}
