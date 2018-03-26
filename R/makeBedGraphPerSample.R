removeDuplicates <- function(ga) {
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

makeFragCoverageBedGraph <- function(sample_info, cores=1L,
                         norm.factor=rep(1, nrow(sample_info)),
                         remove.duplicate=FALSE,
                         max.fragment=700,
                         min.fragment=0,
                         outDir=".") {
    require(parallel)
    ## sanity check: norm.factor, valid length? 
    if (length(norm.factor) != nrow(sample_info))
        stop("Length of norm.factor must equal to the number of sample.")

    message("Scaling pileup by scale factor")
    mcmapply(makeBedGraphPerSample, sample_info$file_bam,
             norm.factor,
             MoreArgs=list(bigWig=bigWig,
                 remove.duplicate=remove.duplicate,
                 outDir=outDir),
             mc.cores=cores, mc.preschedule=FALSE)

    return(invisible())
}

makeBedGraphPerSample <- function(bam_file, 
                                  norm.factor=NULL, bigWig=TRUE,
                                  remove.duplicate=FALSE,
                                  outDir=".") {
    #' this funciton makes spike-normalized pileup bedGraph (bg) and bigWig (bw)
    #' files
    require(GenomicAlignments)
    require(rtracklayer)
    
    message(bam_file)

    #' checking inputs
    if (!file.exists(outDir)) stop(outDir, " does not exists.")
    if (!is.null(scaling_factor)) {
        if (scaling_factor > 1 | scaling_factor < 0) 
            stop("Scaling factor must be NULL or a numerical within the range of (0, 1). Setting scaling factor to NULL")
    }

    #' getting GAlignmentPairs
    flag <- scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE,
                        isPaired=TRUE, isProperPair=TRUE, isDuplicate=FALSE)
    what <- c("qname", "flag", "qwidth", "isize")
    #' Get gap aligment reads and convert to insert coverage
    param <- ScanBamParam(flag = flag, what = what, tag = "XS")
    ga <- readGAlignmentPairs(bam_file, param=param)

    #' remove duplicated fraqment
    if (remove.duplicate) 
        gr <- removeDuplicates(ga)
    else 
        gr <- granges(ga)

    prefix <- sub(".bam", "", basename(bam_file))
    if (!is.null(scaling_factor)) {
        cov <- GenomicAlignments::coverage(gr, weight=scaling_factor)
        prefix <- paste0(prefix, "_scalled", format(scaling_factor, digit=2))
    }
    else cov <- GenomicAlignments::coverage(gr)
    
    #' output bed files
    output <- file.path(outDir, paste0(prefix, ".bdg"))
    rtracklayer::export(cov , con=output, format="bedGraph")
    message("Exported ", output)
    #' output bigWig files
    #if (bigWig) {
    #    output <- file.path(outDir, paste0(prefix, ".bw"))
    #    rtracklayer::export(cov, con=output, format="BigWig")
    #    message("Exported ", output) 
    #} 
    return(invisible()) 
}
