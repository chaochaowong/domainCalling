.RangedSEToBed <- function(rse, file) {
    gr <- granges(rse)
    gr$name <- rownames(rse)
    
    rtracklayer::export(gr, file, format="BED")
}

.RangedSEToXLS <- function(rse, file=NULL) {
    gr <- granges(rse)
    cnt <- assays(rse)[["counts"]]
    cpm <- assays(rse)[["logCPM"]]
    colnames(cpm) <- paste0("logCPM_", colnames(cpm))
    gr$name <- rownames(rse)
    gr$domain_width <- width(rse)
    mcols(gr) <- append(mcols(gr),
                        cbind(as(cnt, "DataFrame"), as(cpm, "DataFrame")))
    if (is.null(file))
        return(gr)

    if (!is.null(file)) {
        write.csv(as.data.frame(gr), file)
        return(invisible())
    }
}

.filterByAvgCount <- function(data, threshold=5L) {
    dge <- csaw::asDGEList(data)
    abundances <- edgeR::aveLogCPM(dge)
    keep <- abundances > edgeR::aveLogCPM(threshold, lib.size=mean(data$totals))
    data <- data[keep]
}

getLogCPM <- function(data) {
    dge <- csaw::asDGEList(data)
    res <- edgeR::cpm(dge, lib.size=data$totals, log=TRUE)
    res
}

.filterByGenomeEnrichment <- function(data, param, spikeNorm) {
    #' find non-specific enrichemnt throughout the genome
    bin.size <- 50000L
    binned <- windowCounts(data$bam.files, bin=TRUE, width=bin.size, param=param)
    if (spikeNorm)
        binned <- .normCntsBySpikeFactor(binned, data$norm.factors)
    filter.stat <- filterWindows(data, background=binned, type="global")
    keep <- filter.stat$filter > log2(3)
    data <- data[keep]
}

.addColData <- function(data, sample_info) {
    #' add spike factor rownames
    colnames(data) <- rownames(sample_info$sample_name)
    colData(data) <- append(colData(data), as(sample_info, "DataFrame"))
    data
}

getSpikeNormFactor <- function(sample_info, fragment){
    #' use edgeR and csaw to determine the spike factor
    param <- readParam(pe="both", max.frag=700)
    lib_size <- sample_info$lib_size
    binned <- windowCounts(sample_info$spike_bam, bin=TRUE,
                           width=2000, param=param,
                           ext=fragment)
    spike_data <- assays(binned)[[1]]
    spike_factor <- normOffsets(spike_data, lib.size=lib_size, type="scaling")
    names(spike_factor) <- sample_info$sample_name
    spike_factor
}

addSpikeFactor <- function(si, multiplying_factor=NA) {
    #' add spike_factor_ratio and spike_factor_edgeR columns to sample_info (si)
    if (is.na(multiplying_factor))
        si$multiplying_factor <- median(si$spike_count)
    
    si$spike_factor_ratio <- si$spike_count / median(si$spike_count)
    
    fragment <- list(init.ext=si$frag_length,
                     final.txt=as.integer(mean(si$frag_length)))
    si$spike_factor_edgeR <-
        domainCalling:::getSpikeNormFactor(si, fragment)
    si
}

estimateSpikeFactor <- function(si, method=c("csaw", "simple_ratio"), multiplying_factor=NA) {
    #' sanitized check si, must be the regular si structure (spike bam? lib_size?)
    method <- match.arg(method, c("csaw", "simple_ratio"))
    if (method == "simple_ratio") {
        if (is.na(multiplying_factor)) {
            md <- median(si$spike_count)
            message("estimateSpikeFactor: mutliplying factor not provided, deafult to median of spike count ", md) 
            si$multiplying_factor <- md
        }
        spike_factor <- si$spike_count / median(si$spike_count)
        names(spike_factor) <- rownames(si)
    }
    
    if (method == "csaw") {
    fragment <- list(init.ext=si$frag_length,
                     final.txt=as.integer(mean(si$frag_length)))
    spike_factor <-
        domainCalling:::getSpikeNormFactor(si, fragment)
    }
    
    spike_factor
}

saveDomains <- function(domains, destination=".", prefix=NULL) {
    if (is.null(prefix) | is.na(prefix))
        prefix <- "domains"
    
    system(paste0("mkdir -p ", destination))
    save(domains, file=file.path(destination, paste0(prefix, ".rda")))
    file <- file.path(destination, paste0(prefix, ".bed"))
    .RangedSEToBed(domains, file=file)
    file <- file.path(destination, paste0(prefix, ".csv"))
    .RangedSEToXLS(domains, file=file)
    return(invisible())
}

csawDomainCalling <- function(sample_info,
                              background_idx=NULL,
                              paired.end=TRUE,
                              minq=0,
                              discard=GRanges(),
                              threshold=5L,
                              window.spacing=100, 
                              window.width=window.spacing,
                              spacing=100, 
                              dedup=FALSE,
                              restrict=NULL,
                              norm.factors=rep(1, nrow(sample_info)),
                              filter.By.GlobalEnrichment=FALSE,
                              min.width=300,
                              max.width=30000,
                              window.gapwidth=2000,
                              cores=1L,
                              destination=".",
                              prefix=NULL) {

    #' use sliding window but find significant changes in binding patterns
    #' find coverage of enriched genomic regions
        
    #' sanity check: sample_info
    #' need sanity on the input parameters
    #checkSampleInfoInput(sample_info, initialize=FALSE)

    bam_files <- sample_info$file_bam
    fragment <- list(init.ext=sample_info$frag_length,
                     final.txt=as.integer(mean(sample_info$frag_length)))
    require(csaw)
    require(BiocParallel)
       
    bpparam <- SnowParam(worker=cores)
    param <- readParam(pe="both", max.frag=700, dedup=dedup,
                       restrict=restrict, minq=minq,
                       discard=discard,
                       BPPARAM=bpparam)
    #' fix for single-end
    if (!paired.end) {
        param <- reform(param, pe="none")
        fragment <- mean(sample_info$read_length)
    }

    #' sanity check on the validity of norm.factors
    #if (!spikeNorm) norm.factors <- rep(1, nrow(sample_info))

    print(param)
    
    #' (1) NO background; use count based filter (count > threshold)
    if (is.null(background_idx)) {
        data <- windowCounts(bam_files, param=param,
                             ext=fragment,
                             filter=2*length(bam_files),
                             width=window.width,
                             spacing=window.spacing)
        data <- .addColData(data, sample_info)
        colnames(data) <- data$sample_name
        data$norm.factors <- norm.factors
        
        if (filter.By.GlobalEnrichment) 
            data <- .filterByGenomeEnrichment(data, param, spikeNorm)
        #' filter window whose average count is above the threshold
        data <- .filterByAvgCount(data, threshold)
    }
    
    #' filter by negative control: non-specific enrichemnt throughout the genome
    if (!is.null(background_idx)) {
        logCPM.threshold <- log2(3)
        counts <- windowCounts(bam_files, param=param,
                               ext=fragment, filter=2*length(bam_files),
                               width=window.width,
                               spacing=window.spacing)
        counts <- .addColData(counts, sample_info)
        counts$norm.factors <- norm.factors
        
        counts$background <- FALSE
        counts$background[background_idx] <- TRUE

        #' filter non-specific enrichment (estimating the background abundance)
        if (filter.By.GlobalEnrichment) 
            counts <- .filterByGenomeEnrichment(counts, param, spikeNorm)

        data <- counts[, -background_idx]
        background <- counts[, background_idx]

        #' composition bias? use larger bins, a way to normalize, need investigation
        in.binned <- windowCounts(bam_files, bin=TRUE, width=8000L,
                                  param=param)
        in.binned$norm.factors <- norm.factors

        data.binned <- in.binned[, -background_idx]
        background.binned <- in.binned[, background_idx]

        #' The filter statistics are defined as the difference between data and
        #' background abundances - the log-fold increase in the counts
        filter.stat <- filterWindows(data=data, background=background,
                                     type="control", prior.count=5)

        keep <- filter.stat$filter > logCPM.threshold
        data <- data[keep]
        #' filter by average count size (average over the non_background samples)
        data <- .filterByAvgCount(data, threshold)
    }

    #' cluster windows into regions
    merged <- mergeWindows(data, tol=window.gapwidth, max.width=max.width)
    
    length(merged$region)
    summary(width(merged$region))
    
    merged_counts <- regionCounts(bam_files, merged$region, ext=fragment,
                                  param=param)
    #' filter by min.width
    keep <- width(merged_counts) > min.width
    merged_counts <- merged_counts[keep]
    
    merged_counts <- .addColData(merged_counts, sample_info=sample_info)
    colnames(merged_counts) <- merged_counts$sample_name
    rownames(merged_counts) <- paste0("domain_", 1:nrow(merged_counts))
    merged_counts$norm.factors <- norm.factors
    assays(merged_counts)[["logCPM"]] <-  getLogCPM(merged_counts)
    
    saveDomains(merged_counts, destination=destination, prefix=prefix)
    
    return(merged_counts)
}


domainCallingFromBed <- function(bed_file, format="bedGraph",
                               threshold=10L, standard.chromosome=TRUE,
                               species="Homo_sapiens",
                               min.range=180L, max.range=25000,
                               gapwidth=800L,
                               outDir=".", prefix=NULL) {
    ## if broad, set min.gapwidth=30 and min.range 100 to 200
    ## if narrow, min.gapwidth=
    ## names(genomeStyles()) for species

    require(rtracklayer)
    if (!file.exists(outDir)) {
        step("Output directory ", outDir, " does not exist.")
    }

    message("Importing ", bed_file)
    gr <- import(bed_file, format=format)
    if (standard.chromosome)
        gr <- keepStandardChromosomes(gr, species=species,
                                      pruning.mode="coarse")
    #' filter
    filt <- gr[gr$score > threshold]
    merge <- reduce(filt, drop.empty.ranges=TRUE,
                    min.gapwidth=gapwidth, ignore.strand=TRUE)
    keep <- width(merge) > min.range & width(merge) < max.range
    peaks <- merge[keep]
    peaks$name <- paste0("peaks_", c(1:length(peaks)))
    #' add score later (number of fragment hits)
    #' export bed files
    if (is.null(prefix))
        prefix <- sub(".bdg", "", basename(bed_file))

    file <- file.path(outDir, paste0(prefix, "_peaks_thres", threshold, ".bed"))
    export(peaks, con=file, format="BED")
    file <- file.path(outDir, paste0(prefix, "_peaks_thres", threshold, ".csv"))
    write.csv(as.data.frame(peaks), file=file)
    message("Exporting ", file)
    
    return(peaks) ## GRanges
}
