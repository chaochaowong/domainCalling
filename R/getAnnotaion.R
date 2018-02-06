getAnnotationFromGR <- function(peaks.gr, TxDb, EnsDb=NULL, annoDb=NULL, rmsk=NULL,
                           annotate.repeat=FALSE) {
    require(ChIPseeker)
    if (is.null(peaks.gr$name))
        peaks.gr$name <- names(peaks.gr)
    #names(peaks.gr) <- peaks.gr$name
    
    peakAnno <- annotatePeak(peaks.gr, tssRegion = c(-3000, 3000),
                             genomicAnnotationPriority=c("Exon", "Intron",
                                 "Promoter", "5UTR", "3UTR",
                                 "Downstream", "Intergenic"),
                             overlap="all",
                             TxDb = TxDb, annoDb = annoDb)
    print(peakAnno)

    #' convert to GRanges
    #' (1) adding symbols from EnsDb
    genes <- genes(EnsDb)
    tmp <- as.GRanges(peakAnno)
    tmp$symbol <- genes[tmp$geneId]$gene_name
    tmp$regions <- sapply(strsplit(tmp$annotation, " (", fixed=TRUE), "[[", 1)
    
    df <- mcols(tmp)
    rownames(df) <- tmp$name
    peaks.mcols <- mcols(peaks.gr)
    rownames(peaks.mcols) <- peaks.gr$name

    #' add df to mcols(peaks.gr), one by one
    for (i in 2:ncol(df)) {
        peaks.mcols[, colnames(df)[i]] <- NA
        peaks.mcols[rownames(df), i] <- df[, i]
    }
    mcols(peaks.gr) <- peaks.mcols
    
    if (annotate.repeat) {
        anno_rmsk <- SeqPipelineTools:::olRMSK(peaks.gr, rmsk)
        mcols(peaks.gr) <- append(mcols(peaks.gr), anno_rmsk)
    }

    return(list(peaks.gr=peaks.gr, peakAnno=peakAnno))
    
}

getAnnotationFromGR.2 <- function(rse, TxDb, annoDb=NULL, rmsk=NULL,
                                  destination=NULL, prefix=NULL) {
    require(ChIPseeker)
    peaks.gr <- granges(rse)
    #' clean up the mcols
    mcols
    if (is.null(peaks.gr$name))
        peaks.gr$name <- names(peaks.gr)

    peakAnno <- annotatePeak(peaks.gr, tssRegion = c(-3000, 3000),
                             genomicAnnotationPriority=c("Exon", "Intron",
                                 "Promoter", "5UTR", "3UTR",
                                 "Downstream", "Intergenic"),
                             overlap="all",
                             TxDb = TxDb, annoDb = annoDb)
    print(peakAnno)

    #' make sure
    tmp <- as.GRanges(peakAnno)
    tmp$regions <- sapply(strsplit(tmp$annotation, " (", fixed=TRUE), "[[", 1)
    
    df <- mcols(tmp)
    rownames(df) <- tmp$name
    peaks.mcols <- mcols(peaks.gr)
    rownames(peaks.mcols) <- peaks.gr$name

    #' add df to mcols(peaks.gr), one by one
    for (i in 2:ncol(df)) {
        peaks.mcols[, colnames(df)[i]] <- NA
        peaks.mcols[rownames(df), i] <- df[, i]
    }
    mcols(peaks.gr) <- peaks.mcols
    
    #' save annotation and counts to csv file
    if (!is.null(destination)) {
        if (is.null(prefix)) prefix <- "domains_annotation"
        
        gr <- .RangedSEToXLS(rse)
        mcols(gr) <- append(mcols(peaks.gr), mcols(gr)[, -1])
        write.csv(as.data.frame(gr),
                  file=file.path(destination, paste0(prefix, ".csv")))
        #' how about figure for the annotation frequency?
        #file <- file.path(destination, paste0(prefix, "_AnnoBar.pdf"))
        #pdf(file, width=6, height=2)
        #plotAnnoBar(peakAnno)
        #dev.off()
    }
    
    return(list(peaks.gr=peaks.gr, peakAnno=peakAnno))
    
}
