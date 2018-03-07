getDomainsFPKM <-  function(rse,
    features = rowRanges(rse),
    counts = assays(rse)$counts,
    lib_size = colData(rse)$lib_size)
{
    n_features <- length(features)
    n_samples <- ncol(rse)
    E <- matrix(NA_real_, nrow = n_features, ncol = n_samples)
    L <- rep(width(features), n_samples)
    E <- matrix(L - 1, nrow = n_features, ncol = n_samples)
    X <- counts / (E * 1e-3)
    X <- sweep(X, 2, lib_size * 1e-6, FUN = "/")
    return(X)
}

getScalingFactorPerSample <- function(spike_file, multiple.factor=10000) {
    message(spike_file)
    flag <- scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE)
    what <- c("qname", "flag", "qwidth", "isize")
    param <- ScanBamParam(flag = flag, what = what, tag = "XS")
    cnt <- countBam(spike_file, param=param)
    print(cnt)
    multiple.factor/(cnt$records[1]/2)
}

bedGraphToBigWig <- function(bdg_file, outDir=".") {
    if (!file.exists(outDir)) stop("Output directory does not exist.")
    
    require(rtracklayer)
    gr <- import(bdg_file)
    prefix <- sub(".bedGraph", "", basename(bdg_file))
    file <- file.path(outDir, paste0(prefix, ".bw"))
    export(gr, con=file, format="BigWig")              
    message(prefix)
}





alternativeAnnotation <- function(peaks.gr, txdb, orgDb) {
    genes <- genes(txdb)
    genes$TSS <- start(genes)
    neg <- strand(genes)=="-"
    genes[neg]$TSS <- end(genes[neg])
    ov <- findOverlaps(peaks.gr, genes, ignore.strand=TRUE)
    split_ov <- split(ov, queryHits(ov))
    names(split_ov) <- NULL
    which <- elementNROWS(split_ov) > 1 # which one has more than one candenate
    keep_nearest <- lapply(split_ov[which], function(x) {
        tss <- genes[subjectHits(x)]$TSS
        which.peak <- peaks.gr[queryHits(x)[1]]
        speak <- start(which.peak)
        keep <- which.min(abs(speak-tss))[1]
        x[keep]
    })
    tmp <- as.data.frame(do.call(c, keep_nearest))
    clean_ov <- rbind(as.data.frame(unlist(split_ov[!which])), tmp)
    clean_ov <- clean_ov[order(clean_ov$queryHits), ]
    rownames(clean_ov) <- clean_ov$queryHits
    keys <- genes$gene_id[clean_ov$subjectHits]
    clean_ov$gene_id <- keys
    #' choose first gene is there are some overlapping    
    clean_ov$sym <- mapIds(orgDb, keys=keys, keytype="ENSEMBL", column="SYMBOL",
                   multiVals="first")
    #' append the overlapping genes
    peaks.gr$overlapSYMBOL  <- NA
    peaks.gr$overlapGeneId  <- NA
    peaks.gr$overalpRegions <- NA
    
    peaks.gr$overlapSYMBOL[clean_ov$queryHits] <- clean_ov$sym
    peaks.gr$overlapGeneId[clean_ov$queryHits] <- clean_ov$gene_id
    peaks.gr$overalpRegions[clean_ov$queryHits] <- "gene body"
    
    peaks.gr$combinedSYMBOL  <- peaks.gr$overlapSYMBOL
    i <- is.na(peaks.gr$combinedSYMBOL)
    peaks.gr$combinedSYMBOL[i]  <- peaks.gr$SYMBOL[i]
    peaks.gr$combinedAnnotation[i] <- peaks.gr$annotation[i]
    
    peaks.gr
}


sanitizeAnnotation <- function(peaks.gr, EnsDb) {
    stop(is.null(peaks.gr$annotation))
    #' Exon, Intron, Downstream (<=3kb) and distal intergenic
    tmp <- sapply(strsplit(peaks.gr$annotation, " (", fixed=TRUE), "[[", 1)
    peaks.gr$simplified.anno <- tmp

    if (!is.null(EnsDb)) {
        genes <- genes(EnsDb)
    }
    
    peaks.gr
}

