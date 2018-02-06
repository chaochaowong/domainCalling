getCountBamPerSample <- function(file_name) {
    flag <- scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE)
    what <- c("qname", "flag", "qwidth", "isize")
    param <- ScanBamParam(flag = flag, what = what, tag = "XS")
    cnt <- countBam(file_name, param=param)
    cnt$records[1]/2
}

getBamInfoPerSample <- function (file_bam)
{   ## this script is copied from SGSeq
    require(Rsamtools)
    if (is(file_bam, "BamFile")) {
        file_tmp <- file_bam
    }
    else {
        file_tmp <- BamFile(file_bam)
    }

    flag <- scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE)
    what <- c("qname", "flag", "qwidth", "isize")
    param <- ScanBamParam(flag = flag, what = what, tag = "XS")
    bam <- scanBam(file = file_tmp, param = param)[[1]]
    read_length <- median(bam$qwidth, na.rm = TRUE)
    isize <- bam$isize
    frag_length <- median(isize[which(isize > 0)], na.rm = TRUE)
    x <- data.frame(read_length = read_length, 
                    frag_length = frag_length, stringsAsFactors = FALSE)
    x$lib_size <- length(unique(bam$qname))
    
    message("Complete ", file_bam, ".")
    return(x)
}

checkSampleInfoInput <- function(sample_info, initialized=TRUE) {
    #' similar to SGSeq:::checkSampleInfo

    isOr <-   function (object, class2) {
        any(sapply(class2, is, object = object))
    }

    col_type <- list(sample_name = "character",
                     file_bam = c("character", "BamFileList"),
                     spike_bam= c("character",  "BamFileList"))
    if (!is(sample_info, "data.frame") && !is(sample_info, "DataFrame")) {
        err <- paste("sample_info must be a data frame with columns", 
            paste(names(col_type), collapse = ", "))
        stop(err, call. = FALSE)
    }

    missing <- !names(col_type) %in% names(sample_info)
    if (any(missing)) {
        err <- paste("sample_info is missing column(s)", paste(names(col_type)[missing], 
            collapse = ", "))
        stop(err, call. = FALSE)
    }

    invalid <- !mapply(isOr, sample_info[names(col_type)], col_type)
    if (any(invalid)) {
        err <- paste("sample_info column(s)", paste(names(col_type)[invalid], 
            collapse = ", "), "\n", "must be of type", paste(unstrsplit(col_type[invalid], 
            "/"), collapse = ", "))
        stop(err, call. = FALSE)
    }

    flag <- file.exists(sample_info$file_bam)
    if (!all(flag)) {
        err <- paste0("Bam file not exists: ", sample_info$file_bam[!flag], "/")
        stop(err, call.=FALSE)
    }

    #flag <- file.exists(sample_info$spike_bam)
    #if (!all(flag)) {
    #    err <- paste0("Spike bam file not exist: ", sample_info$spike_bam[!flag], "/")
    #    stop(err, call.=FALSE)
    #}

    if (!initialized) {
        #' check existance of spike_factor
        flag <- duplicated(sample_info$sample_name)
        if (any(flag)) {
            err <- paste0("sample names are duplicated: ", sample_info$sample_name[flag])
            stop(err, call.=FALSE)
        }
        
        flag <- "frag_length" %in% colnames(sample_info)
        if (!any(flag)) {
            err <- paste0("Must provide the mean fragment length.")
            stop(err, call.=FALSE)
        }
    }
}

getSampleInfo <- function(sample_info, cores=2L) {
    #' si must be a DataFrame or data.frame instance and have mandated
    #' columns sample_name, file_bam and spike_bam 

    #' return with lib_size, spike_counts, spike_percentage, spike_factor
    #' and scaling_factor

    #' sanity check (sample_name, file_bam, spike_bam)
    require(parallel)
    require(Rsamtools)
    
    si <- sample_info

    message("Checking input.")
    checkSampleInfoInput(si, initialize=TRUE)

    message("Getting sample information.")
    info <- mclapply(si$file_bam, getBamInfoPerSample,
                     mc.cores=cores, mc.preschedule=FALSE)
    info <- do.call(rbind, info)
    spike_count <- sapply(si$spike_bam, getCountBamPerSample)
    info$spike_count <- unlist(spike_count)

    si <- cbind(si, info)
    si$spike_percentage <- si$spike_count/si$lib_size
         
    return(si)
}
