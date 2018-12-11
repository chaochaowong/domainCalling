calcDomainFPKM <- function(rse,
    features = rowRanges(rse),
    counts = assays(rse)$counts,
    lib_size = colData(rse)$lib_size)
{   #' Each domain must only have one range and therefore rowRanges must be in the form of
    #' GRanges
    if (class(features) == "GRangesList")
        stop("The rowRanges of the domain must be in the form of GRanges")

    n_features <- length(features)
    n_samples <- ncol(rse)
    E <- matrix(NA_real_, nrow = n_features, ncol = n_samples)
    L <- rep(width(features), n_samples)
    E <- matrix(L - 1, nrow = n_features, ncol = n_samples)
    X <- counts / (E * 1e-3)
    X <- sweep(X, 2, lib_size * 1e-6, FUN = "/")
    return(X)
}
