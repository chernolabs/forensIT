#' @title Convert a pedigree to a paramlink object
#' @description Convert a pedigree to a paramlink object
#' @param x pedigree
#' @param verbose print progress
#' @import pedtools
#' @return paramlink object
#' @export
convertPed <- function(x, verbose = FALSE) { #nolint
  if (!requireNamespace("paramlink", quietly = TRUE))
    stop("Package 'paramlink' is not installed")
  famid <- 1
  mlist <- x$MARKERS
  x$MARKERS <- NULL
  p <- cbind(famid, as.matrix(x), 1)
  colnames(p) <- c("FAMID", "ID", "FID", "MID", "SEX", "AFF")
  y <- paramlink::linkdat(p, verbose = verbose)
  if (!is.null(mlist)) {
    mlist <- lapply(mlist, function(m) {
      attributes(m) <-
        list(dim = dim(m),
             name = pedtools::name(m),
             chrom = if(pedtools::isXmarker(m)) 23 else pedtools::chrom(m),
             pos = pedtools::posMb(m),
             nalleles = pedtools::nAlleles(m),
             alleles = pedtools::alleles(m),
             afreq = as.vector(pedtools::afreq(m)),
             missing = 0,
             mutmat = pedtools::mutmod(m),
             class = "marker")
      m
    })
    y <- paramlink::setMarkers(y, mlist)
  }
  y
}