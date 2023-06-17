#' distKL: KL distribution obtained for specific relative contributor
#'
#' @param ped Reference pedigree. It could be an input from read_fam() function or a pedigree built with pedtools.
#' @param frequency Allele frequency database.
#' @param relative Selected relative.
#' @param numsims Number of simulated genotypes.
#' @param missing Missing person
#' @param cores Enables parallelization.
#' @param frequency Allele frequency database.
#' @return An object of class data.frame with KLs.
#' @export
#' @importFrom pedprobr oneMarkerDistribution
#' @import forrel
#' @importFrom mispitools getfreqs
#'
#' @examples
#' library(forrel)
#' x = linearPed(2)
#' plot(x)
#' x = setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#' x = profileSim(x, N = 1, ids = 2)
#' distKL(ped = x, missing = 5, relative = 1, cores = 10, frequency = NorwegianFrequencies[1:5], numsims = 5)
distKL <- function(ped, missing, relative, frequency, numsims = 100, cores = 1) {
  peds = forrel::profileSim(ped, numsims, ids = relative, numCores = cores)
  out <- list()
  for (j in 1:length(peds)) {
    dat <- perMarkerKLs(peds[[j]], MP = missing, frequency)
    out1 <- c(sum(unlist(dat$KLpopped)), sum(unlist(dat$KLpedpop)))
    names(out1) <- c("KLpopped", "KLpedpop")
    out[[j]] <- out1
  }
  result <- as.data.frame(do.call(rbind, out))
  names(result) <-  c("KLpopped", "KLpedpop")
  return(result)}
