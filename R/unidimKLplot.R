#' unidimKLplot: KL distributions presented in the same units (Log10(LR))
#'
#' @param res output from distKL function.
#' @return A scatterplot.
#' @export
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' library(forrel)
#' x = linearPed(2)
#' plot(x)
#' x = setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#' x = profileSim(x, N = 1, ids = 2)
#' res <- distKL(ped = x, missing = 5, relative = 1, cores = 10, frequency = NorwegianFrequencies[1:5], numsims = 5)
#' unidimKLplot(res)
unidimKLplot <- function(res) {
res2 <- mutate(res, KLpopped = - KLpopped)
res2 <- gather(res2)
res2
ggplot(data=res2, aes(x=value, group=key, fill=key)) +
  geom_density(adjust=1.5, alpha=.4) +
  theme_ipsum()}
