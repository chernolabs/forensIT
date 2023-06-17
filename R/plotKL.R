#' plotKL: scatter plot of obtained KL values from distKL
#'
#' @param res output from distKL function.
#' @return A scatterplot.
#' @export
#' @import ggplot2
#' @import dplyr
#'
#' @examples
#' library(forrel)
#' x = linearPed(2)
#' plot(x)
#' x = setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#' x = profileSim(x, N = 1, ids = 2)
#' res <- distKL(ped = x, missing = 5, relative = 1, cores = 10, frequency = NorwegianFrequencies[1:5], numsims = 5)
#' plotKL(res)
plotKL <- function(res) {
  ggplot(res, aes(x=KLpopped, y=KLpedpop)) + 
  geom_point(size=2, alpha = 0.5) + theme(text = element_text(size = 20)) +
  xlab("KL population to ped") + ylab("KL ped to population") +
  scale_x_continuous(limits = c(0, max(res$KLpopped+0.1)),
                       expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, max(res$KLpedpop+0.1)),
                     expand = c(0, 0))}
