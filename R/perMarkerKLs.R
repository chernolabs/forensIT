#' perMarkerKLs: Obtain KL per maker for a specific pedigree
#'
#' @param ped Reference pedigree. It could be an input from read_fam() function or a pedigree built with pedtools.
#' @param frequency Allele frequency database.
#' @return An object of class data.frame with KLs.
#' @export
#' @import pedprobr
#' @import forrel
#' @import mispitools
#' @import dplyr
#' @import radiant.data
#'
#' @examples
#' library(forrel)
#' x = linearPed(2)
#' plot(x)
#' x = setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#' x = profileSim(x, N = 1, ids = 2)
#' perMarkerKLs(x,  NorwegianFrequencies[1:5])


perMarkerKLs <- function(ped, frequency) {
KLpedpop <- list()
KLpopped <- list()

for (i in 1:length(ped$MARKERS)) {

#i = 3
#ped <- x
#frequency <- NorwegianFrequencies
df <- as.data.frame(pedprobr::oneMarkerDistribution(ped, 5,i))
names(df) <- "CPT"
df <- df %>% rownames_to_column(var = "Genotype")
df <- df %>%
  dplyr::mutate(Allele1 = sapply(strsplit(as.character(Genotype), "/"), `[`, 1),
         Allele2 = sapply(strsplit(as.character(Genotype), "/"), `[`, 2))

pop <- as.data.frame(frequency[i])
pop <- pop %>% rownames_to_column(var = "Allele")
names(pop) <- c("Allele","freq")

df <- df %>% mutate(RPT = ifelse(pop$freq[match(Allele1, pop$Allele)] == pop$freq[match(Allele2, pop$Allele)], 
                                 pop$freq[match(Allele1, pop$Allele)] * pop$freq[match(Allele2, pop$Allele)], 
                                 2 * pop$freq[match(Allele1, pop$Allele)] * pop$freq[match(Allele2, pop$Allele)]))
df <- replace(df, is.na(df) | df == 0, 1e-20)

KLpedpop[[i]]<-sum(df$CPT*(log10(df$CPT) - log10(df$RPT)))
KLpopped[[i]]<-sum(df$RPT*(log10(df$RPT) - log10(df$CPT)))
}
markName <- names(frequency)
data <- cbind(markName, as.data.frame(cbind(KLpopped)),as.data.frame(cbind(KLpedpop)))
return(data)}
