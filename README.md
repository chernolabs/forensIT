#  forensIT: Information Theory Tools for Forensic Analysis

The 'forensIT' package is a comprehensive statistical toolkit tailored for handling missing person cases. 
By leveraging information theory metrics, it enables accurate assessment of kinship, particularly when limited genetic evidence is available. 
With a focus on optimizing statistical power, 'forensIT' empowers investigators to effectively prioritize family members, enhancing the reliability and efficiency of missing person investigations. 

If you want to install forensicIT, please use:

``` r
install.packages("devtools")
library(devtools)
install_github("marsicoFL/forensIT")
library(forensIT)
```

Using other packages
``` r
library(forrel)
library(mispitools)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyr)
library(pedprobr)
```

KL (dits units) per marker could be obtained with the following commands:
``` r
freqs = getfreqs(Argentina[1:15])
x = linearPed(2)
plot(x)
x = setMarkers(x, locusAttributes = freqs)
x = profileSim(x, N = 1, ids = 2)
perMarkerKLs(x,  freqs)
``` 

![](Ped.png)<!-- -->



KL (dits) distributions for a specific relative to be potentially incorporated to the pedigree could be obtained as follows:
``` r
y = linearPed(2)
x = setMarkers(y, locusAttributes = freqs)
res <- distKL(ped = x, missing = 5, relative = 1, cores = 10, frequency = freqs, numsims = 5)
res
``` 
    ##     KLpopped  KLpedpop
    ##    1.1708052 1.6985469
    ##    1.0524863 1.1335089
    ##    1.0738553 1.1462661
    ##    0.9431456 0.9980392
    ##    1.4592184 2.0374789
    ##    0.7972688 0.8059407
    ##    1.0690941 1.0794834
    ##    1.1146577 1.2404207
    ##    0.8978531 1.0942142
    ##    1.0621151 1.2631875



A scatterplot could be obtained with the following code:
``` r
plotKL(res)
``` 
![](distKL.png)<!-- -->

KL distributions presented in the same units (Log10(LR)). Note that this implies plotting -KL(pop to ped)
``` r
unidimKLplot(res)
``` 
![](uniKL.png)<!-- -->
