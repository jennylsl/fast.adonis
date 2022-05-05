## fast.adonis
fast.adonis is a computationally efficient nonparametric multivariate analysis of microbiome data for large-scale studies.

This code generates results consistent with adonis/adonis2 but much faster. For complex sampling studies, fast.adonis integrates sampling weight algebraically without replicating samples to mimic the source population; thus, analysis can be completed much faster without requiring a large amount of memory.

## packages required
R package "vegan" 
## Installation ##


### devtools ###


From an interactive R session:

```{r, eval=FALSE}
library(devtools)
install_github("jennylsl/fast.adonis")
library(fast.adonis)

```
