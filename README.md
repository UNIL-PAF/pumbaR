# pumbaR
MaxQuant provides us with a file *proteinGroups.txt* which contains information about the theoretical masses of proteins and how abbundant peptides from those proteins are found in slices from one dimensional gels. This information can be used to compute a approximation of which masses correspond to each slice.

## Installation

Install *pumbaR* from source in R using *devtools*:

```R
library(devtools)
install_github("UNIL-PAF/pumbaR", build_vignettes = TRUE)
```

## Usage

Look at the vignettes for more information about usage:

```R
browseVignettes(package="pumbaR")
```

