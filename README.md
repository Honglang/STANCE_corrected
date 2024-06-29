# STANCE: Spatial Transcriptomics ANalysis of genes with Cell-type-specific Expression
STANCE is a unified statistical model to detect cell type-specific spatially variable genes (ctSVGs) in spatial transcriptomics.

## Installation
Please run the following codes in R to install STANCE package from GitHub.
```
if (!require("devtools", quietly = TRUE)){
  install.packages("devtools")
}
devtools::install_github("Cui-STT-Lab/STANCE")
```

STANCE relies on some functions in [SPARK](https://xzhoulab.github.io/SPARK/), which can be installed by running the following codes in R:
```
if (!require("devtools", quietly = TRUE)){
  install.packages("devtools")
}
devtools::install_github("xzhoulab/SPARK")
library(SPARK)
```

## Detect cell type-specific spatially variable genes (ctSVGs) with STANCE
The best vignette for getting started with STANCE is [Tutorial](https://github.com/Cui-STT-Lab/STANCE/blob/master/vignettes/tutorial.html).

