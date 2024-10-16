# STANCE: Spatial Transcriptomics ANalysis of genes with Cell-type-specific Expression
STANCE is a unified statistical model to detect cell type-specific spatially variable genes (ctSVGs) in spatial transcriptomics. By integrating gene expression, spatial location, and cell type composition through a linear mixed-effect model, STANCE enables the identification of both SVGs and ctSVGs in an initial stage, followed by a second stage test dedicated to ctSVG detection. Its design ensures robustness in complex scenarios and the results are spatial rotation invariant.

## Installation
Please run the following codes in R to install STANCE package from GitHub.
```
if (!require("devtools", quietly = TRUE)){
  install.packages("devtools")
}
devtools::install_github("Cui-STT-Lab/STANCE")
```

STANCE relies on several packages, including gaston, KRLS, [SPARK](https://xzhoulab.github.io/SPARK/), among others, most of which are installed automatically. If automatic installation fails, you can manually install them in R by running the following codes:
```
install.packages("gaston")
install.packages("KRLS")
if (!require("devtools", quietly = TRUE)){
  install.packages("devtools")
}
devtools::install_github("xzhoulab/SPARK")
library(SPARK)
```

## Detect cell type-specific spatially variable genes (ctSVGs) with STANCE
The best vignette for getting started with STANCE is [Tutorial](https://haroldsu.github.io/STANCE/tutorial.html).

