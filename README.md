exp2flux : convert gene EXPression data to FBA FLUXes
======
The **exp2flux** package was designed as a tool to convert EXPression data to FBA FLUXes.
Install:
--------
This package required R version 2.10 or higher. If you are using an older version of R you will be prompted to upgrade when you try to install the package.

For install the latest stable version this package directly from GitHub:
```
# Install 'devtools' R package
install.packages("devtools")

# Install 'minval' package
devtools::install_github("gibbslab/exp2flux")
library("exp2flux")
```

Available functions:
-------------------
|Function | Description |
|:--------|:------------|
|exp2flux|Convert Gene Expression Data to FBA fluxes|
|fluxDifferences|Report the fold change of fluxes between two models|

Citation
--------
Daniel Osorio, Kelly Botero, Janneth Gonzalez and Andres Pinzon-Velasco (2016). **exp2flux: convert EXPression data to FBA FLUXes**. R package version 0.1.
