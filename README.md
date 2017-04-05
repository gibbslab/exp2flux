exp2flux : Convert gene expression data to FBA flux boundaries.
======

The **exp2flux** package was designed as a tool to convert expression data to FBA flux boundaries.

Install:
--------
This package required R version 2.10 or higher. If you are using an older version of R you will be prompted to upgrade when you try to install the package.

For install the latest stable version this package directly from GitHub:
```{r}
# Install 'minval' package
devtools::install_github("gibbslab/exp2flux")
library("exp2flux")
```

Available functions:
-------------------
|Function | Description |
|:--------|:------------|
|exp2flux|Convert gene expression data to FBA flux boundaries|
|getDifferences||
|getFluxes||
|plotDifferences||

Available datasets:
-------------------
|Function | Description |
|:--------|:------------|
|antibioticsExpression||
|mutantsExpression||

Citation
--------
Daniel Osorio, Kelly Botero, Janneth Gonzalez and Andres Pinzon (2017). **exp2flux: convert EXPression data to FBA FLUXes**. R package version 0.2.
