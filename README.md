exp2flux : Convert gene expression data to FBA flux boundaries.
======
The computational simulations of metabolism through constraint-based modeling approaches can help to predict the metabolic phenotype of an organism in response to different stimuli or treatments. To recreate specific metabolic phenotypes and enhance the model predictive accuracy, several methods for the integration of transcriptomics data into constraint-based models have been proposed. The majority of available implemented methods are based on the discretization of data to incorporate constraints into the metabolic models via boolean logic representation, which reduces the accurate of physiological representations. The implemented methods for gene-expression data integration as continuous values are very few. The **exp2flux** package was designed as a tool to incorporate in a continuous way the gene-expression data as FBA flux boundaries in a metabolic model. Also, a function to calculate the differences between model fluxes in different metabolic scenarios was included. This is an implementation of the algorithm described by Lee *et al.* (2012) in  DOI: 10.1186/1752-0509-6-73

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
