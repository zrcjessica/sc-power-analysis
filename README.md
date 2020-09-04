# sc-power-analysis
A Snakemake pipeline for scRNA-seq power analyses. 

# Introduction
This pipeline estimates power to detect differentially expressed genes from scRNA-seq data at various numbers of samples/cells and effect sizes. It models mean and dispersion parameters on a Negative Binomial distribution from reference data and simulates representative UMI count data from these estimates. The user can define various simulation parameters, such as the number of differentially expressed genes, effect sizes, number of samples, and number of cells. 

# Download
Clone the repository to your desired destination.
```Shell
https://github.com/zrcjessica/sc-power-analysis.git
```
