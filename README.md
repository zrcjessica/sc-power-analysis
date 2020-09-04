# sc-power-analysis
A snakemake pipeline for scRNA-seq power analyses. 

# Introduction
This pipeline estimates power to detect differentially expressed genes from scRNA-seq data at various numbers of samples/cells and effect sizes. It models mean and dispersion parameters on a Negative Binomial distribution from reference data and simulates representative UMI count data from these estimates. The user can define various simulation parameters, such as the number of differentially expressed genes, effect sizes, number of samples, and number of cells. 

## Parameter estimation on reference data
Given a reference count matrix, we fit the observed counts for each gene to a Negative Binomial distribution and estimates the mean and dispersion parameters. Our likelihood function includes a scaling factor for the mean equivalent to the sequencing depth of each cell. 

## Data simulation
For each replicate, we generate a simulated count matrix reflecting a total number of genes and cells specified in the inputs. We simulate data for a control condition and a treatment condition. To begin, we randomly sample from the sequencing depths observed in the reference data to generate a vector of scaling factors equivalent in size to the number of cells specified for the simulatoin.
For each gene, we randomy sample a mean and dispersion parameter (with replacement) from the estimated parameters in the previous step and randomly sample from the Negative Binomial distribution 

![workflow](workflow.png)

# Download
Clone the repository to your desired destination.
```Shell
git clone https://github.com/zrcjessica/sc-power-analysis.git
```
# Dependencies
This pipeline was written with snakemake version 5.20.1. We recommend that you follow these instructions for [snakemake installation via conda](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda). This will install snakemake in an isolated environment, which is the recommended best practice for avoiding package conflicts. This allows you to [use the `--use-conda` flag](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management) when running the pipeline, which will automatically install the package specified in the `environment.yml` file. If you do not install snakemake in an isolated environment, then make sure you have installed all of the conda packages listed in `environment.yml`. 

# Execution

All input parameters are defined in the `config.yml` file and must be edited accordingly to suit your own data. 

## Inputs

### `data`
The `data` variable in the `config.yml` file contains the names of all samples for which you have count data available. Each sample name also stores the paths to the count matrix files under `counts` variable. 

This pipeline is designed to take UMI count matrices in compressed Matrix Market format (`*.mtx.gz`). We suggest using the `matrix.mtx.gz` file associated with the [filtered feature-barcode matrix](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices) from [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count).

Example:
```YAML
data:
  SampleA:
    counts: data/SampleA/outs/filtered_feature_bc_matrix/matrix.mtx.gz
  Sample B:
    counts: data/Sample B/outs/filtered_feature_bc_matrix/matrix.mtx.gz
  Sample C:
    counts: data/Sample C/outs/filtered_feature_bc_matrix/matrix.mtx.gz
```

### `sample`
Include the name(s) of the sample you want to use as reference panel(s) for your power analyses. 
