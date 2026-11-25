# LDpred2-LDSC

**R-based statistical genetics workflows for analysing GWAS summary statistics, estimating SNP heritability and polygenic architecture, and investigating genetic correlations across complex neurological disorders.**

## Overview

Genome-wide association studies (GWAS) provide large-scale information about genetic variants associated with complex diseases. However, understanding the underlying genetic architecture requires statistical methods that account for linkage disequilibrium (LD), polygenicity and shared genetic effects across traits.

This repository contains computational workflows implementing LDpred2, LD Score regression (LDSC) and GenomicSEM to investigate the genetic architecture and shared genetic components of complex neurological disorders.

The repository comprises two complementary analytical approaches:

1. **LDpred2:** Bayesian modelling of polygenic architecture and estimation of posterior SNP effect sizes using GWAS summary statistics and HapMap3+ LD reference data.
2. **Multivariate LDSC / GenomicSEM:** Estimation of SNP heritability, genetic covariance and genetic correlations across neurological disorders using harmonised GWAS summary statistics.

The workflows were developed and applied to large-scale neurological disease GWAS datasets, including Alzheimer's disease, amyotrophic lateral sclerosis (ALS), Lewy body dementia (LBD) and Parkinson's disease (PD).

## Repository structure

```text
LDpred2-LDSC/
|
|-- ldpred2_Hapmap3+.R
|-- ldpred2_Hapmap3+_paths.R
|-- gsem_LDSC.R
|-- README.md
```

## Workflow 1: LDpred2 – Bayesian polygenic modelling

**Scripts:** `ldpred2_Hapmap3+.R`, `ldpred2_Hapmap3+_paths.R`

### Objective

Estimate SNP heritability, polygenic architecture and posterior SNP effect sizes from GWAS summary statistics using LDpred2.

The workflow integrates HapMap3+ reference data and precomputed chromosome-specific LD correlation matrices to perform Bayesian modelling across genome-wide genetic variants.

### Analytical workflow

**1. GWAS summary-statistics preparation and harmonisation**

- Imports GWAS summary statistics and HapMap3+ reference SNP information.
- Matches genomic variants between the GWAS dataset and the LD reference using chromosome positions and allele information.
- Filters and harmonises variants to ensure compatible effect alleles, effect estimates, standard errors and sample sizes.
- Retains the variants available in the reference dataset for downstream analysis.

**2. Linkage disequilibrium matrix construction**

- Loads precomputed chromosome-specific LD correlation matrices.
- Extracts the relevant SNP subsets based on the harmonised GWAS variants.
- Combines chromosome-specific matrices into a genome-wide sparse file-backed big matrix (SFBM).
- Uses disk-backed matrix storage to handle large LD reference datasets more efficiently.

**3. SNP heritability estimation and LDpred2-inf**

- Applies LD Score regression through `snp_ldsc()` to estimate SNP heritability.
- Uses the estimated heritability and genome-wide LD information to calculate posterior SNP effects under the infinitesimal model.

**4. LDpred2-auto**

- Implements Bayesian estimation of genome-wide SNP effects and polygenic architecture.
- Runs multiple MCMC/Gibbs chains with 50 alternative initial values for the proportion of causal variants.
- Evaluates chain behaviour using diagnostic quantities derived from the sampling process.
- Filters chains according to the distribution of their correlation-scale estimates.
- Combines retained chains to derive posterior SNP effect estimates.

**5. Posterior inference and diagnostics**

The main implementation includes an additional LDpred2-auto inference stage using 500 burn-in iterations and 500 retained iterations.

It extracts posterior estimates for:
- SNP heritability (h²).
- Proportion of causal variants (p).
- Alpha (α), describing the relationship between allele frequency and effect-size variance.

Posterior medians and 95% credible intervals are calculated from the retained chains.

### Technologies and implementation

- R
- bigsnpr
- bigstatsr
- bigreadr
- data.table
- dplyr
- ggplot2
- Matrix

### Outputs

The main implementation generates:

- Chain convergence and parameter-trajectory visualisation (`ldpred2_chainpath_plot.png`).
- Posterior summaries of SNP heritability, polygenicity and alpha, printed during execution.
- Posterior SNP effect estimates available within the R session.

The alternative `_paths.R` implementation additionally:

- Visualises the trajectories of six selected MCMC/Gibbs chains.
- Identifies and exports the ten variants with the largest absolute estimated effect sizes.
- Generates `ldpred2_selected_chain_paths.png`.
- Exports `top10_snps_beta_auto.csv`.

The final SNP effect vectors are not currently saved as a dedicated output file by the main script.

## Workflow 2: Multivariate LDSC and GenomicSEM

**Script:** `gsem_LDSC.R`

### Objective

Investigate the shared genetic architecture of neurological disorders by estimating SNP heritability, genetic covariance and genetic correlations across multiple GWAS datasets.

The workflow uses multivariate LD Score regression through the GenomicSEM framework.

### Analytical workflow

**1. GWAS preprocessing and variant harmonisation**

- Imports GWAS summary statistics from Alzheimer's disease, ALS, LBD and PD studies.
- Processes heterogeneous summary-statistics formats.
- Matches genomic variants against a common reference dataset using chromosome, position and allele information.
- Annotates compatible variants with reference SNP identifiers and standardises summary-statistics columns.

**2. Quality control and summary-statistics harmonisation**

- Applies GenomicSEM's `munge()` procedure to standardise and QC GWAS summary statistics.
- Uses a common HapMap3 SNP reference.
- Applies information-score (INFO) filtering and minor allele frequency (MAF) filtering.
- Prepares harmonised summary statistics for multivariate LD Score regression.

The current implementation specifies INFO filtering at 0.9 and MAF filtering at 0.01.

**3. Multivariate LD Score regression**

- Applies `ldsc()` from GenomicSEM to multiple harmonised GWAS datasets.
- Uses chromosome-specific European LD scores and regression-weight references.
- Estimates SNP heritability for each trait.
- Calculates genetic covariance between pairs of neurological disorders.
- Derives the genetic correlation matrix from the estimated covariance matrix.

**4. Statistical result extraction**

- Extracts the genetic covariance matrix and its sampling covariance matrix.
- Calculates SNP heritability and associated standard errors and Z-statistics.
- Estimates pairwise genetic correlations.
- Structures the results into tables for downstream statistical analysis and interpretation.

### Technologies and implementation

- R
- GenomicSEM
- LDSC
- bigsnpr
- bigstatsr
- bigreadr
- data.table
- dplyr
- tidyverse
- Matrix
- ggplot2
- reshape2
- knitr

### Outputs

The workflow generates:

| Output file | Description |
|---|---|
| `LDSC_heritabilities.csv` | SNP heritability estimates, standard errors and Z-statistics for each neurological disorder |
| `LDSC_genetic_correlations.csv` | Pairwise genetic correlation estimates between neurological disorders |

## Input datasets and reference resources

The analyses require external GWAS datasets and reference files, which are not included in this repository.

### LDpred2 workflow

Required resources include:

- GWAS summary statistics containing variant identifiers, chromosome positions, alleles, effect estimates, standard errors and sample-size information.
- HapMap3+ SNP map (`map_hm3_plus.rds`).
- Precomputed chromosome-specific LD correlation matrices (`LD_with_blocks_chr*.rds`).

### Multivariate LDSC / GenomicSEM workflow

Required resources include:

- GWAS summary statistics for Alzheimer's disease, ALS, LBD and PD.
- A common variant reference dataset, including a 1000 Genomes Phase 3 reference `.bim` file.
- HapMap3 SNP reference list (`w_hm3.snplist`).
- Chromosome-specific European LD-score reference files.
- Corresponding LD regression-weight reference files.
- Appropriate sample and population prevalence estimates for the case-control traits.

**Important:** GWAS datasets must be harmonised consistently with respect to genome assembly, chromosome coordinates, effect alleles and reference ancestry before performing genetic correlation analyses.

## Requirements and computational environment

### Software

- R
- A compatible R package environment.
- Access to the required GWAS summary statistics and LD reference resources.

### R packages

The scripts use the following R packages:

| Package | Purpose |
|---|---|
| bigsnpr | SNP matching, LD handling and LDpred2 modelling |
| bigstatsr | Large-scale and file-backed matrix processing |
| bigreadr | Efficient reading of large datasets |
| data.table | Efficient tabular data processing |
| dplyr / tidyverse | Data manipulation and harmonisation |
| Matrix | Sparse matrix operations |
| ggplot2 | Statistical visualisation |
| GenomicSEM | Multivariate LD Score regression |
| reshape2 | Result transformation |
| knitr | Table formatting |
| cowplot | Plot composition in the alternative LDpred2 workflow |

The exact R and package versions have not yet been pinned in this repository.

## Running the workflows

**Configuration note:** The current scripts contain dataset-specific file paths and reference directories. Update these paths and confirm that all external input files are available before execution. The GenomicSEM script also contains trait-specific sample and population prevalence values that must be reviewed when changing datasets.

### 1. LDpred2 workflow

Configure the paths to:

- GWAS summary-statistics files.
- HapMap3+ SNP map.
- Chromosome-specific LD correlation matrices.
- Temporary and output directories.

Run the main implementation:

```bash
Rscript "ldpred2_Hapmap3+.R"
```

Alternatively, run the chain-diagnostic and SNP-ranking implementation:

```bash
Rscript "ldpred2_Hapmap3+_paths.R"
```

The scripts are designed to use parallel processing. The number of computational cores is currently defined as 22 and should be adjusted according to available resources.

### 2. Multivariate LDSC / GenomicSEM workflow

Configure the paths to the four GWAS summary-statistics datasets, common SNP reference, HapMap3 list and European LD-score resources.

Verify the trait-specific prevalence parameters before execution.

Run:

```bash
Rscript gsem_LDSC.R
```

The analysis exports SNP heritability and genetic correlation tables as CSV files.

## Quality control and analytical considerations

The workflows incorporate several analytical QC procedures.

**LDpred2:**
- Variant matching and allele harmonisation against HapMap3+.
- Construction of compatible chromosome-specific LD matrices.
- Evaluation of MCMC/Gibbs-chain behaviour.
- Filtering of unsuitable chains before combining posterior estimates.
- Visualisation of parameter trajectories.

**GenomicSEM / LDSC:**
- Variant matching against a common reference.
- Summary-statistics harmonisation.
- INFO and MAF filtering.
- Use of common LD-score reference resources.
- Extraction of genetic covariance and sampling covariance estimates.

The interpretation of estimates depends on the quality of GWAS summary statistics, reference compatibility, ancestry matching and appropriate prevalence assumptions.

## Limitations and reproducibility

This repository contains research-oriented statistical analysis scripts rather than a packaged software application.

Important considerations:

- Input datasets and large reference resources are not distributed with the repository.
- File paths are currently dataset- and environment-specific and require adaptation.
- The GenomicSEM implementation uses European LD-score reference resources.
- Sample and population prevalence parameters are specified for particular datasets and should not be transferred uncritically to other GWAS.
- The alternative LDpred2 scripts contain overlapping analytical steps and should be consolidated or documented carefully.
- Automated unit tests, integration tests and continuous integration workflows are not currently included.
- R package versions and a complete computational environment are not currently pinned.

The LDpred2 outputs represent statistical estimates of SNP effects and polygenic architecture. Independent validation is required before interpreting derived polygenic scores as clinically predictive models.

## References and useful resources

- Privé F. et al. (2020). *Efficient Bayesian mixed-model analysis increases association power in large cohorts.* Bioinformatics. [LDpred2 publication](https://doi.org/10.1093/bioinformatics/btaa1029).
- [LDpred2 implementation and documentation](https://github.com/privefl/bigsnpr).
- [GenomicSEM framework](https://github.com/GenomicSEM/GenomicSEM).
- [HapMap3+ LD reference dataset](https://figshare.com/articles/dataset/LD_reference_for_HapMap3_/21305061).

---

**Author:** Diego Núñez Martínez

**Repository:** [LDpred2-LDSC](https://github.com/diegonez1912/LDpred2-LDSC)
