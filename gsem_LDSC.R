
## gSEM HDL tutorial with HapMap3+ ## 

library(GenomicSEM)
library(data.table)
library(dplyr)
library(ggplot2)
library(Matrix)
library(bigsnpr)
library(bigstatsr)
library(tidyverse)
library(reshape2)
library(knitr)

# Prevent nested parallelism
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.setenv(BLAS_NUM_THREADS = 1)
NCORES <- 22 # One per each chromosome


# FIRST STEP --> DOWNLOAD SUMMARY STATISTICS AN PREPARE THEM:

# 1.1 Load summary statistics:
ad_initial <- bigreadr::fread2("/AD_datasets/PGCALZ2sumstatsExcluding23andMe.txt")
str(ad_initial)
als_initial <- bigreadr::fread2("AS_datasets/GCST90027164_buildGRCh37.tsv")
str(als_initial)
lbd_intitial <- bigreadr::fread2("/LBD_datasets/GCST90001390_buildGRCh38.tsv")
str(lbd_intitial)
pd_initial <- bigreadr::fread2("/PD_datasets/nallsEtAl2019_excluding23andMe_allVariants.tab")
str(pd_initial)

# 1.2 Identify rsIDs for each SNP based on chromosome, position, and alleles:
ref <- read_table(
  "/gSEM_HDL/1000G_phase3_common_norel.bim",
  col_names = c("CHR", "SNP", "CM", "BP", "A1", "A2")) %>%
  select(CHR, BP, A1, A2, SNP)

ad_sumstats <- ad_initial %>%
  rename(CHR = chr, BP = PosGRCh37, A1 = testedAllele, A2 = otherAllele)
ad_annotated <- ad_sumstats %>%
  left_join(ref , by = c("CHR", "BP", "A1", "A2")) %>%
  drop_na()
fwrite(ad_annotated, "sumstats_AD_annot.txt", sep = "\t")

pd_sumstats <- pd_initial %>%
  separate(SNP, into = c("CHR", "BP"), sep = ":", remove = FALSE) %>%
  mutate(
    CHR = as.integer(gsub("chr", "", CHR)),
    BP  = as.integer(BP),
    N   = N_cases + N_controls
  ) %>%
  rename(A1 = A1, A2 = A2, BETA = b, SE = se, P = p)
pd_annot <- pd_sumstats %>%
  left_join(ref , by = c("CHR", "BP", "A1", "A2")) %>%
  drop_na()
pd_annot_fixed <- pd_annot %>%
  rename(SNP = SNP.y) %>%
  select(SNP, CHR, BP, A1, A2, freq, BETA, SE, P, N)
fwrite(pd_annot_fixed, "sumstats_PD_annot.txt", sep = "\t")

# SECOND STEP --> RUN munge() TO PREPARE SUMMARY STATISTICS FOR gSEM HDL ANALYSIS:
# This function will harmonize the summary statistics with the SNP map, and perform quality control.
# It produces a standardized and QC'ed version of the summary statistics. 

# 2.1 Multiple-trait munge in ONE CALL
files_munge <- c("/sumstats_AD_annot.txt", 
                 "/ALS_datasets/GCST90027164_buildGRCh37.tsv",
                 "/LBD_datasets/LBD_GRCh37_ready.tsv",
                 "/PD_datasets/sumstats_PD_annot.txt")
hm3_file <- "/w_hm3.snplist"
trait.names <- c("AD", "ALS", "LBD", "PD")

# 2.2 Munge summary statistics:
munge(
  files = files_munge,
  hm3 = hm3_file,
  trait.names = trait.names,
  info.filter = 0.9,
  maf.filter = 0.01,
  parallel = TRUE,
  cores = NCORES,
  log.name = "/ldsc_APOE/munge"
  )    


# THIRD STEP --> RUN ldsc() TO GENETIC CORRELATIONS AND HERITABILITIES:
# This function will compute the multivariate LD Score Regression.
# It computes: SNP-heritability for each trait, genetic covariance between pairs, and genetic correlation matrix.
ldsc_trait <- ldsc(
  c("AD.sumstats.gz", "ALS.sumstats.gz", "LBD.sumstats.gz", "PD.sumstats.gz"), # input GWAS results
  sample.prev = c(0.158, 0.197, 0.392, 0.026), # For case-control GWAS --> proportion of cases in the sample
  population.prev = c(0.05, 0.0002, 0.01, 0.01), # For case control GWAS --> True prevalence in population
  trait.names = c("AD", "ALS", "LBD", "PD"),
  ld = "/gSEM_HDL/eur_w_ld_chr/", # Folder with LD scores for each SNP, split by chromosomes.
  wld = "/gSEM_HDL/eur_w_ld_chr/", # Folder with regression weights used by LDSC, split by chromosomes.
)
str(ldsc_trait)


### FOURTH STEP --> Extract results from GenomicSEM LDSC object

# S = genetic covariance matrix
S <- ldsc_trait$S

# V = sampling covariance matrix of S
V <- ldsc_trait$V

# 4.1 SNP heritability:
# Heritability is the diagonal of S
h2 <- diag(S)

# V is 10x10; extract variance of diagonal S elements
# indices of diagonal elements in vectorized S:
idx_h2 <- c(1, 3, 6, 10)

h2_se <- sqrt(diag(V)[idx_h2])

h2_table <- data.frame(
  Trait = colnames(S),
  h2 = h2,
  SE = h2_se,
  Z = h2 / h2_se
)

# 4.2 Genetic correlations:
# Correlation matrix
cor_mat <- cov2cor(S)

# Convert matrix to a long table
rg_table <- melt(cor_mat, varnames = c("Trait1", "Trait2"), value.name = "rg")

# Remove diagonal rows (rg = 1)
rg_table <- rg_table[rg_table$Trait1 != rg_table$Trait2, ]


# 4.3 Save results:
write.csv(h2_table, "LDSC_heritabilities.csv", row.names = FALSE)
write.csv(rg_table, "LDSC_genetic_correlations.csv", row.names = FALSE)


# 4.4 Nicely formatted tables:
kable(h2_table, caption = "LDSC SNP-Heritability (Covariance-Derived)")
kable(rg_table, caption = "LDSC Genetic Correlations (Covariance-Derived)")


cat("Analysis completed successfully!\n")
