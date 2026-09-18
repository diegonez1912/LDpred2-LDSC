
## LDpred2-auto Script with HapMap3+ ## 

library(bigsnpr)
library(bigstatsr)
library(ggplot2)
library(dplyr)
library(data.table)
library(Matrix)
library(cowplot) # for combining plots

# Prevent nested parallelism
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.setenv(BLAS_NUM_THREADS = 1)
NCORES <- 22 # One per each chromosome


# FIRST STEP --> DOWNLOAD DATA AND SUMMARY STATISTICS

# 1.1 Load HapMap3+ SNP map:
SNP_map <- readRDS("/home/diegoonez/project/ldpred2/map_hm3_plus.rds")
SNP_map$pos_initial <- SNP_map$pos
SNP_map$pos <- SNP_map$pos_hg38
str(SNP_map)

# 1.2 Read external summary statistics --> A table with one row per SNP, and all of the characteristics of the SNP as columns.
sumstats <- bigreadr::fread2("/home/diegoonez/project/ldpred2/pap2025/GCST90565439.ldpred2.tsv")
sumstats$pos <- sumstats$base_pair_location
str(sumstats)


# SECOND STEP --> MATCH VARIANTS BETWEEN HAPMAP AND SUMMARY STATISTICS:

# The genotype data provides SNP positions (pos, chr) and alleles (a0 and a1).
# GWAS summary statistics provides effect sizes (beta), standard errors (beta_se) and sample sizes (n_eff).
# To apply GWAS effect sizes to your genotypes, you must make sure the SNPs are the same and the effect allele matches.
# These variables are used to match variants between the two data frames (we match by rsIDs instead of positions).
# snp_match() --> We get a cleaned version of summary statistics that are aligned with SNPs that are present in both datasets.
	# Variants to be matched --> Number of rows in summary statistics (number of SNPs).
	# Ambiguous SNPs --> SNPs where the alleles could pair between them (A/T and C/G).
	# Variants matched --> Number of succesful matches.
	# Variants flipped --> Same allels but on different strands in each dataset (after flipping, they are the same allele).
	# Variants reversed --> Same alleles but in the reverse order (e.g. summary has A1=A, A0=G while map has A1=G, A0=A).

# 2.1 Make sure that the required columns exist:
# Keep only SNPs present in HapMap3+
sumstats$n_eff <- as.numeric(sumstats$N)
sumstats$chr <- as.integer(sumstats$chr)
sumstats$pos <- as.integer(sumstats$pos)

# 2.2 Match GWAS summary statistics with SNPs from HapMap3+:
df_beta <- snp_match(sumstats, SNP_map, join_by_pos = TRUE)
cat("Number of SNPs:", nrow(df_beta), "\n")
str(df_beta)

# 2.3 Make sure that effect sizes, SEs and n_eff are numerical variables:
df_beta$beta <- as.numeric(df_beta$beta)
df_beta$beta_se <- as.numeric(df_beta$beta_se)
df_beta$n_eff <- as.numeric(df_beta$n_eff)


# THIRD STEP --> COMPUTE SFBM CORRELATION MATRIX:

# In this model, we run all SNPs together across all chromosomes, as LD (linkage disequilibrium) is a genome-wide phenomenon.
	# Thus, we build 22 giant correlation (LD) matrices, one per each chromosome.
	# They are stored in a Sparse Filebacked Big Matrix (SFBM) → a disk-based, memory-efficient representation of correlations.
	# Take into account that, with HapMap3+, we have a LD reference that has already been filtered, QC has been applied and LD values have been precomputed.
	# Thus, we have all the results (LD values) in ".rds" sparse matrices (SFBM).
	# Thus, we just need to load the files.

# 3.1 Define the path to the directory where we can find the LD references:
ldref_dir <- "/home/diegoonez/project/ldpred2/hapmap3_plus_ldref"

# 3.2 Create temporary file for SFBM:
tmp <- tempfile(tmpdir = "tmp-data")

# 3.3 Load correlation matrices chromosome by chromosome and build SFBM:
	# We get the indices for each chromosome in df_beta.
	# We load the correlation matrix for each chromosome.
	# We get the indices in the correlation matrix of each chromosome that correspond to the indices of SNPs in df_beta.
	# We subset the correlation matrix so that it only includes our SNPs.
	# We finally add the correlation matrix of each chromosome to SFBM, appending columns and rows in each chr.
for (chr in 1:22) {
	
	cat("Processing chromosome", chr, "...\n")
	ind.chr <- which(df_beta$chr == chr)
	if (length(ind.chr) == 0) {
		cat("No SNPs for chromosome", chr, "\n")
		next
	}
	
	corr_file <- file.path(ldref_dir, paste0("LD_with_blocks_chr", chr, ".rds"))
	if (!file.exists(corr_file)) {
		cat("File not found:", corr_file, "\n")
		next
	}
	corr_chr <- readRDS(corr_file)
	
	snps_chr_full <- which(SNP_map$chr == chr)
	ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
	ind.chr3 <- match(ind.chr2, snps_chr_full)
	
	corr_subset <- corr_chr[ind.chr3, ind.chr3]
	
	if (chr == 1) {
		corr <- as_SFBM(corr_subset, tmp, compact = TRUE)
	} else {
		corr$add_columns(corr_subset, nrow(corr))
	}
}

# 3.4 Check the structure of "corr":
class(corr)
str(corr)


# FOURTH STEP --> RUN LDpred2-inf TO ESTIMATE THE HERITABILITY (h2) OF THOSE SNPs:

# The idea behind it is to use the LD scores from the correlation, and we assume that all SNPs are causal (p=1).
# Then, we compute the posterior effect sizes under the infinitesimal prior and multiply them by the genotype matrix.
# By doing this, we get Polygenic Risk Scores (PRS) for the test set.

# 4.1 SNP-heritability estimation from LD score regression
# Use precomputed LD scores from SNP_map instead of computing them
ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(SNP_map),
								chi2 = (beta / beta_se)^2, 
								sample_size = n_eff, 
								ncores = NCORES))

ldsc_h2_est <- ldsc[["h2"]]
cat("Estimated heritability (h2):", ldsc_h2_est, "\n")

# 4.8 Run LDpred2-inf with the filtered data to get posterior effect sizes for each important SNP:
	# corr --> Genome-wide SFBM correlation matrix (the one with LD values).
	# df_beta --> It contains summary statistics for the important SNPs.
	# h2 = ldsc_h2_est --> Supplies the SNP heritability estimate for the infinitesimal prior.
	# Here we can find per-allele effect sizes ("weights") estimated from GWAS summary statistics. 
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = ldsc_h2_est)


# FIFTH STEP --> LDpred2-auto

# LDpred2-auto does not need a validation set as it directly infers values for h2 and p.
# It does this by running multiple MCMC/Gibbs chains from different starting points and then combining the good chains.

# 5.1 Set up before the run:
	# coef_shrink --> It is is a regulation coefficient used to reduce the impact of the LD signals.
	# It gets lower (up to 0.4) when the LD reference doesn’t perfectly match your GWAS data.
coef_shrink <- 0.95
set.seed(1)

# 5.2 Running many automatic chains --> Each chain runs an MCMC-like procedure that alternates between sampling SNP effect sizes and updating hyper-parameters (h2,p,...).
	# corr and df_beta --> LD reference and the matched summary-stats we already prepared.
	# h2_init --> Starting value for SNP heritability (here the LD Score estimate, step 4.1).
	# vec_p_init --> Vector of initial values for the fraction of causal variants p.
	# allow_jump_sign = FALSE --> Prevents abrupt sign flips in effect estimates during sampling, which improves robustness.
	# shrink_corr = coef_shrink --> Uses the shrinkage/regularization described above.
	# The function returns a list of chain results, and each chain might have hundreds or thousands of iterations. 
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = ldsc_h2_est,
								vec_p_init = seq_log(1e-4, 0.2, length.out = 50), ncores = NCORES,
								allow_jump_sign = FALSE, shrink_corr = coef_shrink)

# 5.3 Quick check of what is inside each chain:
	# beta_est --> Final estimated SNP effects/weights for that chain (length = #SNPs).
		# This defines the "weight" that each SNP has in creating a specific phenotype. 
		# This estimate is the per-allele effect --> how much the phenotype changes when you add one more copy of the allele.
		# 0 = no alleles --> effect = beta_est·0 // 2 = 2 alleles (max) --> effect = beta_est·2
	# postp_est --> Posterior inclusion probabilities (probability each SNP is causal).
	# corr_est --> Effect estimates expressed in the correlation scale used internally by LDpred2 (useful for QC).
	# sample_beta --> Sparse matrix of sampled betas during the run (the MCMC samples).
	# path_p_est --> Trajectory of p across iterations (how the chain moved from initial p to its final estimate).
	# path_h2_est --> Trajectory of h2 across iterations.
	# path_alpha_est --> Trajectory of alpha across iterations (it models how effect-size variance depends on allele frequency).
	# h2_est, p_est, alpha_est --> Final estimates for hyper-parameters.
	# h2_init, p_init --> Initial values used by that chain.
str(multi_auto[[1]], max.level = 1)

# 5.4 Filtering bad chains --> It is a QC measure based on the scale of a chain’s corr_est vector:
	# `range` should be between 0 and 2.
	# diff(range(auto$corr_est)) = max(corr_est) − min(corr_est) for that chain.
	# Chains that have a very small spread or with a extreme scale are defined as outliers. 
	# We only keep the chains whose range is greater than the 95% of all chain ranges (range close to the top of the whole distribution). 
range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
keep <- which(range > (0.8 * quantile(range, 0.8, na.rm = TRUE)))
cat("Keeping", length(keep), "out of", length(multi_auto), "chains\n")

# 5.5 Combine the chains that passed the filtering and produce final predictions for posterior effect sizes:
	# beta_auto = Final combined vector of SNP weights (averaged across good chains, and one per SNP).
		# rowMeans() --> It makes the combined vector be an average across important chains (in this case, rows).
	# Here we can find per-allele effect sizes ("weights") estimated from GWAS summary statistics. 
beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))


# SIXTH STEP --> Plot paths of selected chains:

chains_to_plot <- c(1, 10, 20, 30, 40, 50)

get_chain_paths <- function(multi_auto, chains) {
  df_list <- lapply(chains, function(i) {
    auto <- multi_auto[[i]]
    data.frame(
      Iteration = seq_along(auto$path_p_est),
      p = auto$path_p_est,
      h2 = auto$path_h2_est,
      Chain = paste0("Chain_", i)
    )
  })
  do.call(rbind, df_list)
}

paths_df <- get_chain_paths(multi_auto, chains_to_plot)

p_p <- ggplot(paths_df, aes(x = Iteration, y = p, color = Chain)) +
  geom_line(size = 0.8) +
  scale_y_log10() +
  labs(y = "p (causal variants)", x = "Iteration") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        panel.grid.major = element_line(color = "gray80"),
        panel.grid.minor = element_line(color = "gray90"))

p_h2 <- ggplot(paths_df, aes(x = Iteration, y = h2, color = Chain)) +
  geom_line(size = 0.8) +
  labs(y = "h2 (SNP heritability)", x = "Iteration") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        panel.grid.major = element_line(color = "gray80"),
        panel.grid.minor = element_line(color = "gray90"))

combined_plot <- plot_grid(p_p, p_h2, ncol = 1, align = 'v')
ggsave("ldpred2_selected_chain_paths.png", combined_plot, width = 10, height = 8, dpi = 300)  


# SEVENTH STEP --> Select Top 10 SNPs by effect size:

top_snps <- df_beta %>%
  mutate(effect_size = beta_auto) %>%
  arrange(desc(abs(effect_size))) %>%
  slice_head(n = 10) %>%
  select(rsid, chr, pos, a0, a1, effect_size)

print(top_snps)
write.csv(top_snps, "top10_snps_beta_auto.csv", row.names = FALSE)

# Cleanup
file.remove(paste0(tmp, ".sbk"))

cat("Analysis completed successfully!\n")

