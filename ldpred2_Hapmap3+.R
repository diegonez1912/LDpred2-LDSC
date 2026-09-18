
## LDpred2-auto Script with HapMap3+ ## 

library(bigsnpr)
library(bigstatsr)
library(ggplot2)
library(dplyr)
library(data.table)
library(Matrix)

# Prevent nested parallelism
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.setenv(BLAS_NUM_THREADS = 1)
NCORES <- 22 # One per each chromosome


# FIRST STEP --> DOWNLOAD DATA AND SUMMARY STATISTICS

# 1.1 Load HapMap3+ SNP map:
SNP_map <- readRDS("/ldpred2/map_hm3_plus.rds")
SNP_map$pos_initial <- SNP_map$pos
SNP_map$pos <- SNP_map$pos_hg38
str(SNP_map)

# 1.2 Read external summary statistics --> A table with one row per SNP, and all of the characteristics of the SNP as columns.
sumstats <- bigreadr::fread2("/ldpred2/GCST90565439.ldpred2.tsv")
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
ldref_dir <- "/ldpred2/hapmap3_plus_ldref"

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

# 5.4 Plot the path of the chains --> In this case, we only plot the path of the first chain:
auto <- multi_auto[[1]]
plot_grid(
  qplot(y = auto$path_p_est) + 
    theme_bigstatsr() + 
    geom_hline(yintercept = auto$p_est, col = "blue") +
    scale_y_log10() +
    labs(y = "p"),
  qplot(y = auto$path_h2_est) + 
    theme_bigstatsr() + 
    geom_hline(yintercept = auto$h2_est, col = "blue") +
    labs(y = "h2"),
  ncol = 1, align = "hv"
) -> p_path

ggsave("ldpred2_chainpath_plot.png", plot = p_path, width = 8, height = 6, dpi = 300)

# 5.5 Filtering bad chains --> It is a QC measure based on the scale of a chain’s corr_est vector:
	# `range` should be between 0 and 2.
	# diff(range(auto$corr_est)) = max(corr_est) − min(corr_est) for that chain.
	# Chains that have a very small spread or with a extreme scale are defined as outliers. 
	# We only keep the chains whose range is greater than the 95% of all chain ranges (range close to the top of the whole distribution). 
range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
keep <- which(range > (0.8 * quantile(range, 0.8, na.rm = TRUE)))
cat("Keeping", length(keep), "out of", length(multi_auto), "chains\n")

# 5.6 Combine the chains that passed the filtering and produce final predictions for posterior effect sizes:
	# beta_auto = Final combined vector of SNP weights (averaged across good chains, and one per SNP).
		# rowMeans() --> It makes the combined vector be an average across important chains (in this case, rows).
	# Here we can find per-allele effect sizes ("weights") estimated from GWAS summary statistics. 
beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))


# SIXTH STEP --> Inference with LDpred2-auto

# We estimate unknown parameters and their uncertainty from the data (GWAS summary statistics + LD reference).
# We do this by using the sampler built into snp_ldpred2_auto, and we produce posterior samples for:
	# h2 --> SNP heritability (how much variance is attributable to all SNPs).
	# p --> Fraction (proportion) of SNPs that are causal.
	# α --> It controls whether rare or common SNPs tend to have bigger effects.
	# per-SNP effect vectors (many sampled beta draws).
# From those posterior samples we:
	# Summarize each parameter (median, 95% interval).
	# Compute derived quantities (predicted R², posterior inclusion probabilities per SNP).
	# Decide which chains / samples are trustworty.

# 6.1 Set up before the run:
        # coef_shrink --> It is is a regulation coefficient used to reduce the impact of the LD signals.
	# It gets lower (up to 0.4) when the LD reference doesn’t perfectly match your GWAS data.
coef_shrink <- 0.95
set.seed(1)

# 6.2 Running many chains --> Is the same process and code as in step 6.2, but with some minor changes:
	# burn_in = 500 --> Initial iterations to discard (warm-up).
	# num_iter = 500 --> Number of iterations kept after burn-in (we collect 500 post-burn samples per chain).
	# report_step = 20 --> Progress printout frequency.
multi_auto_inf <- snp_ldpred2_auto(
  corr, df_beta, h2_init = ldsc_h2_est,
  vec_p_init = seq_log(1e-4, 0.2, length.out = 50), ncores = NCORES,
  burn_in = 500, num_iter = 500, report_step = 20,
  allow_jump_sign = FALSE, shrink_corr = coef_shrink)

# 6.3 Filtering bad chains --> Is the same process and code as in step 5.5:
range_inf <- sapply(multi_auto_inf, function(auto) diff(range(auto$corr_est)))
keep_inf <- which(range_inf > (0.8 * quantile(range, 0.8, na.rm = TRUE)))

# 6.4 Calculation of unknown parameters --> We follow the same procedure for h2, p and alpha:
	# For each kept chain we extract the last num_iter entries of the h2 trace (these are the post-burn-in samples). 
	# sapply() --> It returns a matrix (columns = chains, rows = samples).
	# quantile(...) --> It coerces values to a single vector (all samples from all kept chains) and computes median and CI.

all_h2 <- sapply(multi_auto_inf[keep_inf], function(auto) tail(auto$path_h2_est, 500))
h2_quantiles <- quantile(all_h2, c(0.5, 0.025, 0.975))
cat("h2 posterior: median =", h2_quantiles[1], 
    "95% CI = [", h2_quantiles[2], ",", h2_quantiles[3], "]\n")

all_p <- sapply(multi_auto_inf[keep_inf], function(auto) tail(auto$path_p_est, 500))
p_quantiles <- quantile(all_p, c(0.5, 0.025, 0.975))
cat("p posterior: median =", p_quantiles[1], 
    "95% CI = [", p_quantiles[2], ",", p_quantiles[3], "]\n")

all_alpha <- sapply(multi_auto_inf[keep_inf], function(auto) tail(auto$path_alpha_est, 500))
alpha_quantiles <- quantile(all_alpha, c(0.5, 0.025, 0.975))
cat("alpha posterior: median =", alpha_quantiles[1], 
    "95% CI = [", alpha_quantiles[2], ",", alpha_quantiles[3], "]\n")

# Cleanup
file.remove(paste0(tmp, ".sbk"))

cat("Analysis completed successfully!\n")

