
################################################################################
# 0. Set working directory
################################################################################

setwd("~/seg_epigenomics_meth/seg_data/methylation_analysis/DSS")

################################################################################
# 1. Install required packages
################################################################################

# Install BiocManager if it is not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required Bioconductor packages
bioc_packages <- c("DSS", "bsseq", "GenomicRanges")

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

# Install required CRAN packages
cran_packages <- c("dplyr", "ggplot2")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

################################################################################
# 2. Load libraries
################################################################################

library(DSS)
library(bsseq)
library(dplyr)
library(GenomicRanges)
library(ggplot2)

################################################################################
# 3. Read in methylation data
################################################################################

dat1.1 <- read.table("danio_4hpf_rep1.DSS.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat1.2 <- read.table("danio_4hpf_rep2.DSS.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat2.1 <- read.table("danio_36hpf_rep1.DSS.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat2.2 <- read.table("danio_36hpf_rep2.DSS.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

################################################################################
# 4. Create the BSseq object
################################################################################

BSobj <- makeBSseqData(
  list(dat1.1, dat1.2, dat2.1, dat2.2),
  c("4hpf_rep1", "4hpf_rep2", "36hpf_rep1", "36hpf_rep2")
)

# The input tables are no longer needed after creating the BSseq object
rm(dat1.1, dat1.2, dat2.1, dat2.2)
gc()

# Save the BSseq object so that the analysis does not have to be repeated
#saveRDS(BSobj, file = "BSobj.rds")

# Load a previously saved BSseq object
#BSobj <- readRDS("BSobj.rds")

# Inspect the object
#BSobj
#head(BSobj)

# BSobj is a BSseq object containing the methylation data
# from all four samples.
#
# It stores the genomic position of each CpG together with
# the methylated read counts and the total number of reads
# covering that position in each sample.
#
# Printing BSobj gives an overview of the object, including
# the number of CpG sites and the number of samples.
#
# head(BSobj) shows the first few genomic positions
# and the corresponding methylation data stored in the object.

################################################################################
# 5. Test for differential methylation
################################################################################

# A DML (Differentially Methylated Locus) is a CpG site that shows
# a statistically significant difference in methylation between two groups.
#
# In our case, we compare:
# Group 1: 4 hpf
# Group 2: 36 hpf

# DMLtest performs the statistical analysis at the CpG level.
# It estimates methylation levels and dispersion at each CpG site
# and then performs a statistical test to assess whether methylation
# differs between the two groups.

# Methylation levels at nearby CpG sites are often spatially correlated.
# Smoothing can use information from nearby CpG sites to improve the
# estimation of methylation levels.
#
# For WGBS data, CpG sites are relatively dense, so smoothing can be useful.
# For sparse data such as RRBS, smoothing may provide less benefit.

dmlTest.sm <- DMLtest(
  BSobj,
  group1=c("4hpf_rep1", "4hpf_rep2"),
  group2=c("36hpf_rep1", "36hpf_rep2"),
  smoothing=TRUE,
  ncores=6
)

# Inspect the results
head(dmlTest.sm)

# dmlTest.sm contains the statistical results from DMLtest
# for the comparison between 4 hpf and 36 hpf.
#
# Each row corresponds to a CpG site tested for differential methylation.
#
# The results include information about the estimated methylation
# difference between the two groups and the statistical evidence
# for this difference.
#
# At this stage, we are looking at the statistical results for the
# tested CpG sites. We have not yet selected individual significant DMLs.
#
# head(dmlTest.sm) shows the first few rows of the statistical results.
#
# These results will be used in the next step to identify
# differentially methylated regions (DMRs).

# Save the results so that the statistical test does not have to be repeated
# saveRDS(dmlTest.sm, file = "dmlTest_sm.rds")

# Load previously saved DMLtest results
# dmlTest.sm <- readRDS("dmlTest_sm.rds")

################################################################################
# 6. DMR detection using callDMR
################################################################################

# DMRs (Differentially Methylated Regions) are genomic regions containing
# multiple CpG sites that show a consistent difference in methylation
# between the two groups.

# callDMR identifies regions based on the results from DMLtest.

# Main parameters:
# delta: minimum difference in methylation between the two groups
# p.threshold: statistical significance threshold
# minlen: minimum length of the DMR
# dis.merge: maximum distance between nearby DMRs that can be merged
# pct.sig: minimum percentage of significant CpG sites within the region

dmrs <- callDMR(
  dmlTest.sm,
  delta=0.3,
  p.threshold=0.05,
  minlen=50,
  dis.merge=100,
  pct.sig=0.5
)

# Save the DMR results so that the analysis does not have to be repeated
#saveRDS(dmrs, file = "dmrs.rds")

# Load previously saved DMR results
#dmrs <- readRDS("dmrs.rds")

# Inspect the results
head(dmrs)
nrow(dmrs)

# Each row represents one DMR, i.e. a genomic region containing
# multiple CpG sites with a consistent difference in methylation
# between 4 hpf and 36 hpf.
#
# The table contains information about the genomic coordinates
# of each DMR, the number of CpG sites within the region,
# the methylation difference between the two groups,
# and the statistical significance of the region.
#
# head(dmrs) shows the first DMRs returned by callDMR().
# These DMRs are shown in the order in which they appear in the results,
# and not necessarily in genomic order.
#
# nrow(dmrs) shows the total number of DMRs identified.

################################################################################
# 7. Generate BED files for hypermethylated and hypomethylated DMRs
################################################################################

# Remove DMRs with missing methylation difference
dmrs <- dmrs[!is.na(dmrs$diff.Methy), ]

# DSS calculates diff.Methy as:
# mean methylation in group 1 - mean methylation in group 2.
#
# Here:
# group 1 = 4 hpf
# group 2 = 36 hpf
#
# Therefore:
# diff.Methy > 0: higher methylation at 4 hpf
#                 -> hypomethylated at 36 hpf
#
# diff.Methy < 0: lower methylation at 4 hpf
#                 -> hypermethylated at 36 hpf

dmrs$direction <- ifelse(dmrs$diff.Methy > 0,
                         "36hpf_hypomethylated",
                         "36hpf_hypermethylated")

make_DMR_beds <- function(dmrs, prefix){
  
  # Create a copy for BED output
  dmrs_bed <- dmrs
  
  # Convert coordinates to 0-based BED format
  dmrs_bed$start <- dmrs_bed$start - 1
  
  # BED: hypermethylated at 36 hpf compared with 4 hpf
  bed_hyper <- dmrs_bed[dmrs_bed$direction == "36hpf_hypermethylated",
                        c("chr", "start", "end")]
  
  write.table(
    bed_hyper,
    file = paste0(prefix, "_36hpf_hypermethylated_vs_4hpf.bed"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  # BED: hypomethylated at 36 hpf compared with 4 hpf
  bed_hypo <- dmrs_bed[dmrs_bed$direction == "36hpf_hypomethylated",
                       c("chr", "start", "end")]
  
  write.table(
    bed_hypo,
    file = paste0(prefix, "_36hpf_hypomethylated_vs_4hpf.bed"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

make_DMR_beds(dmrs, "DMRs")


################################################################################
# 8. Calculate methylation levels for each DMR
################################################################################

# Get methylated reads and total coverage
M <- getBSseq(BSobj, type="M")
Cov <- getBSseq(BSobj, type="Cov")

# Merge biological replicates within each condition
M_4hpf <- M[, "4hpf_rep1"] + M[, "4hpf_rep2"]
Cov_4hpf <- Cov[, "4hpf_rep1"] + Cov[, "4hpf_rep2"]

M_36hpf <- M[, "36hpf_rep1"] + M[, "36hpf_rep2"]
Cov_36hpf <- Cov[, "36hpf_rep1"] + Cov[, "36hpf_rep2"]

# Get CpG positions
cpg_gr <- rowRanges(BSobj)

# Calculate methylation for each DMR 
# 
# For each DMR, we sum the methylated reads and total coverage across all CpGs 
# within the region. 
# 
# DMR methylation = total methylated reads / total coverage
dmr_methylation <- lapply(seq_len(nrow(dmrs)), function(i) {
  
  idx <- as.character(seqnames(cpg_gr)) == dmrs$chr[i] &
    start(cpg_gr) >= dmrs$start[i] &
    start(cpg_gr) <= dmrs$end[i]
  
  data.frame(
    DMR = paste0("DMR_", i),
    direction = dmrs$direction[i],
    methylation_4hpf = sum(M_4hpf[idx], na.rm = TRUE) /
      sum(Cov_4hpf[idx], na.rm = TRUE),
    methylation_36hpf = sum(M_36hpf[idx], na.rm = TRUE) /
      sum(Cov_36hpf[idx], na.rm = TRUE)
  )
}) %>%
  bind_rows()
# Remove intermediate objects
rm(M, Cov, M_4hpf, Cov_4hpf, M_36hpf, Cov_36hpf, cpg_gr)
gc()


################################################################################
# 9. Prepare data for plotting
################################################################################

plot_df <- rbind(
  data.frame(
    DMR = dmr_methylation$DMR,
    direction = dmr_methylation$direction,
    condition = "4 hpf",
    methylation = dmr_methylation$methylation_4hpf
  ),
  data.frame(
    DMR = dmr_methylation$DMR,
    direction = dmr_methylation$direction,
    condition = "36 hpf",
    methylation = dmr_methylation$methylation_36hpf
  )
)

# Remove missing values
plot_df <- plot_df[!is.na(plot_df$methylation), ]


################################################################################
# 10. Plot methylation levels
################################################################################

ggplot(plot_df, aes(x=condition, y=methylation)) +
  
  geom_boxplot(
    width=0.6,
    outlier.shape=NA
  ) +
  
  geom_jitter(
    width=0.15,
    size=1,
    alpha=0.4
  ) +
  
  facet_wrap(~direction) +
  
  labs(
    x=NULL,
    y="DNA methylation (mCG/CG)"
  ) +
  
  theme_classic(base_size=12) +
  
  theme(
    strip.text=element_text(face="bold")
  )

# Each point represents one DMR. 
# The boxplots show the distribution of methylation levels across the identified 
# DMRs in each condition. 
# 
# The two panels separate DMRs that are hypermethylated or hypomethylated at 
# 36 hpf.



