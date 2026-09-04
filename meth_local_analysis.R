library(methylKit)

# 1. Set the working directory to your local project folder
setwd("/path/to/your/working_directory")

# 2. Load each sample individually
sample1 <- methRead(
  location = "extracted_output/SRR8208203_4hpf_rep1_R1_bismark_bt2_pe.deduplicated.bismark.cov.gz",
  sample.id = "4hpf_rep1",
  assembly = "danRer11",
  treatment = 0,
  context = "CpG",
  pipeline = "bismarkCoverage"
)

sample2 <- methRead(
  location = "extracted_output/SRR8208211-13_36hpf_rep1_R1_bismark_bt2_pe.deduplicated.bismark.cov.gz",
  sample.id = "36hpf_rep1",
  assembly = "danRer11",
  treatment = 1,
  context = "CpG",
  pipeline = "bismarkCoverage"
)

sample3 <- methRead(
  location = "extracted_output/SRR8208230-33_4hpf_rep2_R1_bismark_bt2_pe.deduplicated.bismark.cov.gz",
  sample.id = "4hpf_rep2",
  assembly = "danRer11",
  treatment = 0,
  context = "CpG",
  pipeline = "bismarkCoverage"
)

sample4 <- methRead(
  location = "extracted_output/SRR8208242-45_36hpf_rep2_R1_bismark_bt2_pe.deduplicated.bismark.cov.gz",
  sample.id = "36hpf_rep2",
  assembly = "danRer11",
  treatment = 1,
  context = "CpG",
  pipeline = "bismarkCoverage"
)

# 3. Combine individual samples into a methylRawList object
myobj <- new("methylRawList", list(sample1, sample2, sample3, sample4))

# 4. Explore descriptive statistics (methylation percentage and coverage)
getMethylationStats(myobj[[1]], plot = TRUE, both.strands = FALSE)
getCoverageStats(myobj[[1]], plot = TRUE, both.strands = FALSE)

# 5. Filter samples based on read coverage
filtered.myobj <- filterByCoverage(myobj, lo.count = 10, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)

# 6. Merge samples to find bases covered across all samples
meth <- unite(filtered.myobj, destrand = FALSE)

# 7. Sample correlation and hierarchical clustering quality control
getCorrelation(meth, plot = TRUE)
clusterSamples(meth, dist = "euclidean", method = "ward.D2", plot = TRUE)

# 8. Principal Component Analysis (PCA)
PCASamples(meth, screeplot = TRUE)
PCASamples(meth)

# 9. Calculate differential methylation
my.diffMeth <- calculateDiffMeth(meth)

# 10. Get differentially methylated regions (DMRs/bases) with specific thresholds
# (e.g., methylation difference > 25% and q-value < 0.01)
my.diffMeth.per <- getMethylDiff(my.diffMeth, difference = 25, qvalue = 0.01)

# 11. Plot methylation statistics and distribution of DMRs
getMethylDiffStats(my.diffMeth.per, plot = TRUE)
