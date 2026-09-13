################################################################################
# 0. Set working directory
################################################################################

setwd("/path/methylation_analysis/rstudio_analysis")

################################################################################
# 1. Install required packages
################################################################################

# Install BiocManager if it is not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required Bioconductor packages
packages <- c("methylKit", "GenomicRanges")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg, update = FALSE)
  }
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

library(methylKit)
library(GenomicRanges)
library(dplyr)
library(ggplot2)

################################################################################
# 3. Read in methylation data
################################################################################

dat1.1 <- read.table("danio_4hpf_rep1_DSS.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat1.2 <- read.table("danio_4hpf_rep2_DSS.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat2.1 <- read.table("danio_36hpf_rep1_DSS.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat2.2 <- read.table("danio_36hpf_rep2_DSS.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)


unique(dat1.1$chr)
unique(dat1.2$chr)
unique(dat2.1$chr)
unique(dat2.2$chr)

################################################################################
# 4. Bisulfite conversion efficiency using lambda DNA
################################################################################

lambda_data <- bind_rows(
  dat1.1 %>% filter(chr=="lambda") %>% mutate(sample="4hpf_rep1"),
  dat1.2 %>% filter(chr=="lambda") %>% mutate(sample="4hpf_rep2"),
  dat2.1 %>% filter(chr=="lambda") %>% mutate(sample="36hpf_rep1"),
  dat2.2 %>% filter(chr=="lambda") %>% mutate(sample="36hpf_rep2")
)

lambda_conversion <- lambda_data %>% group_by(sample) %>% summarise(conversion_efficiency=1-sum(X)/sum(N))

lambda_conversion

################################################################################
# 5. Bisulfite conversion efficiency using mitochondrial DNA
################################################################################

mt_data <- bind_rows(
  dat1.1 %>% filter(chr=="MT") %>% mutate(sample="4hpf_rep1"),
  dat1.2 %>% filter(chr=="MT") %>% mutate(sample="4hpf_rep2"),
  dat2.1 %>% filter(chr=="MT") %>% mutate(sample="36hpf_rep1"),
  dat2.2 %>% filter(chr=="MT") %>% mutate(sample="36hpf_rep2")
)

mt_conversion <- mt_data %>% group_by(sample) %>% summarise(conversion_efficiency=1-sum(X)/sum(N))

mt_conversion

################################################################################
# 6. Create methylKit objects
################################################################################

myobj <- methRead(
  list(
    "danio_4hpf_rep1_methylKit.txt",
    "danio_4hpf_rep2_methylKit.txt",
    "danio_36hpf_rep1_methylKit.txt",
    "danio_36hpf_rep2_methylKit.txt"
  ),
  sample.id = list(
    "4hpf_rep1",
    "4hpf_rep2",
    "36hpf_rep1",
    "36hpf_rep2"
  ),
  assembly = "danio",
  treatment = c(0, 0, 1, 1),
  context = "CpG"
)

################################################################################
# 7. Sample correlation 
################################################################################

target_seq <- "4"

myobj_target <- new("methylRawList", 
                    lapply(myobj, function(x) {
                      df <- getData(x)
                      x[df$chr == target_seq, ]
                    }),
                    treatment=myobj@treatment)

meth_target <- unite(myobj_target, destrand=FALSE)

getCorrelation(meth_target, method="pearson", plot=TRUE)

################################################################################
# 8. Sample clustering
################################################################################

clusterSamples(meth_target,
               dist="correlation",
               method="ward.D2",
               plot=TRUE)

################################################################################
# 9. Principal component analysis
################################################################################

PCASamples(
  meth_target,
  screeplot = FALSE,
  adj.lim = c(0.5, 1.5)
)
