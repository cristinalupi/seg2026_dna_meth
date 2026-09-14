# seg2026_dna_meth


## 1. File Organization & Dataset

The raw sequencing data for the *Danio rerio* (zebrafish) samples comprise paired-end reads across two developmental stages (4hpf and 36hpf) with two biological replicates each. The raw data can be retrieved from the NCBI SRA database using the following accessions:

* [SRX5027434](https://www.ncbi.nlm.nih.gov/sra?term=SRX5027434) (4hpf, rep 1)
* [SRX5027442](https://www.ncbi.nlm.nih.gov/sra?term=SRX5027442) (4hpf, rep 2)
* [SRX5027437](https://www.ncbi.nlm.nih.gov/sra?term=SRX5027437) (36hpf, rep 1)
* [SRX5027445](https://www.ncbi.nlm.nih.gov/sra?term=SRX5027445) (36hpf, rep 2)

Library preparation was performed using the Accel-NGS Methyl-Seq DNA Library Kit (Swift Biosciences).

```text
/data/classes/seg_epig2026/data/methylation_data/
├── danio_36hpf_rep1_R1.fastq.gz
├── danio_36hpf_rep1_R2.fastq.gz
├── danio_36hpf_rep2_R1.fastq.gz
├── danio_36hpf_rep2_R2.fastq.gz
├── danio_4hpf_rep1_R1.fastq.gz
├── danio_4hpf_rep1_R2.fastq.gz
├── danio_4hpf_rep2_R1.fastq.gz
├── danio_4hpf_rep2_R2.fastq.gz
├── GRCz11.fa
└── lambda.fasta
```
### 1.1. Create a Symbolic Link to the Input Data
```bash
mkdir -p methylation_analysis/input_data
ln -s /data/classes/seg_epig2026/data/meth/methylation_analysis/* methylation_analysis/input_data/
```
---

## 2. Quality Control & Hard Trimming

* **Quality Control:** General quality trimming has already been performed. 
* **Hard Trimming:** Focuses exclusively on library-specific requirements to eliminate sequencing artifacts and bisulfite conversion biases from read extremities.

```bash
mkdir -p fastq_trimmed

for sample in danio_36hpf_rep1 danio_36hpf_rep2 danio_4hpf_rep1 danio_4hpf_rep2; do
  fastp -i input_data/${sample}_R1.fastq.gz \
        -I input_data/${sample}_R2.fastq.gz \
        -o fastq_trimmed/${sample}_R1_trimmed.fastq.gz \
        -O fastq_trimmed/${sample}_R2_trimmed.fastq.gz \
        --trim_front1 10 \
        --trim_front2 15 \
        --trim_tail1 10 \
        --trim_tail2 10
done
```
- [Bismark Library Type Documentation](https://felixkrueger.github.io/Bismark/usage/library-types/)
- [fastp Documentation](https://github.com/opengene/fastp)
## 3. Genome Preparation

### 3.1. Create the Bismark genome directory

Create the directory required by Bismark.

```bash
mkdir -p genome_dir
```

### 3.2. Combine the reference genomes

Combine the zebrafish reference genome with the lambda sequence and save the resulting reference genome in the Bismark genome directory.

```bash
cat input_data/GRCz11.fa input_data/lambda.fasta > genome_dir/GRCz11_lambda.fa
```

### 3.3. Prepare the genome for bisulfite alignment

Prepare the reference genome and generate the required bisulfite-converted genomes.
- [Bismark Genome Preparation Documentation](https://felixkrueger.github.io/Bismark/usage/genome-preparation/)
```bash
bismark_genome_preparation genome_dir/
```


## 4. Alignment

Once the genome is prepared and the reads are trimmed, the next step is to align the paired-end WGBS reads against the bisulfite-converted reference genome using `bismark` in a loop over the trimmed files.
- [Bismark Alignment Documentation](https://felixkrueger.github.io/Bismark/usage/alignment/) 

```bash
mkdir -p bismark_aligned

for fwd in fastq_trimmed/*_R1_trimmed.fastq.gz; do
    rev="${fwd/_R1_trimmed.fastq.gz/_R2_trimmed.fastq.gz}"
    base=$(basename "$fwd" "_R1_trimmed.fastq.gz")

    bismark \
        --genome genome_dir/ \
        -1 "$fwd" \
        -2 "$rev" \
        --output_dir bismark_aligned \
        --parallel 4 \
        --bowtie2 \
        --local \
        -X 2000
done
```
Parameters:
- `--bowtie2`: uses Bowtie 2 for read alignment.
- `--local`: performs local alignment, allowing reads to align without requiring the entire read to match the reference.
- `-X 2000`: sets the maximum allowed insert size for paired-end reads to 2000 bp.
- `--parallel 4`: runs four parallel alignment processes.
- `--output_dir`: specifies the directory where Bismark output files are saved.

## 5. Deduplication
- [Bismark Library Type Documentation](https://felixkrueger.github.io/Bismark/usage/library-types/)
- [Bismark Deduplication Documentation](https://felixkrueger.github.io/Bismark/usage/deduplication/)

```bash
mkdir -p bismark_deduplicated

for bam in bismark_aligned/*_pe.bam; do
    deduplicate_bismark \
        --bam \
        --paired \
        --output_dir bismark_deduplicated \
        "$bam" &
done
wait
```
## 6. Methylation Extraction
- [MethylDackel Documentation](https://github.com/dpryan79/methyldackel)
### 6.1 Sort bam file
```bash
mkdir -p bismark_deduplicated/sorted

for bam in bismark_deduplicated/*_pe.deduplicated.bam; do
    base=$(basename "$bam" .bam)

    samtools sort \
        -o "bismark_deduplicated/sorted/${base}.sorted.bam" \
        "$bam" &
done
wait
```
### 6.2 Index sorted bam files
```bash
for bam in bismark_deduplicated/sorted/*.sorted.bam; do
    samtools index "$bam" &
done
wait
```
### 6.3 Index genome file
```bash
samtools faidx genome_dir/GRCz11_lambda.fa
```

### 6.4 M-bias Analysis

```bash
mkdir -p methylation_bias

for bam in bismark_deduplicated/sorted/*_pe.deduplicated.sorted.bam; do
    prefix=$(basename "$bam" .bam)

    MethylDackel mbias \
        -@ 8 \
        genome_dir/GRCz11_lambda.fa \
        "$bam" \
        methylation_bias/"$prefix"
done

```
### 6.5 Define M-bias filtering parameters
```bash
declare -A MBIAS=(
    [danio_4hpf_rep1]="--OT 0,0,0,0 --OB 0,0,0,0"
    [danio_4hpf_rep2]="--OT 0,0,0,0 --OB 0,0,0,0"
    [danio_36hpf_rep1]="--OT 0,0,0,0 --OB 0,0,0,0"
    [danio_36hpf_rep2]="--OT 0,0,0,0 --OB 0,0,0,0"
)
```

### 6.6 Extract methylation information

```bash
mkdir -p methylation_output

for sample in danio_4hpf_rep1 danio_4hpf_rep2 danio_36hpf_rep1 danio_36hpf_rep2; do
    MethylDackel extract \
        genome_dir/GRCz11_lambda.fa \
        bismark_deduplicated/sorted/${sample}_R1_trimmed_bismark_bt2_pe.deduplicated.sorted.bam \
        --minOppositeDepth 10 \
        --maxVariantFrac 0.5 \
        ${MBIAS[$sample]} \
        --mergeContext \
        -p 8 \
        -o methylation_output/${sample}
done
```

Parameters:
* `--minOppositeDepth 10`: requires a minimum of 10 reads covering the opposite strand.
* `--maxVariantFrac 0.5`: excludes positions where more than 50% of the reads support a non-reference base, helping to remove potential genetic variants.
* `--OT`, `--OB`: specify the positions to exclude from each read according to the M-bias analysis. The values are sample-specific and should be filled in using the results obtained from the M-bias plots.
* `--mergeContext`: merges methylation calls from both strands at CpG sites.
* `-p 8`: uses 8 threads for the extraction.
* `-o`: specifies the output prefix for each sample.

## 7. Generate DSS files
The MethylDackel methylation output is converted into the format required by DSS for downstream differential methylation analysis.
- [DSS documentation](https://www.bioconductor.org/packages/release/bioc/vignettes/DSS/inst/doc/DSS.html#3_Using_DSS_for_BS-seq_differential_methylation_analysis)
### DSS and methylKit input preparation

**Input:** MethylDackel `*.bedGraph` files containing the following information:

| Column | Description |
|---|---|
| `chr` | Chromosome |
| `start` | Start coordinate |
| `end` | End coordinate |
| `methylation` | Methylation percentage |
| `methylated` | Number of methylated reads |
| `unmethylated` | Number of unmethylated reads |

The same MethylDackel `*.bedGraph` files are used as input for both **DSS** and **methylKit**. The files are converted into the respective formats required by each downstream analysis.

#### DSS format

**Transformation:**

- `N` = methylated reads + unmethylated reads
- `X` = methylated reads

**Output:** One DSS-formatted file per sample:

| Column | Description |
|---|---|
| `chr` | Chromosome |
| `pos` | Genomic position |
| `N` | Total number of reads |
| `X` | Number of methylated reads |

#### methylKit format

**Transformation:**

- `chrBase` = chromosome + `"."` + genomic position
- `chr` = chromosome
- `base` = genomic position
- `strand` = `"+"`
- `coverage` = methylated reads + unmethylated reads
- `freqC` = 100 × methylated reads / coverage
- `freqT` = 100 × unmethylated reads / coverage

**Output:** One methylKit-formatted file per sample:

| Column | Description |
|---|---|
| `chrBase` | Unique chromosome-position identifier |
| `chr` | Chromosome |
| `base` | Genomic position |
| `strand` | Strand, set to `+` |
| `coverage` | Total number of reads |
| `freqC` | Percentage of methylated reads |
| `freqT` | Percentage of unmethylated reads |

```bash
mkdir -p rstudio_analysis

for sample in danio_4hpf_rep1 danio_4hpf_rep2 danio_36hpf_rep1 danio_36hpf_rep2; do

    input="methylation_output/${sample}_CpG.bedGraph"

    # methylKit
    awk 'BEGIN {
        OFS="\t";
        print "chrBase","chr","base","strand","coverage","freqC","freqT"
    }
    /^track/ {next}
    {
        chr=$1
        pos=$2
        N=$5+$6
        X=$5

        print chr "." pos, chr, pos, "+", N, 100*X/N, 100*(N-X)/N
    }' "$input" > "rstudio_analysis/${sample}_methylKit.txt"


    # DSS
    awk 'BEGIN {
        OFS="\t";
        print "chr","pos","N","X"
    }
    /^track/ {next}
    {
        chr=$1
        pos=$2
        N=$5+$6
        X=$5

        print chr, pos, N, X
    }' "$input" > "rstudio_analysis/${sample}_DSS.txt"

done
```







