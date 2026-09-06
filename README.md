# seg2026_dna_meth


## 1. File Organization & Dataset

The raw sequencing data for the *Danio rerio* (zebrafish) samples comprise paired-end reads across two developmental stages (4hpf and 36hpf) with two biological replicates each. The raw data can be retrieved from the NCBI SRA database using the following accessions:

* [SRX5027445](https://www.ncbi.nlm.nih.gov/sra?term=SRX5027445) (4hpf, rep 1)
* [SRX5027437](https://www.ncbi.nlm.nih.gov/sra?term=SRX5027437) (4hpf, rep 2)
* [SRX5027442](https://www.ncbi.nlm.nih.gov/sra?term=SRX5027442) (36hpf, rep 1)
* [SRX5027434](https://www.ncbi.nlm.nih.gov/sra?term=SRX5027434) (36hpf, rep 2)

Library preparation was performed using the Accel-NGS Methyl-Seq DNA Library Kit (Swift Biosciences).

---

## 2. Quality Control & Hard Trimming

* **Quality Control:** General quality trimming has already been performed. 
* **Hard Trimming:** Focuses exclusively on library-specific requirements to eliminate sequencing artifacts and bisulfite conversion biases from read extremities.

### Execution

```bash
mkdir -p fastq_trimmed

for sample in danio_36hpf_rep1 danio_36hpf_rep2 danio_4hpf_rep1 danio_4hpf_rep2; do
  fastp -i ${sample}_R1.fastq.gz \
        -I ${sample}_R2.fastq.gz \
        -o fastq_trimmed/trimmed_${sample}_R1.fastq.gz \
        -O fastq_trimmed/trimmed_${sample}_R2.fastq.gz \
        --trim_front1 10 \
        --trim_front2 15 \
        --trim_tail1 10 \
        --trim_tail2 10 \
        -w 12
done
```
### Procedure

1. **Output Directory Creation:** Creates the `fastq_trimmed/` directory if it does not already exist.
2. **Batch Iteration:** Loops through each of the four biological replicates (`danio_36hpf_rep1`, `danio_36hpf_rep2`, `danio_4hpf_rep1`, `danio_4hpf_rep2`).
3. **Hard Trimming Execution:** Applies `fastp` with 12 parallel threads to remove the library-specific fixed bases at the 5' and 3' extremities (`--trim_front1 10`, `--trim_front2 15`, `--trim_tail1 10`, `--trim_tail2 10`), saving the output files to the target directory.

## 3. Genome Preparation

Before mapping bisulfite-seq reads, the reference genome (`GRCz11.fa`) must be prepared using Bismark to generate two copies of the genome with converted cytosines (C-to-T and G-to-A).

For detailed usage, refer to the [Bismark Genome Preparation Documentation](https://felixkrueger.github.io/Bismark/usage/genome-preparation/).

### Execution

```bash
bismark_genome_preparation genome_dir/
```
## 4. Alignment

Once the genome is prepared and the reads are trimmed, the next step is to align the paired-end WGBS reads against the bisulfite-converted reference genome using `bismark` in a loop over the trimmed files.

### Script Execution

```bash
mkdir -p bismark_aligned

for sample in danio_36hpf_rep1 danio_36hpf_rep2 danio_4hpf_rep1 danio_4hpf_rep2; do
  bismark --genome genome_dir/ \
        -1 fastq_WGBS/fastq_trimmed/trimmed_${sample}_R1.fastq.gz \
        -2 fastq_WGBS/fastq_trimmed/trimmed_${sample}_R2.fastq.gz \
        --parallel 8 \
        --output_dir bismark_aligned
done
```
## 5. Deduplication

To remove PCR duplicates generated during library amplification, use `deduplicate_bismark` on the paired-end alignment files.

### Script Execution

```bash
for bam in *_pe.bam; do
  deduplicate_bismark --bam --paired "$bam" &
done
wait
```
## 6. Methylation Extraction

After deduplication, the methylation state of each cytosine is extracted using `bismark_methylation_extractor`.

### Script Execution

```bash
mkdir -p methylation_output

for bam in *_pe.deduplicated.bam; do
  bismark_methylation_extractor --paired-end --comprehensive --bedGraph \
        --output_dir methylation_output \
        "$bam" &
done
wait
```
scp seg2user@158.42.124.228:/remote/path/file.ext .
