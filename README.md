# seg2026_dna_meth

## Input data
### WGBS dataset:
https://www.ncbi.nlm.nih.gov/sra?term=SRX5027445
```
bsseq-analysis/
├── input_fastq_directory/
│   ├── danio_4hpf_rep1_R1.fastq
│   ├── danio_4hpf_rep1_R2.fastq
│   ├── danio_4hpf_rep2_R1.fastq
│   ├── danio_4hpf_rep2_R2.fastq
│   ├── danio_36hpf_rep1_R1.fastq
│   ├── danio_36hpf_rep1_R2.fastq
│   ├── danio_36hpf_rep2_R1.fastq​
│   └── danio_36hpf_rep2_R2.fastq
└── input_genome_directory/
```
## 1. Fastq quality control
Link to the tool page: https://github.com/s-andrews/fastqc
```
# 1.1. Create the output_dir
mkdir fastqc_out
# 1.2. Quality control
## General command
fastqc -t8 -o <output_dir> <file.fastq>​
## Automated command
fastqc –t8 -o fastqc_out *fastq​
```
## 2. Short read processing
Link to the tool page: https://github.com/opengene/fastp
```
# 2.1. Create the output_dir
mkdir fastp_out
# 2.2. Short read processing
INPUT_DIR="$1"
THREADS="${2:-8}"

OUTPUT_DIR="${INPUT_DIR}/trimmed_fastq"

mkdir -p "$OUTPUT_DIR"

shopt -s nullglob

for R1 in "${INPUT_DIR}"/*_R1.fastq.gz
do
    R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
    SAMPLE=$(basename "$R1" _R1.fastq.gz)

    [ -f "$R2" ] || { echo "Missing pair for $R1"; exit 1; }

    fastp \
        -i "$R1" -I "$R2" \
        -o "${OUTPUT_DIR}/${SAMPLE}_R1.trimmed.fastq.gz" \
        -O "${OUTPUT_DIR}/${SAMPLE}_R2.trimmed.fastq.gz" \
        --detect_adapter_for_pe \
        --cut_front --cut_front_window_size 5 --cut_front_mean_quality 20 \
        --cut_tail --cut_tail_window_size 5 --cut_tail_mean_quality 20 \
        --trim_poly_g \
        --trim_poly_x \
        --length_required 25 \
        --thread <threads>\
        --html "${OUTPUT_DIR}/${SAMPLE}_fastp.html"
done
```
## 3. Bismark alignment
Link to the tool page: https://felixkrueger.github.io/Bismark/



`
bismark_methylation_extractor --paired-end --comprehensive --bedGraph --output_dir methylation_output/extracted/ deduplicated_bam/archivo_dedup.bam
``

