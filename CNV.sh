#!/bin/bash

set -e
set -u

start_time=$(date +%Y%m%d_%H%M%S)

WORK_DIR="reference"
NORMAL_BAM_DIR="gastric_work/normal/alignData"
GENOME_FILES_DIR="$WORK_DIR/genome_files"
FASTA="hg38.fa"
REFFLAT="refFlat.txt"
ACCESS_BED="$GENOME_FILES_DIR/access-1kb.hg38.bed"

mkdir -p $WORK_DIR
LOG_DIR="$WORK_DIR/logs_${start_time}"
mkdir -p $LOG_DIR

NORMAL_BAMS=(
    "$NORMAL_BAM_DIR/G005-2N.sort.rmdup.bam"
    "$NORMAL_BAM_DIR/G008-2N.sort.rmdup.bam"
    "$NORMAL_BAM_DIR/G1-2N.sort.rmdup.bam"
    "$NORMAL_BAM_DIR/G136-2N.sort.rmdup.bam"
    "$NORMAL_BAM_DIR/G295-2N.sort.rmdup.bam"
    "$NORMAL_BAM_DIR/G3-4N.sort.rmdup.bam"
    "$NORMAL_BAM_DIR/G41N.sort.rmdup.bam"
    "$NORMAL_BAM_DIR/G42N.sort.rmdup.bam"
    "$NORMAL_BAM_DIR/G43N.sort.rmdup.bam"
    "$NORMAL_BAM_DIR/G44N.sort.rmdup.bam"
    "$NORMAL_BAM_DIR/G45N.sort.rmdup.bam"
    "$NORMAL_BAM_DIR/G46N.sort.rmdup.bam"
    "$NORMAL_BAM_DIR/G47N.sort.rmdup.bam"
    "$NORMAL_BAM_DIR/G48N.sort.rmdup.bam"
    "$NORMAL_BAM_DIR/G659-2N.sort.rmdup.bam"
    "$NORMAL_BAM_DIR/G744-2N.sort.rmdup.bam"
    "$NORMAL_BAM_DIR/G754-1N.sort.rmdup.bam"
    "$NORMAL_BAM_DIR/G847-2N.sort.rmdup.bam"
    "$NORMAL_BAM_DIR/G849-2N.sort.rmdup.bam"
    "$NORMAL_BAM_DIR/G967-2N.sort.rmdup.bam"
    "$NORMAL_BAM_DIR/G968-2N.sort.rmdup.bam"
)

processed_count=0
failed_samples=()

for bam in "${NORMAL_BAMS[@]}"; do
    sample_name=$(basename "$bam" .sort.rmdup.bam)
    processed_count=$((processed_count + 1))

    if [ -f "$WORK_DIR/${sample_name}.targetcoverage.cnn" ]; then
        mv "$WORK_DIR/${sample_name}.targetcoverage.cnn" "$WORK_DIR/${sample_name}.targetcoverage.cnn.bak"
    fi

    if cnvkit.py coverage "$bam" "$ACCESS_BED" \
        -o "$WORK_DIR/${sample_name}.targetcoverage.cnn" \
        -p 2; then
    else
        failed_samples+=("$sample_name")
    fi
done

if [ ${#failed_samples[@]} -gt 0 ]; then
    printf '%s\n' "${failed_samples[@]}"
    exit 1
fi

missing_files=0
for bam in "${NORMAL_BAMS[@]}"; do
    sample_name=$(basename "$bam" .sort.rmdup.bam)
    if [ ! -f "$WORK_DIR/${sample_name}.targetcoverage.cnn" ]; then
        missing_files=1
    fi
done


COVERAGE_FILES=()
for bam in "${NORMAL_BAMS[@]}"; do
    sample_name=$(basename "$bam" .sort.rmdup.bam)
    COVERAGE_FILES+=("$WORK_DIR/${sample_name}.targetcoverage.cnn")
done


cnvkit.py reference "${COVERAGE_FILES[@]}" \
    --fasta "$FASTA" \
    -o "$WORK_DIR/reference.cnn" \
    2>&1 | tee "$LOG_DIR/reference_build.log"

if [ -f "$WORK_DIR/reference.cnn" ]; then

    cnvkit.py metrics "$WORK_DIR/reference.cnn" > "$WORK_DIR/reference_metrics.txt"

    cat "$WORK_DIR/reference_metrics.txt"

    cnvkit.py diagram "$WORK_DIR/reference.cnn" \
        -o "$WORK_DIR/reference_diagram.pdf" \
        2>&1 | tee "$LOG_DIR/diagram.log"

    md5sum "$WORK_DIR/reference.cnn" > "$WORK_DIR/reference.cnn.md5"
else
    exit 1
fi

EOF



TUMOR_BAM="${sample}.sort.rmdup.bam"
REFERENCE_CNN="/reference.cnn"
OUT_DIR="result/${sample}"
OUTPUT_BASENAME=$(basename "${TUMOR_BAM}" .bam)

cnvkit.py batch "$TUMOR_BAM" \
    -r "$REFERENCE_CNN" \
    --output-dir "$OUT_DIR" \
    -p 8

echo "CNVkit batch 分析完成。"


CNR_FILE="$OUT_DIR/${OUTPUT_BASENAME}.cnr"
CNS_FILE="$OUT_DIR/${OUTPUT_BASENAME}.cns"


cnvkit.py scatter "$CNR_FILE" \
    -s "$CNS_FILE" \
    -o "$OUT_DIR/${OUTPUT_BASENAME}_scatter.pdf"


cnvkit.py diagram "$CNR_FILE" \
    -s "$CNS_FILE" \
    -o "$OUT_DIR/${OUTPUT_BASENAME}_diagram.pdf"

