#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Activate conda environment
# -----------------------------
cd /Users/rus/Library/CloudStorage/Dropbox/SEPAL_AI/4_Huntington_NSC
source /Users/rus/miniconda3/etc/profile.d/conda.sh
conda activate rna_tools
# -----------------------------------------------
# CONFIGURATION
# -----------------------------------------------
THREADS=16
BASE_DIR="/Users/rus/Library/CloudStorage/Dropbox/SEPAL_AI/4_Huntington_NSC"
INPUT_DIR="${BASE_DIR}/fastq_downsampled"
ALIGN_DIR="${BASE_DIR}/aligned"
LOG_DIR="${BASE_DIR}/logs"
COUNT_DIR="${BASE_DIR}/counts"
GENOME_INDEX_DIR="${BASE_DIR}/hisat2_index"
GTF_FILE="${BASE_DIR}/gencode.v43.annotation.gtf"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz"

mkdir -p "$ALIGN_DIR" "$LOG_DIR" "$COUNT_DIR" "$GENOME_INDEX_DIR"

# -----------------------------------------------
# STEP 1: Download HISAT2 index for GRCh38 (human)
# -----------------------------------------------
cd "$GENOME_INDEX_DIR"
if [ ! -d "grch38" ]; then
  echo "📥 Downloading HISAT2 index for human (GRCh38)..."
  curl -O ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz
  tar -xzf grch38.tar.gz
else
  echo "✅ HISAT2 index already exists. Skipping download."
fi
cd "$BASE_DIR"

# -----------------------------------------------
# STEP 2: Download GTF annotation file
# -----------------------------------------------
if [ ! -f "$GTF_FILE" ]; then
  echo "📥 Downloading GENCODE annotation (v43)..."
  curl -o "${GTF_FILE}.gz" "$GTF_URL"
  gunzip -f "${GTF_FILE}.gz"
else
  echo "✅ GTF annotation already exists. Skipping download."
fi

# -----------------------------------------------
# STEP 3: Align with HISAT2
# -----------------------------------------------
echo "🎯 Starting alignment..."
SAMPLES=($(ls ${INPUT_DIR}/*_1.subsampled.fastq.gz | sed 's/.*\///; s/_1.subsampled.fastq.gz//'))

for sample in "${SAMPLES[@]}"; do
  echo "🔧 Aligning $sample..."
  hisat2 -p $THREADS \
    -x $GENOME_INDEX_DIR/grch38/genome \
    -1 ${INPUT_DIR}/${sample}_1.subsampled.fastq.gz \
    -2 ${INPUT_DIR}/${sample}_2.subsampled.fastq.gz \
    2> ${LOG_DIR}/${sample}_hisat2.log | \
    samtools sort -@ $THREADS -o ${ALIGN_DIR}/${sample}.bam
  samtools index ${ALIGN_DIR}/${sample}.bam
  echo "✅ $sample alignment done."
done

# -----------------------------------------------
# STEP 4: Count with featureCounts
# -----------------------------------------------
echo "🔢 Running featureCounts..."
featureCounts -T $THREADS -p -t exon -g gene_name \
  -a "$GTF_FILE" \
  -o ${COUNT_DIR}/raw_counts_gene_sym.txt ${ALIGN_DIR}/*.bam \
  &> ${LOG_DIR}/featureCounts_gene_sym.log

# -----------------------------------------------
# STEP 5: Clean counts matrix
# -----------------------------------------------
echo "🧹 Cleaning counts table..."
{ printf "GeneSymbol\t"; head -n 2 ${COUNT_DIR}/raw_counts_gene_sym.txt | tail -n 1 | cut -f7-; } > ${COUNT_DIR}/final_counts_symbols.tsv
tail -n +3 ${COUNT_DIR}/raw_counts_gene_sym.txt | \
  awk -v OFS="\t" '{ out=$1; for(i=7;i<=NF;i++) out=out OFS $i; print out }' >> ${COUNT_DIR}/final_counts_symbols.tsv
sed -i '' '1 s|aligned/||g; 1 s|\.bam||g' ${COUNT_DIR}/final_counts_symbols.tsv

echo "✅ RNA-seq pipeline complete."
head ${COUNT_DIR}/final_counts_symbols.tsv | column -t
