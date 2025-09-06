#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Output folder setup
# -----------------------------
cd /Users/rus/Library/CloudStorage/Dropbox/SEPAL_AI/4_Huntington_NSC
mkdir -p fastq1
cd fastq1

# -----------------------------
# HTT-KO and Control (IC1) SRR runs
# -----------------------------
SRRS=(
  "SRR29493482"  # HTT mutation (109CAG)
  "SRR29493483"  # HTT mutation (109CAG)
  "SRR29493484"  # HTT mutation (109CAG)
  "SRR29493485"  # HTT mutation (109CAG)
  "SRR29493486"  # IC1 (Control)
  "SRR29493487"  # IC1 (Control)
  "SRR29493488"  # IC1 (Control)
)


# Number of parallel jobs
THREADS=4

# -----------------------------
# Download function
# -----------------------------
download_srr () {
  SRR=$1
  echo "🔍 Resolving $SRR via ENA API..."

  # Get FASTQ FTP links
  URLS=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR}&result=read_run&fields=fastq_ftp&format=tsv" \
    | tail -n +2 | tr ';' '\n' | sed 's|^|ftp://|')

  if [[ -z "$URLS" ]]; then
    echo "❌ No FASTQ URLs found for $SRR"
    return 1
  fi

  for URL in $URLS; do
    FILE=$(basename "$URL")
    if [[ -f "$FILE" ]]; then
      echo "⚠️  Skipping existing $FILE"
    else
      echo "➡️  Downloading $FILE ..."
      wget -c "$URL"
    fi
  done

  echo "✅ Done: $SRR"
}

# Export the function for xargs parallel execution
export -f download_srr

# -----------------------------
# Run downloads in parallel
# -----------------------------
echo "🚀 Starting parallel downloads for HTT-KO and IC1..."

printf "%s\n" "${SRRS[@]}" | xargs -n 1 -P "$THREADS" -I {} bash -c 'download_srr "$@"' _ {}

echo "🎉 All FASTQ files downloaded to ./fastq1/"
