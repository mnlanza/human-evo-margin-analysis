#!/bin/bash

# Exit on error, undefined var, or pipeline fail
set -euo pipefail

# === Usage Check ===
if [ $# -lt 6 ]; then
  echo "Usage: $0 POS GENE_START GENE_END SEQ_ID AID INPUT_FASTA [left_margin=2000] [right_margin=1000]"
  exit 1
fi

# === Positional and Optional Arguments ===
pos="$1"
gene_start="$2"
gene_end="$3"
seq_id="$4"
aid="$5"
input_fasta="$6"
left_margin="${7:-2000}"
right_margin="${8:-1000}"

# Check required files exist
if [ ! -f "$input_fasta" ]; then
  echo "Error: Input FASTA file $input_fasta not found"
  exit 1
fi

if [ ! -f "input/codon_table" ]; then
  echo "Error: input/codon_table not found"
  exit 1
fi

# Create required directories
mkdir -p output jobs figures

# Define job ID first (moved up)
job_id="$(echo "$aid" | tr '[:upper:]' '[:lower:]')"
job_version="$(echo "${seq_id//_/-}-${pos}" | tr '[:upper:]' '[:lower:]')"

# Create aid-specific directories
mkdir -p "output/${job_id}" "figures/${job_id}"

# File naming
fasta_out="output/${job_id}/query_${seq_id}_${pos}.fasta"
codon_out="output/${job_id}/query_${seq_id}_${pos}.tab"
job_dir="jobs/${job_id}-${job_version}"
compare_out="output/${job_id}/compare_strands_${seq_id}_${pos}.tab"

# Create job directory before running anything
mkdir -p "$job_dir"


# generate all codon variants at position 83
# change seq-id for each population
python3 scripts/generate_codons_variants.py \
  --fasta "$input_fasta" \
  --codon-table input/codon_table \
  --aa-coord "$pos" \
  --seq-id "$seq_id" \
  --output-fasta "$fasta_out" \
  --output-codon-table "$codon_out" \
  --left-margin "$left_margin" \
  --right-margin "$right_margin" \
  --gene-start "$gene_start" \
  --gene-end "$gene_end"

# submit job
evo_gcp submit --job "$job_id" \
  --output_type summary_only \
  --input_fasta "$(pwd)/$fasta_out" \
  --job_version "$job_version" \
  --wait

# download results
evo_gcp download --job "$job_id" \
  --job_version "$job_version" \
  --jobs_dir "$(pwd)/jobs" \
  --output_type summary_only

# Create strand comparison table
Rscript -e "
source('scripts/create_strand_table.r')
create_strand_table(
  ifn_tsv='${job_dir}/output/input_summary.txt',
  ofn='$compare_out')
"

# Plot strand comparison
Rscript -e "
source('scripts/plot_strand_scatter.r')
plot_strand_scatter(
  ifn_tab='$compare_out',
  ifn_codon='$codon_out',
  title='${seq_id}_pos_${pos}',
  fdir='figures/${job_id}')
"
