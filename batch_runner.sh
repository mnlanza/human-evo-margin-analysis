#!/bin/bash
set -euo pipefail

# Create all required directories upfront
mkdir -p input output jobs figures

# Shared parameters
left_margin=2000
right_margin=1000
input_fasta="input/human_contigs_src.fasta"

# Read TSV file line by line, skipping header
while IFS=$'\t' read -r aid gene contig start end strand flipped src_codon tgt_codon mut_pos || [ -n "$contig" ]; do
  # Skip header line
  if [ "$contig" = "contig" ]; then
    continue
  fi

  echo "Running $contig with center position $mut_pos and gene range $start-$end (aid: $aid)"

  # Generate range from mut_pos - 2 to mut_pos + 2
  for offset in {-2..2}; do
    pos=$((mut_pos + offset))
    echo "  → Running position $pos"
    ./scripts/FR_all_pipeline.sh "$pos" "$start" "$end" "$contig" "$aid" "$input_fasta" "$left_margin" "$right_margin"
  done
done < "input/updated_data.tsv"