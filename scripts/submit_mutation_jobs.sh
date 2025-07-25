#!/bin/bash
set -euo pipefail

# Create all required directories upfront
mkdir -p input output jobs figures

# Create aid-specific directories from TSV
while IFS=$'\t' read -r aid _ _ _ _ _ _ _ _ _ || [ -n "$aid" ]; do
  if [ "$aid" = "aid" ]; then continue; fi
  aid_lower="$(echo "$aid" | tr '[:upper:]' '[:lower:]')"
  mkdir -p "output/${aid_lower}" "figures/${aid_lower}"
done < "input/updated_data.tsv"

# Shared parameters
left_margin=2000
right_margin=1000
input_fasta="input/human_contigs_src.fasta"

# Process each mutation and submit jobs
while IFS=$'\t' read -r aid gene contig start end strand flipped src_codon tgt_codon mut_pos || [ -n "$contig" ]; do
  # Skip header line
  if [ "$contig" = "contig" ]; then
    continue
  fi

  echo "Submitting jobs for $contig with center position $mut_pos and gene range $start-$end (aid: $aid)"

  # Generate range from mut_pos - 2 to mut_pos + 2
  for offset in {-2..2}; do
    pos=$((mut_pos + offset))
    echo "  → Submitting position $pos"
    
    # Define variables needed for job submission
    job_id="$(echo "$aid" | tr '[:upper:]' '[:lower:]')"
    job_version="$(echo "${contig//_/-}-${pos}" | tr '[:upper:]' '[:lower:]')"
    fasta_out="output/${job_id}/query_${contig}_${pos}.fasta"
    codon_out="output/${job_id}/query_${contig}_${pos}.tab"
    job_dir="jobs/${job_id}-${job_version}"
    
    # Create job directory
    mkdir -p "$job_dir"

    # Generate codon variants
    python3 scripts/generate_codons_variants.py \
      --fasta "$input_fasta" \
      --codon-table input/codon_table \
      --aa-coord "$pos" \
      --seq-id "$contig" \
      --output-fasta "$fasta_out" \
      --output-codon-table "$codon_out" \
      --left-margin "$left_margin" \
      --right-margin "$right_margin" \
      --gene-start "$start" \
      --gene-end "$end"

    # Submit job without wait flag
    evo_gcp submit --job "$job_id" \
      --output_type summary_only \
      --input_fasta "$(pwd)/$fasta_out" \
      --job_version "$job_version"
  done
done < "input/updated_data.tsv"

echo "All jobs submitted successfully!" 