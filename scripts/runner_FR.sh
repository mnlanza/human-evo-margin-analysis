#!/bin/bash

# Usage:
#   ./scripts/format_human_data.sh
#
# Description:
#   This script processes human data by running process_human_data.py with preset arguments.
#   It will generate two output files:
#   1. TSV file: input/updated_data.tsv (1 based indexing)
#   2. FASTA file: input/human_contigs.fasta
#
# Note:
#   Input directories should be named with _in suffix (e.g., BAA_in, BAD_in, etc.)
#
# To run:
#   1. Make sure you're in the project root directory
#   2. Make the script executable: chmod +x scripts/format_human_data.sh
#   3. Run the script: ./scripts/format_human_data.sh

mkdir -p input

echo "=== Generating source FASTA file ==="
# output-tsv and fasta commands are optional
python scripts/process_human_data.py \
  --output-tsv input/updated_data.tsv \
  --write-fasta \
  --fasta-path input/human_contigs_src.fasta

echo -e "\n=== Generating mutated FASTA file ==="
# for generating fasta files with reference and mutations
python scripts/FR_format_human_data.py \
  --write-fasta \
  --fasta-path input/mut_human_src.fasta
