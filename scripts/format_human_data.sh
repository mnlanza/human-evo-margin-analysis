#!/bin/bash

# Usage:
#   ./scripts/format_human_data.sh
#
# Description:
#   This script processes human data by running process_human_data.py with preset arguments.
#   It will generate two output files:
#   1. TSV file: output/updated_data.tsv
#   2. FASTA file: output/human_contigs.fasta
#
# To run:
#   1. Make sure you're in the project root directory
#   2. Make the script executable: chmod +x scripts/format_human_data.sh
#   3. Run the script: ./scripts/format_human_data.sh

# output-tsv and fasta commands are optional
python scripts/process_human_data.py \
  --output-tsv input/updated_data.tsv \
  --write-fasta \
  --fasta-path input/human_src_FR.fasta
