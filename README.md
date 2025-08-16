## margin_human

### Overview
End-to-end utilities for creating per-gene variant FASTA queries, submitting them to EVO 2 (via `evo_gcp`), and analyzing plus vs minus strand likelihoods for the mutation position and up to 2 positions upstream and downstream of the mutation. Includes helpers to validate gene/contig coordinates and to plot strand preference.

### Requirements
- Python 3.9+
  - Packages: `pandas`, `biopython`
- R 4.0+
  - Packages: `ggplot2`, `ggrepel`, `reticulate`, `seqinr`
  - Python package `numpy` available to R via `reticulate`
- CLI: `evo_gcp` (available in your PATH)

### Key inputs
- `input/variants.txt`: Variant and taxonomy table (includes `species`)
- `input/genes.txt`: Gene coordinates and metadata per subject (`aid`)
- `input/amplified_genes.txt`: Subject-to-gene mapping used to select targets
- `input/codon_table`: Tab-delimited `codon\taa` mapping used for translation checks (codon table 11)
- `input/<AID>_in/contigs.fa`: Source contigs per subject (e.g., `BAA_in/contigs.fa`)

### Scripts
- `scripts/process_human_data.py`
  - Joins `variants.txt`, `genes.txt`, and `amplified_genes.txt` to a per-target table; validates contig orientation and mutation positions; optionally writes a multi-FASTA of contigs.
  - Outputs: updated dataframe (TSV optional) and, if requested, FASTA via `generate_fasta_v1.py`.
- `scripts/FR_format_human_data.py`
  - Variant of the above with slightly different mutation/contig handling; useful if you need both forward and mutated sequences written.
- `scripts/generate_fasta_v1.py`
  - Helpers: `generate_fasta`, `generate_mult_fasta`, and `generate_all_83_variants` (utility to enumerate codon substitutions at a position).
- `scripts/generate_codons_variants.py`
  - Given a contig (`--fasta`, `--seq-id`) and gene region (`--gene-start`, `--gene-end`), writes a FASTA containing all 64 codon variants at an amino-acid coordinate (`--aa-coord`). Adds configurable left/right margins around the gene. Also writes a small TSV of the reference codon for highlighting in future plots.
  - Inputs: `input/codon_table`, contig FASTA. Outputs: `output/<aid>/query_<seq>_<pos>.fasta` and `.../query_<seq>_<pos>.tab`.
- `scripts/analyze_codon_position.sh`
  - Orchestrates a full run for one target: generate variant FASTA -> submit to `evo_gcp` -> download results -> build P vs M table -> plot.
  - Usage: `$ scripts/analyze_codon_position.sh POS GENE_START GENE_END SEQ_ID AID INPUT_FASTA [left_margin] [right_margin]`
  - Defaults: `left_margin=2000`, `right_margin=1000` (nucleotides)
  - Produces: `jobs/<aid>-<version>/output/input_summary.txt`, `output/<aid>/compare_strands_<seq>_<pos>.tab`, and `figures/<aid>/PvM_<AID>_pos_<POS>.pdf`.
- `run_mutation_analysis.sh`
  - Batch driver over `input/updated_data.tsv`: for each row, runs positions `mut_pos-2..mut_pos+2` via `scripts/analyze_codon_position.sh`.
  - Filters: optional `AID` and optional `POSITION` (requires `AID` if given).
  - Usage:
    - All aids and positions: `./run_mutation_analysis.sh`
    - All positions for one aid: `./run_mutation_analysis.sh BAA`
    - One position for one aid: `./run_mutation_analysis.sh BAA 154`
  - Notes: uses default margins `left_margin=2000` and `right_margin=1000`; edit the variables at the top of the script to change. Requires `input/updated_data.tsv` and `input/human_contigs_src.fasta` (see Workflow step 2).
- `scripts/create_strand_table.r`
  - Reads model `input_summary.txt` and creates a table of plus vs minus log-likelihood for each variant ID.
- `scripts/plot_strand_scatter.r`
  - Creates a labeled scatterplot of plus vs minus strand log-likelihoods; highlights the reference codon.
- `scripts/format_human_data.sh`
  - Thin wrapper to run `process_human_data.py` with defaults; demonstrates writing TSV and FASTA.
- `scripts/utils.r`
  - R helpers for reading logits from Evo2 via `reticulate` and lightweight plotting utilities.
- `testing/open_variants_table.R`
  - Tiny helper to preview `input/variants.txt` as a CSV in a spreadsheet app.

### Typical workflow
1) Prepare inputs
- Ensure `input/genes.txt`, `input/variants.txt`, `input/amplified_genes.txt`, `input/codon_table`, and `input/<AID>_in/contigs.fa` exist.

2) Validate gene/contig context (optional, recommended)
```bash
python scripts/process_human_data.py \
  --output-tsv input/updated_data.tsv \
  --write-fasta \
  --fasta-path input/human_contigs_src.fasta
```
- Confirms coordinates/orientation and prints source vs target codon at the mutation site.

3) Generate variants, run model, and analyze one target
```bash
# Example values
POS=152
GENE_START=41459
GENE_END=44200
SEQ_ID=k77_114799
AID=BAA
INPUT_FASTA=input/human_contigs_src.fasta

./scripts/analyze_codon_position.sh \
  "$POS" "$GENE_START" "$GENE_END" "$SEQ_ID" "$AID" "$INPUT_FASTA"
```
- Outputs:
  - `output/<aid>/query_<seq>_<pos>.fasta` and `.../query_<seq>_<pos>.tab`
  - `jobs/<aid>-<version>/output/input_summary.txt`
  - `output/<aid>/compare_strands_<seq>_<pos>.tab`
  - `figures/<aid>/PvM_<AID>_pos_<POS>.pdf`

3b) Batch run across aids and nearby positions
```bash
# All aids and +/-2 around each mutation
./run_mutation_analysis.sh

# All positions for one aid
./run_mutation_analysis.sh BAA

# Only position 154 for one aid
./run_mutation_analysis.sh BAA 154
```
- Consumes `input/updated_data.tsv` and uses `input/human_contigs_src.fasta`.
- Produces the same outputs as step 3, for each aid/position.

4) Inspect results
- Open the PDF under `figures/<aid>/` for plus vs minus strand comparison.
- Open the comparison table: `output/<aid>/compare_strands_<seq>_<pos>.tab`.

### Margins and plots
- Left/right margins
  - Defaults: `left_margin=2000` nt, `right_margin=1000` nt.
  - Purpose: provide genomic context flanking the gene so EVO 2 evaluates variants within realistic sequence neighborhoods.
  - How to change: pass as optional args to `scripts/analyze_codon_position.sh` or edit the variables in `run_mutation_analysis.sh`.
- Plot interpretation (`PvM_<AID>_pos_<POS>.pdf`)
  - Axes: X = total_log_likelihood on the plus strand; Y = total_log_likelihood on the minus strand.
  - Points: one per codon variant at the target amino-acid position; labels show `AA_CODON`.
  - Colors: reference codon highlighted in red; all others in blue.
  - Aspect: square (equal scales). Points near the diagonal indicate no strand preference; right of the diagonal favors plus; above the diagonal favors minus.

### Quick helpers
- Open `variants.txt` as a spreadsheet:
```bash
Rscript testing/open_variants_table.R
```
- View amplified genes mapping in-editor: open `input/amplified_genes.md`.
