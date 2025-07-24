#!/usr/bin/env python3

import os
import sys
import argparse
from pathlib import Path
# Now import the required packages
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from generate_fasta_v1 import generate_mult_fasta

def extract_gene_id(full_name):
    """Extract gene ID from full name (e.g., 'BAA_g34760_1' -> 'g34760_1')"""
    return full_name.split('_', 1)[1]

def get_amplified_genes():
    """Read amplified_genes.txt and return a dictionary with aid and gene info"""
    df = pd.read_csv('input/amplified_genes.txt', sep='\t')
    # Create a dictionary mapping aid (Subject) to gene info
    result = {}
    for _, row in df.iterrows():
        aid = row['Subject']
        full_gene = row['gene']
        gene_id = extract_gene_id(full_gene)
        result[aid] = {
            'full_gene': full_gene,
            'gene_id': gene_id
        }
    return result

def read_input_files():
    """Read variants and genes files"""
    variants = pd.read_csv("input/variants.txt", sep='\t')
    genes = pd.read_csv("input/genes.txt", sep='\t')
    return variants, genes

def make_df():
    """Create a new dataframe with gene information by joining relevant tables"""
    variants, genes = read_input_files()
    amplified_genes = get_amplified_genes()
    collected_data = []
    # print(variants.head())
    # print(genes.head())
    for aid, info in amplified_genes.items():
        gene_name = info['full_gene']
        formatted_name = extract_gene_id(gene_name)
        gene_rows = genes[(genes['aid'] == aid) & (genes['gene'] == formatted_name)]
        variants_rows = variants[(variants['aid'] == aid) & (variants['gene'] == gene_name)]
        gene_row_data = {
            'aid': aid,
            'gene': formatted_name,
            'contig': gene_rows['contig'].values[0],
            'start': gene_rows['start'].values[0],
            'end': gene_rows['end'].values[0],
            'strand': gene_rows['strand'].values[0],
            'flipped': variants_rows['flipped'].values[0],
            'src_codon': variants_rows['src_codon'].values[0],
            'tgt_codon': variants_rows['tgt_codon'].values[0],
            'mut_pos': variants_rows['start_dist'].values[0]}

        collected_data.append(gene_row_data)
    final_df = pd.DataFrame(collected_data)
    return final_df
        

def get_contig(fasta_file: str, contig: str) -> str:
    """Extract gene sequence from FASTA file based on contig and coordinates
    
    Args:
        fasta_file: Path to the FASTA file
        contig: Contig ID to search for
        start_pos: Start position in the contig
        end_pos: End position in the contig
    
    Returns:
        str: The extracted DNA sequence
    """
    # Read the FASTA file
    fasta_obj = SeqIO.parse(fasta_file, "fasta")
    
    # Search for and pull out contig sequence from FASTA
    contig_seq = None
    for record in fasta_obj:
        if record.id.startswith(contig):
            contig_seq = str(record.seq)
            break
    
    if contig_seq is None:
        raise ValueError(f"Contig {contig} not found in {fasta_file}")
    
    return contig_seq

# RC if needed and sets the source codon at the mutated position 
# Note contig_mut_pos is the position of the mutated codon in the contig, not the position of the mutated base
# This code might only work for mutation at the second position of the codon (works for given 8 human genes)
def check_contigs(contig_df: pd.DataFrame, write_fasta: bool = False, output_path: str = None):
    contigs = []
    for index, row in contig_df.iterrows():
        aid = row['aid']
        input_fasta = f"input/{aid}_in/contigs.fa"
        contig_seq = get_contig(input_fasta, row['contig'])
        ogg_start = row['start']
        ogg_end = row['end']
        mutg_pos = row['mut_pos']
        if row['strand'] == '-':
            contig_seq_obj = Seq(contig_seq)
            rev_contig_seq = contig_seq_obj.reverse_complement()
            gene_start = len(contig_seq) - ogg_end
            contig_mut_pos = gene_start + mutg_pos - 1
            contig_df.loc[index, 'start'] = gene_start + 1 # using base 1 indexing
            contig_df.loc[index, 'end'] = len(contig_seq) - ogg_start + 1 # using base 1 indexing
            working_contig_seq = str(rev_contig_seq)
        elif row['strand'] == '+':  
            contig_mut_pos = ogg_start + mutg_pos - 2
            working_contig_seq = contig_seq # already a string
        else:
            raise ValueError(f"Invalid strand: {row['strand']}")
        contig_df.loc[index, 'mut_pos'] = ((mutg_pos + 2 )// 3)
        desired_codon = row['src_codon']
        codon = working_contig_seq[contig_mut_pos:contig_mut_pos+3]
        # using source codon if flipped is true
        if row['flipped']:
            desired_codon = row['tgt_codon']
            changed_contig_seq = working_contig_seq[:contig_mut_pos] + row['src_codon'] + working_contig_seq[contig_mut_pos+3:]
        else:
            changed_contig_seq = working_contig_seq
        # check if the codon is correctly mutated to source
        if codon != desired_codon:
            print(f"Codon {codon} is incorrect")
            print(f"Expected: {desired_codon}, Actual: {codon}")
            
        # Append original contig
        contigs.append((changed_contig_seq, row['contig']))
        ####### Appending mutated contig ######### (commented out if want only forward strand)
        mut_contig_seq = changed_contig_seq[:contig_mut_pos] + row['tgt_codon'] + working_contig_seq[contig_mut_pos+3:]
        # print(f"Mutated codon: {mut_contig_seq[contig_mut_pos:contig_mut_pos+3]}")
        # print(f"Source codon: {changed_contig_seq[contig_mut_pos:contig_mut_pos+3]}")
        contigs.append((mut_contig_seq, f"{row['contig']}_mut"))
    if write_fasta:
        generate_mult_fasta(contigs, output_path = output_path)
    return(contig_df)

def main():
    parser = argparse.ArgumentParser(description='Process human data and optionally save outputs')
    parser.add_argument('--output-tsv', '-o', 
        help='Output TSV file path. If provided, will save the updated dataframe to this file.')
    parser.add_argument('--write-fasta', '-w', action='store_true',
        help='If set, will write the FASTA output')
    parser.add_argument('--fasta-path', '-f',
        help='Path for FASTA output (required if --write-fasta is set)')
    
    args = parser.parse_args()

    if args.write_fasta and not args.fasta_path:
        parser.error("--fasta-path is required when --write-fasta is set")
    # print(f'Original dataframe: {make_df()}')
    updated_df = check_contigs(make_df(), write_fasta=args.write_fasta,
                             output_path=args.fasta_path)
    # print(f'Updated dataframe: {updated_df}')
    if args.output_tsv:
        updated_df.to_csv(args.output_tsv, sep='\t', index=False)
        print(f"\nSaved updated dataframe to: {args.output_tsv}")
    
    if args.write_fasta:
        print(f"\nSaved FASTA sequences to: {args.fasta_path}\n")

if __name__ == "__main__":
    main()



