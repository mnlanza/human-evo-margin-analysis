import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import itertools

def generate_fasta(seq: str, output_path: str, seq_id: str, description: str = ""):
    """
    Generate a single sequence FASTA file.
    
    Args:
        seq: The sequence string
        output_path: Full path where to save the FASTA file (e.g., "path/to/output/file.fasta")
        seq_id: Identifier for the sequence
        description: Optional description for the sequence
    """
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    record = SeqRecord(Seq(seq), id=seq_id, description=description)
    SeqIO.write(record, output_path, "fasta")
    print(f"FASTA written to {output_path}")


def generate_mult_fasta(seqs, output_path: str):
    """
    Generate a multi-sequence FASTA file.
    
    Args:
        seqs: list of tuples (sequence_str, seq_id, [description])
        output_path: Full path where to save the FASTA file (e.g., "path/to/output/file.fasta")
    """
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    records = []
    for entry in seqs:
        if len(entry) == 3:
            seq_str, seq_id, description = entry
        elif len(entry) == 2:
            seq_str, seq_id = entry
            description = ""
        else:
            raise ValueError("Each item in `seqs` must be a tuple of (sequence_str, seq_id, [description])")

        record = SeqRecord(Seq(seq_str), id=seq_id, description=description)
        records.append(record)

    SeqIO.write(records, output_path, "fasta")
    print(f"Wrote {len(records)} sequences to {output_path}")

    """
    Seq formatting for generate_mult_fasta(): 
    
    seqs = [
        ("ATGCGT", "seq1", "example sequence 1"),
        ("ATGAAA", "seq2", "example sequence 2"),
        ("ATGTTT", "seq3", "example sequence 3"),
    ]
    """

def generate_all_83_variants(seq: str, output_path: str, codon_index: int = 83):
    """
    Generate all possible codon variants at a specified position (default: 83) in the sequence.
    
    Args:
        seq: The reference sequence string
        output_path: Full path where to save the FASTA file (e.g., "path/to/output/file.fasta")
        codon_index: The 1-based position of the codon to vary (default: 83)
    
    Returns:
        None. Writes a FASTA file containing all 64 codon variants.
    """
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Calculate codon position (0-based)
    start = (codon_index - 1) * 3
    end = start + 3
    
    # Generate all 64 possible codons
    nucleotides = ['A', 'C', 'G', 'T']
    codons = [''.join(c) for c in itertools.product(nucleotides, repeat=3)]
    
    # Create records with each codon substituted at the specified position
    records = []
    for codon in codons:
        modified = seq[:start] + codon + seq[end:]
        record = SeqRecord(
            Seq(modified), 
            id=f"{codon_index}_{codon}",
            description=""  # Empty description to match original format
        )
        records.append(record)
    
    # Write all variants to FASTA file
    SeqIO.write(records, output_path, "fasta")
    print(f"Wrote {len(records)} codon variants to {output_path}")



