#!/usr/bin/env python3
"""
fasta_stats.py - Analyze FASTA files for sequence length and GC content statistics

Usage:
    fasta-stats input.fasta [-T] [-summarize] [-o output.txt]

Description:
    Compute sequence length and GC content from a FASTA file, with optional summarization.
    Can output in human-readable format or machine-readable tab-separated format.
    Provides comprehensive statistics including N50, N90 values and top longest sequences.
"""

import argparse

# Try to import BioPython, provide helpful error message if not available
try:
    from Bio import SeqIO
except ImportError:
    print("Error: BioPython is required for this tool.")
    print("Install it with: conda install biopython")
    exit(1)

def calculate_gc_content(sequence):
    """
    Calculate GC content percentage for a given DNA/RNA sequence.
    
    Parameters:
    - sequence: DNA/RNA sequence string
    
    Returns:
    - Float representing GC content percentage (0-100)
    """
    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0

def calculate_nx(lengths, x):
    """
    Calculate Nx statistic (e.g., N50, N90) for a list of sequence lengths.
    
    Parameters:
    - lengths: List of sequence lengths (sorted in descending order)
    - x: Percentage threshold (e.g., 50 for N50, 90 for N90)
    
    Returns:
    - Length value at which x% of total bases are contained in sequences of this length or longer
    """
    threshold = sum(lengths) * (x / 100.0)
    cumulative_length = 0
    for length in lengths:
        cumulative_length += length
        if cumulative_length >= threshold:
            return length
    return 0

def process_fasta(fasta_file, machine_readable, summarize, output_file):
    """
    Process FASTA file and calculate statistics for each sequence.
    
    Parameters:
    - fasta_file: Path to input FASTA file
    - machine_readable: Boolean, output in tab-separated format if True
    - summarize: Boolean, output summary statistics instead of per-sequence data
    - output_file: Path to output file (None to print to stdout)
    """
    results = []
    
    # Parse FASTA file and calculate stats for each sequence
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_length = len(record.seq)
        gc_content = calculate_gc_content(record.seq)
        results.append((record.id, seq_length, gc_content))
    
    output_lines = []
    
    if summarize:
        # Generate summary statistics
        lengths = sorted([length for _, length, _ in results], reverse=True)
        n10 = calculate_nx(lengths, 10)
        n20 = calculate_nx(lengths, 20)
        n30 = calculate_nx(lengths, 30)
        n50 = calculate_nx(lengths, 50)
        n70 = calculate_nx(lengths, 70)
        n80 = calculate_nx(lengths, 80)
        n90 = calculate_nx(lengths, 90)
        top_10 = lengths[:10]
        
        output_lines.append("Summary Statistics:")
        output_lines.append(f"Total sequences: {len(results)}")
        output_lines.append(f"Total bases: {sum(lengths):,}")
        output_lines.append(f"Average length: {sum(lengths) / len(lengths):.2f}")
        output_lines.append(f"Median length: {lengths[len(lengths)//2]:,}")
        output_lines.append("")
        output_lines.append("Nx Statistics:")
        output_lines.append(f"N10: {n10:,}\nN20: {n20:,}\nN30: {n30:,}\nN50: {n50:,}")
        output_lines.append(f"N70: {n70:,}\nN80: {n80:,}\nN90: {n90:,}")
        output_lines.append("")
        output_lines.append("Top 10 sequence lengths: " + ", ".join(f"{x:,}" for x in top_10))
    else:
        # Generate per-sequence output
        if machine_readable:
            # Tab-separated output for machine processing
            output_lines.append("ID\tLength\tGC_Content")
            for seq_id, length, gc in results:
                output_lines.append(f"{seq_id}\t{length}\t{gc:.2f}")
        else:
            # Human-readable formatted output
            output_lines.append(f"{'ID':<20}{'Length':<10}{'GC_Content (%)':<15}")
            output_lines.append("-" * 45)
            for seq_id, length, gc in results:
                output_lines.append(f"{seq_id:<20}{length:<10}{gc:.2f}")
    
    # Write output to file or print to stdout
    if output_file:
        with open(output_file, "w") as f:
            f.write("\n".join(output_lines) + "\n")
        print(f"Results written to {output_file}")
    else:
        print("\n".join(output_lines))

def main():
    """
    Main function to handle command-line arguments and execute FASTA analysis.
    """
    parser = argparse.ArgumentParser(
        description="Compute sequence length and GC content from a FASTA file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  fasta-stats sequences.fasta                    # Human-readable output
  fasta-stats sequences.fasta -T                 # Machine-readable output  
  fasta-stats sequences.fasta -summarize         # Summary statistics only
  fasta-stats sequences.fasta -o results.txt     # Save to file
  fasta-stats genome.fasta -T -summarize -o stats.tsv  # Combined options
        """
    )
    parser.add_argument("fasta", 
                       help="[REQUIRED] Input FASTA file path")
    parser.add_argument("-T", action="store_true", 
                       help="[OPTIONAL] Enable machine-readable output (tab-separated)")
    parser.add_argument("-summarize", action="store_true", 
                       help="[OPTIONAL] Output Nx statistics (N10-N90) and top 10 longest sequences")
    parser.add_argument("-o", "--output", 
                       help="[OPTIONAL] Save results to a file instead of printing to stdout")
    
    args = parser.parse_args()
    
    # Execute FASTA processing with parsed arguments
    process_fasta(args.fasta, args.T, args.summarize, args.output)

if __name__ == "__main__":
    main()