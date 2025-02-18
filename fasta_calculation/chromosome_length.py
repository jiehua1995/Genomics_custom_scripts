# -*- coding: utf-8 -*-
# Copyright (c) 2025 Jie Hua<Jie.Hua@lmu.de>
# Created Date: Tuesday, February 18th 2025, 11:55:48 am
# Compute sequence length and GC content from a FASTA file, with optional summarization.


import argparse
from Bio import SeqIO

def calculate_gc_content(sequence):
    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0

def calculate_nx(lengths, x):
    threshold = sum(lengths) * (x / 100.0)
    cumulative_length = 0
    for length in lengths:
        cumulative_length += length
        if cumulative_length >= threshold:
            return length
    return 0

def process_fasta(fasta_file, machine_readable, summarize, output_file):
    results = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_length = len(record.seq)
        gc_content = calculate_gc_content(record.seq)
        results.append((record.id, seq_length, gc_content))
    
    output_lines = []
    if summarize:
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
        output_lines.append(f"N10: {n10}\nN20: {n20}\nN30: {n30}\nN50: {n50}\nN70: {n70}\nN80: {n80}\nN90: {n90}")
        output_lines.append("Top 10 sequence lengths: " + ", ".join(map(str, top_10)))
    else:
        if machine_readable:
            output_lines.append("ID\tLength\tGC_Content")
            for seq_id, length, gc in results:
                output_lines.append(f"{seq_id}\t{length}\t{gc:.2f}")
        else:
            output_lines.append(f"{'ID':<20}{'Length':<10}{'GC_Content (%)':<15}")
            output_lines.append("-" * 45)
            for seq_id, length, gc in results:
                output_lines.append(f"{seq_id:<20}{length:<10}{gc:.2f}")
    
    if output_file:
        with open(output_file, "w") as f:
            f.write("\n".join(output_lines) + "\n")
    else:
        print("\n".join(output_lines))

def main():
    parser = argparse.ArgumentParser(description="Compute sequence length and GC content from a FASTA file.")
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument("-T", action="store_true", help="Enable machine-readable output (tab-separated)")
    parser.add_argument("-summarize", action="store_true", help="Output N10, N20, N30, N50, N70, N80, N90 statistics and top 10 longest sequences")
    parser.add_argument("-o", "--output", help="Save results to a file instead of printing to stdout")
    args = parser.parse_args()
    
    process_fasta(args.fasta, args.T, args.summarize, args.output)

if __name__ == "__main__":
    main()
