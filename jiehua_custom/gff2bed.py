#!/usr/bin/env python3
"""
gff2bed.py - Convert GFF file to BED file (only chrom/start/end)

Usage:
    gff2bed -i input.gff -o output.bed

Description:
    - GFF files with lines starting with '#' will be skipped
    - GFF format is 1-based inclusive, BED format is 0-based start, end-exclusive
    - The output BED file will contain only three columns: chrom, start, end
"""

import pandas as pd
import argparse
import sys
import os

# try to use tqdm for nicer progress bars, but fall back to simple text progress
try:
    from tqdm import tqdm
    HAS_TQDM = True
except Exception:
    HAS_TQDM = False

def gff_to_bed(input_file, output_file, chunksize=100000):
    """
    Convert GFF to BED keeping only chrom, start, end.
    Adjust start coordinate to 0-based for BED.

    This implementation reads the input in chunks so we can show a progress
    bar. If `tqdm` is installed it will be used; otherwise a simple text
    progress display is printed.
    """
    # Count non-comment lines to estimate total (used for progress)
    total = 0
    try:
        with open(input_file, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line and not line.startswith("#") and line.strip():
                    total += 1
    except Exception as e:
        print(f"Error opening input file for counting lines: {e}")
        sys.exit(1)

    # Use pandas chunked reader so we can process large files and update progress
    try:
        reader = pd.read_csv(
            input_file,
            sep="\t",
            comment="#",
            header=None,
            usecols=[0, 3, 4],
            dtype={0: str, 3: int, 4: int},
            chunksize=chunksize,
            iterator=True,
        )
    except Exception as e:
        print(f"Error reading GFF file: {e}")
        sys.exit(1)

    processed = 0
    first_chunk = True
    pbar = None
    if HAS_TQDM:
        pbar = tqdm(total=total, unit="line", desc="Converting", ncols=80)

    try:
        for chunk in reader:
            # Convert GFF 1-based start to BED 0-based start
            chunk[3] = chunk[3] - 1

            # write first chunk with mode 'w', then append
            mode = "w" if first_chunk else "a"
            chunk.to_csv(output_file, sep="\t", header=False, index=False, mode=mode)
            first_chunk = False

            processed += len(chunk)
            if pbar:
                pbar.update(len(chunk))
            else:
                # simple textual progress
                if total:
                    percent = processed / total * 100
                    print(f"Processed {processed}/{total} ({percent:.1f}%)", end="\r")
                else:
                    print(f"Processed {processed} lines", end="\r")
    except Exception as e:
        if pbar:
            pbar.close()
        print(f"\nError processing GFF file: {e}")
        sys.exit(1)

    if pbar:
        pbar.close()
    else:
        print()

    print(f"BED file saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Convert GFF to BED (chrom/start/end only, 0-based coordinates)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  gff2bed -i annotations.gff -o annotations.bed
  gff2bed -i genome.gff3 -o genome_features.bed
        """
    )
    parser.add_argument("-i", "--input", required=True, 
                       help="[REQUIRED] Input GFF file path")
    parser.add_argument("-o", "--output", required=True, 
                       help="[REQUIRED] Output BED file path")
    args = parser.parse_args()

    gff_to_bed(args.input, args.output)

if __name__ == "__main__":
    main()