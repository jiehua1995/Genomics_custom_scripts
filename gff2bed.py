#!/usr/bin/env python3
"""
gff2bed.py - Convert GFF file to BED file (only chrom/start/end)

Usage:
    python gff2bed.py -i input.gff -o output.bed

Description:
    - GFF files with lines starting with '#' will be skipped
    - GFF format is 1-based inclusive, BED format is 0-based start, end-exclusive
    - The output BED file will contain only three columns: chrom, start, end
"""

import pandas as pd
import argparse
import sys

def gff_to_bed(input_file, output_file):
    """
    Convert GFF to BED keeping only chrom, start, end.
    Adjust start coordinate to 0-based for BED.
    """
    try:
        gff = pd.read_csv(
            input_file,
            sep="\t",
            comment="#",       # skip comment lines
            header=None,
            usecols=[0, 3, 4], # columns: chrom, start, end
            dtype={0: str, 3: int, 4: int}
        )
    except Exception as e:
        print(f"Error reading GFF file: {e}")
        sys.exit(1)

    # Convert GFF 1-based start to BED 0-based start
    gff[3] = gff[3] - 1

    # Save to BED file
    try:
        gff.to_csv(output_file, sep="\t", header=False, index=False)
    except Exception as e:
        print(f"Error writing BED file: {e}")
        sys.exit(1)

    print(f"BED file saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Convert GFF to BED (chrom/start/end only, 0-based coordinates)"
    )
    parser.add_argument("-i", "--input", required=True, help="Input GFF file")
    parser.add_argument("-o", "--output", required=True, help="Output BED file")
    args = parser.parse_args()

    gff_to_bed(args.input, args.output)

if __name__ == "__main__":
    main()
