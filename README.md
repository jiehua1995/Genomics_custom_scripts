# jiehua_custom

A collection of bioinformatics command-line utilities for common genomics tasks.

## Overview

This package provides three standalone command-line tools that can be installed via conda or pip:

- **`gff2bed`** - Convert GFF files to BED format with progress tracking
- **`file-checksum`** - Calculate file hashes for data integrity verification  
- **`fasta-stats`** - Analyze FASTA files for sequence length and GC content statistics

Each tool is designed to be independent and can handle large files efficiently with progress indicators.

## Installation

### From GitHub (Recommended)
```bash
# Clone the repository
git clone https://github.com/jiehua1995/Genomics_custom_scripts.git
cd Genomics_custom_scripts

# Install dependencies
pip install -r requirements.txt

# Install the package in development mode
pip install -e .
```

### Alternative: Direct pip install from GitHub
```bash
pip install git+https://github.com/jiehua1995/Genomics_custom_scripts.git
```

### Requirements
- **Python 3.8+**
- **pandas** (≥1.3.0) - for GFF/BED file processing
- **biopython** (≥1.78) - for FASTA file parsing  
- **tqdm** (≥4.60.0) - for progress bars

## Tools

### gff2bed
Convert GFF files to BED format with coordinate adjustment (1-based to 0-based).

**Usage:**
```bash
gff2bed -i input.gff -o output.bed
```

**Parameters:**
- `-i, --input` (REQUIRED): Input GFF file path
- `-o, --output` (REQUIRED): Output BED file path

**Features:**
- Handles large files with chunked processing and progress tracking
- Automatically skips comment lines (starting with #)
- Converts GFF 1-based coordinates to BED 0-based format
- Shows progress bar if tqdm is available

**Example:**
```bash
# Convert annotations from GFF to BED format
gff2bed -i genome_annotations.gff -o genome_annotations.bed
```

### file-checksum  
Calculate hash values for all files in a directory recursively.

**Usage:**
```bash
file-checksum -p /path/to/directory -o checksums.txt -a sha256
```

**Parameters:**
- `-p, --path` (REQUIRED): Directory path to scan for files (recursive)
- `-o, --output` (REQUIRED): Output file to save hash results
- `-a, --algorithm` (OPTIONAL): Hash algorithm to use
  - **Default:** `sha256`
  - **Options:** `md5`, `sha1`, `sha256`, `sha224`, `sha512`

**Features:**
- Outputs tab-separated file with relative file paths and hash values
- Progress tracking for large directories (shows files processed)
- Recursive directory scanning
- Uses relative paths to avoid absolute path exposure

**Examples:**
```bash
# Calculate SHA256 hashes for all files in a directory
file-checksum -p ./raw_data -o data_integrity.txt

# Use MD5 algorithm instead
file-checksum -p ./sequences -o md5_checksums.txt -a md5

# Verify large genomics dataset integrity
file-checksum -p /path/to/genome_data -o genome_checksums.txt -a sha512
```

### fasta-stats
Analyze FASTA files for sequence statistics and quality metrics.

**Usage:**
```bash
# Basic per-sequence statistics (human-readable)
fasta-stats sequences.fasta

# Machine-readable tab-separated output
fasta-stats sequences.fasta -T

# Summary statistics only (N50, N90, etc.)
fasta-stats sequences.fasta -summarize

# Save results to file
fasta-stats sequences.fasta -o results.txt
```

**Parameters:**
- `fasta` (REQUIRED): Input FASTA file path
- `-T` (OPTIONAL): Enable machine-readable tab-separated output format
- `-summarize` (OPTIONAL): Output summary statistics instead of per-sequence data
- `-o, --output` (OPTIONAL): Output file path (default: print to stdout)

**Output Formats:**
- **Default**: Human-readable table with sequence ID, length, and GC content
- **Machine-readable (-T)**: Tab-separated values for easy parsing
- **Summary (-summarize)**: Nx statistics (N10, N20, N30, N50, N70, N80, N90) and top 10 longest sequences

**Examples:**
```bash
# Basic analysis - view sequence stats in terminal
fasta-stats genome.fasta

# Generate machine-readable output for downstream analysis
fasta-stats contigs.fasta -T -o contig_stats.tsv

# Get assembly quality metrics
fasta-stats assembly.fasta -summarize -o assembly_summary.txt

# Combined: detailed stats saved to file
fasta-stats sequences.fasta -T -o detailed_stats.tsv
```

**Summary Statistics Include:**
- Total sequences and bases
- Average and median sequence lengths  
- N10, N20, N30, N50, N70, N80, N90 values
- Top 10 longest sequences

## Dependencies

- **Python 3.8+**
- **pandas** (≥1.3.0) - for GFF/BED file processing
- **biopython** (≥1.78) - for FASTA file parsing  
- **tqdm** (≥4.60.0) - for progress bars (optional but recommended)

## Examples

### Complete workflow example
```bash
# 1. Calculate checksums for raw data integrity
file-checksum -p ./raw_data -o data_checksums.txt

# 2. Analyze FASTA file quality
fasta-stats genome.fasta -summarize -o genome_stats.txt

# 3. Convert annotation from GFF to BED
gff2bed -i annotations.gff -o annotations.bed
```

## Development

### Running tests
```bash
pip install -e ".[dev]"
pytest tests/
```

### Code formatting
```bash
black jiehua_custom/
flake8 jiehua_custom/
```

## License

Apache License 2.0 - see LICENSE file for details.

## Author

**Jie Hua**  
Email: Jie.Hua@lmu.de  
GitHub: [@jiehua1995](https://github.com/jiehua1995)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Changelog

### v1.0.0
- Initial release with three core tools
- gff2bed: GFF to BED conversion with progress tracking
- file-checksum: Directory hash calculation with multiple algorithms
- fasta-stats: Comprehensive FASTA sequence analysis
