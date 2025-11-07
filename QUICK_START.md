# Quick Start Guide for jiehua_custom

This guide shows you how to quickly install and use the jiehua_custom bioinformatics toolkit.

## Quick Installation

```bash
# Clone and install from GitHub
git clone https://github.com/jiehua1995/Genomics_custom_scripts.git
cd Genomics_custom_scripts
pip install -r requirements.txt
pip install -e .
```

## Quick Usage Examples

### 1. Get Help Information
```bash
# Show all available tools
jiehua_custom --help

# Get help for specific tools  
gff2bed --help
file-checksum --help
fasta-stats --help
```

### 2. Convert GFF to BED
```bash
# REQUIRED parameters: -i (input) and -o (output)
gff2bed -i annotations.gff -o annotations.bed
```

### 3. Calculate File Checksums
```bash
# REQUIRED: -p (path) and -o (output)
# OPTIONAL: -a (algorithm, default: sha256)
file-checksum -p ./data_directory -o checksums.txt

# Using different hash algorithm
file-checksum -p ./data_directory -o checksums.txt -a md5
```

### 4. Analyze FASTA Files
```bash
# REQUIRED: fasta file path
# All other parameters are OPTIONAL

# Basic analysis (human-readable)
fasta-stats genome.fasta

# Machine-readable output
fasta-stats genome.fasta -T

# Summary statistics only
fasta-stats genome.fasta -summarize  

# Save to file
fasta-stats genome.fasta -o results.txt
```

## Parameter Reference

| Tool | Required Parameters | Optional Parameters |
|------|-------------------|-------------------|
| **gff2bed** | `-i` (input GFF), `-o` (output BED) | None |
| **file-checksum** | `-p` (directory path), `-o` (output file) | `-a` (algorithm: md5/sha1/sha256/sha224/sha512) |
| **fasta-stats** | `fasta` (input file) | `-T` (machine-readable), `-summarize`, `-o` (output file) |

## License

Apache License 2.0