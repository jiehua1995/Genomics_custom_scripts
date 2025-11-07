#!/usr/bin/env python3
"""
Setup script for jiehua_custom package

This setup.py configures the package installation and creates CLI entry points
for each bioinformatics tool so they can be run directly from the command line
after installation.
"""

from setuptools import setup, find_packages
import os

# Read version from __init__.py
def get_version():
    """Extract version from package __init__.py"""
    with open(os.path.join("jiehua_custom", "__init__.py")) as f:
        for line in f:
            if line.startswith("__version__"):
                return line.split("=")[1].strip().strip('"').strip("'")
    return "1.0.0"

# Read long description from README
def get_long_description():
    """Read long description from README.md if available"""
    try:
        with open("README.md", "r", encoding="utf-8") as f:
            return f.read()
    except FileNotFoundError:
        return "A collection of bioinformatics command-line utilities"

setup(
    # Package metadata
    name="jiehua_custom",
    version=get_version(),
    author="Jie Hua",
    author_email="Jie.Hua@lmu.de",
    description="A collection of bioinformatics command-line utilities",
    long_description=get_long_description(),
    long_description_content_type="text/markdown",
    url="https://github.com/jiehua1995/Genomics_custom_scripts",
    
    # Package configuration
    packages=find_packages(),
    include_package_data=True,
    
    # Python version requirement
    python_requires=">=3.8",
    
    # Dependencies - these will be installed automatically
    install_requires=[
        "pandas>=1.3.0",      # Required for gff2bed (GFF/BED file processing)
        "biopython>=1.78",    # Required for fasta_stats (FASTA file parsing)
        "tqdm>=4.60.0",       # Optional but recommended for progress bars
    ],
    
    # Optional dependencies that enhance functionality
    extras_require={
        "dev": [
            "pytest>=6.0",
            "black",
            "flake8",
        ],
    },
    
    # Entry points - these create the CLI commands
    # Each entry point maps a command name to a module:function
    entry_points={
        "console_scripts": [
            # Main package entry point - shows help and lists all tools
            "jiehua_custom=jiehua_custom.main:main",
            
            # gff2bed command - converts GFF files to BED format
            "gff2bed=jiehua_custom.gff2bed:main",
            
            # file-checksum command - calculates file hashes for integrity verification
            "file-checksum=jiehua_custom.file_checksum:main",
            
            # fasta-stats command - analyzes FASTA files for sequence statistics
            "fasta-stats=jiehua_custom.fasta_stats:main",
        ],
    },
    
    # Classification metadata for PyPI
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    
    # Keywords for package discovery
    keywords="bioinformatics genomics fasta gff bed hash checksum",
)